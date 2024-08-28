
import click
import pandas as pd
import os
from multiprocessing import Pool
import pyBigWig
from hatchet.utils.combine_counts import compute_baf_task_multi
import numpy as np
from datetime import datetime
import tqdm
import pysam 

def chr2key(chr):
    if chr.startswith('chr'):
        tkn = chr[3:]
    if tkn == 'X':
        return 23
    elif tkn == 'Y':
        return 24
    else:
        return int(tkn)

def adaptive_bin_ont(bams, bigwigs, cent_starts, cent_ends, 
                     interval=10000, min_reads=1000, threads=1):
    """
        Construct adaptive bins across chromosome arms for ONT data
        params:
        -reads: DataFrame containing individual reads 
            with at least "chr", "haplogroup" (phase set), and "haplotype" columns 
        -bigwig_filename: path to deeptools.bamCoverage results
        -cent_starts: pandas series containing centromere start locations
        -cent_ends: pandas series containing centromere end locations
        -interval (int): number of bases between candidate thresholds to consider
        -min_reads (int): minimum total number of reads per bin
    """

    chromosomes = sorted(cent_starts.index, key = chr2key)
    
    # remove sex chromosomes since HATCHet2 doesn't work on them
    chromosomes = [ch for ch in chromosomes if not (ch.endswith('X') or ch.endswith('Y'))]


    bws = [pyBigWig.open(bw) for bw in bigwigs]
    chromosome_ends = {ch:max([bw.chroms()[str(ch)] for bw in bws]) for ch in chromosomes}

    adaptive_bins = []
    all_params = []

    for ch in chromosomes:
        if ch.endswith('X') or ch.endswith('Y'):
            continue
        
        cent_start = cent_starts.loc[ch, 'start']
        cent_end = cent_ends.loc[ch, 'end'] + 1
        chr_end = chromosome_ends[ch]
            
        for arm_start, arm_end in [(1, cent_start), (cent_end, chr_end)]:
            if arm_end <= arm_start:
                assert ch in set(['13', '14', '15', '21', '22',
                                  'chr13', 'chr14', 'chr15', 'chr21', 'chr22'])
                continue
            all_params.append((bams, bigwigs, ch, arm_start, arm_end, interval, min_reads))

    if threads > 1:
        with Pool(threads) as p:
            result = p.map(arm_bin_wrapper, all_params)
    else:
        result = [arm_bin_wrapper(p) for p in all_params]
    [adaptive_bins.extend(l) for l in result]
    print([len(a) for a in result])
    
    df = pd.DataFrame(adaptive_bins)
    df['RD'] = df.bamcov / df.bamcov.mean()

    print(datetime.now(), 'Done binning')
    return df

def call_em(records, chromosome, start, end, polling_interval=int(1e4), 
            min_counts_per_sample=1):
    r = []

    for i, sample_df in enumerate(records):
        sample_df = sample_df[(sample_df.chr == chromosome) & (
            (sample_df.start >= start) & (sample_df.start < end) | 
            (sample_df.end >= start) & (sample_df.end < end) | 
            (sample_df.start < start) & (sample_df.end >= end))].copy()
        for ps, df in sample_df.groupby('haplogroup'):
            pos = start
            
            while pos < end:
                my_reads = df[(df.start <= pos) & (df.end > pos)]
                tot = len(my_reads)
                alt = len(my_reads[my_reads.haplotype == 1])
                ref = tot - alt
                r.append({'#CHR':chromosome, 'POS':f'{ps}_{pos}', 
                          'SAMPLE':i, 'REF':ref, 'ALT':alt, 'TOTAL':tot})
                
                pos += polling_interval
    input_df = pd.DataFrame(r)
    counts_by_sample = input_df.pivot(index='SAMPLE', columns='POS', values='TOTAL')
    valid_pos = counts_by_sample.columns[
    np.where(np.all(counts_by_sample.values >= min_counts_per_sample, axis = 0))[0]]
    input_df = input_df[input_df.POS.isin(valid_pos)].sort_values(
        by=['POS', 'SAMPLE']).reset_index(drop=True)
    b = compute_baf_task_multi(input_df, 0, 0, 0)
    return counts_by_sample.shape[1], b

def arm_bin_wrapper(params):
    return adaptive_bin_arm_ont(*params)

def adaptive_bin_arm_ont(bams, bigwigs, chromosome, arm_start, arm_end, interval=10000, min_reads=1000):
    """
        Construct adaptive bins across chromosome arms for ONT data
        params:
        -reads: DataFrame containing individual reads 
            with at least "chr", "haplogroup" (phase set), and "haplotype" columns 
        -bigwig_filename: path to deeptools.bamCoverage results
        -chromosome: chromosome to bin
        -arm_start: start position of chromosome arm
        -arm_end: end position of chromosome arm
        -interval (int): number of bases between candidate thresholds to consider
        -min_reads (int): minimum total number of reads per bin
    """

    # load haplotagged reads
    records = []
    for i, bam in enumerate(sorted(bams)):
        print('loading sample', bam)
        df = load_reads(bam, chromosome, arm_start, arm_end)
        df['startbin'] = (df.start // interval).astype(int)
        df['endbin'] = (df.end // interval).astype(int)
        df['sample'] = f'sample{i}'
        records.append(df)
    n_samples = len(bams)

    bin_min = min([df.startbin.min() for df in records])
    bin_max = max([df.endbin.max() for df in records])

    # count the number of reads that overlap each interval
    interval_counts = np.zeros((n_samples, bin_max - bin_min + 1), int)

    for i, df in enumerate(records):
        for _, r in tqdm.tqdm(df.iterrows()):
            interval_counts[i, r.startbin - bin_min:r.endbin - bin_min + 1] += 1

    # count othe number of reads that start in each interval
    interval_starts = np.zeros((n_samples, bin_max - bin_min + 1), int)
    for i, df in enumerate(records):
        vc = df.startbin.value_counts()
        interval_starts[i, vc.index - bin_min] = vc.values

    # load coverage ratios
    bws = [pyBigWig.open(bigwig) for bigwig in bigwigs]
    
    arm_bins = []

    # interval indices
    my_start = int(arm_start // interval)
    first_start = my_start

    curr_end = my_start + 1
    curr_reads = interval_counts[:, 0].copy()

    while curr_end < bin_max:    
        if np.all(curr_reads > min_reads):
            bin_start = my_start * interval if my_start != first_start else arm_start
            bin_end = curr_end * interval
            
            # draw bin boundary and add record
            n_phase_groups, em_result = call_em(records, chromosome, bin_start, bin_end)
            for s in range(n_samples):
                try:
                    _, _, mhf, alpha, beta = em_result[s]
                except KeyError as e:
                    print(em_result.keys())
                    raise e
                arm_bins.append({'chr':chromosome, 'start':bin_start, 'end':bin_end, 'sample':f'sample{s}',
                                'bamcov':bws[s].stats(chromosome, bin_start, bin_end)[0],
                                'n_phase_groups':n_phase_groups, 'mhf':mhf, 'alpha':alpha, 'beta':beta})

            my_start = curr_end
            curr_end = my_start + 1
            curr_reads = interval_counts[:, my_start - bin_min].copy()

        else:
            # add new reads to this bin
            curr_reads += interval_starts[:, curr_end - bin_min]
            curr_end += 1     

    if curr_end * interval < arm_end:
        if len(arm_bins) > 0:
            # remove previous bin boundary and recompute
            penultimate_bin = arm_bins[-1]
            arm_bins = arm_bins[:-1 * n_samples] 
            
            bin_start = penultimate_bin['start'] 
        else:
            bin_start = arm_start

        bin_end = arm_end

        # draw bin boundary and add record
        n_phase_groups, em_result = call_em(records, chromosome, bin_start, bin_end)
        for s in range(n_samples):
            try:
                _, _, mhf, alpha, beta = em_result[s]
            except KeyError as e:
                print(em_result.keys())
                raise e
            arm_bins.append({'chr':chromosome, 'start':bin_start, 'end':bin_end, 'sample':f'sample{s}',
                            'bamcov':bws[s].stats(chromosome, bin_start, bin_end)[0],
                            'n_phase_groups':n_phase_groups, 'mhf':mhf, 'alpha':alpha, 'beta':beta})

    return arm_bins

def load_reads(bamfile, chromosome, start, end):
    records = []

    samfile = pysam.AlignmentFile(bamfile, "rb")
    for line in tqdm.tqdm(samfile.fetch(chromosome, start, end)):
        d = line.to_dict()
        if int(d['map_quality']) > 30:
            if any([a.startswith('PS') for a in d['tags']]):
                records.append([d['name'], d['ref_name'], d['ref_pos'], d['map_quality'],
                                len(d['seq']), d['tags']])
    assert all([x[-1][-1].startswith('PS:i:') for x in records])
    haplotags = [x[-1][-1][5:] for x in records]
    haplotype = [x[-1][-3][5:] for x in records]

    df = pd.DataFrame([a[:-1] for a in records],
                columns = ['read_name', 'chr', 'start', 'qual', 'length'])
    df.start = df.start.astype(int)
    df.qual = df.qual.astype(int)
    df['end'] = df.start + df.length
    df['haplogroup'] = haplotags
    df.haplogroup = df.haplogroup.astype(int)
    df['haplotype'] = haplotype
    df.haplotype = df.haplotype.astype(int)
    return df

@click.command()
@click.option('--bams', help='Haplotagged BAM files (colon-separated)')
@click.option('--bigwigs', help='BigWig files containing bamCoverage results from deepTools (colon_separated)')
@click.option('--outfile', help='output filename')
@click.option('--centromere_starts', help='CSV file containing centromere start positions')
@click.option('--centromere_ends', help='CSV file containing centromere end positions')
@click.option('--mtr', help="Minimum number of total reads per bin", type=int, default=1000)
@click.option('--interval', help="Number of bases between candidate bin thresholds", type=int, default=10000)
@click.option('--threads', help="Number of worker threads", type=int, default=1)
def main(bams, bigwigs, outfile, centromere_starts, centromere_ends, mtr, interval, threads):
    bams = bams.split(':')
    bigwigs = bigwigs.split(':')

    assert len(bams) == len(bigwigs), 'number of BAM files and number of bigwig files must match - 1 per sample'
    for bam in bams:
        assert os.path.exists(bam), "missing haplotagged BAM file"
    for bigwig in bigwigs:
        assert os.path.exists(bigwig), "missing bamCoverage bigwig file"
    assert mtr > 0, 'Min total reads must be positive'
    assert interval > 0, 'Interval between candidate bin thresholds must be positive'
    
    centromere_starts = pd.read_csv(centromere_starts, dtype={'chr':'str'}).set_index('chr')
    centromere_ends = pd.read_csv(centromere_ends, dtype={'chr':'str'}).set_index('chr')

    result = adaptive_bin_ont(bams, bigwigs, centromere_starts, centromere_ends,
                              interval=interval, min_reads=mtr, threads=threads)
    result = result.rename(columns={'chr':'#CHR', 'start':'START', 'end':'END', 'sample':'SAMPLE',
                       'bamcov':'COV', 'alpha':'ALPHA', 'beta':'BETA', 'mhf':'BAF',
                       'n_phase_groups':'#SNPS'})
    result = result[['#CHR', 'START', 'END', 'SAMPLE', 'RD', '#SNPS', 'COV', 'ALPHA', 'BETA', 'BAF']]
    result.to_csv(outfile, index=False, sep='\t')

if __name__ == '__main__':
    main()