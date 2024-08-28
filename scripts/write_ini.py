
import click
import pandas as pd

@click.command()
@click.option('--work_dir', help='Working directory for HATCHet2 run')
@click.option('--tumor_1bed', help='Path to pre-clustering table')
@click.option('--ini_filename', help='Filename for output hatchet.ini files')
@click.option('--processes', help="Number of concurrent processes", default=1)
@click.option('--min_clones', help="Minimum number of clones for solver", default=2)
@click.option('--max_clones', help="Maximum number of clones for solver", default=4)
@click.option('--maxcn_diploid', help="Minimum copy number for diploid solutions", default=6)
@click.option('--maxcn_tetraploid', help="Maximum number of clones for tetraploid solutions", default=14)
def main(work_dir, tumor_1bed, ini_filename, min_clones, max_clones, maxcn_diploid, maxcn_tetraploid,
         processes):
    assert min_clones <= max_clones, (min_clones, max_clones)
    assert min_clones > 0, min_clones
    assert maxcn_diploid >= 2, maxcn_diploid
    assert maxcn_tetraploid >= 4, maxcn_tetraploid
    df = pd.read_table(tumor_1bed)    
    
    with open(ini_filename, 'w') as f:
        f.write('[run]\n')
        f.write('download_panel=False\n')
        f.write('count_reads=False\n')
        f.write('genotype_snps=False\n')
        f.write('phase_snps = False\n')
        f.write('fixed_width = False\n')
        f.write('count_alleles=False\n')
        f.write('combine_counts=False\n')
        f.write('cluster_bins=True\n')
        f.write('loc_clust = True\n')
        f.write('plot_bins=True\n')
        f.write('compute_cn=True\n')
        f.write('plot_cn=True\n\n')

        f.write('reference=/data1/shahs3/users/myersm2/reference/GRCh37-lite.fa\n')
        f.write('reference_version=hg19\n')
        f.write(f'processes={processes}\n')
        f.write(f'samples={" ".join(sorted(df.SAMPLE.unique()))}\n')
        f.write(f'output={work_dir}\n\n')
        
        f.write('[cluster_bins]\n')
        f.write('minK=5\n')
        f.write('maxK=20\n')
        f.write('diploidbaf=0.04\n')
        f.write('tau=1e-12\n\n')

        f.write('[compute_cn]\n')
        f.write('solver=cpp\n')
        f.write(f'diploidcmax={maxcn_diploid}\n')
        f.write(f'tetraploidcmax={maxcn_tetraploid}\n')
        f.write(f'clones={min_clones+1},{max_clones+1}\n')
        
if __name__ == '__main__':
    main()