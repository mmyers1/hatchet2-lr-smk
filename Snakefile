configfile: "config.yaml"
import pandas as pd

data_dir = config['data_dir']
results_dir = config['results_dir']
repo_dir =  config['repo_dir']
ref_genome =  config['ref_genome']

cohort =  pd.read_csv(config['cohort'])
centromere_starts_file = config['centromere_starts']
centromere_ends_file = config['centromere_ends']

min_reads_per_bin = config['min_reads_per_bin']
bin_interval = config['bin_interval']
threads_per_job = config['threads_per_job']
min_clones = config['min_clones']
max_clones = config['max_clones']
maxcn_diploid = config['maxcn_diploid']
maxcn_tetraploid = config['maxcn_tetraploid']
min_mapq = config['min_mapq']


patients = sorted(cohort.patient_id.values)
patients.remove('P-0063815')

def get_orig_vcf(wildcards):
    return cohort.set_index('patient_id').loc[wildcards.patient, 'phased_vcf']

def get_haplotag_inputs(wildcards):
    bams = cohort.set_index('patient_id').loc[wildcards.patient, 'ont_tumor_bam'].split(':')
    if len(bams) == 1:
        my_bam = bams[0]
    else:
        my_bam = bams[int(wildcards.sample[6:])]
    return {
        'bam':my_bam,
        'vcf':os.path.join(data_dir, wildcards.patient, 'phased.vcf.gz'),
        'tbi':os.path.join(data_dir, wildcards.patient, 'phased.vcf.gz.tbi')
    }

def get_prepare_hatchet_inputs(wildcards):
    n_bams = len(cohort.set_index('patient_id').loc[wildcards.patient, 'ont_tumor_bam'].split(':'))
    bams = [os.path.join(data_dir, f'{wildcards.patient}/{wildcards.patient}_sample{s}_haplotagged_filtered.bam') for s in range(n_bams)]
    bais = [os.path.join(data_dir, f'{wildcards.patient}/{wildcards.patient}_sample{s}_haplotagged_filtered.bam.bai') for s in range(n_bams)]
    bigwigs = [os.path.join(data_dir, f'{wildcards.patient}/{wildcards.patient}_sample{s}_haplotagged.bw') for s in range(n_bams)]
    return {
        'bams':bams,
        'bais':bais,
        'bigwigs':bigwigs
    }

def get_bams_string(wildcards):
    n_bams = len(cohort.set_index('patient_id').loc[wildcards.patient, 'ont_tumor_bam'].split(':'))
    return ':'.join(
        [os.path.join(data_dir, f'{wildcards.patient}/{wildcards.patient}_sample{s}_haplotagged_filtered.bam') for s in range(n_bams)])
    
def get_bigwigs_string(wildcards):
    n_bams = len(cohort.set_index('patient_id').loc[wildcards.patient, 'ont_tumor_bam'].split(':'))
    return ':'.join(
        [os.path.join(data_dir, f'{wildcards.patient}/{wildcards.patient}_sample{s}_haplotagged.bw') for s in range(n_bams)])
    

rule all:
    input:
        #expand(os.path.join(results_dir,  'hatchet2', '{patient}', 'bb', 'bulk.bb'), patient=patients)
        expand(os.path.join(results_dir,  'hatchet2', '{patient}', 'results', 'best.bbc.ucn'), patient=patients)

rule compress_vcf:
    input: get_orig_vcf
    output: 
        os.path.join(data_dir, '{patient}', 'phased.vcf.gz'),
    resources:
        time="0:59:59",
        mem_mb=8000
    shell: """
    bgzip {input} -c > {output}
    """

rule index_vcf:
    input: '{filename}.vcf.gz'
    output: '{filename}.vcf.gz.tbi'
    resources:
        time="1:59:59",
        mem_mb=32000
    shell: """
    tabix {input}
    """

rule sort_index_bam:
    input: '{bamfile}'
    output: '{bamfile}.bai'
    resources:
        time="5:59:59",
        mem_mb=32000
    shell: """
    samtools sort {input}
    samtools index {input}
    """

rule haplotag:
    input: 
        unpack(get_haplotag_inputs)
    threads: 16
    resources:
        time="11:59:59",
        mem_mb=32000
    output:
        os.path.join(data_dir, '{patient}/{patient}_{sample}_haplotagged.bam')
    params:
        reference=ref_genome
    shell:"""
    whatshap haplotag --ignore-read-groups -o {output} --reference {params.reference} {input.vcf} {input.bam} --output-threads=16
    """

rule bamcov:
    input: 
        bam=os.path.join(data_dir, '{patient}/{patient}_{sample}_haplotagged.bam'),
        bai=os.path.join(data_dir, '{patient}/{patient}_{sample}_haplotagged.bam.bai'),
    output: 
        bigwig=os.path.join(data_dir, '{patient}/{patient}_{sample}_haplotagged.bw')
    resources:
        time="11:59:59",
        mem_mb=64000
    shell:"""
    bamCoverage --bam {input.bam} -o {output.bigwig} --binSize 1000 --ignoreForNormalization chrX chrY
    """

rule filter_bam:
    input: 
        bam=os.path.join(data_dir, '{patient}/{patient}_{sample}_haplotagged.bam'),
        bai=os.path.join(data_dir, '{patient}/{patient}_{sample}_haplotagged.bam.bai'),
    output: 
        os.path.join(data_dir, '{patient}/{patient}_{sample}_haplotagged_filtered.bam'),
    threads: 16
    resources:
        time="11:59:59",
        mem_mb=64000
    shell:"""
    samtools view {input.bam} -q {min_mapq} -d PS -b -@ 16 > {output}
    """


rule prepare_hatchet_input:
    input:
        unpack(get_prepare_hatchet_inputs),
        centromere_starts=centromere_starts_file,
        centromere_ends=centromere_ends_file,
    output:
        os.path.join(results_dir, 'hatchet2', '{patient}', 'bb', 'bulk.bb'),
    params:
        bam_str=get_bams_string,
        bigwig_str=get_bigwigs_string,
        mtr=min_reads_per_bin,
        interval=bin_interval,
        outdir=os.path.join(results_dir, 'hatchet2', '{patient}', 'bb'),
        threads=threads_per_job
    threads: threads_per_job
    resources:
        time="5:59:59",
        mem_mb=32000,
    shell: """
    mkdir -p {params.outdir}
    python scripts/construct_hatchet_bins.py --bams {params.bam_str} --bigwigs {params.bigwig_str} --outfile {output} \
    --centromere_starts {input.centromere_starts} --centromere_ends {input.centromere_ends} \
    --mtr {params.mtr} --interval {params.interval} --threads {params.threads}
    """

rule write_ini:
    input:
        os.path.join(results_dir, 'hatchet2', '{patient}', 'bb', 'bulk.bb'),
    output:
        os.path.join(results_dir, "hatchet2", '{patient}',  "{patient}.ini")
    params:
        min_clones=min_clones,
        max_clones=max_clones,
        maxcn_diploid=maxcn_diploid,
        maxcn_tetraploid=maxcn_tetraploid,
        outdir=os.path.join(results_dir, 'hatchet2', '{patient}'),
        processes=threads_per_job
    resources:
        time="0:10:00",
        mem_mb=1000,
    shell: """
    mkdir -p {params.outdir}
    python scripts/write_ini.py --tumor_1bed {input} --work_dir {params.outdir} --ini_filename {output} \
    --min_clones {params.min_clones} --max_clones {params.max_clones} --maxcn_diploid {params.maxcn_diploid} \
    --maxcn_tetraploid {params.maxcn_tetraploid} --processes {params.processes}
    """

rule run_hatchet:
    input:
        bb=os.path.join(results_dir, 'hatchet2', '{patient}', 'bb', 'bulk.bb'),
        ini=os.path.join(results_dir, "hatchet2", '{patient}', "{patient}.ini")
    resources:
        time="95:59:59",
        mem_mb=64000,
    threads: 16
    output:
        os.path.join(results_dir, "hatchet2", '{patient}', "results", "best.bbc.ucn") #TODO: list more outputs
    shell: """
    module load gurobi
    hatchet run {input.ini}
    """  

    