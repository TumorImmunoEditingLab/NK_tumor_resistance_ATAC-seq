export username=$(whoami)
NXF_SINGULARITY_CACHEDIR=/home/${username}/bioinf_isilon/core_bioinformatics_unit/Internal/pipelines/atacseq/atacseq-2.0/singularity-images 
export NXF_SINGULARITY_CACHEDIR=$NXF_SINGULARITY_CACHEDIR

# design file
input_design="atacseq_EXP9_5_6_7_design.csv"  #for atac-seq-2.0 the format is sample, fastq_1, fastq_2, replicate

# reference files
ref_fasta=/home/${username}/bioinf_isilon/core_bioinformatics_unit/Internal/references/Mus_musculus_GRCm38_v102/Mus_musculus.GRCm38.dna.primary_assembly.fa
ref_bwa_index=/home/${username}/bioinf_isilon/core_bioinformatics_unit/Internal/references/Mus_musculus_GRCm38_v102/
ref_gtf=/home/${username}/bioinf_isilon/core_bioinformatics_unit/Internal/references/Mus_musculus_GRCm38_v102/Mus_musculus.GRCm38.102.gtf
ref_bed=/home/${username}/bioinf_isilon/core_bioinformatics_unit/Internal/references/Mus_musculus_GRCm38_v102/Mus_musculus.GRCm38.102.bed
ref_tss_bed=/home/${username}/bioinf_isilon/core_bioinformatics_unit/Internal/references/Mus_musculus_GRCm38_v102/Mus_musculus.GRCm38.102.tss.bed
ref_blacklist=/home/${username}/bioinf_isilon/core_bioinformatics_unit/Internal/references/Mus_musculus_GRCm38_v102/mm10.blacklist.bed

# Effective genome size for GRCm38	2308125349
# https://deeptools.readthedocs.io/en/develop/content/feature/effectiveGenomeSize.html
# --max_cpus; --max_memory; --max_time is "Maximum number of CPUs, memory, time that can be requested for any single job."! not max system cpus (as in local executor)

nextflow run nf-core/atacseq -r 2.0 \
    -c /data_synology_rg3/aleksandr_b/ab3/nfcr-atacseq-2.0_run_full_batch/biohazard_low.config \
    --outdir /data_synology_rg3/aleksandr_b/ab3/nfcr-atacseq-2.0_run_full_batch/results/ \
    -profile singularity \
    --max_cpus 12 \
    --max_memory '10.GB' \
    --max_time '240.h' \
    --macs_gsize 2308125349 \
    --input ${input_design} \
    --seq_center 'CeMM-BSFandCCRI' \
    --fasta ${ref_fasta} \
    --bwa_index ${ref_bwa_index} \
    --gtf ${ref_gtf} \
    --gene_bed ${ref_bed} \
    --tss_bed ${ref_tss_bed} \
    --blacklist ${ref_blacklist} \
    --mito_name 'MT' \
    --igenomes_ignore \
    --skip_deseq2_qc > atacseq_run_20230223.log 2>&1



