params.seq_stats = "/home/bkutambe/data/miseq_runs/Fastq/seq_stats.csv"
seq_stats = file(params.seq_stats)
params.root_dir = "/home/bkutambe/data/miseq_runs/Fastq"
params.ref = "/home/bkutambe/salmonella/references/2022.11.29/AE006468.2.fasta"
ref = file(params.ref)

Channel
    .fromFilePairs("${params.root_dir}/{CNS8I6,CWM10U,CXL13B,CXL19X,CXM115,CXM141,CXM1BG,CXM1KM,CXM220,CXM25G,CXM2ER,CXM2J7,CXM2LS,CXT10N,CDNC1A,CRA14D,CWN10N,CXL13F,CXL1B4,CXM13T,CXM183,CXM1DY,CXM1L0,CXM23N,CXM25R,CXM2EU,CXM2JU,CXMZEV,CXZL1BP,CNS1G6,CRB13A,CXL12D,CXL14M,CXL1G8,CXM13V,CXM18D,CXM1II,CXM1N2,CXM243,CXM25U,CXM2FN,CXM2JW,CXMZKZ,CNS4TI,CRB13B,CXL12S,CXL19F,CXL1H6,CXM140,CXM18T,CXM1IU,CXM1RR,CXM24H,CXM2EI,CXM2I2,CXM2KQ,CXMZLW}/*_bbduk_{1,2}.fastq.gz", flat: true)
    .ifEmpty { exit 1, "Read pairs could not be found: " }
    .set { read_pairs_phenix }


process calculateStats {
    publishDir "${params.root_dir}", mode: "copy"

    input:
    tuple val(dataset_id), file(forward), file(reverse)
    file ref

    output:
    stdout

    script:
    """
    if [ ! -f ${seq_stats} ]; then
        echo "Sample ID,Average sequencing depth,SNPs,Indels,Ns,%reads mapped" > ${seq_stats}
    fi

    average_depth=\$(samtools depth ${params.root_dir}/${dataset_id}/phenix_bbduk/${dataset_id}.bam | awk '{sum+=\$3} END {print sum/NR}')
    snp_count=\$(grep -vc "^#" ${params.root_dir}/${dataset_id}/phenix_bbduk/${dataset_id}.filtered.vcf || echo 0)
    indel_count=\$(grep -vc "^#" ${params.root_dir}/${dataset_id}/phenix_bbduk/${dataset_id}.filtered.vcf | grep "INDEL" || echo 0)
    N_count=\$(samtools view -F 4 ${params.root_dir}/${dataset_id}/phenix_bbduk/${dataset_id}_all.fasta | grep -o N | wc -l)
    reads_cov=\$(samtools flagstat ${params.root_dir}/${dataset_id}/phenix_bbduk/${dataset_id}.bam | awk '/mapped \\(/ {mapped = \$1} /in total/ {total = \$1} END {printf "%.2f", (mapped / total) * 100}')

    echo "${dataset_id},\${average_depth},\${snp_count},\${indel_count},\${N_count},\${reads_cov}" >> ${seq_stats}
    """
}

workflow {
    calculateStats(read_pairs_phenix, ref)
}
