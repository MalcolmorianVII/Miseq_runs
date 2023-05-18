nextflow.enable.dsl=2

// General configuration variables
params.root_dir = "/home/bkutambe/data/miseq_runs/Fastq"
params.ref = "/home/bkutambe/salmonella/references/2022.11.29/AE006468.2.fasta"
ref = file(params.ref)

params.output_dir_name = 'phenix_bbduk'
phenix_config = file("/home/bkutambe/salmonella/references/2022.11.29/phenix_config.yml")

params.help = false

Channel
    .fromFilePairs("${params.root_dir}/{CNS8I6,CWM10U,CXL13B,CXL19X,CXM115,CXM141,CXM1BG,CXM1KM,CXM220,CXM25G,CXM2ER,CXM2J7,CXM2LS,CXT10N,CDNC1A,CRA14D,CWN10N,CXL13F,CXL1B4,CXM13T,CXM183,CXM1DY,CXM1L0,CXM23N,CXM25,CXM2EU,CXM2JU,CXMZEV,CXZL1BP,CNS1G6,CRB13A,CXL12D,CXL14M,CXL1G8,CXM13V,CXM18D,CXM1II,CXM1N2,CXM243,CXM25U,CXM2FN,CXM2JW,CXMZKZ,CNS4TI,CRB13B,CXL12S,CXL19F,CXL1H6,CXM140,CXM18T,CXM1IU,CXM1RR,CXM24H,CXM2EI,CXM2I2,CXM2KQ,CXMZLW}/*_bbduk_{1,2}.fastq.gz", flat: true)
    .ifEmpty{exit 1, "Read pairs could not be found: "}
    .set {read_pairs_phenix}

process phenixing {
    // errorStrategy 'ignore'
    conda = "/home/bkutambe/miniconda3/envs/phenix"
    tag { dataset_id }
    publishDir "${params.root_dir}/${dataset_id}", mode: "copy"

    input:
    tuple val(dataset_id), file(forward), file(reverse) 
    file(phenix_config)
    file(ref)

    output:
    file "${params.output_dir_name}/${dataset_id}*"

    """
    phenix.py prepare_reference -r $ref \
    --mapper bwa \
    --variant gatk
    
    phenix.py run_snp_pipeline \
    -r1 $forward \
    -r2 $reverse \
    -r ${ref} \
    -c ${phenix_config} \
    --keep-temp \
    --sample-name ${dataset_id} \
    -o ${params.output_dir_name}

    phenix.py vcf2fasta \
    -i ${params.output_dir_name}/${dataset_id}.filtered.vcf \
    -o ${params.output_dir_name}/${dataset_id}_all.fasta \
    --reference ${params.ref} \
    
    /home/bkutambe/miniconda3/bin/seqkit grep -n -i -v -p "reference" ${params.output_dir_name}/${dataset_id}_all.fasta > ${params.output_dir_name}/${dataset_id}.fasta
    rm ${params.output_dir_name}/${dataset_id}_all.fasta
    gzip ${params.output_dir_name}/${dataset_id}.vcf
 
    """
}

workflow {
	phenixing(read_pairs_phenix,phenix_config,ref)
}

workflow.onComplete {
    log.info "Nextflow Version: $workflow.nextflow.version"
    log.info "Command Line:     $workflow.commandLine"
    log.info "Container:        $workflow.container"
    log.info "Duration:     $workflow.duration"
}
