
process prepareReference {
    label "mbull"

    publishDir "${params.ref_dir}", mode: "copy"

    input:
    file(ref_dir)

    output:
    file "${FASTA}.*"

    script:
    """
    phenix.py prepare_reference -r ${ref_dir}/${FASTA} \
    --mapper bwa \
    --variant gatk
    """
}

process runSNPPipeline {
    label "mbull"

    tag { dataset_id }
    publishDir "${params.root_dir}/${dataset_id}", mode: "copy"

    input:
    tuple val(dataset_id), file(forward), file(reverse)
    file(phenix_config)
    file(ref)

    output:
    file "${params.output_dir_name}/${dataset_id}.*"

    script:
    """
    phenix.py run_snp_pipeline \
    -r1 $forward \
    -r2 $reverse \
    -r ${ref} \
    -c ${phenix_config} \
    --keep-temp \
    --sample-name ${dataset_id} \
    -o ${params.output_dir_name}
    """
}

process convertVCFtoFASTA {
    label "mbull"

    tag { dataset_id }
    publishDir "${params.root_dir}/${dataset_id}", mode: "copy"

    input:
    tuple val(dataset_id), file(vcf)

    output:
    file "${params.output_dir_name}/${dataset_id}.*"

    script:
    """
    phenix.py vcf2fasta \
    -i ${vcf} \
    -o ${params.output_dir_name}/${dataset_id}_all.fasta \
    --reference ${params.ref}
    """
}

process cleanUpFiles {

    tag { dataset_id }
    publishDir "${params.root_dir}/${dataset_id}", mode: "copy"

    input:
    tuple val(dataset_id), file(fasta), file(vcf)

    script:
    """
    seqkit grep -n -i -v -p "reference" ${fasta} > ${params.output_dir_name}/${dataset_id}.fasta
    rm ${fasta}
    gzip ${vcf}
    """
}

workflow {
    Channel
        .fromFilePairs("${params.root_dir}/{CNS8I6,CWM10U,CXL13B,CXL19X,CXM115,CXM141,CXM1BG,CXM1KM,CXM220,CXM25G,CXM2ER,CXM2J7,CXM2LS,CXT10N,CDNC1A,CRA14D,CWN10N,CXL13F,CXL1B4,CXM13T,CXM183,CXM1DY,CXM1L0,CXM23N,CXM25,CXM2EU,CXM2JU,CXMZEV,CXZL1BP,CNS1G6,CRB13A,CXL12D,CXL14M,CXL1G8,CXM13V,CXM18D,CXM1II,CXM1N2,CXM243,CXM25U,CXM2FN,CXM2JW,CXMZKZ,CNS4TI,CRB13B,CXL12S,CXL19F,CXL1H6,CXM140,CXM18T,CXM1IU,CXM1RR,CXM24H,CXM2EI,CXM2I2,CXM2KQ,CXMZLW}/*_bbduk_{1,2}.fastq.gz", flat: true)
        .ifEmpty { exit 1, "Read pairs could not be found: " }
        .set { read_pairs_phenix }

    prepareReference(read_pairs_phenix, params.ref)
    runSNPPipeline(read_pairs_phenix, phenix_config, params.ref)
    convertVCFtoFASTA(runSNPPipeline.out, params.ref)
    cleanUpFiles(convertVCFtoFASTA.out)
}

workflow.onComplete {
    log.info "Nextflow Version: $workflow.nextflow.version"
    log.info "Command Line:     $workflow.commandLine"
    log.info "Container:        $workflow.container"
    log.info "Duration:         $workflow.duration"
}
