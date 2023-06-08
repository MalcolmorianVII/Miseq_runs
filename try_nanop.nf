nextflow.enable.dsl=2

workflow {
    //Get long read
    reads_ch = Channel.of('BZD62J','BZD6D3','BZD2E5','CDN5H1','BZD50K','CDN33K','BZD5DD','BZD6GK','BZD2ZC','BZD44J','BZD5JZ','CDN5ZB')
    //Run flye
    DISPLAY(reads_ch)
    second(DISPLAY.out) | view
}


process DISPLAY {
    fair true
    tag "Just list cool stuff"
    // publishDir "${params.fastq_pass}/${dataset_id}",mode:"copy",overwrite: false

    input:
    val(dataset_id)

    output:
    path "${dataset_id}.txt"

    script:
    """
    echo ${dataset_id} > ${dataset_id}.txt 
    """
}

process second {

    input:
    file sample_txt

    output:
    stdout 

    script:
    """
    cat ${sample_txt}  
    """
}

