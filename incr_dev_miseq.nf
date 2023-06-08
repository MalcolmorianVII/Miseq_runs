nextflow.enable.dsl=2

process batchDir {
    // Mv batchdir from where it is received from the miseq directory to process dir
    input:
    path sourceDir
    path targetDir

    output:
    path targetDir 

    script:
    """
    mv  ${sourceDir}/${BATCH} ${targetDir}
    """
}

process fofn {
    publishDir params.todo,mode:"copy"

    input:
    path directory

    output:
    file 'fofn.txt'

    script:
    """
    find ${directory}/${BATCH} -type f -name '*.fastq.gz' | sed -e 's#.*/##' -e 's/_.*//' | sort -u > fofn.txt
    """
}

process order_fastqs {
    // Organize fastqs per sample directory & rename them????
    publishDir "${params.des}/${BATCH}",mode:"copy"
    input:
    file fofn
    path directory

    output:
    path "C*"  

    script:
    """
    while read -r line
    do
        mkdir -p \$line
        mv ${directory}/${BATCH}/\${line}*.fastq.gz \$line
    done < ${fofn}
    """
}

process validate {
    // Create two files i.e sufficient.txt & insufficient.txt
    publishDir params.todo,mode:"copy", overwrite: false

    input:
    path sample
    // file directory

    output:
    file 'filtered_fofn.txt'
    // stdout

    script:
    """
    find ${sample}/* -type f -name "*.fastq.gz" -size -100c -exec dirname {} \\; | sort -u >> filtered_fofn.txt
    """
}

workflow {
    // source_ch = Channel.fromPath('/Users/malcolmorian/Documents/Bioinformatics/Projects2023/source_data/2023.05.24_RUN')
    // dest_ch = Channel.fromPath('/Users/malcolmorian/Documents/Bioinformatics/Projects2023/dest_data')
    source_dir = '/Users/malcolmorian/Documents/Bioinformatics/Projects2023/source_data'
    dest_dir = '/Users/malcolmorian/Documents/Bioinformatics/Projects2023/dest_data'
    batchDir(source_dir,dest_dir)
    fofn(batchDir.out)
    order_fastqs(fofn.out,batchDir.out)
    def samples = order_fastqs.out.flatten()
    validate(samples) 
   
}