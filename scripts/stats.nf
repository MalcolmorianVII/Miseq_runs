params {
    alignment_file = "/path/to/alignment.bam"
    reference_file = "/home/bkutambe/salmonella/references/2022.11.29/AE006468.2.fasta"
    variant_file = "/path/to/variants.vcf"
}

process calculateSequencingDepth {
    input:
    file alignment from alignment_file

    output:
    stdout result

    script:
    """
    samtools depth $alignment | awk '{sum+=$3} END {print "Average sequencing depth: ", sum/NR}'
    """
}

process countSNPs {
    input:
    file variants from variant_file

    output:
    stdout result

    script:
    """
    grep -vc "^#" $variants
    """
}

process countIndels {
    input:
    file variants from variant_file

    output:
    stdout result

    script:
    """
    grep -vc "^#" $variants | grep "INDEL"
    """
}

process countNs {
    input:
    file reference from reference_file

    output:
    stdout result

    script:
    """
    samtools faidx $reference | awk '{sum+=$2} END {print "Number of Ns in reference: ", sum}'
    """
}

process calculateProportionMapped {
    input:
    file alignment from alignment_file

    output:
    stdout result

    script:
    """
    samtools flagstat $alignment | grep "mapped (" | awk '{print "Proportion of reads mapped: ", $1}'
    """
}

workflow {
    sequencingDepth = calculateSequencingDepth(alignment_file)
    snpCount = countSNPs(variant_file)
    indelCount = countIndels(variant_file)
    nsCount = countNs(reference_file)
    proportionMapped = calculateProportionMapped(alignment_file)

    sequencingDepth.out.stdout.into { depth }
    snpCount.out.stdout.into { snps }
    indelCount.out.stdout.into { indels }
    nsCount.out.stdout.into { ns }
    proportionMapped.out.stdout.into { mapped }

    output:
    file "sequencing_depth.txt" into sequencing_depth_stats
    file "snps.txt" into snp_stats
    file "indels.txt" into indel_stats
    file "ns.txt" into ns_stats
    file "proportion_mapped.txt" into mapped_stats

    depth.into(sequencing_depth_stats)
    snps.into(snp_stats)
    indels.into(indel_stats)
    ns.into(ns_stats)
    mapped.into(mapped_stats)
}
