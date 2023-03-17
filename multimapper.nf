#!/usr/bin/env nextflow
nextflow.enable.dsl=2


def validate_choice_param(param_option, param, choices) {
    param_name = (param_option - "--").replaceAll("_", " ")
    if (param) {
        if (!choices.any { it.contains(param.toString()) }) {
            log.error("Please specify the ${param_name} using the ${param_option} option. Possibilities are ${choices}.")
            return 1
        }
    } else {
        log.error("Please specify the ${param_name} using the ${param_option} option")
        return 1
    }
    return 0
}

def validate_parameters() {
    def errors = 0

    errors += validate_choice_param("--mapper", params.mapper, ["smalt", "bwa-mem2, bowtie2"])

    if (errors > 0) {
        log.error(String.format("%d errors detected", errors))
        exit 1
    }
}

def build_tool_options(tool_box) {
    // accepts parameter list of command line alterations for the tool used
    // initiate parameter list that are in use
    def tool_belt = []

    for( tool : tool_box ) {
        tool_setting = tool.split(' ');
        // tools are set up with '-flag input', split by gap and take second item
        tool_status = tool_setting[1]

        // tools are marked 'unused' if they are not in use by default, take only those in use
        if( tool_status != 'unused' ) {
            tool_belt.add(tool)
    }
}
    // return a string of all used flags in the form '-flag one -flag2 two'
    equipped_tools = tool_belt.join(' ')
    return equipped_tools
}

process QC {

    input:
    tuple val(id), val(seq1), val(seq2)

    output:
    tuple val(id), val(seq1), val(seq2) //TODO: add filter for qc failures

    script:
    if (params.qc == 'on') {
    """
    trim_galore --paired --illumina --fastqc -q 20 ${seq1} ${seq2}
    """
    }
}

process SMALT_INDEXING {
    input:
    path(ref)
    output:
    tuple val(index), path("${index}.smi"), path("${index}.sma")
    script:
    index = ref.getBaseName()
    """
    smalt index -k $params.smltk -s $params.smlts ${index} ${ref}
    """
}

process SMALT_MAPPING {
  input:
    tuple val(id), path(read1), path(read2)
    tuple val(index), path(smi), path(sma)

  output:
    path('*.sam')

  script:
  options = build_tool_options(params.smlt)
   """
   smalt map $options $index $read1 $read2 > ${id}.sam
   """
}

process BWA_INDEXING {
    input:
    path(ref)
    output:
    tuple val(index), path("${index}.*")
    script:
    index = ref.getBaseName()
    """
    bwa-mem2 index $ref -p ${index}
    """
}

process BWA_MAPPING {
  input:
    tuple val(id), path(read1), path(read2)
    tuple val(index), path(x)

  output:
    path('*.sam')

  script:
      """
      bwa-mem2 mem ${index} $read1 $read2 > out.sam
      """
}

process BOWTIE_INDEXING {
    input:
    path(ref)
    output:
    tuple val(index), path("${index}.*")
    script:
    index = ref.getBaseName()
    """
    bowtie2-build $ref ${index}
    """
}

process BOWTIE_MAPPING {
  input:
    tuple val(id), path(read1), path(read2)
    tuple val(index), path(x)

  output:
    path('*.sam')

  script:
      """
      bowtie2 -x ${index} -1 $read1 -2 $read2
      """
}

process BAMIFY {
    input:
        tuple val(id), path(read1), path(read2)
    output:
        tuple val(id), path(read_out)
    script:
    read_out = '${id}.sam'
    """
    samtools view -b -h -O BAM -@ 2 -o ${id}.bam ${id}.sam
    """
}

process VCF_SORT {
    input:
        path(vcf)
        path(ref)
        tuple val(id), val(seq1), val(seq2)

    output:
        path("*.bam")

    script:
    """
    sambamba sort ${id}
    sambamba markdup ${ref} ${id} > out.bam
    """

}

process SH_TIDY {
    input:
        path(bam)
        path(ref)
    output:
        path(endvcf)

    script:
    endvcf = "var.vcf"
    """
    samtools sort $bam > ${bam}_tmp
	picard MarkDuplicates INPUT=${bam}_tmp OUTPUT=${bam}_tmp2 METRICS_FILE=${bam}_metrics.txt
    gatk -I ${bam}_tmp2  -R ${ref} -T RealignerTargetCreator -o ${bam}
	gatk -I ${bam}_tmp2  -R ${ref} -T IndelRealigner --filter_bases_not_stored -targetIntervals ${bam}.intervals -o ${bam}
    """
}

process SH_SNP {
    input:
        path(ref)
    output:
        path(endvcf)

    script:
    endvcf = "var.vcf"
    """
    join_dna_files_with_indels.py -r ${ref} -o {aln} -t "+tmpname+"_mfas.txt
    'summarise_snps.py -g -w -r '+options.ref.split("/")[-1].split(".")[0]+' -o '+options.output+' -i '+options.output+'.aln'

    """
}

process FREEBAYES {
    input:
        path(bam)
        path(ref)
    output:
        path(endvcf)

    script:
    bamname = bam.getBaseName()
    endvcf = "{bamname}.vcf"
    """
    freebayes -f ${ref} ${bam} > {bamname}.vcf
    """
}

process TREEBUILD {
    input:
        path(bam)
        path(ref)

    """
    RAxML-NG -f ${ref} ${bam} > tree.vcf
    """

}
workflow SNP_VARIANTS {
    take:
    sams
    ref

    main:
    if (params.snp == 'sh16') {
    SH_SNP(ref)
    SH_TIDY(ref, SH_SNP.out)
    vcfs = SH_SNP.out
    }
    if (params.snp == 'freebayes') {
    FREEBAYES(sams, ref)
    vcfs = FREEBAYES.out
    }

    emit:
    vcfs
}

workflow MAPPING {
    take:
    fastqc
    ref

    main:
    if (params.mapper == 'smalt') {
    SMALT_INDEXING(ref)
    SMALT_MAPPING(fastqc, SMALT_INDEXING.out)
    sams = SMALT_MAPPING.out
    }
    if (params.mapper == 'bwa-mem2') {
    BWA_INDEXING(ref)
    BWA_MAPPING(fastqc, BWA_INDEXING.out)
    sams = BWA_MAPPING.out
    }
    if (params.mapper == 'bowtie2') {
    BOWTIE_INDEXING(ref)
    BOWTIE_MAPPING(fastqc, BWA_INDEXING.out)
    sams = BWA_MAPPING.out
    }

    emit:
    sams
}

workflow {
    if (params.help) {
        printHelp()
        exit 0
    }

    validate_parameters()

    sequences = Channel.fromPath(params.in + "/*.gz").map{[it.getSimpleName().replaceFirst(/_[12FRL]/, ""), it]}.groupTuple(size:2).map{it.flatten()}
    reference = Channel.fromPath(params.ref)
    sequences.view()

    QC(sequences)

    MAPPING(QC.out, reference)

//     VCF_SORT(MAPPING.out.vcfs, reference, QC.out)

    SNP_VARIANTS(MAPPING.out.sams, reference)

//     TREEBUILD(VCF_SORT.out.collect(), reference)
//     COLLECT_SNP_DATA(VCF_SORT.out.collect())
}
