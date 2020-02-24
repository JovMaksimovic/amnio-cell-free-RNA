load '/group/bioi1/jovanam/scripts/bpipe-pipelines/software_paths.groovy'
load '/group/bioi1/jovanam/scripts/bpipe-pipelines/rnaseq_stages.groovy'
load '/group/bioi1/jovanam/scripts/bpipe-pipelines/common_stages.groovy'

//STAR parameters
TRIM = "" //this is required for no trimming
GEN_DIR = "\$GENOMES/hg38/star_2_5/" //directory of the genome index to be used for mapping
GEN_FASTA = "\$GENOMES/hg38/fasta/hg38.fa" //genome fasta file
GTF = "\$GENOMES/hg38/star/hg38_GENCODEV20_Comp.gtf" //GTF file (used by STAR for junction mapping)
SAF = "\$GENOMES/hg38/saf/hg38_GENCODEV20_Comp.saf" //SAF file (used by featureCounts for read counting)
SJBOHANG = 149 //Length of genomic sequence around annotated junctions; max(ReadLength)-1; default: 100
OUTSJCOLL = 5000000 //Default: 1000000
SJDBINSNSJ = 2000000 //--limitSjdbInsertNsj default: 1000000

//SALMON parameters
SALMONINDEX = "/usr/local/transcriptomes/hg38/indexes/salmon/salmon-v0.14.1/hg38_gencode"

//Other parameters
NTHREADS = 8
MULTIQCDIR = "."

//Paired-end rna-seq pipeline with salmon quantification
run {
    ~".*(CMV[0-9]+|Corriel|NTC-2).*.fastq.gz" * [
      ~".*_R(1|2).fastq.gz" * [ concat_fastq ] ] +
    "%.fastq.gz" * [ fastqc ] + star_map_1pass_PE +
    "%_R*.fastq.gz" * [ salmon_quant_PE, star_map_2pass_PE ] + count_reads_RNA +
    "%.bam" * [ sort_bam + index_bam ] + multiqc
}

