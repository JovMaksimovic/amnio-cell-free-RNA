//load file with paths to necessary software
load '/oshlack_lab/jovana.maksimovic/scripts/bpipe-pipelines/software_paths.groovy'
//load bpipe stages required for this pipeline
load '/oshlack_lab/jovana.maksimovic/scripts/bpipe-pipelines/rnaseq_stages.groovy'
load '/oshlack_lab/jovana.maksimovic/scripts/bpipe-pipelines/common_stages.groovy'
load '/oshlack_lab/jovana.maksimovic/scripts/bpipe-pipelines/dnaseq_stages.groovy'

//path to adapter file for adapter type used in this experiment
ADAPTERS = "/config/binaries/trimmomatic/0.39/adapters/TruSeq3-PE-2.fa"
//trimmomatic parameters
LEADING = 3
TRAILING = 3
MINLEN = 35
CLIPSTRING = "ILLUMINACLIP:${ADAPTERS}:2:30:10:2:keepBothReads"

// directory of the genome index to be used for mapping
GEN_DIR = "/oshlack_lab/shared/genomes/hg38/gencode_v34/indices/star/star_v2.7.5b/"
// genome fasta file
GEN_FASTA = "/oshlack_lab/shared/genomes/hg38/gencode_v34/GRCh38.primary_assembly.genome.fa"
// GTF file (used by STAR for junction mapping)
GTF = "/oshlack_lab/shared/genomes/hg38/gencode_v34/gencode.v34.annotation.gtf"
FLATGTF = "/home/jmaksimovic/team-folder/projects/MCRI/lisa.hui/amnio-cell-free-RNA/data/star-genome-analysis/gencode.v34.annotation.fc.gtf"

//STAR: length of genomic sequence around annotated junction; max(ReadLength)-1; default: 100
SJBOHANG = 149
OUTSJCOLL = 5000000 //Default: 1000000
SJDBINSNSJ = 2000000 //--limitSjdbInsertNsj default: 1000000
TRIM = "trim."

//SALMONINDEX = "/oshlack_lab/shared/genomes/hg38/gencode_v34/indices/salmon/salmon_v1.3.0"

//number of treads per sample for multi-threaded tools
NTHREADS = 8 //no. threads to use for mapping, then counting

REMOVE_DUPLICATES = "TRUE"

MULTIQCOUTDIR = "."
MULTIQCANALYSISDIR = "."
MULTIQCIGNORESTR = "--ignore old/ --ignore test/ --ignore mapped/1pass/"

INDEXOUTDIR = "deduped"
STRAND = "1"

// gene identifier attribute for featureCounts
IDENT = "gene_id"

//Paired-end rna-seq pipeline
run {
  ~".*(CMV[0-9]+|Corriel|NTC-2).*.fastq.gz" * [
    ~".*_R(1|2).fastq.gz" * [ concat_fastq ] ] +
    "%.fastq.gz" * [ fastqc ] +
    "%_R*.fastq.gz" * [ trim_PE ] + star_map_1pass_PE +
    [ "%.trim.R*.fastq.gz" * [ star_map_2pass_PE +
                                 sort_bam +
                                 dedup +
                                 index_bam ] + count_reads_RNA_PE ],
      "%.trim.unpaired.R*.fastq.gz" * [ star_map_2pass_SE +
                                          sort_bam +
                                          dedup +
                                          index_bam ] + count_reads_RNA_SE ] + multiqc
}
