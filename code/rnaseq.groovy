//load file with paths to necessary software
load '/oshlack_lab/jovana.maksimovic/scripts/bpipe-pipelines/software_paths.groovy'
//load bpipe stages required for this pipeline
load '/bpipe_stages.groovy'

FASTP = "/home/jmaksimovic/fastp-0.20.1/fastp"
TRIM_FRONT1 = 5
TRIM_TAIL2 = 5

// directory of the genome index to be used for mapping
//GEN_DIR = "/oshlack_lab/shared/genomes/hg38/gencode_v34/indices/star/star_v2.7.5b/"
// genome fasta file
//GEN_FASTA = "/oshlack_lab/shared/genomes/hg38/gencode_v34/GRCh38.primary_assembly.genome.fa"
// GTF file (used by STAR for junction mapping)
//GTF = "/oshlack_lab/shared/genomes/hg38/gencode_v34/gencode.v34.annotation.gtf"

// directory of the genome index to be used for mapping
GEN_DIR = "/oshlack_lab/shared/genomes_mcri/hg38/star_v2.7.5b/"
// genome fasta file
GEN_FASTA = "/oshlack_lab/shared/genomes_mcri/hg38/fasta/hg38.fa"
// GTF file (used by STAR for junction mapping)
GTF = "/oshlack_lab/shared/genomes_mcri/hg38/gtf/gencode.v34.annotation.gtf"

//STAR: length of genomic sequence around annotated junction; max(ReadLength)-1; default: 100
SJBOHANG = 149
OUTSJCOLL = 5000000 //Default: 1000000
SJDBINSNSJ = 2000000 //--limitSjdbInsertNsj default: 1000000
PEOVERLAPMIN = 10 //min overlap between PE reads to trigger merge and remap

SALMONINDEX = "/oshlack_lab/shared/genomes/hg38/gencode_v34/indices/salmon/salmon_v1.3.0"

//number of treads per sample for multi-threaded tools
NTHREADS = 8 //no. threads to use

MULTIQCOUTDIR = "."
MULTIQCDIR = "."

run {
	~".*(CMV[0-9]+|Corriel|NTC-2).*.fastq.gz" * [
      	~".*_R(1|2).fastq.gz" * [ concat_fastq ] ] +
	"%_R*.fastq.gz" * [ fastp ] + star_map_1pass_PE +
	//"%.trim.R*.fastq.gz" * [ star_map_2pass_PE + samtools_index, salmon_quant ] + multiqc
	"%.trim.R*.fastq.gz" * [ star_map_2pass_PE, salmon_quant ] + multiqc
}

