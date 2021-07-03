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

SALMONINDEX = "/oshlack_lab/shared/genomes/hg38/gencode_v34/indices/salmon/salmon_v1.3.0"
EQ = "--dumpEq"

//number of treads per sample for multi-threaded tools
NTHREADS = 12 //no. threads to use for mapping, then counting

REMOVE_DUPLICATES = "TRUE"

MULTIQCOUTDIR = "."
MULTIQCANALYSISDIR = "."
MULTIQCIGNORESTR = "--ignore .bpipe"

run {
	~".*(CMV[0-9]+|Corriel|NTC-2).*.fastq.gz" * [
      	~".*_R(1|2).fastq.gz" * [ concat_fastq ] ] +
	"%.fastq.gz" * [ fastqc ] +
	"%_R*.fastq.gz" * [ trim_PE ] +
	"%.fastq.gz" * [ fastqc ] +
  "%.trim.R*.fastq.gz" * [ run_salmon ] + multiqc
}
