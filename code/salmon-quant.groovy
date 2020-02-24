//load file with paths to necessary software (meerkat paths)
load '/group/bioi1/jovanam/scripts/bpipe-pipelines/software_paths.groovy'
load '/group/bioi1/jovanam/scripts/bpipe-pipelines/common_stages.groovy'

SALMONINDEX = "/usr/local/transcriptomes/hg38/indexes/salmon/salmon-v0.14.1/hg38_gencode"
NTHREADS = 8

salmon_quant_PE = {

  doc "Quasi-map and quantify paired-end RNAseq reads using Salmon"
  def sample_name = branch.name
  output.dir = "quants/${sample_name}_quant/"

produce("quant.sf"){
  exec """
        $SALMON quant -i $SALMONINDEX -l A \
        -1 $input1.gz \
        -2 $input2.gz \
        -p $NTHREADS \
        --seqBias \
        --gcBias \
        --validateMappings \
        -o ${output.dir}
    ""","salmon"
}
}

run {
     "%_R*.fastq.gz" * [ salmon_quant_PE ]
}
