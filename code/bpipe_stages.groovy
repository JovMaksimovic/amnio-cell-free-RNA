// File preparation
concat_fastq = {
    // Concatenate fastq files from the same sample
    doc "Concatenate fastq files from the same run"
    output.dir = "catfastq"

    produce(output.gz.prefix.prefix.prefix + ".fastq.gz"){
        exec """
                cat $inputs.gz > $output.gz
        ""","tiny"
    }
}

fastp = {
    // fastqc quality control
    doc "Quality control using FASTP"
    output.dir = "fastp"
    def sample_name = branch.name

    produce(sample_name + '.trim.R1.fastq.gz',
            sample_name + '.trim.R2.fastq.gz') {
        exec """
              $FASTP -i $input1.gz
              -I $input2.gz
              -o $output1.gz
              -O $output2.gz
              -y --detect_adapter_for_pe
              --trim_front1 $TRIM_FRONT1
              --trim_tail2 $TRIM_TAIL2
              -h ${sample_name}.fastp.html
              -j ${sample_name}.fastp.json
    """
    }
}

// Mapping
star_genome_gen = {
    //Generate STAR genome index
    doc "Generate STAR genome index"

    def gen_dir = "${GEN_DIR}";

    produce(gen_dir + "Genome") {
            exec """module load star/2.7.5b;

            STAR --runMode genomeGenerate
		            --runThreadN $NTHREADS
		            --genomeDir $GEN_DIR
		            --genomeFastaFiles $GEN_FASTA
		            --sjdbGTFfile $GTF
		            --sjdbOverhang $SJBOHANG
            ""","stargenind"
    }

    forward inputs
}

//STAR --runMode genomeGenerate --runThreadN 8 --genomeDir star_v2.7.5b --genomeFastaFiles fasta/hg38.fa --sjdbGTFfile gtf/gencode.v34.annotation.gtf --sjdbOverhang 149

star_map_1pass_PE = {
  //Map paired-end reads using the STAR aligner: 1st pass
  String inString = "${inputs}";
  def inList = inString.split(' ').collect{it as String}
  inList.removeAll { it.toLowerCase().contains('unpaired') }

  String R1 = inList[0];
  String R2 = inList[1];

  for (int i = 2; i < inList.size(); i++) {
        if (i % 2 == 0){
                R1 = R1+","+inList[i];
        } else {
                R2 = R2+","+inList[i];
        }
  }

  String files = R1+" "+R2;
  println files;

  doc "Map paired-end reads using the STAR aligner: 1st pass"
  output.dir = "mapped/1pass"

  from("*.fastq.gz") {
    produce ("SJ.out.tab") {
      exec """
        $STAR --genomeDir $GEN_DIR
        --readFilesIn ${files}
        --readFilesCommand zcat
        --outSAMtype None
        --limitOutSJcollapsed $OUTSJCOLL
        --peOverlapNbasesMin $PEOVERLAPMIN
        --runThreadN $NTHREADS
        --sjdbGTFfile $GTF
        --outFileNamePrefix ${output.dir}/
      ""","star1pass"
    }
  }
}

star_map_2pass_PE = {
  //Map paired-end reads using the STAR aligner: 2nd pass
  doc "Map paired-end reads using the STAR aligner: 2nd pass"
  output.dir = "mapped"

   def sample_name = branch.name.substring(0, branch.name.length() - 2)

   produce(sample_name + "Aligned.sortedByCoord.out.bam") {
    exec """
        $STAR --genomeDir $GEN_DIR
        --readFilesIn $input1.gz $input2.gz
	      --sjdbFileChrStartEnd ${output.dir}/1pass/SJ.out.tab
	      --sjdbOverhang $SJBOHANG
        --outFileNamePrefix ${output.dir}/$sample_name
        --readFilesCommand zcat
        --limitSjdbInsertNsj $SJDBINSNSJ
        --peOverlapNbasesMin $PEOVERLAPMIN
        --outSAMtype BAM SortedByCoordinate
        --quantMode TranscriptomeSAM GeneCounts
        --outReadsUnmapped Fastx
        --runThreadN $NTHREADS
    ""","star2pass"
  }
}

samtools_index = {
    // index bam files
    doc "Index bam files"
    output.dir = "mapped"

    transform("bam") to ("bam.bai") {
        exec """
		      $SAMTOOLS index $input.bam
	      ""","srtindex"
    }
    forward input
}

salmon_quant = {
    doc "Quasi-map and quantify paired-end RNAseq reads using Salmon"
    output.dir = "quants/" + branch.name.substring(0, branch.name.length() - 2) + "/"

    def reads = '-r ' + inputs
    if (reads.size() > 1) {
        reads = '-1 ' + inputs[0] + ' -2 ' + inputs[1]
    }

    produce('quant.sf') {
        exec """
        $SALMON quant --seqBias --gcBias --dumpEq
        --index $SALMONINDEX
        -l A $reads
        -p $NTHREADS
        -o ${output.dir}
        """, "salmon"
    }
}

multiqc = {
   // summarise statistics from all tools using multiqc
   doc "Pipeline summary report by MultiQC"

   exec """
   module unload `module -t list 2>&1 | xargs`;
   module load multiqc/1.8;

	 multiqc -f -o $MULTIQCOUTDIR
	 --ignore .bpipe/
	 --ignore-samples 1pass
	 --ignore ignore-overlap-mapping/
	 $MULTIQCDIR
   """, "multiqc"
}
