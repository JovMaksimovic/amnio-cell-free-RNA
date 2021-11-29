FASTQC = "fastqc"
TRIMMOMATIC = "java -jar /config/binaries/trimmomatic/0.39/trimmomatic.jar"
STAR = "STAR"
FEATURECOUNTS = "featureCounts"
MARK_DUPLICATES = "java -jar /config/binaries/picard/2.17.3/picard.jar MarkDuplicates"
SAMTOOLS = "samtools"

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

trim_SE = {
    // trim single-end reads using trimmomatic
    doc "Trim poor quility bases and/or adapter sequences from reads using Trimmomatic"
    output.dir = "trimmed"

    produce("trim.fastq.gz"){
        exec """
            $TRIMMOMATIC SE -threads $NTHREADS $input.gz $output.gz $CLIPSTRING
            TRAILING:$TRAILING LEADING:$LEADING MINLEN:$MINLEN
        ""","trimmomatic"
    }
}

trim_PE = {
   // trim paired-end reads using trimmomatic
   doc "Trim poor quility bases and/or adapter sequences from reads using Trimmomatic"
   output.dir="trimmed"
   def sample_name = branch.name

   //println inputs;

   produce(sample_name+'.trim.R1.fastq.gz',
           sample_name+'.trim.unpaired.R1.fastq.gz',
           sample_name+'.trim.R2.fastq.gz',
           sample_name+'.trim.unpaired.R2.fastq.gz') {
        exec """
                $TRIMMOMATIC PE -threads $NTHREADS
                $input1.gz $input2.gz
                $output1.gz ${output2.prefix}.gz
                $output3.gz ${output4.prefix}.gz $CLIPSTRING
                LEADING:$LEADING TRAILING:$TRAILING MINLEN:$MINLEN
        ""","trimmomatic"
   }
}

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
        --runThreadN $NTHREADS
        --sjdbGTFfile $GTF
        --outReadsUnmapped Fastx
        --outFileNamePrefix ${output.dir}/
      ""","star1pass"
    }
  }
}

star_map_2pass_SE = {
    //Map single-end reads using the STAR aligner; 2nd pass
    String files = "${inputs}";
    files = files.replaceAll(' ',',');

    doc "Map single-end reads using the STAR aligner: 2nd pass"
    output.dir = "mapped"
    def sample_name = branch.name

   produce(sample_name+".SE.Aligned.out.bam") {
        exec """
	      	$STAR --genomeDir $GEN_DIR
	      	--readFilesIn ${files}
	      	--outFileNamePrefix ${output.prefix.prefix.prefix}.
		      --sjdbFileChrStartEnd ${output.dir}/1pass/SJ.out.tab
		      --sjdbOverhang $SJBOHANG
	      	--readFilesCommand zcat
	      	--limitSjdbInsertNsj $SJDBINSNSJ
	      	--outSAMtype BAM Unsorted
	      	--runThreadN $NTHREADS
      	""","star2pass"
   }
}

star_map_2pass_PE = {
  //Map paired-end reads using the STAR aligner: 2nd pass
  doc "Map paired-end reads using the STAR aligner: 2nd pass"
  output.dir = "mapped"

   def sample_name = branch.name

   produce(sample_name+".PE.Aligned.out.bam") {
    exec """
        $STAR --genomeDir $GEN_DIR
        --readFilesIn $input1.gz $input2.gz
	      --sjdbFileChrStartEnd ${output.dir}/1pass/SJ.out.tab
	      --sjdbOverhang $SJBOHANG
        --outFileNamePrefix ${output.prefix.prefix.prefix}.
        --readFilesCommand zcat
        --limitSjdbInsertNsj $SJDBINSNSJ
        --outSAMtype BAM Unsorted
        --runThreadN $NTHREADS
    ""","star2pass"
  }
}

count_reads_RNA_PE = {
    //Count reads across features using RNA data
    doc "Count reads across features from RNA data"
    output.dir = "counts-pe"

    produce("counts.txt") {
        exec """
            $FEATURECOUNTS
            -p
            -M
            --primary
            -t exon
            -g $IDENT
            -T $NTHREADS
            -F GTF
            -a $GTF
            -s $STRAND
            -o $output $inputs.bam
        ""","count"
    }
}

count_reads_RNA_SE = {
    //Count reads across features using RNA data
    doc "Count reads across features from RNA data"
    output.dir = "counts-se"

    produce("counts.txt") {
        exec """
            $FEATURECOUNTS
            -M
            --primary
            -t exon
            -g $IDENT
            -T $NTHREADS
            -F GTF
            -a $GTF
            -s $STRAND
            -o $output $inputs.bam
        ""","count"
    }
}

sort_bam = {
    // sort bam files
    doc "Sort bam files"
    output.dir = "sorted"

    filter("srt") {
        exec """
		$SAMTOOLS sort -m 1G -@ $NTHREADS -o $output $input
	""","srtindex"
    }
}

index_bam = {
    // index bam files
    doc "Index bam files"
    output.dir = "$INDEXOUTDIR"

    transform("bam") to ("bam.bai") {
        exec """
		      $SAMTOOLS index $input.bam
	      ""","srtindex"
    }
    forward input
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
   """, "multiqc"
}
