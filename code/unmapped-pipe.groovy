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

  from("*.mate*") {
    produce ("SJ.out.tab") {
      exec """
        module load star/2.7.5b;

        STAR --genomeDir $GEN_DIR
        --readFilesIn ${files}
        --outSAMtype None
        --limitOutSJcollapsed $OUTSJCOLL
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

   def sample_name = branch.name

   produce(sample_name + "Aligned.sortedByCoord.out.bam") {
    exec """
        module load star/2.7.5b;

        STAR --genomeDir $GEN_DIR
        --readFilesIn $input1 $input2
	      --sjdbFileChrStartEnd ${output.dir}/1pass/SJ.out.tab
	      --sjdbOverhang $SJBOHANG
        --outFileNamePrefix ${output.dir}/$sample_name
        --limitSjdbInsertNsj $SJDBINSNSJ
        --outSAMtype BAM SortedByCoordinate
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
        module load samtools/1.11;

		      samtools index $input.bam
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
	 --ignore ignore-overlap-mapping/
	 $MULTIQCDIR
   """, "multiqc"
}

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

//number of treads per sample for multi-threaded tools
NTHREADS = 8 //no. threads to use

MULTIQCOUTDIR = "."
MULTIQCDIR = "."

run {
	star_genome_gen + star_map_1pass_PE +
	"%_OvationSoloRNA_L000Unmapped.out.mate*" * [ star_map_2pass_PE + samtools_index ] + multiqc
}

