GFF = "/home/jmaksimovic/team-folder/projects/MCRI/lisa.hui/amnio-cell-free-RNA/data/star-genome-analysis/gencode.v34.annotation.DEXSeq.chr.gff"

count_exons_RNA_PE = {
    //Count reads across features using RNA data
    doc "Count reads across features from RNA data"
    output.dir = "count-exons-pe"

    exec """
      module load miniconda3;
      python /home/jmaksimovic/.local/share/rstudio/library/DEXSeq/python_scripts/dexseq_count.py -p yes -f bam $GFF $input.bam $output.txt
    """, "count_exons"
}

count_exons_RNA_SE = {
    //Count reads across features using RNA data
    doc "Count reads across features from RNA data"
    output.dir = "count-exons-se"

    exec """
      module load miniconda3;
      python /home/jmaksimovic/.local/share/rstudio/library/DEXSeq/python_scripts/dexseq_count.py -f bam $GFF $input.bam $output.txt
    """, "count_exons"
}

run {
   [ "%PE.Aligned.out.srt.dedup.bam" * [ count_exons_RNA_PE ],
   "%SE.Aligned.out.srt.dedup.bam" * [ count_exons_RNA_SE ] ]
}
