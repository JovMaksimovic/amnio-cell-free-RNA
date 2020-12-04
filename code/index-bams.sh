#! /bin/bash
for i in *.bam
do
  echo "Sorting ${i}"
  samtools index ${i}
done
