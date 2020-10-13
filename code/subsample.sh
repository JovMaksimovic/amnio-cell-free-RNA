#!/bin/bash

module load seqtk/1.0;

for f in catfastq/*CMV3*.fastq.gz
do
samp=`basename ${f}`
echo $samp
seqtk sample -s42 $f 10000 | gzip > ./test/catfastq/${samp}

done

exit
