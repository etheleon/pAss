#!/usr/bin/env sh

#Run this inside the /pAss/example folder
#clear up the rubbish files
rm -rf {testData, testOut}
cp -r $PWD/data $PWD/testData
mkdir testOut

../maxDiversity \
    --format \
    --outputDIR=$PWD/testOut \
    --refseqKO=$PWD/testData/refSeqProtDB \
    --contigs=$PWD/testData/contigs \
