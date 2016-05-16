#!/usr/bin/env sh

#cpanm local::lib
#Run this inside the /pAss/example folder

#clear up the rubbish files
#rm -rf testData
#rm -rf testOut
#cp -r $PWD/data $PWD/testData
#cp -r $PWD/out $PWD/testOut

#download MEGAN (need to test)
#wget http://ab.inf.uni-tuebingen.de/data/software/megan5/download/MEGAN_unix_5_11_3.sh
#bash MEGAN_unix_5_11_3.sh < megan_install

#run
../maxDiversity \
    --format \
    --outputDIR=$PWD/testOut \
    --refseqKO=$PWD/testData/refSeqProtDB \
    --contigs=$PWD/testData/contigs \
    --megan=$PWD/local/bin/MEGAN \
    --meganLicense=$PWD/MEGAN5-academic-license.txt
