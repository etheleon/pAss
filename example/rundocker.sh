cp -r ./SingleCopyGene ./SCG
docker run --rm \
	-v `pwd`/SCG/data/konr:/data/refSeqProtDB \
	-v `pwd`/SCG/data/newbler:/data/contigs \
	-v `pwd`/SCG/out:/data/out \
	-v `pwd`/SCG/misc:/data/misc \
	pass \
	/tmp/pAss/maxDiversity --outputDIR /data/out --format --threads 8 --refseqKO /data/refSeqProtDB  --contigs /data/contigs  --megan /usr/local/bin/MEGAN --meganLicense /data/misc/MEGAN5-academic-license.txt 

docker run --rm \
	-v `pwd`/SCG/data/konr:/data/refSeqProtDB \
	-v `pwd`/SCG/data/newbler:/data/contigs \
	-v `pwd`/SCG/out:/data/out \
	-v `pwd`/SCG/misc:/data/misc \
	pass \
	chown -R `id -u`:`id -g` /data
