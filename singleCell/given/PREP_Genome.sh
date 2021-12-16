### CREATE GENOME

source $CODEBASE/tfcf/setup.sh

cd $HOME/omicstmp

mkdir newGenome/
oldGenomePath="$GFS/RESOURCES/Genomes/refdata-gex-mm10-2020-A/"

cp $oldGenome/fasta/genome.fa newGenome/genome.fa
cp $oldGenome/genes/genes.gtf newGenome/genes.gtf

x="GFP"
for x in GFP BFP; do
	pathFA=$CODEBASE/tfcf/metadata/${x}.fa

	numberBases=$(cat $pathFA | grep -v "^>" | tr -d "\n" | wc -c)
	gtf="$x\tunknown\texon\t1\t${numberBases}\t.\t+\t.\tgene_id xxx${x}xxx; transcript_id xxx${x}xxx; gene_name xxx${x}xxx; gene_biotype xxxprotein_codingxxx;"
	gtf=$(echo $gtf | sed 's/xxx/"/g')

	echo -e $gtf > ${x}.gtf

	cat ${x}.gtf >> newGenome/genes.gtf
	fold -w 60 $pathFA >> newGenome/genome.fa
	
	echo -e "" >> newGenome/genome.fa

	rm ${x}.gtf
done


grep ">" newGenome/genome.fa
tail -5 newGenome/genes.gtf
tail -200 newGenome/genome.fa


~/code/cellranger-6.0.1/cellranger mkref --genome=newGenomeExtended --fasta=newGenome/genome.fa --genes=newGenome/genes.gtf