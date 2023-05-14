nthreads=8
juicerJ="/home/jmendietaes/programas/juicer_v1.6/juicer_tools_1.22.01.jar"

tadbit export --workdir ${finalOut} --resolution 5000 --format hic \
        --output ${finalOut}/${filename}_5Kb.hic --cpus ${nthreads} \
        --chr_name ${chromCheck} --juicerjar ${juicerJ}


java -Xmx32g -jar ${juicerJ} pre -j ${nthreads} hic_export_3c9f06620b.tsv ${finalOut}/${filename}_5Kb.hic hic_3c9f06620b.chrom.sizes

calderR="/home/jmendietaes/programas/miniconda3/envs/calder2/bin/R"




conda activate python3
nthreads=16
juicerJ="/home/julen/programas/Juicer/juicer_tools_1.22.01.jar"
resol=50000
outpath="/media/julen/Elements/pruebas"
filename="DM_HiC_normal"
inBam="/scratch/julen/HiC/bamfiles/intersection_94e4921e80.bam"
chromCheck=$(for n in {1..19}; do echo chr${n}; done)

resolN=$(echo "print(${resol}//1000)" | python)
tadbit export --workdir ${outpath} --bam ${inBam} --resolution ${resol} --format hic \
        --output ${outpath}/${filename}_${resolN}Kb.hic --cpus ${nthreads} \
        --chr_name ${chromCheck} --juicerjar ${juicerJ}