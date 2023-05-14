module load SAMtools/1.12-GCC-10.2.0
module load BWA/0.7.17-foss-2018b
module load BEDTools/2.27.1-foss-2018b

sbatch --job-name=MAPS --cpus-per-task=8 \
            --mem=30G --time=20:00:00 \
            -p short run_pipeline_trial1_Mye.sh