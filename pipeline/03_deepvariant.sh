#!/bin/bash -l
#SBATCH --partition=gpu
#SBATCH --gres=gpu:a100:2
#SBATCH --mem=100g
#SBATCH --time=0-12:00:00
#SBATCH --mail-type=ALL
#SBATCH --job-name="pb-pipeline"
#SBATCH --output=log/%x_%A_%a.log
##################################################################


## Get configuration file
if [ -f config.txt ]; then
    source config.txt
else
    echo "You're missing the config.txt file"
    exit
fi

# Load modules and environments
module load singularity/3.9.3
module load parabricks/4.2.1-1
module load bwa/0.7.17
module load samtools/1.19.2

# Job variables
N=${SLURM_ARRAY_TASK_ID}
NUM_CORES=2

# Define directories
RAW_DIR="varseq_data"
DEEPVARIANT_DIR="varseq_results"
FASTA="${RAW_DIR}/tair10.fasta"
SIF="/opt/linux/rocky/8.x/x86_64/pkgs/parabricks/4.2.1-1/parabricks_4.2.1-1.sif"
METADATA="metadata.csv"

# Running workflow
echo "## Starting PB Pipeline Var-seq Workflow ##########"

# Reading the metadata CSV file
sed -n ${N}p ${METADATA} | while IFS="," read SAMPLENAME FASTQ1 FASTQ2
do
    SAMPLENAME=$(echo ${SAMPLENAME} | xargs)  # Remove any leading/trailing spaces or tabs
    echo "Processing sample = ${SAMPLENAME}"

    # Generate output directories
    DIR=${DEEPVARIANT_DIR}/${SAMPLENAME}
    [ ! -d "$DIR" ] && mkdir -p "${DIR}"


        # Run DeepVariant using the sorted BAM file
        singularity run --nv ${SIF} \
            pbrun deepvariant_germline \
            --in-fq ${RAW_DIR}/${FASTQ1} ${RAW_DIR}/${FASTQ2} \
            --num-gpus ${NUM_CORES} \
        	--ref ${FASTA} \
        	--gvcf \
       	 	--out-duplicate-metrics ${DEEPVARIANT_DIR}/${SAMPLENAME}/${SAMPLENAME}.${INDEX_NAME}.metrics.txt \
        	--out-bam ${DEEPVARIANT_DIR}/${SAMPLENAME}/${SAMPLENAME}.${INDEX_NAME}.bam \
        	--out-variants ${DEEPVARIANT_DIR}/${SAMPLENAME}/${SAMPLENAME}.${INDEX_NAME}.g.vcf.gz 
done
