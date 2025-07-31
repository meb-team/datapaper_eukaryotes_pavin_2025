#!/bin/bash
#SBATCH --output=log/slurm-%A_%a.out
#SBATCH --error=log/slurm-%A_%a.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=48G
#SBATCH --tmp=600G
#SBATCH --time=12:00:00

# Setup the environment
module purge
module load gcc/8.1.0
module load bwa/0.7.17 samtools/1.16.1
module load rclone sratoolkit/3.0.0
module load conda/23.3.1

echo "START job:" $(date)
echo "HOSTNAME:" $(hostname)

function usage {
    echo -e "sbatch ThisScript.sh SRRidFile\n"
    echo "Provides the IDs in a file: 'SRRidFile'. One ID per line"
    exit 2
}

# Read argument passed to the command-line. There must be one, a file
# with accession IDs. I cannot check the validity of accessions, but I can check
# presence of an argument and test if it is a file!
if ! [[ $# -eq 1 ]]; then
    # Not enough or too much argument
    usage
fi
if ! [[ -f $1 ]]; then
    # the argument is not a file
    usage
fi

# Setup the working and temp directory, then move in it
WORKDIR=/storage/scratch/$USER/${SLURM_JOB_ID}_${SLURM_ARRAY_TASK_ID}
TMPDIR=$WORKDIR/tmp
mkdir -p $WORKDIR
mkdir -p $TMPDIR
export TMPDIR
export TEMPDIR=$TMPDIR
export TMP=$TMPDIR
export TEMP=$TMPDIR

echo "move in $WORKDIR"
cd $WORKDIR

# Extract current sample
### Get the $SLURM_ARRAY_TASK_ID Ist line from file $1
ACC=$(sed -n "$SLURM_ARRAY_TASK_ID p" $1)

# Create result and logs directory and move in it
RESDIR=$WORKDIR/$ACC
mkdir $RESDIR
mkdir $WORKDIR/log
mkdir $WORKDIR/flagstat
mkdir $WORKDIR/coverage

# Download the reads
## Note: prefetch as a max limit of 20GB by default, "--max-size 30G" can leverage it
echo "Get $ACC via SRAToolkit 'prefetch'" $(date)
prefetch --max-size 65G -O sra_data $ACC

echo "Extract the reads $ACC via SRAToolkit 'fastq-dump'" $(date)
## It seems important to customise the defline as with the default, 'bwa mem'
## raises the error "paired reads have different names". I ended with the
## conclusion that .1 and .2 are considered as part of the reads' name, while
## /1 and /2 are recognised properly.
fastq-dump -v --skip-technical --readids --dumpbase --split-files --clip \
    --defline-seq '@$ac.$sn/$ri' --defline-qual '+$ac.$sn/$ri' sra_data/$ACC
# Clean .sra
rm -r sra_data/

# Clean reads - I have FASTP installed in a Conda env
#####
#
# MAKE SURE THE CONDA ENVIRONMENT WITH FASTP IS ACCESSIBLE !!
#
#####
conda activate $HOME/conda/envs/fastp
echo "Clean $ACC reads " $(date)
fastp --in1 ${ACC}_1.fastq --in2 ${ACC}_2.fastq \
    --out1 ${ACC}_1.ok.fq --out2 ${ACC}_2.ok.fq \
    --detect_adapter_for_pe --n_base_limit 0 --average_qual 30 \
    --thread $SLURM_CPUS_PER_TASK -j $WORKDIR/log/$ACC.json \
    -h $WORKDIR/log/$ACC.html

# Clean raw reads
rm ${ACC}_1.fastq ${ACC}_2.fastq

# COPY the index of SAGs and MAGs
#####
#
# MAKE SUR THE BWA INDEX OF yOUR DATA IS IN YOUR WORKING DIRECTORY!
#
#####
INDEX="MAGs_SAGs_euk.fa.gz"

# Mapping
echo -e "Start mapping $ACC" $(date)
bwa mem -t $SLURM_CPUS_PER_TASK $INDEX ${ACC}_1.ok.fq ${ACC}_2.ok.fq |
    samtools view -bS -F4 -q30 - |
    samtools sort -@$SLURM_CPUS_PER_TASK -m 2G - >$ACC/$ACC.bam

echo "Flagstat initial BAM file"
samtools flagstat $ACC/${ACC}.bam >$WORKDIR/flagstat/${ACC}.flagstat

# Filtering and coverage - WITH A SINGULARITY CONTAINER
#####
#
# MAKE SURE THE SINGULARITY CONTAINER IS IN YOUR WORKING DIRECTORY !
#   ==> msamtools_coverm.sif
#
#####
MP=/data # Mount point within the singularity container

# Filter
## Note for the usage of Singularity container:
## I have to use the mounting point to READ a file WITHIN the container,
## but I can redirect the result printed in stdout directly into a file!
##
## Like I do not use special filter like '--best-hit', I there is no need to
## sort the reads by name first!
echo "Filtering $ACC with mSAMtools (-l 50 -p 95 -z 80)" $(date)
singularity exec -B $WORKDIR:$MP $WORKDIR/msamtools_coverm.sif \
    msamtools filter -b -l 50 -p 95 -z 80 ${MP}/$ACC/${ACC}.bam \
    >$ACC/${ACC}.filter.bam

echo "Flagstat after the filter" $(date)
samtools flagstat $ACC/${ACC}.filter.bam >$WORKDIR/flagstat/${ACC}.filter.flagstat

echo "Index filtered BAM" $(date)
samtools index $ACC/${ACC}.filter.bam

echo "Compute coverage from the filtered BAM" $(date)
samtools coverage -o $WORKDIR/coverage/${ACC}.filter.cov.tsv $ACC/${ACC}.filter.bam

echo "GZip the 'coverage' file"
gzip $WORKDIR/coverage/${ACC}.filter.cov.tsv

echo "END computing:" $(date)

#### Save the results
echo "Save:"
echo " - FILTERED BAM file"
echo " - Fastp log file"
echo " - Flagstat files"
echo " - Coverage file"

#####
#
# DO NOT FORGET TO COPY THE RESULTS IN A SAFE PLACE !
#
#####

# Wait and clean
echo "wait 10 s..."
sleep 10

conda deactivate

rm -r $WORKDIR

echo "END job:" $(date)
