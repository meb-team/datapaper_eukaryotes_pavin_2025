#!/bin/bash
#SBATCH --job-name=back_OL
#SBATCH --output=log/slurm-%A-%a.out
#SBATCH --error=log/slurm-%A-%a.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=24G
#SBATCH --array=0-15%6

###########
# There are 16 samples from the condition OXIC + LARGE fraction
###########

# Setup the environment
module purge
module load gcc/8.1.0
module load bwa/0.7.17 samtools/1.9
module load rclone

#### test this ####
unset http_proxy
unset https_proxy
###################

echo "START job:" $(date)
echo "HOSTNAME:" $(hostname)

# Setup the working directory directly in RAM and move in it
WORKDIR=/storage/scratch/$USER/${SLURM_JOB_ID}_${SLURM_ARRAY_TASK_ID}
mkdir -p $WORKDIR
cd $WORKDIR

echo "move in $WORKDIR"

#### The assembly index for BWA ####
### Do not forget the file ".amb" ".ann" ".bwt" ".pac" ".sa" !
ASM_FASTA="assembly.2500.fa.gz"
INDEX=${WORKDIR}/$ASM_FASTA

#### Get all samples for the current condition
SAMPLES=(AA AE AI AM AR AT AU AY BC BD BJ BK BR BS BY BZ)

#### The current sample!
CUR_SAMPLE=${SAMPLES[${SLURM_ARRAY_TASK_ID}]}
SM=CIM_${CUR_SAMPLE}

### Setup the directory with the results and MOVE in it
RESDIR=$WORKDIR/$SM
mkdir $RESDIR

echo "Move in $RESDIR"
cd $RESDIR

#### Copy all the reads
echo "Copy the reads in $RESDIR"

#####
#
#
# Here copy ALL fastq files for your sample
#
# Multiple libraries of the same sample can be used
#
#
#####

#### For the samples that have two or more libraries, ie multiple R1 and R2 files
NUM_FILE=$(ls ${SM}*.fastq.gz | wc -l)

#### Samples with mutiples pairs of read
if [[ $NUM_FILE -gt 2 ]]; then
    echo "Several pairs of read for sample ${CUR_SAMPLE}"
    echo "Map each pair separatelly, merge the results later"

    # Get all R1 ### THIS MIUST REQUIRE A BIT OF ADAPTATION TO FIT YOUR SAMPLES
    R1s=($(ls $RESDIR/CIM_${CUR_SAMPLE}*_?_1_*.fastq.gz))

    #### Process each pair, one at a time
    for read1 in ${R1s[@]}; do
        # Get the file for the read2

        ####
        #
        # The step bellow may require a bit of adaptation!
        #
        ####

        read2=$(echo $read1 | sed 's/\(CIM_[A-Z]\{6\}_[0-9]_\)1_/\12_/')

        # Get Rample Name and Read group to track the reads after the merge
        ID=$(basename $read1 | cut -f 1 -d ".")
        RG=$(echo "@RG\tID:${ID}\tSM:${SM}")

        # Mapping
        echo -e "Mapping files\n - ${read1}\n - ${read2}"
        bwa mem -R $RG -t $SLURM_CPUS_PER_TASK -o $ID.sam \
            $INDEX $read1 $read2

        # Flagstats + SAM to BAM
        echo "Flagstat of the raw mapping, then SAM to BAM"
        samtools flagstat ${ID}.sam >${ID}.flagstat

        samtools view -bS -F4 -q10 -o ${ID}.sam.raw.bam ${ID}.sam

        # Delete SAMfile
        echo "Unsorted BAM file for $SM -- $ID is ready, remove ${ID}.sam"
        rm ${ID}.sam
    done

    echo "Concatenate the unsorted BAM files"
    # Concat the BAMfiles
    samtools cat -o ${SM}.merge.bam *.sam.raw.bam

    echo "Sort by name (samtools sort -n) the merged BAM file"
    # Sort by name
    samtools sort -@4 -m 2G -n -o ${SM}.namesorted.bam ${SM}.merge.bam

    echo "Flagstat after the Samtools view"
    samtools flagstat ${SM}.namesorted.bam >${SM}.namesorted.flagstat

    # Remove
    echo "Clean sub-BAM and unsorted BAM files "
    rm ${SM}.merge.bam
    rm *.sam.raw.bam

#### Samples with a uniq pair of reads
else
    echo "1 pair for sample ${CUR_SAMPLE}"
    #### Get files of reads ### THIS MIUST REQUIRE A BIT OF ADAPTATION TO FIT YOUR SAMPLES
    read1=$(ls $RESDIR/CIM_${CUR_SAMPLE}*_?_1_*.fastq.gz)
    read2=$(ls $RESDIR/CIM_${CUR_SAMPLE}*_?_2_*.fastq.gz)

    ID=$(basename $read1 | cut -f 1 -d ".")
    RG=$(echo "@RG\tID:${ID}\tSM:${SM}")

    # Mapping
    echo -e "Mapping files\n - ${read1}\n - ${read2}"
    bwa mem -R $RG -t $SLURM_CPUS_PER_TASK -o ${SM}.sam $INDEX $read1 $read2

    # SAM to BAM and sort
    echo "Flagstat of the raw mapping"
    samtools flagstat ${SM}.sam >${SM}.flagstat

    echo "SAM to BAM and Sort by name (samtools sort -n) the merged BAM file"
    samtools view -bS -F4 -q10 -o ${SM}.sam.raw.bam ${SM}.sam
    samtools sort -@4 -m 2G -n -o ${SM}.namesorted.bam \
        ${SM}.sam.raw.bam

    echo "Flagstat after the Samtools view"
    samtools flagstat ${SM}.namesorted.bam >${SM}.namesorted.flagstat

    # Remove
    echo "Clean SAM and unsorted BAM files"
    rm ${SM}.sam
    rm ${SM}.sam.raw.bam
fi

#### SETUP the filtering and Coverage

### Mounting point within the singularity container
MP=/data

### Filter
###### Note for the usage of Singularity container:
###### I have to use the mounting point to READ a file WITHIN the container,
###### but I can redirect the result printed in stdout directly into a file!

echo "Filter with mSAMtools, min. of 95% identity and alignment of 75 nucl. min."
singularity exec -B $RESDIR:$MP $RESDIR/msamtools_coverm.sif \
    msamtools filter -b -l 75 -p 95 --besthit \
    ${MP}/${SM}.namesorted.bam >${SM}.namesorted.filter.bam

echo 'Sort by position (default) the BAM filtered and index'
samtools sort -@4 -m 2G -o ${SM}.namesorted.filter.sort.bam \
    ${SM}.namesorted.filter.bam
samtools index ${SM}.namesorted.filter.sort.bam

echo "Flagstat after the Msamtools"
samtools flagstat ${SM}.namesorted.filter.sort.bam >${SM}.namesorted.filter.sort.flagstat

echo "END computing:" $(date)

#### Move the results in the Bucket
echo "Save the results"
echo "Files that will be saved:"
ls -lh $RESDIR/*.flagstat
ls -lh $RESDIR/${SM}.namesorted.filter.sort.bam*

###
#
#
# MOVE THE RESULTS WHERE YOU WISH
#
#
###

echo "wait 10 s..."
sleep 10

rm -r $WORKDIR

echo "END job:" $(date)
