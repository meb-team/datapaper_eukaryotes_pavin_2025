#!/bin/bash
#SBATCH --job-name=testMAG
#SBATCH --output=log/slurm-%j.out
#SBATCH --error=log/slurm-%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=64

# Setup the environment
module purge
module load gcc/8.1.0
module load rclone

unset http_proxy
unset https_proxy

echo "START job:" $(date)
echo "HOSTNAME:" $(hostname)

# General setup
### WorkDir
## if you have enough RAM, like 1TB, work in /dev/shm, it is much faster
WORKDIR=/dev/shm/$USER/${SLURM_JOB_ID}
## Otherviwe work in a scratch dir
#WORKDIR=/storage/scratch/$USER/${SLURM_JOB_ID}

mkdir -p $WORKDIR
echo "move in $WORKDIR"
cd $WORKDIR

# Prepare the output
RESDIR=00_mmseqs_identify_proteins
SEQS=seqs
mkdir $RESDIR $SEQS

### Get the fasta for the MAGs
echo "Copy contigs of MAGs in 'seqs/'"
#####
#
# Copy your MAG sequences in the directory "seqs"
#
#####

### Singularity container
MMimg=mmseqs-avx2.sif
MP=/data

#####
#
# Do not forget to copy the singularity container here
#
#####

## Get the database for MMSeqs
echo "Copy query database (METdb + UniRef90)" $(date)
qDB=Uniref90andMETdbDB
MEM=300 # in GB
TMPDIR=$(basename $(mktemp -u -p .))
mkdir -p $WORKDIR/$TMPDIR $WORKDIR/$MMdir

#####
#
# Do not forget to copy the database here singularity container here
#
#####

echo "Start MMseqs:" $(date)
### Create the target database:
singularity exec -B $WORKDIR:$MP $MMimg mmseqs \
    createdb $SEQS/*.fa $MP/MAGsDB

### MMseqs search similar to BLASTX (prot vs trans nucc DB),
### with sensitivity increment: 3 steps from 1 to 5 (=1, 3, 5)
singularity exec -B $WORKDIR:$MP $MMimg mmseqs \
    search $MP/$qDB $MP/MAGsDB $MP/MAGs_alnDB $MP/$TMPDIR \
    --start-sens 1 --sens-steps 3 -s 5 --split-memory-limit ${MEM}G \
    --threads $SLURM_CPUS_PER_TASK -e 0.00001 --min-length 50 --translate 1

### Extract the alignments
singularity exec -B $WORKDIR:$MP $MMimg mmseqs convertalis \
    $MP/$qDB $MP/MAGsDB $MP/MAGs_alnDB $MP/$MMdir/MAGs_aln.m8 \
    --threads $SLURM_CPUS_PER_TASK

### Split the results for all MAGs
##### CHECK THAT THE COMMAND BELLOW GET THE MAG IDs!!
MAG_LIST=($(find $SEQS -type 'f' -name '*.fa' | rev | cut -d "/" -f -1 |
    rev | cut -d '-' -f 2))

for MAG in ${MAG_LIST[@]}; do
    grep $MAG MAGs_aln.m8 | awk -F"\t" '{ if ($12 >= 200) print }' |
        cut -f 1 | uniq | sort | uniq >$RESDIR/${MAG}_${qDB}.list

    singularity exec -B $WORKDIR:$MP $MMimg mmseqs createsubdb \
        $MP/$RESDIR/${MAG}_${qDB}.list $MP/$qDB \
        $MP/${MAG}_${qDB}_filteredDB --id-mode 1

    singularity exec -B $WORKDIR:$MP $MMimg mmseqs convert2fasta \
        $MP/${MAG}_${qDB}_filteredDB $MP/$RESDIR/${MAG}_${qDB}_filteredseqs.fa

    rm ${MAG}_${qDB}_filteredDB*

done

echo "END computing:" $(date)

# Save the results
echo "Send results to ...."

echo "wait 10 s..."
sleep 10

# Last clean
rm -r $WORKDIR

echo "END job:" $(date)
