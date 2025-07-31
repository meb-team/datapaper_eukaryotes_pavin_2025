#!/bin/bash
echo "START :" $(date)
# General setup

## $1=path/to/workdir
## $2=Name of the MAG
## $3=Number of CPUs
## $4=MINCONTIG length
## $5=Augustus species
## $6=Fasta of the Proteins

if [ $# -ne 6 ]; then
    echo "This script needs 6 positionnal arguments:"
    echo -e "- \$1 = path/to/workdir , where you want to work"
    echo -e "- \$2 = the MAG Fasta file, must end by '.fa'"
    echo -e "- \$3 = Number of CPUs"
    echo -e "- \$4 = MINCONTIG length, see the GeneMark documentation"
    echo -e "- \$5 = Augustus species, see the Augustus documentation"
    echo -e "- \$6 = Proteins hints, as fasta"

    exit 1
fi

### ------- ### ------- ### ------- ### ------- ### ------- ### ------- ###
## Setup
WORKDIR=$1
CPU=$3
MINCTG=$4
SPECIES=$5
PROTS=$6

RESDIR=gene_call-${MAG}
mkdir $WORKDIR/$RESDIR

### Singularity image
PREDimg=augustus_gmove_genemark.sif
MP=/data
#### Get the license files from the container and add it in your $HOME
if [[ ! -f $HOME/.gm_key ]]; then
    singularity exec -B $(pwd):$MP $PREDimg \
        cp -r $MP/../tools/gm_key_64 $MP/gm_key
    mv gm_key $HOME/.gm_key
fi

### MAG name
MAG_FASTA=$2
MAG=$(basename $2 | sed "s/\.fa$//")

echo "Run the gene prediction for the MAG $MAG"

### Setup a tempdir, especially for "parallel"
TMPDIR=$(basename $(mktemp -u -p .))
mkdir $TMPDIR

echo "### Initial setup ok, result directory: $RESDIR"

### ------- ### ------- ### ------- ### ------- ### ------- ### ------- ###
# 2. Sample the proteins ==> reduce CPU time and disk usages
### Compute the size of the MAG: rm deflines + rm newlines + count characters
echo "Sub-sample proteins:" $(date)
MAG_SIZE=$(grep -v "^>" $MAG_FASTA | tr -d "\n" | wc -c)

subsample_prots() {
    local NUMPROT=$1
    ## Subsample protein, "-s 12546" is a seed to set the random generator
    singularity exec -B $(pwd):$MP $PREDimg seqtk sample \
        -s 12546 $MP/$PROTS $NUMPROT >sub_selected_proteins.fa
    PROTS=sub_selected_proteins.fa
    ## Keep their IDS
    grep '^>' $PROTS | tr -d '>' | sed "s/ .\+//" \
        >$RESDIR/sub_selected_proteins-ids.txt
}

## Subsample the umber of proteins to align according to the MAG size
### 1M proteins if MAG < 30Mb; 2M prots MAG between 30 and 60 Mb; 3M prot >60Mb
if [[ $MAG_SIZE -le 30000000 ]]; then
    echo '### Subsampling 1 000 000 proteins'
    subsample_prots 1000000
elif [[ $MAG_SIZE -le 60000000 ]]; then
    echo '### Subsampling 2 000 000 proteins'
    subsample_prots 2000000
else
    echo '### Subsampling 3 000 000 proteins'
    subsample_prots 3000000
fi

gzip $RESDIR/sub_selected_proteins-ids.txt

### ------- ### ------- ### ------- ### ------- ### ------- ### ------- ###
# 3. Generate hints for Augustus: MiniProt
echo "Start Miniprot" $(date)
MINIdir=$RESDIR/01_miniprot
HINTS=$MINIdir/${MAG}_Uni90METdb_selected_proteins.gff
mkdir -p $WORKDIR/$MINIdir

singularity exec -B $(pwd):$MP $PREDimg miniprot \
    -G10k -t $CPU -F69 --gff $MP/$MAG_FASTA $MP/$PROTS \
    >$HINTS

### Save space: removes in-place lines with alignment as PAF
echo "### Clean and prepare protein2genome GFF file for Augustus"
sed -i "/^##PAF/d" $HINTS

### Format miniprot outputs for Augustus
### "stop_codon" => "stop" / "CDS" => "CDSpart" / "source=P"
awk '{FS="\t";OFS=FS}{if ($3 != "mRNA") print $0}' <$HINTS |
    awk '{FS="\t"; OFS=FS}{ if ($3 == "CDS") $3 = "CDSpart"; print $0}' |
    grep -v "^#" >tempfile1
sort -t "$(printf '\t')" -k1,1 -k4,4n <tempfile1 >tempfile2

sed -e "s/stop_codon/stop/" -e "s/$/;source=P/" -e "1i ##gff-version 3" \
    <tempfile2 >${HINTS//.gff/}.augustus.gff
rm tempfile1 tempfile2

### ------- ### ------- ### ------- ### ------- ### ------- ### ------- ###
# 4. Augustus
## 4.1 First prediction with pre-built model
### Copy the config files from the image, to access them
echo "Start Augustus 1:" $(date)

singularity exec -B $(pwd):$MP $PREDimg \
    cp -r $MP/../opt/augustus-3.4.0/config $MP/augustus_config

AUGdir1=$RESDIR/02_augustus_ab_initio
mkdir -p $AUGdir1

### Split the assembly to run in parallel, bricks of ~1Mb
mkdir split

MAGsplit=$(basename ${MAG_FASTA//.fa/}.split)

singularity exec -B $(pwd):$MP $PREDimg splitMfasta.pl \
    $MP/$MAG_FASTA --outputpath=$MP/split --minsize=1000000

numSplits=$(ls split/$MAGsplit.*.fa | wc -l)

### Split the hints too ==> Reduce the amount of RAM used by 'parallel'
mkdir split_hints

### BEWARE, the range exapansion {1..$numSplits} does not works with all
### versions of BASH :(
for ((i = 1; i <= $numSplits; i++)); do
    grep "^>" split/$MAGsplit.$i.fa | sed 's/^>//' >part.$i.lst
    grep -w -f part.$i.lst ${HINTS//.gff/}.augustus.gff \
        >split_hints/hints.split.$i.gff
    rm part.$i.lst
done

### Adjust the number of threads.
### RUN ONLY 8 PARALLEL JOBS! I REACHED OOM ERROR WITH HIGHER NUMBER
### Check num of parallel job compared to: 1) CPU per taks, 2) Number of splits
declare -i AUG_THREADS=8
if [[ $AUG_THREADS -gt $CPU ]]; then
    $AUG_THREADS=$CPU
fi

if [[ $AUG_THREADS -gt $numSplits ]]; then
    AUG_THREADS=$numSplits
fi

echo "### Augustus 1, $AUG_THREADS Augustus jobs in parallel, for a total of" \
    "$numSplits jobs"

### Run Augustus in parallel
seq 1 $numSplits | parallel --tmpdir $TMPDIR -j $AUG_THREADS \
    "singularity exec -B $(pwd):$MP $PREDimg augustus \
    --species=$SPECIES --uniqueGeneId=true --progress=false --gff3=true \
    --outfile=$MP/$AUGdir1/$MAG.predi.01.{}.gff --genemodel=partial \
    --hintsfile=$MP/split_hints/hints.split.{}.gff \
    --extrinsicCfgFile=extrinsic.MP.cfg \
    --protein=off --UTR=off --AUGUSTUS_CONFIG_PATH=$MP/augustus_config \
    --softmasking=False --stopCodonExcludedFromCDS=false \
    $MP/split/$MAGsplit.{}.fa 2>augustus.{}.err"

### Merge the predictions
for ((i = 1; i <= $numSplits; i++)); do
    cat $AUGdir1/${MAG}.predi.01.$i.gff >>augustus_not_merged.gff
done

echo "### Augustus1 is done. Merge and clean"
singularity exec -B $(pwd):$MP $PREDimg join_aug_pred.pl <augustus_not_merged.gff \
    >$AUGdir1/${MAG}.predi.01.gff

# Clean
rm augustus.*err augustus_not_merged.gff $AUGdir1/${MAG}.predi.01.*.gff

## 4.2 Second prediction: generate a custom model
echo "Start Augustus 2:" $(date)

AUGdir2=$RESDIR/03_augustus_custom_model
mkdir -p $AUGdir2

### Get the length of flanking region
singularity exec -B $(pwd):$MP $PREDimg computeFlankingRegion.pl \
    $MP/$AUGdir1/${MAG}.predi.01.gff >tempfile 2>fuzzy.stderr

flank_region=$(cat tempfile | grep "flanking_DNA" | cut -d ":" -f 2 |
    sed "s/^ //" | cut -d " " -f 1)

singularity exec -B $(pwd):$MP $PREDimg gff2gbSmallDNA.pl \
    $MP/$AUGdir1/${MAG}.predi.01.gff $MP/$MAG_FASTA $flank_region \
    $MP/$AUGdir2/${MAG}.hints.02.gb

### Generate files for the new species
singularity exec -B $(pwd):$MP $PREDimg new_species.pl \
    --species=${MAG} --AUGUSTUS_CONFIG_PATH=$MP/augustus_config

### Train the custom gene model
echo "### Augustus 2: initial training of the new model"
singularity exec -B $(pwd):$MP $PREDimg etraining \
    --species=${MAG} --AUGUSTUS_CONFIG_PATH=$MP/augustus_config \
    $MP/$AUGdir2/${MAG}.hints.02.gb &>etraining.out

### Filter-out genes that raised an error
grep "in-frame stop codon" etraining.out | grep -v "^ExonModel" |
    cut -d " " -f7 | sed "s/:$//" >bad_genes.not_sorted
grep "exon doesn't end in stop codon" etraining.out |
    cut -d " " -f7 | sed "s/:$//" >>bad_genes.not_sorted
grep "coding length not a multiple of 3." etraining.out |
    cut -d " " -f7 | sed "s/:$//" >>bad_genes.not_sorted
sort -u bad_genes.not_sorted >$AUGdir2/${MAG}.badgenes.txt

echo "### Augustus 2: clean the hints from error-raising genes"
if [[ $(cat $AUGdir2/${MAG}.badgenes.txt | wc -l) -gt 0 ]]; then
    singularity exec -B $(pwd):$MP $PREDimg filterGenes.pl \
        $MP/$AUGdir2/${MAG}.badgenes.txt \
        $MP/$AUGdir2/${MAG}.hints.02.gb >ok.gb

    mv ok.gb $AUGdir2/${MAG}.hints.02.gb
fi

### Train the custom model without the "bad genes"
echo "### Augustus 2: training of the new model"
singularity exec -B $(pwd):$MP $PREDimg etraining \
    --species=${MAG} --AUGUSTUS_CONFIG_PATH=$MP/augustus_config \
    $MP/$AUGdir2/${MAG}.hints.02.gb &>etraining.out

rm tempfile fuzzy.stderr etraining.out bad_genes.not_sorted

### Save the current parameters in the result directory
cp -r augustus_config/species/${MAG} $AUGdir2/${MAG}_augustus_training

### Run the second Augustus prediction with the custom model
echo "### Augustus 2: gene prediction"
seq 1 $numSplits | parallel --tmpdir $TMPDIR -j $AUG_THREADS \
    "singularity exec -B $(pwd):$MP $PREDimg augustus \
    --species=$MAG --uniqueGeneId=true --progress=false --gff3=true \
    --outfile=$MP/$AUGdir2/$MAG.predi.02.{}.gff --genemodel=partial \
    --protein=on --codingseq=on --softmasking=False --UTR=off \
    --hintsfile=$MP/split_hints/hints.split.{}.gff \
    --extrinsicCfgFile=extrinsic.MP.cfg \
    --AUGUSTUS_CONFIG_PATH=$MP/augustus_config \
    $MP/split/$MAGsplit.{}.fa 2>augustus.{}.err"

echo "### Augustus1 is done. Merge and clean"
if [[ -f augustus_not_merged.gff ]]; then
    rm augustus_not_merged.gff
fi

### Merge the predictions
for ((i = 1; i <= $numSplits; i++)); do
    cat $AUGdir2/${MAG}.predi.02.$i.gff >>augustus_not_merged.gff
done

singularity exec -B $(pwd):$MP $PREDimg join_aug_pred.pl \
    <augustus_not_merged.gff >$AUGdir2/${MAG}.predi.02.gff

### Extract the sequences for the predicted genes: DNA and Amino-acid
singularity exec -B $(pwd):$MP $PREDimg getAnnoFasta.pl \
    $MP/$AUGdir2/${MAG}.predi.02.gff

### Clean
rm augustus.*.err augustus_not_merged.gff $AUGdir2/${MAG}.predi.02.*.gff
rm -r augustus_config split split_hints

### ------- ### ------- ### ------- ### ------- ### ------- ### ------- ###
# 5. ProtHint
echo "Start ProtHint:" $(date)
echo "### Run GenMark-ES for ProtHint"

PH=$RESDIR/04_prothint
mkdir -p $PH
#### To store the GenMark-ES results for ProtHint
mkdir -p $PH/GM_ES

singularity exec -B $(pwd):$MP $PREDimg gmes_petap.pl \
    --ES --cores $CPU --v --work_dir $MP/$PH/GM_ES \
    --min_contig $MINCTG --seq $MP/$MAG_FASTA

echo "### Run ProtHint"
singularity exec -B $(pwd):$MP $PREDimg prothint.py \
    --workdir $MP/$PH --threads $CPU \
    --geneMarkGtf $MP/$PH/GM_ES/genemark.gtf \
    $MP/$MAG_FASTA $MP/$PROTS
echo "### ProtHint done"

### ------- ### ------- ### ------- ### ------- ### ------- ### ------- ###
# 6. GeneMark
echo "Start GeneMark EP+:" $(date)

GM=$RESDIR/05_genemark_ep
mkdir -p $GM

## case where there are few evidences: DO NOT RUN in --EP mode!
num_evidences=$(grep -c "$(printf '\tIntron\t')" $PH/prothint.gff)
if [[ $num_evidences -ge 100 ]]; then
    echo "Run GeneMark WITH evidences"
    singularity exec -B $(pwd):$MP $PREDimg gmes_petap.pl \
        --cores $CPU --v --format GFF3 --work_dir $MP/$GM \
        --min_contig $MINCTG --EP $MP/$PH/prothint.gff \
        --evidence $MP/$PH/evidence.gff --seq $MP/$MAG_FASTA
else
    echo "Run GeneMark WITHOUT evidences"
    singularity exec -B $(pwd):$MP $PREDimg gmes_petap.pl \
        --ES --cores $CPU --v --format GFF3 --work_dir $MP/$GM \
        --min_contig $MINCTG --seq $MP/$MAG_FASTA
fi

# Clean GeneMark and ProtHint directories
rm -r $PH/diamond $PH/Spaln/ $PH/nuc.fasta $PH/seed_proteins.faa
rm -r $PH/GM_ES/data $PH/GM_ES/run $PH/GM_ES/output
rm -r $GM/run $GM/data $GM/output

### ------- ### ------- ### ------- ### ------- ### ------- ### ------- ###
# 7. Final Clean
echo "End - Last cleaning:" $(date)

## Miniprot hints: remove file used by Augustus, Gzip the original file
gzip $HINTS
rm $MINIdir/*.augustus.gff

## Gzip other files
gzip $AUGdir1/${MAG}.*
gzip $AUGdir2/${MAG}.*

echo "END computing:" $(date)

echo "wait 10 s..."
sleep 10

echo "END job:" $(date)
