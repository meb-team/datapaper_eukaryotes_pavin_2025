#!/bin/bash

########### Requirements
# Python3 and the ETE3 toolkit packages must be installed
# Also the tools seqtk and seqkit must be in your $PATH!

if [ ! $1 ]; then
	echo "Need the path to blastn results"
	exit
fi

if [ ! $2 ]; then
	echo "Need the sample name"
	exit
fi

# Activate Conda environment
CONDA_BASE=$(conda info --base)
source $CONDA_BASE/etc/profile.d/conda.sh
conda activate ete3

# Start
cat ${1%/}/result.blastn* >blastn_${2}

python parse_blastn_taxonomyID.py --input blastn_${2}

cut -f 1-3 <blastn_${2}.taxo >blastn_${2}.taxo.cut

grep -u 'Bacteria	' blastn_${2}.taxo.cut >blastn_${2}.taxo.cut.bacteria
grep -u 'Eukaryota	' blastn_${2}.taxo.cut >blastn_${2}.taxo.cut.eukaryota
grep -u 'other	' blastn_${2}.taxo.cut >blastn_${2}.taxo.cut.other
grep -u 'Homo' blastn_${2}.taxo.cut.eukaryota >blastn_${2}.taxo.cut.homo
grep -u 'Human' blastn_${2}.taxo.cut.eukaryota >blastn_${2}.taxo.cut.Hum

cat blastn_${2}.taxo.cut.homo blastn_${2}.taxo.cut.Hum >blastn_${2}.taxo.cut.human

grep -v 'Homo' blastn_${2}.taxo.cut.eukaryota >blastn_${2}.taxo.cut.eukaryota.exclude.homo
grep -v 'Human' blastn_${2}.taxo.cut.eukaryota.exclude.homo >blastn_${2}.taxo.cut.eukaryota.exclude.homo.human

cut -f 1-2 <blastn_${2}.taxo.cut.bacteria >blastn_${2}.taxo.cut.bacteria.cut2
uniq blastn_${2}.taxo.cut.bacteria.cut2 >blastn_${2}.taxo.cut.bacteria.cut2.uniq

cut -f 1-2 <blastn_${2}.taxo.cut.other >blastn_${2}.taxo.cut.other.cut2
uniq blastn_${2}.taxo.cut.other.cut2 >blastn_${2}.taxo.cut.other.cut2.uniq

cut -f 1-2 <blastn_${2}.taxo.cut.eukaryota.exclude.homo.human >blastn_${2}.taxo.cut.eukaryota.exclude.homo.human.cut2
uniq blastn_${2}.taxo.cut.eukaryota.exclude.homo.human.cut2 >blastn_${2}.taxo.cut.eukaryota.exclude.homo.human.cut2.uniq

cut -f 1-2 <blastn_${2}.taxo.cut.human >blastn_${2}.taxo.cut.human.cut2
uniq blastn_${2}.taxo.cut.human.cut2 >blastn_${2}.taxo.cut.human.cut2.uniq

# Write sequences as FASTA
seqtk subseq ${1%/}/contigs_${2}_sup_300.fasta blastn_${2}.taxo.cut.bacteria.cut2.uniq >${1%/}/liste_bacteria.fasta

seqtk subseq ${1%/}/contigs_${2}_sup_300.fasta blastn_${2}.taxo.cut.other.cut2.uniq >${1%/}/liste_other.fasta

seqtk subseq ${1%/}/contigs_${2}_sup_300.fasta blastn_${2}.taxo.cut.eukaryota.exclude.homo.human.cut2.uniq >${1%/}/liste_eukaryota.fasta

seqtk subseq ${1%/}/contigs_${2}_sup_300.fasta blastn_${2}.taxo.cut.human.cut2.uniq >${1%/}/liste_hum

# Clean human sequences
grep -u '>' ${1%/}/liste_hum >${1%/}/liste_hum2
sed -i -e 's/>//g' ${1%/}/liste_hum2
uniq ${1%/}/liste_hum2 >${1%/}/liste_hum3

seqtk subseq ${1%/}/contigs_${2}_sup_300.fasta ${1%/}/liste_hum3 >${1%/}/liste_human.fasta

rm ${1%/}/liste_hum
rm ${1%/}/liste_hum2
rm ${1%/}/liste_hum3

# Write Unknown sequences
cat ${1%/}/liste_bacteria.fasta ${1%/}/liste_eukaryota.fasta ${1%/}/liste_other.fasta ${1%/}/liste_human.fasta >${1%/}/homology
grep -u '>' ${1%/}/homology >${1%/}/homology2
sed -i -e 's/>//g' ${1%/}/homology2
uniq ${1%/}/homology2 >${1%/}/homology3

seqtk subseq ${1%/}/contigs_${2}_sup_300.fasta ${1%/}/homology3 >${1%/}/homology4

grep -u '>' ${1%/}/homology4 >${1%/}/homology5
sed -i -e 's/>//g' ${1%/}/homology5
uniq ${1%/}/homology5 >${1%/}/homology6

seqtk subseq ${1%/}/contigs_${2}_sup_300.fasta ${1%/}/homology6 >${1%/}/homology7

seqkit sort ${1%/}/homology7 >${1%/}/homology8

grep -vxFf ${1%/}/homology8 ${1%/}/contigs_${2}_sup_300.fasta >${1%/}/liste_unknown.fasta

rm blastn_${2}.taxo.cut.homo
rm blastn_${2}.taxo.cut.Hum
mv blastn_${2}.taxo* ${1}/
rm blastn_${2}
rm ${1%/}/homology*
rm ${1%/}/blastn_${2}*

echo "$2 done"
