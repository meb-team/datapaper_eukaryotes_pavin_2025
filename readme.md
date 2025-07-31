# Script and code snippets

Content :

- Singularity containers recipes
- Metagenomic co-assembly pre-processing
- MAG gene prediction
- Mapping public metagenomes

## Singularity containers recipes

To build a container, use the following command

```bash
sudo singularity build IMAGE.sif IMAGE.def
```

All recipes, "_.def_" files, are present in the directory
`singularity_containers`:

- `msamtools_coverm.def` mostly used to clean mapping
- `mmseqs_avx2.def`: latest version of MMseqs, for "modern" CPUs
- `mmseqs.def` : latest version of MMseqs, for CPUs that do not support _avx2_
  instructions
- `augustus_gmove_genemark.def` : Augustus and GeneMark, **requires extra steps, see bellow**
-

### Add GeneMark in the dedicated container

This tool requires a license file. It is free to get one but it is so anoying...

1. visit [http://exon.gatech.edu](http://exon.gatech.edu/genemark/license_download.cgi)
2. SELECT "_GeneMark-ES/ET/EP+_", and "_ LINUX 64 kernel 3.10 - 5_"
3. fill the form at the bottom of the page. Tips, you can use a random name,
   affiliation and a disposable e-mail address

Make sure you have files `gm_key_64.gz` and `gmes_linux_64_4.tar.gz` next to
the Singularity `augustus_gmove_genemark.def` file to built the image.

## Metagenomic co-assembly pre-processing

Reference of tools used below:

- [Anvi'o v8](https://github.com/merenlab/anvio) => Install it with _Conda_
  as all required tools are installed in the same environment
- [Seqtk](https://github.com/lh3/seqtk/)
- [Metabat2](https://bitbucket.org/berkeleylab/metabat/src/master/)

The process described bellow was performed on all co-assemblies

```bash
CPUS=8 # a number of CPU to run steps that uses multi-threading
ASM=assembly.fa.gz # Your metagenomic co-assembly here, MUST ENDs by ".fa.gz"

# Drop contigs lower thant 2.5 kb
seqtk seq -L 2500 $ASM > ${ASM//.fa.gz}.2500.fa
gzip -k ${ASM//.fa.gz}.2500.fa  # keep original file for the next step

# Prepare t(he contig database
mkdir contigs_db

ctg=contigs_db/${ASM//.fa.gz}.2500.db

anvi-gen-contigs-database -f ${ASM//.fa.gz}.2500.fa \
    -T 4 -n ${ASM//.fa.gz}.2500 -o $ctg

rm ${ASM//.fa.gz}.2500.fa # This was was only used to generate anvi'o contigs db

# Search single copy core genes sets for bacteria, archaea and eukaryotes
anvi-run-hmms -c $ctg -I Bacteria_71,Protista_83,Archaea_76 -T $CPUS
```

Then the backmapping, _i.e._ mapping of the samples that were co-assembled back
on the assembly.

```bash
## Prepare the index of the assembly
bwa index assembly.2500.fa.gz
```

An example of script that was used for this step. It was written to run on a HPC
system controlled by _SLURM_ :

`assembly_binning_refinement/backmapping_example.sbatch`

Then after the backmapping of each samples on the assembly, _anvi'o_ has to
"profile" each sample and then merge these profiles. Let's assume you put all
_BAM_ files in the directory `mapping/`:

```bash
mkdir single_profiles/

for sample in mapping/*.bam
do
    NAME=$(basename $fsample | sed "s/\.bam$//") # get the sample name
    # profile it
    anvi-profile -c $ctg -T $CPUS -i $sample -o single_profiles/${NAME} -S $NAME
done

# Merge individual profiles in a single "merged" profile
mkdir merged_profile

# Linkage will not be performed as there are >25 000 splits (contigs)
anvi-merge -c $ctg -S MyAssemby -o merged_profile \
    single_profiles/*/PROFILE.db

## Uncomment the next line to remove single profiles and save disk space
#rm -r single_profiles/*
```

Everything is ready to run the binning. _Metabat2_ must be installed on your machine :

```bash
# This is a wrapper for Matabat2
anvi-cluster-contigs -p merged_profile/PROFILE.db -c $ctg \
    -C METABAT2 -T $CPUS --driver metabat2 --just-do-it

# This will summarise the binning.
anvi-estimate-genome-completeness -c $ctg -p merged_profile/PROFILE.db \
    -C METABAT2 -o binning_metabat2.txt
```

The files above gives a summary of each bin, and also estimate the to which
domain it belongs. This is results is based on single copy core genes sets searched
earlier. To refine an interesting bin, via the _anvi'o_ interactive interface,
use

```bash
# here we want to visualise the bin "METABAT__152"
anvi-refine -c $ctg -p merged_profile/PROFILE.db -C METABAT2 \
    -b METABAT__152
```

For more details on this step, please visit [anvio.org](https://anvio.org) and/or
read the this [tutorial about `anvi-refine`](https://merenlab.org/2015/05/11/anvi-refine/).

## MAG gene prediction

This requires two steps. First search protein hints from the public databases to
help the gene prediction tools, and then in a second time run the gene prediction.

### Identify proteins from close organisms, in the databases

#### Get the data from [METdb](https://metdb.sb-roscoff.fr/)

The data are structured by organisms and seems hard to extract alone. The script
`mag_gene_prediction/download_metdb.sh` will do the job, except for the following
organisms for which the protein were downloaded manualy :

- _Fabrice-METDB_00194_
- _MMETSP-METDB_00388_
- _MMETSP-metdb_00357_ (no protein prediction)
- _MMETSP-METDB_00415_
- _MMETSP-METDB_00015_

It is important to know that:

- _METDB-00195_ and _METDB-00336_ have redundant protein names, this must be fixed
- _RCC-METDB_00336_ and _RCC-METDB_00485_ comes from the same organisms, and like
  above, some proteins have redundant name, which must be fixed.

#### Prepare the `mmseqs` protein database

Create the database for METdb data :

```bash
mmseqs createdb data/*/*.fa.gz database/METdbDB
```

Then download and merge with UniRef90 :

```bash
# download the database
mmseqs databases UniRef90 UniRef90DB tmp

# concatenate UniRef90 + METdb
mmseqs concatdbs database/METdbDB uniref90DB \
    database/Uniref90andMETdbDB
mmseqs concatdbs database/METdbDB_h uniref90DB_h \
    database/Uniref90andMETdbDB_h

# Clean UniRef90 to save disk space, if you want
# rm uniref90DB*
```

#### Align proteins _vs_ MAGs

Second step, identify candidate protein that could serve as protein hint.

All steps are described in the script `mag_gene_prediction/identify_candidate_proteins.sbatch`

This script returns a fasta file with proteins that will serve as hints for the
gene prediction

#### Gene prediction : _Augustus_ and _GeneMark_

Run the script `mag_gene_prediction/run_single_MAG_gene_prediction.sh`.
The Singularity image `augustus_gmove_genemark.sif` must be in the same directory
as the script! This script get the GeneMark license key from the Singularity
container and copy it in your home directory as `$HOME/.gm_key` if it is not
present.

The parameters _MINCTG_ (`$4`) and _SPECIES_ (`$5`) are respectivelly use by
GeneMark and Augustus. Please check their documentation for more information

All results are stored in the directory `gene_call-<NameofTheMAG>`.

## Mapping public metagenomes

### The mapping

Example of script, `mapping_public_metagenomes/mapping_public_metagenome.sbatch`.
It was designed to run on a HPC cluster with the SLURM scheduler, more precisely
with a job array. It takes as input a file with a list of SRA run identifiers,
one per line.
The script requires a _Conda_ environment called _fastp_, which contains the tool
[fastp](https://github.com/OpenGene/fastp/), and the Singularity container
`msamtools_coverm.sif` (_cf the Singularity container section_).

### Compute MAG/SAG detection

This script transforms a number of read per contig into a number of read per MAG.
It also compute the MAG detection (breadth of coverage).  
All informations used come from the `samtools coverage` command ran on a _BAM_ file.

This scripts needs :

- a two-columns tab file that link a contig name to its MAG, (with no header!!) e.g.

```
MAG_AS_0007_000000000003	MAG_AS_0007
MAG_AS_0007_000000000005	MAG_AS_0007
MAG_AS_0005_000000004006	MAG_AS_0005
MAG_AS_0005_000000004007	MAG_AS_0005
MAG_AS_0005_000000004008	MAG_AS_0005
MAG_OS_0050_000000004639	MAG_OS_0050
MAG_OS_0050_000000004640	MAG_OS_0050
```

- a two-columns tab file with the length of each MAG, as follow (with headers!)

```
MAG	length
mag000	5487655
mag005	25687963
```

Outputs two files :

- `detection_per_MAG.tsv`
- `sum_read_per_MAG.tsv`

### UMAP

This is an example script. It requires a table with the number of reads mapped
par MAG per metagenome, the MAG taxonomy, the MAG detection per metagenome,
and a metadata table of the metagenomes.
