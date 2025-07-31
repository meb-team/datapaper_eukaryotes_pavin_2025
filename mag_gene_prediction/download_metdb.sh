#!/bin/bash

# This scripts download proteins from the METdb resource,

URL="https://metdb.sb-roscoff.fr/metdb/download/download"
DONE="done"
if [ -f "index.html" ]; then
    rm index.html
fi

echo "Discovering $URL"
wget --quiet $URL/

for d in $(cat index.html | grep "metdb" | cut -d '"' -f 2 | sed "s@/@@"); do
    # In theory, once we have this direectory, the files we want are stored
    # in the subdirectory "Annotations"
    #
    # I checked and we can also extract a source and a metdb_id from the name in $d
    # The source is "mmetsp" or "rcc" and the id has the format "metdb_\d{5}"
    sourceID=$(echo $d | cut -f 1 -d "-")
    metdbID=$(echo $d | cut -f 2 -d "-")

    # Create the output directories for the current dataset
    if [[ ! -d data/$sourceID ]]; then
        mkdir -p data/$sourceID
    fi

    ISDONE=$(grep "$d" $DONE)
    if [ -z "$ISDONE" ]; then
        echo ""
        echo "Fetching $URL/$d"
        echo "##########################################################"
        wget --quiet $URL/$d/Annotations/

        # Only two files are of interest. To save time, energy and disk usage
        # let's only download these two files:
        # ===> *transdecoder.pep
        # ===> *annotations_reduced
        fasta=$(grep "transdecoder.pep" index.html.1 | cut -d '"' -f 2)
        annot=$(grep "annotations_reduced" index.html.1 | cut -d '"' -f 2)

        echo "Downloading $URL/$d/Annotations/$fasta"
        wget --quiet $URL/$d/Annotations/$fasta

        echo "Downloading $URL/$d/Annotations/$annot"
        wget --quiet $URL/$d/Annotations/$annot

        echo "Parse the downloaded files"
        # Select only the proteins that have an annotation in the
        # "annotations_reduced" file
        sed "1d" $annot | cut -f 1 | uniq | sort | uniq >list_prot
        # Clean the header in the fasta file
        sed -i "s/ .\+$//" $fasta
        # Select the proteins
        seqtk subseq $fasta list_prot | sed "s/*$//" \
            >data/$sourceID/${sourceID}-${metdbID}-selected_proteins.fa

        # Move the results in the appropriate directory
        mv $annot data/$sourceID/${sourceID}-${metdbID}-selected_annot.tsv
        # Clean
        rm list_prot $fasta index.html.1

        echo "$d" >>$DONE

        # Wait 2 secondes, to avoid being black-listed by the remote server :)
        sleep 5
    else
        echo "$URL/$d already done :-)"
    fi

done

rm index.html*
