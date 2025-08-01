#!/usr/bin/env python3

import sys
import argparse
from ete3 import NCBITaxa

# Script input
parser = argparse.ArgumentParser(description="Parse NCBI taxID values")
parser.add_argument(
    "--input",
    help="tab-delimited BLAST result, with TAXIDS in THIRTEENTH (13th) column, "
    'like from the command `blastn -outfmt "6 std staxid"`',
    metavar="",
)
parser.add_argument(
    "--update",
    help="Update the ETE3 toolkit taxonomy database",
    action="store_true",
    default=False,
)

args = parser.parse_args()

# Print the help when zero argument
if len(sys.argv) == 1:  # In the case where nothing is provided
    parser.print_usage(file=sys.stderr)
    sys.exit(1)

# Load the ETE3 toolkit taxnomoy data
ncbi = NCBITaxa()

# Check update
if args.update:
    ncbi.update_taxonomy_database()
# Only update?
if not args.input:
    # This may be useless, but let's keep it to stay safe
    print("You only ask for update OR did not provide the '--input' " + "argument, bye")
    sys.exit(0)

# Read the BLAST file
with open(args.input, "r") as fi:
    # Create a result file
    with open(args.input + ".taxo", "w") as fo:

        # Read the results line by line
        for line in fi.readlines():
            taxids_t = line.split("\t")[12]  # Get the taxID from the thirteen field

            # Split TaxIDs if multiple (eg use of "staxids")
            if ";" in taxids_t:
                taxids = [taxids_t.split(";")]
            else:
                taxids = [taxids_t]

            # Find the lineages,
            lineages = []  # Complete lineage from the root to the "TaxID"
            for tax in taxids:
                # If "staxid" was used, a single value is present
                lineages.append(ncbi.get_lineage(int(tax)))

            # Make Search for certain TaxID to classify the results
            target = "other"  # Default value
            for lineage in lineages:
                val = [str(x) for x in lineage]

                if "2759" in val:
                    target = "Eukaryota"
                elif "2" in val:
                    target = "Bacteria"
                elif "1940" in val:
                    target = "Archaea"

            # Print the result by replacing the "taxID" with a value choosen
            # above
            print(
                line.split("\t")[0],
                target,
                "\t".join(line.rstrip().split("\t")[1:]),
                sep="\t",
                file=fo,
            )
