# -*- coding: utf-8 -*-
"""
Created on Fri Apr 15 22:38:27 2022

@author: khahuynhnov
"""

########################################################################
# For the purposes of this workshop, don't worry about how this works, #
# a key point here is if you run this you'll all get the same file.    #
########################################################################

# Biopython's module to access the NCBI Entrez Programming Utilities
from Bio import Entrez

# The NCBI likes to know who is using their services in case of problems,
Entrez.email = "your.name.here@example.org"

accession = "CP053601"

print("Fetching %s from NCBI..." % accession)

# Return type "gbwithparts" matches "GenBank (full)" on the website
fetch_handle = Entrez.efetch("nuccore", id=accession, rettype="gbwithparts", retmode="text")

# Open an output file, and write all the data from the NCBI to it
with open(accession + ".gbk", "w") as output_handle:
    output_handle.write(fetch_handle.read())

print("Saved %s.gbk" % accession)

# Biopython's SeqIO module handles sequence input/output
from Bio import SeqIO

def get_cds_feature_with_qualifier_value(seq_record, name, value):
    """Function to look for tRNA feature by annotation value in sequence record.

    """
    # Loop over the features
    for feature in genome_record.features:
        if feature.type == "tRNA" and value in feature.qualifiers.get(name, []):
            return feature
    # Could not find it
    return None
########################################################################
genome_record = SeqIO.read("CP053601.gbk", "genbank")
cds_feature = get_cds_feature_with_qualifier_value(genome_record, "locus_tag", "HO397_01015")
print(cds_feature)
########################################################################

print(cds_feature.location)
gene_sequence = cds_feature.extract(genome_record.seq)
print("tRNA nucleotide sequence:")
print(gene_sequence)
print("Start codon is %s" % gene_sequence[:3])  # Python's way to get first three letters
print("Stop codon is %s" % gene_sequence[-3:])  # Python trick for last three letters

########################################################################

old_tags = ["HO397_01015","HO397_01020","HO397_01035","HO397_01085","HO397_01240","HO397_02495","HO397_03140","HO397_03145",
            "HO397_03150","HO397_03155","HO397_03160","HO397_03165","HO397_03170","HO397_03555","HO397_03560",
            "HO397_03565","HO397_03570","HO397_03575","HO397_03580","HO397_03585","HO397_04435","HO397_04875","HO397_05170",
            "HO397_06075","HO397_06085",
            "HO397_08210","HO397_08215","HO397_09450","HO397_09455","HO397_09460","HO397_09525","HO397_09535",
            "HO397_09565","HO397_09575","HO397_09590","HO397_10580","HO397_11350","HO397_11510","HO397_11515",
            "HO397_11535","HO397_11540","HO397_11545","HO397_11550","HO397_12460","HO397_12795","HO397_12800",
            "HO397_12805","HO397_12810","HO397_12815","HO397_13415","HO397_13420","HO397_13425","HO397_13600",
            "HO397_14125","HO397_14785","HO397_15295","HO397_15310","HO397_15770","HO397_15785","HO397_15790",
            "HO397_17160","HO397_17775","HO397_18400","HO397_18415","HO397_18420","HO397_18590","HO397_18595",
            "HO397_18600","HO397_18605","HO397_18950","HO397_18955","HO397_19535","HO397_19565","HO397_19570",
            "HO397_19575","HO397_19580","HO397_19730","HO397_20375","HO397_20520","HO397_20525","HO397_20530",
            "HO397_21055","HO397_21630","HO397_21635","HO397_21640"]

with open("nucleotides.fasta", "w") as nt_output:
    for tag in old_tags:
        print("Looking at " + tag)
        cds_feature = get_cds_feature_with_qualifier_value(genome_record, "locus_tag", tag)
        gene_sequence = cds_feature.extract(genome_record.seq)
     
        # Output FASTA records - note \n means insert a new line.
        # This is a little lazy as it won't line wrap the sequence:
        nt_output.write(">%s\n%s\n" % (tag, gene_sequence))

print("Done")


