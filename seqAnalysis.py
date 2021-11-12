# Source of documentation: 
# https://biopython.org/wiki/

# Import BioPython modules
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO

# Get a DNA sequence
sample_dna = Seq("TAGTGTATTGACATGATAGAAGCACTCTACTATATTCTCAATAGGTCCACG")
sample_dna_name = "tatabox"
record = SeqRecord(sample_dna, id="{sample_dna_name}", annotations={"Molecule Type": "DNA"})
# print(sample_dna)
sample_dna_info = dict()
for line in str(record).split('\n'):
    lineOfRecord = str(line).split(":")
    if len(lineOfRecord) == 2:
        attribute = lineOfRecord[0].strip()
        detail = lineOfRecord[1].strip()
        sample_dna_info[attribute] = detail
    elif lineOfRecord[0].startswith("S"):
        unformatted_seq = str(lineOfRecord)
        sample_dna_info["Sequence"] = unformatted_seq[7:-4]
    else:
        # Add annotation to the info dict
        molecule_type = lineOfRecord[0]
        newEntries = molecule_type.split("=")
        print(f"new entries are {newEntries}")
        sample_dna_info[newEntries[0][1:]] = newEntries[1]

# DNA Transcript 
sample_rna = sample_dna.transcribe()
# print(sample_rna)

# Translate into protein
sample_protein = str(sample_rna.translate())
sample_dna_info["Amino acid sequence"] = sample_protein
print(sample_protein)
print(sample_dna_info)
# print(sample_protein)