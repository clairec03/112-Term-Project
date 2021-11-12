# Source of documentation: 
# https://biopython.org/wiki/Alphabet

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
        sample_dna_info[newEntries[0]] = newEntries[1]
print(sample_dna_info)

# DNA Transcript 
sample_rna = sample_dna.transcribe()
# print(sample_rna)

# Translate into protein
sample_protein = sample_rna.translate()
# print(sample_protein)
'''
alanine - ala - A 
arginine - arg - R
asparagine - asn - N 
aspartic acid - asp - D 
cysteine - cys - C
glutamine - gln - Q
glutamic acid - glu - E
glycine - gly - G
histidine - his - H
isoleucine - ile - I
leucine - leu - L
lysine - lys - K
methionine - met - M
phenylalanine - phe - F
proline - pro - P
serine - ser - S
threonine - thr - T 
tryptophan - trp - W
tyrosine - tyr - Y
valine - val - V
'''