from Bio import Entrez
import os

Entrez.email = "taavi.pall@ut.ee"
Entrez.api_key = os.environ.get("NCBI_APIKEY")
handle = Entrez.esearch(db="nucleotide", retmax=1000, term="SARS-CoV-2", idtype="acc")
record = Entrez.read(handle)
handle.close()
acc = record["IdList"]
handle = Entrez.efetch(db="nucleotide", id=",".join(acc), rettype="gb", retmode="text")
with open("data/sequences.gb", "w") as h:
    h.writelines(handle)
handle.close()
