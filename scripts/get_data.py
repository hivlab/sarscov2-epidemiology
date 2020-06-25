from Bio import Entrez
import os


def batch(iterable, n=1):
    l = len(iterable)
    for ndx in range(0, l, n):
        yield iterable[ndx : min(ndx + n, l)]


Entrez.email = snakemake.params["email"]
Entrez.api_key = snakemake.params["api_key"]
Entrez.max_tries = 6
handle = Entrez.esearch(db="nucleotide", retmax=10000, term="SARS-CoV-2", idtype="acc")
record = Entrez.read(handle)
handle.close()
acc = record["IdList"]
with open(snakemake.output[0], "w+") as h:
    for b in batch(acc, 20):
        handle = Entrez.efetch(
            db="nucleotide", id=",".join(b), rettype="gb", retmode="text"
        )
        h.writelines(handle)
        handle.close()
