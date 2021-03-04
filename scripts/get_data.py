from Bio import Entrez
import os
from http.client import IncompleteRead


def batch(iterable, n=1):
    l = len(iterable)
    for ndx in range(0, l, n):
        yield iterable[ndx : min(ndx + n, l)]


term = '("Severe acute respiratory syndrome coronavirus 2"[Organism] OR ("Severe acute respiratory syndrome coronavirus 2"[Organism] OR SARS-CoV-2[All Fields])) AND "Severe acute respiratory syndrome coronavirus 2"[Primary Organism] AND biomol_genomic[PROP]'

Entrez.email = snakemake.params["email"]
Entrez.api_key = snakemake.params["api_key"]
Entrez.max_tries = 6
handle = Entrez.esearch(
    db="nucleotide", retmax=snakemake.params["retmax"], term=term, idtype="acc"
)
record = Entrez.read(handle)
handle.close()
acc = record["IdList"]
with open(snakemake.output[0], "w+") as h:
    for b in batch(acc, snakemake.params["batches"]):
        handle = Entrez.efetch(
            db="nucleotide", id=",".join(b), rettype="gb", retmode="text"
        )
        try:
            h.writelines(handle)
        except IncompleteRead as e:
            print(e)
            continue
        handle.close()
