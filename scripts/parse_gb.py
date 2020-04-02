from Bio import SeqIO
from collections import OrderedDict, Counter
import pandas as pd
from datetime import datetime
import re


def fix_date(x):
    try:
        d = datetime.strptime(x, "%d-%b-%Y")
        fixed = datetime.strftime(d, "%Y-%m-%d")
    except Exception:
        fixed = x
        pass
    return fixed


sars = "Severe acute respiratory syndrome coronavirus 2 (isolate )?"
gen = ", complete genome"


ord_list = []
with open("data/sequences.fasta", "w") as fasta_handle:
    for seq_record in SeqIO.parse("data/sequences.gb", "genbank"):
        if len(seq_record.seq) > 28000:
            q = seq_record.features[0]
            ids = {"accession": seq_record.id, "description": seq_record.description, "length": len(seq_record.seq)}
            qualifiers = {k: v[0] for k, v in q.qualifiers.items()}
            if "isolate" in qualifiers.keys():
                ids.update({"strain": qualifiers["isolate"]})
            if "strain" not in ids.keys():
                ids.update({"strain": ids["accession"]})
            ids.update(qualifiers)
            ord_list.append(ids)
            fasta_handle.write(
                ">{}\n{}\n".format(ids["strain"], seq_record.seq)
            )


# Parsing metadata
df = pd.DataFrame(ord_list, columns=ord_list[0].keys()).set_index("strain", drop=False)
df_renamed = df.rename(
    columns={"organism": "virus", "collection_date": "date"}
)
df_renamed["date"] = df_renamed["date"].apply(lambda x: fix_date(x))
new = df_renamed["country"].str.split(": ?", expand = True)
df_renamed["division"] = new[1]
df_renamed["country"] = new[0]


# Fix missing country
jp = ["SARS-CoV-2/Hu/DP/Kng/19-027", "SARS-CoV-2/Hu/DP/Kng/19-020"]
df_renamed.loc[jp, "country"] = "Japan"


# Writing metadata to file
df_renamed.to_csv("data/metadata.tsv", sep="\t", index=False)
