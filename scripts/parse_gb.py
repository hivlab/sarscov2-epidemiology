from Bio import SeqIO
from collections import OrderedDict, Counter
import pandas as pd
from datetime import datetime


def fix_date(x):
    try:
        d = datetime.strptime(x, "%d-%b-%Y")
        fixed = datetime.strftime(d, "%Y-%m-%d")
    except Exception:
        fixed = x
        pass
    return fixed


ord_list = []
with open("data/sequences.fasta", "w") as fasta_handle:
    for seq_record in SeqIO.parse("data/sequences.gb", "genbank"):
        q = seq_record.features[0]
        ids = {"id": seq_record.id, "description": seq_record.description}
        qualifiers = {k: v[0] for k, v in q.qualifiers.items()}
        ids.update(qualifiers)
        ord_list.append(ids)
        fasta_handle.write(
            ">{} {}\n{}\n".format(seq_record.id, seq_record.description, seq_record.seq)
        )

df = pd.DataFrame(ord_list, columns=ord_list[0].keys())
df_rename1 = df.rename(columns={"strain": "strain2"})
df_renamed = df_rename1.rename(
    columns={"id": "accession", "organism": "virus", "collection_date": "date"}
)
df_renamed["strain"] = (
    df_renamed["description"]
    .str.replace("Severe acute respiratory syndrome coronavirus 2 (isolate )?", "")
    .str.replace(", complete genome", "")
)
df_renamed["date"] = df_renamed["date"].apply(lambda x: fix_date(x))
new = df_renamed["country"].str.split(": ", expand = True)
df_renamed["division"] = new[1]
df_renamed["country"] = new[0]
df_renamed.to_csv("data/metadata.tsv", sep="\t")
