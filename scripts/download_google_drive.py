import gdown
import tempfile
import pandas as pd
from Bio import SeqIO
import io

# Download GISAID cov2020 acknowledgements file from Google drive
excel_url = "https://drive.google.com/uc?id=1g85nEcuiVnmO75Hh8yWAty5uW8P_RSiR"

# Download sars-cov-2 genomic sequences fasta file
fasta_url = "https://drive.google.com/uc?id=1ds-tnGObN7wIRDb12yq1Ey0-Ch8sR5V4"

with tempfile.TemporaryFile() as tmp:
    gdown.download(excel_url, tmp, quiet=False)
    df = pd.read_excel(tmp, skiprows=[0, 1, 3])

df.columns = df.columns.str.strip().str.lower().str.replace(" ", "_")
df = df.set_index("accession_id")
country = ["Estonia", "Finland", "Latvia", "Russia", "Lietuva", "Sweden", "Norway"]
ids = []
for i in country:
    ids.append(df[df["virus_name"].str.contains(i)])

df_country = pd.concat(ids)
df_country = df_country[
    df_country["submitting_lab"]
    != "Charité Universitätsmedizin Berlin, Institute of Virology"
]

with open("data/metadata_gisaid.tsv", "w") as oh:
    df_country.to_csv(oh, sep="\t")

ids_list = df_country.index.tolist()
with tempfile.TemporaryFile() as tmp:
    gdown.download(fasta_url, tmp, quiet=False)
    tmp.seek(0)
    for record in SeqIO.parse(io.StringIO(tmp.read().decode("utf8")), "fasta"):
        try:
            id = record.id.split("|")[1]
        except IndexError:
            id = ""
        if id in ids_list:
            SeqIO.write(record, "data/sequences_gisaid.fasta", "fasta")
