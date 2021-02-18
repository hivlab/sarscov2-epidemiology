import pandas as pd
import requests
import xml.etree.ElementTree as ET
import itertools
import argparse


def chunked_iterable(iterable, size):
    it = iter(iterable)
    while True:
        chunk = tuple(itertools.islice(it, size))
        if not chunk:
            return
        yield chunk


def parse_sample_xml(xml_string):
    root = ET.fromstring(xml_string)
    sample = {}
    for x in root.iter("SAMPLE_SET"):
        samples = x.findall("SAMPLE")
        for s in samples:
            PRIMARY_ID = s.find("IDENTIFIERS/PRIMARY_ID").text
            sample_attrib = s.findall("SAMPLE_ATTRIBUTES/SAMPLE_ATTRIBUTE")
            tags = {}
            for e in sample_attrib:
                tags.update({e.find("TAG").text: e.find("VALUE").text})
            sample.update({PRIMARY_ID: tags})
    return sample


def main(accessions, output, size=200):
    url = "https://www.ebi.ac.uk/ena/browser/api/xml/{accession}"
    headers = {"accept": "application/xml"}
    metadata = {}
    for chunk in chunked_iterable(accessions, size):
        r = requests.get(url.format(accession=",".join(chunk)), headers=headers)
        r.raise_for_status
        metadata.update(parse_sample_xml(r.content))
    return pd.DataFrame.from_dict(metadata, orient="index").to_csv(
        output, index_label="PRIMARY_ID"
    )


if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--infile",
        "-I",
        dest="infile",
        required=True,
        help="Path to tsv file with ENA analysis accessions.",
    )
    parser.add_argument(
        "--outfile",
        "-O",
        dest="outfile",
        required=True,
        help="Path to csv output.",
    )
    args = parser.parse_args()
    acc = pd.read_csv(args.infile, delimiter=",")

    main(
        accessions=acc.PRIMARY_ID,
        output=args.outfile,
    )
