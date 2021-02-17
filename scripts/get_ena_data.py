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


def parse_analysis_xml(xml_string):
    root = ET.fromstring(xml_string)
    sample = {}
    for x in root.iter("ANALYSIS"):
        PRIMARY_ID = x.find("SAMPLE_REF/IDENTIFIERS/PRIMARY_ID").text
        SUBMITTER_ID = x.find("IDENTIFIERS/SUBMITTER_ID").attrib["namespace"]
        SEQUENCE_ASSEMBLY = x.find("ANALYSIS_TYPE/SEQUENCE_ASSEMBLY")
        NAME = SEQUENCE_ASSEMBLY.find("NAME").text
        TYPE = SEQUENCE_ASSEMBLY.find("TYPE").text
        PARTIAL = SEQUENCE_ASSEMBLY.find("PARTIAL").text
        COVERAGE = SEQUENCE_ASSEMBLY.find("COVERAGE").text
        PROGRAM = SEQUENCE_ASSEMBLY.find("PROGRAM").text
        PLATFORM = SEQUENCE_ASSEMBLY.find("PLATFORM").text
        MOL_TYPE = SEQUENCE_ASSEMBLY.find("MOL_TYPE").text
        analysis = {
            "SUBMITTER_ID": SUBMITTER_ID,
            "NAME": NAME,
            "TYPE": TYPE,
            "PARTIAL": PARTIAL,
            "COVERAGE": COVERAGE,
            "PROGRAM": PROGRAM,
            "PLATFORM": PLATFORM,
            "MOL_TYPE": MOL_TYPE,
        }
        FASTA = x.findall("FILES/FILE")
        files = {}
        for e in FASTA:
            files.update(
                {
                    e.attrib["filetype"]: e.attrib["filename"],
                    f"checksum_{e.attrib['filetype']}": e.attrib["checksum"],
                }
            )
        analysis.update(files)
        sample.update({PRIMARY_ID: analysis})
    return sample


def main(accessions, output, size=200):
    url = "https://www.ebi.ac.uk/ena/browser/api/xml/{accession}"
    headers = {"accept": "application/xml"}

    metadata = {}
    for chunk in chunked_iterable(accessions, size):
        r = requests.get(url.format(accession=",".join(chunk)), headers=headers)
        r.raise_for_status
        metadata.update(parse_analysis_xml(r.content))

    pd.DataFrame.from_dict(metadata, orient="index").to_csv(
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
    acc = pd.read_csv(args.infile, delimiter="\t")

    main(
        accessions=acc.analysis_accession,
        output=args.outfile,
    )
