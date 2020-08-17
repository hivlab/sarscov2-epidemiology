import os
from snakemake.logging import logger
import sys
from packaging import version
from snakemake.remote.S3 import RemoteProvider as S3RemoteProvider
S3 = S3RemoteProvider(access_key_id=os.environ.get("AWS_ACCESS_KEY_ID"), secret_access_key=os.environ.get("AWS_SECRET_ACCESS_KEY"))

MIN_AUGUR_VERSION = "7.0.2"

try:
    from augur.__version__ import __version__ as augur_version
except ModuleNotFoundError:
    logger.error("ERROR: Could not find augur. Follow installation instructions at https://nextstrain.org/docs/ and try again.")
    sys.exit(1)

if version.parse(augur_version) < version.parse(MIN_AUGUR_VERSION):
    logger.error("ERROR: Found version '%s' of augur, but version '%s' or greater is required" % (augur_version, MIN_AUGUR_VERSION))
    sys.exit(1)

rule all:
    input:
        auspice_json = "auspice/sarscov2.json"


include_strains = "config/include_strains.txt",
reference = "config/sarscov2_outgroup.gb",
auspice_config = "config/auspice_config.json"
min_length = 25000


# Fetch SARS-CoV2 data from NCBI
rule getdata:
    output:
        "data/sequences.gb"
    params:
        email = "taavi.pall@ut.ee",
        api_key = os.environ.get("NCBI_APIKEY")
    script:
        "scripts/get_data.py"


# Parse downloaded genbank sequences
# Filter by min length
rule parsegb:
    input:
        "data/sequences.gb"
    output:
        fasta = "data/global_sequences.fasta",
        metadata = "data/global_metadata.tsv"
    params:
        min_length = min_length
    script:
        "scripts/parse_gb.py"


# Parse reference genbank to fasta 
rule reference:
    input:
        reference
    output:
        fasta = "data/sarscov2_outgroup.fasta"
    script:
        "scripts/parse_gb.py"


# Merge local data to global data 
rule merge_metadata:
    input:
        our_metadata = S3.remote("sc2-consensus-seqs/our_metadata.tsv"),
        metadata = rules.parsegb.output.metadata
    output:
        metadata = "data/metadata.tsv"
    run:
        import pandas as pd
        md = pd.read_csv(input.metadata, sep = "\t")
        lmd = pd.read_csv(input.our_metadata, sep = "\t")
        concatenated = pd.concat([md, lmd], sort=False)
        concatenated.to_csv(output.metadata, sep = "\t")


rule merge_fasta:
    input:
        our_fasta = S3.remote("sc2-consensus-seqs/our_sequences.fasta"),
        fasta = rules.parsegb.output.fasta
    output:
        fasta = "data/sequences.fasta"
    shell:
        "cat {input.fasta} {input.our_fasta} > {output}"


# Fetch colors for country metadata
rule colors:
    input:
        "data/metadata.tsv"
    output:
        col = "config/colors.tsv",
        loc = "config/lat_longs.tsv"
    script:
        "scripts/country_colors.py" 


# Main nextstrain workflow
rule filter:
    message:
        """
        Filtering to
          - {params.sequences_per_group} sequence(s) per {params.group_by!s}
          - from {params.min_date} onwards
        """
    input:
        sequences = "data/sequences.fasta",
        metadata = "data/metadata.tsv",
        include = include_strains
    output:
        sequences = "results/filtered.fasta"
    params:
        group_by = "country",
        sequences_per_group = 100,
        min_date = 2018,
        min_length = min_length
    shell:
        """
        augur filter \
            --sequences {input.sequences} \
            --metadata {input.metadata} \
            --include {input.include} \
            --output {output.sequences} \
            --group-by {params.group_by} \
            --sequences-per-group {params.sequences_per_group} \
            --min-date {params.min_date} \
            --min-length {params.min_length}
        """


rule align:
    message:
        """
        Aligning sequences to {input.reference}
          - filling gaps with N
        """
    input:
        sequences = rules.filter.output.sequences,
        reference = rules.reference.output.fasta
    output:
        alignment = "results/aligned.fasta"
    threads: 4
    shell:
        """
        augur align \
            --sequences {input.sequences} \
            --reference-sequence {input.reference} \
            --output {output.alignment} \
            --nthreads {threads} \
            --fill-gaps
        """


rule mask:
    message:
        """
        Mask bases in alignment
          - masking {params.mask_from_beginning} from beginning
          - masking {params.mask_from_end} from end
          - masking other sites: {params.mask_sites}
        """
    input:
        alignment = rules.align.output.alignment
    output:
        alignment = "results/masked.fasta"
    log:
        "results/mask.log"
    params:
        mask_from_beginning = 100,
        mask_from_end = 50,
        mask_sites = "13402 18529 24389 24390"
    shell:
        """
        python3 scripts/mask-alignment.py \
            --alignment {input.alignment} \
            --mask-from-beginning {params.mask_from_beginning} \
            --mask-from-end {params.mask_from_end} \
            --mask-sites {params.mask_sites} \
            --output {output.alignment} 2>&1 | tee {log}
        """



rule tree:
    message: "Building tree"
    input:
        alignment = rules.mask.output.alignment
    output:
        tree = "results/tree_raw.nwk"
    threads: 4
    shell:
        """
        augur tree \
            --alignment {input.alignment} \
            --output {output.tree} \
            --nthreads {threads}
        """

# Clock rate and sd: https://www.sciencedirect.com/science/article/pii/S1567134820301829?via%3Dihub
rule refine:
    message:
        """
        Refining tree
          - estimate timetree
          - use {params.coalescent} coalescent timescale
          - estimate {params.date_inference} node dates
          - filter tips more than {params.clock_filter_iqd} IQDs from clock expectation
        """
    input:
        tree = rules.tree.output.tree,
        alignment = rules.mask.output.alignment,
        metadata = "data/metadata.tsv"
    output:
        tree = "results/tree.nwk",
        node_data = "results/branch_lengths.json"
    params:
        root = "best", 
        clock_rate = 0.0006,
        clock_std_dev = 0.00008,
        coalescent = "skyline",
        date_inference = "marginal",
        divergence_unit = "mutations",
        clock_filter_iqd = 4
    shell:
        """
        augur refine \
            --tree {input.tree} \
            --alignment {input.alignment} \
            --metadata {input.metadata} \
            --output-tree {output.tree} \
            --output-node-data {output.node_data} \
            --root {params.root} \
            --timetree \
            --clock-rate {params.clock_rate} \
            --clock-std-dev {params.clock_std_dev} \
            --coalescent {params.coalescent} \
            --date-inference {params.date_inference} \
            --divergence-unit {params.divergence_unit} \
            --date-confidence \
            --no-covariance \
            --clock-filter-iqd {params.clock_filter_iqd}
        """


rule ancestral:
    message: "Reconstructing ancestral sequences and mutations"
    input:
        tree = rules.refine.output.tree,
        alignment = rules.mask.output.alignment
    output:
        node_data = "results/nt_muts.json"
    params:
        inference = "joint"
    shell:
        """
        augur ancestral \
            --tree {input.tree} \
            --alignment {input.alignment} \
            --output-node-data {output.node_data} \
            --inference {params.inference} \
            --infer-ambiguous
        """


rule translate:
    message: "Translating amino acid sequences"
    input:
        tree = rules.refine.output.tree,
        node_data = rules.ancestral.output.node_data,
        reference = reference
    output:
        node_data = "results/aa_muts.json"
    shell:
        """
        augur translate \
            --tree {input.tree} \
            --ancestral-sequences {input.node_data} \
            --reference-sequence {input.reference} \
            --output-node-data {output.node_data} \
        """


rule traits:
    message: "Inferring ancestral traits for {params.columns!s}"
    input:
        tree = rules.refine.output.tree,
        metadata = "data/metadata.tsv"
    output:
        node_data = "results/traits.json",
    params:
        columns = "country"
    shell:
        """
        augur traits \
            --tree {input.tree} \
            --metadata {input.metadata} \
            --output-node-data {output.node_data} \
            --columns {params.columns} \
            --confidence
        """


rule clades:
    message: "Labeling clades as specified in config/clades.tsv"
    input:
        tree = rules.refine.output.tree,
        aa_muts = rules.translate.output.node_data,
        nuc_muts = rules.ancestral.output.node_data,
        clades = "config/clades.tsv"
    output:
        clade_data = "results/clades.json"
    shell:
        """
        augur clades --tree {input.tree} \
            --mutations {input.nuc_muts} {input.aa_muts} \
            --clades {input.clades} \
            --output {output.clade_data}
        """


rule recency:
    message: "Use metadata on submission date to construct submission recency field"
    input:
        metadata = "data/metadata.tsv"
    output:
        "results/recency.json"
    shell:
        """
        python3 scripts/construct-recency-from-submission-date.py \
            --metadata {input.metadata} \
            --output {output}
        """


rule export:
    message: "Exporting data files for auspice"
    input:
        tree = rules.refine.output.tree,
        metadata = "data/metadata.tsv",
        branch_lengths = rules.refine.output.node_data,
        traits = rules.traits.output.node_data,
        nt_muts = rules.ancestral.output.node_data,
        aa_muts = rules.translate.output.node_data,
        colors = rules.colors.output.col,
        lat_longs = rules.colors.output.loc,
        clades = rules.clades.output.clade_data,
        auspice_config = auspice_config,
        recency = rules.recency.output
    output:
        auspice_json = rules.all.input.auspice_json,
    shell:
        """
        augur export v2 \
            --tree {input.tree} \
            --metadata {input.metadata} \
            --node-data {input.branch_lengths} {input.traits} {input.nt_muts} {input.aa_muts} {input.clades} {input.recency} \
            --colors {input.colors} \
            --lat-longs {input.lat_longs} \
            --auspice-config {input.auspice_config} \
            --output {output.auspice_json}
        """

rule clean:
    message: "Removing directories: {params}"
    params:
        "results ",
        "auspice"
    shell:
        "rm -rfv {params}"
