# ruleorder: align>chunk

rule whole_genome:
    input:
        auspice_json = "auspice/Wuhan-BCov_full.json",

rule chunks:
    input:
        auspice_json = expand("auspice/Wuhan-BCov_s{x:05d}.json", x=range(0,30000,5000))

reference = "data/MK062183.gb",
def get_reference(w):
    print(w)
    if w.region=='full':
        return reference
    else:
        return f"results/{w.region}_reference.gb"


rule cat:
    input:
        "data/VIPR_SARS_like.fasta",
        "data/WH-Human_1.fasta"
    output:
        sequences = "results/all_sequences.fasta"
    shell:
        """
        cat {input} > {output.sequences}
        """

rule filter:
    input:
        sequences = rules.cat.output.sequences
    output:
        sequences = "results/filtered.fasta"
    run:
        from Bio import SeqIO
        with open(output.sequences, 'wt') as fh:
            for seq in SeqIO.parse(input.sequences, 'fasta'):
                if any([x in seq.description for x in ['ExoN1', 'FJ88', 'FV53']]):
                    continue
                else:
                    seq.id = seq.id.replace('|NA', '|')
                    seq.name = seq.id
                    seq.description = seq.id
                    SeqIO.write(seq, fh, 'fasta')


rule parse:
    message: "Parsing fasta into sequences and metadata"
    input:
        sequences = rules.filter.output.sequences
    output:
        metadata = "results/metadata_vipr.tsv",
        sequences = "results/sequences_vipr.fasta"
    params:
        fasta_fields = "strain labid x date host country y virus",
    shell:
        """
        augur parse \
            --sequences {input.sequences} \
            --output-metadata {output.metadata} \
            --output-sequences {output.sequences} \
            --fields {params.fasta_fields}
        """

rule parse_gisaid:
    message: "Parsing fasta into sequences and metadata"
    input:
        sequences = "data/gisaid.fasta"
    output:
        metadata =  "results/metadata_gisaid.tsv",
        sequences = "results/sequences_gisaid.fasta"
    params:
        fasta_fields = "strain labid date country division city host submitter",
    shell:
        """
        augur parse \
            --sequences {input.sequences} \
            --output-metadata {output.metadata} \
            --output-sequences {output.sequences} \
            --fields {params.fasta_fields}
        """

rule combine:
    input:
        m_gisaid = rules.parse_gisaid.output.metadata,
        s_gisaid = rules.parse_gisaid.output.sequences,
        m_vipr = rules.parse.output.metadata,
        s_vipr = rules.parse.output.sequences
    output:
        metadata="results/metadata.tsv",
        sequences="results/sequences.tsv"
    run:
        import pandas as pd
        d = pd.concat([pd.read_csv(x, sep='\t') for x in [input.m_gisaid, input.m_vipr]])
        d.to_csv(output.metadata, sep='\t')

        from Bio import SeqIO
        SeqIO.write(list(SeqIO.parse(input.s_gisaid, 'fasta'))\
                   +list(SeqIO.parse(input.s_vipr,'fasta')), output.sequences, 'fasta')


rule align:
    message:
        """
        Aligning sequences to {input.reference}
          - filling gaps with N
        """
    input:
        sequences = rules.combine.output.sequences,
        reference = reference
    output:
        alignment = "results/full_aligned.fasta"
    shell:
        """
        augur align \
            --nthreads 4 \
            --sequences {input.sequences} \
            --reference-sequence {input.reference} \
            --output {output.alignment} \
            --remove-reference
        """

rule chunk:
    input:
        reference = reference,
        alignment = rules.align.output.alignment
    output:
        alignment = "results/{region}_aligned.fasta",
        reference = "results/{region}_reference.gb"
    params:
        size=5000
    run:
        from Bio import SeqIO, AlignIO
        from Bio import SeqFeature

        start = int(wildcards.region[1:])
        stop = start+params.size
        print(input.reference[0])
        ref = SeqIO.read(input.reference[0], format='genbank')
        source = ref.features[0]
        source.location = SeqFeature.FeatureLocation(0,params.size)
        chunk_ref = ref[start:stop]
        if start:
            chunk_ref.features.insert(0,source)
        SeqIO.write(chunk_ref, output.reference, format="genbank")

        aln = AlignIO.read(input.alignment, 'fasta')
        SeqIO.write([s[start:stop] for s in aln], output.alignment, 'fasta')



rule tree:
    message: "Building tree"
    input:
        alignment = "results/{region}_aligned.fasta"
    output:
        tree = "results/{region}_tree_raw.nwk"
    shell:
        """
        augur tree \
            --nthreads 4 \
            --alignment {input.alignment} \
            --output {output.tree}
        """

rule refine:
    input:
        tree = rules.tree.output.tree,
        alignment = "results/{region}_aligned.fasta",
        metadata = rules.combine.output.metadata
    output:
        tree = "results/{region}_tree.nwk",
        node_data = "results/{region}_branch_lengths.json"
    params:
        coalescent = "opt",
        date_inference = "marginal",
        clock_filter_iqd = 4
    shell:
        """
        augur refine \
            --tree {input.tree} \
            --root HI553383 KY352407 \
            --alignment {input.alignment} \
            --metadata {input.metadata} \
            --output-tree {output.tree} \
            --output-node-data {output.node_data} \
            --coalescent {params.coalescent} \
        """

rule ancestral:
    message: "Reconstructing ancestral sequences and mutations"
    input:
        tree = rules.refine.output.tree,
        alignment = "results/{region}_aligned.fasta",
    output:
        node_data = "results/{region}_nt_muts.json"
    params:
        inference = "joint"
    shell:
        """
        augur ancestral \
            --tree {input.tree} \
            --alignment {input.alignment} \
            --output {output.node_data} \
            --inference {params.inference}
        """

rule translate:
    message: "Translating amino acid sequences"
    input:
        tree = rules.refine.output.tree,
        node_data = rules.ancestral.output.node_data,
        reference = get_reference
    output:
        node_data = "results/{region}_aa_muts.json"
    shell:
        """
        augur translate \
            --tree {input.tree} \
            --ancestral-sequences {input.node_data} \
            --reference-sequence {input.reference} \
            --output {output.node_data} \
        """

rule traits:
    message:
        """
        Inferring ancestral traits for {params.columns!s}
          - increase uncertainty of reconstruction by {params.sampling_bias_correction} to partially account for sampling bias
        """
    input:
        tree = rules.refine.output.tree,
        metadata = rules.combine.output.metadata
    output:
        node_data = "results/{region}_traits.json",
    params:
        columns = "country host",
        sampling_bias_correction = 3
    shell:
        """
        augur traits \
            --tree {input.tree} \
            --metadata {input.metadata} \
            --output {output.node_data} \
            --columns {params.columns} \
            --confidence \
            --sampling-bias-correction {params.sampling_bias_correction}
        """

title = "SARS-like Betacoronavirus phylogeny including novel coronavirus from Wuhan"\
        " using data generated by the China CDC, Shanghai Public Health Clinical, "\
        " Institute of Pathogen Biology, Chinese Academy of Medical Sciences & Peking Union Medical College"\
        " Wuhan Institute of Virology, Chinese Academy of Sciences shared via GISAID"
rule export:
    message: "Exporting data files for for auspice"
    input:
        tree = rules.refine.output.tree,
        metadata = rules.combine.output.metadata,
        branch_lengths = rules.refine.output.node_data,
        traits = rules.traits.output.node_data,
        nt_muts = rules.ancestral.output.node_data,
        aa_muts = rules.translate.output.node_data,
    output:
        auspice_json = "auspice/Wuhan-BCov_{region}.json"
    shell:
        """
        augur export v2 \
            --tree {input.tree} \
            --metadata {input.metadata} \
            --color-by-metadata host country \
            --title {title:q} \
            --maintainer "Richard Neher <https://neherlab.org>" \
            --node-data {input.branch_lengths} {input.traits} {input.nt_muts} {input.aa_muts} \
            --output {output.auspice_json}
        """

rule clean:
    message: "Removing directories: {params}"
    params:
        "results ",
        "auspice"
    shell:
        "rm -rfv {params}"
