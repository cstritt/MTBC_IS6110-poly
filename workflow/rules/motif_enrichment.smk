
rule motif_prep:
    """
    Build weighted foreground FASTAs (5' and 3') and a background FASTA.

    Use all insertions with mappable reference positions from detettore6110.
    Join the 5' and 3' cFS of each mapped insertion and cut to the same length, 
    such that the insertion site is in the middle. 

    The foreground sequences are weighted by the number of independent IS births
    at that site (from ASR site_summary).

    The background FASTA is a random sample of reference sequences from the 
    reference genome, with the same length as the sequences in the foreground fasta.

    """
    input:
        reference_insertions = outdir + '/detettore/ALL_reference_insertions.tsv',
        anchors_tsv = outdir + '/detettore/ALL_anchors.tsv',
        births_5      = outdir + '/detettore/AS_5prime.site_summary.tsv',
        births_3      = outdir + '/detettore/AS_3prime.site_summary.tsv',
        anchor_map = outdir + '/detettore/anchor_map.tsv',
        reference     = config['detettore']['reference']
    output:
        fg   = outdir + '/motif/foreground.fasta',
        bg    = outdir + '/motif/background.fasta'
    params:
        min_length = config.get('motif', {}).get('min_length', 50),
        max_length = config.get('motif', {}).get('max_length', 150),
        max_weights = config.get('motif', {}).get('max_weights', 10),
        max_total_seqs = config.get('motif', {}).get('max_total_seqs', 1000),
        seed         = config.get('motif', {}).get('seed',         42)
    log:
        outdir + '/logs/motif_prep.log'
    conda:
        '../envs/motif.yaml'          # biopython, pandas
    script:
        '../scripts/motif_prep.py'


rule motif_streme:
    """Run STREME without background. In this case foreground sequences
    are shuffled test for enrichment."""
    input:
        fg = outdir + '/motif/foreground.fasta'
    output:
        xml = outdir + '/motif/streme/streme.xml',
        txt = outdir + '/motif/streme/streme.txt'
    params:
        outdir   = outdir + '/motif/streme',
        minw     = config.get('motif', {}).get('minw',  4),
        maxw     = config.get('motif', {}).get('maxw',  20),
        thresh   = config.get('motif', {}).get('thresh', 0.05),
        nmotifs  = config.get('motif', {}).get('nmotifs', 10),
        extra    = config.get('motif', {}).get('streme_extra', '')
    log:
        outdir + '/logs/motif_streme_shuffled.log'
    conda:
        '../envs/motif.yaml'          # meme-suite
    shell:
        """
        streme \
            --p   {input.fg}      \
            --oc  {params.outdir} \
            --dna                 \
            --minw  {params.minw} \
            --maxw  {params.maxw} \
            --thresh {params.thresh} \
            --nmotifs {params.nmotifs} \
            {params.extra}        \
        > {log} 2>&1
        """


rule motif_streme_bg:
    """Run STREME on the 5'-anchor side foreground vs. background."""
    input:
        fg      = outdir + '/motif/foreground.fasta',
        bg      = outdir + '/motif/background.fasta'
    output:
        xml = outdir + '/motif/streme_bg/streme.xml',
        txt = outdir + '/motif/streme_bg/streme.txt'
    params:
        outdir   = outdir + '/motif/streme_bg',
        minw     = config.get('motif', {}).get('minw',  4),
        maxw     = config.get('motif', {}).get('maxw',  12),
        thresh   = config.get('motif', {}).get('thresh', 0.05),
        nmotifs  = config.get('motif', {}).get('nmotifs', 10),  
        extra    = config.get('motif', {}).get('streme_extra', '')
    log:
        outdir + '/logs/motif_streme.log'
    conda:
        '../envs/motif.yaml'          # meme-suite
    shell:
        """
        streme \
            --p   {input.fg}      \
            --n   {input.bg}      \
            --oc  {params.outdir} \
            --dna                 \
            --minw  {params.minw} \
            --maxw  {params.maxw} \
            --thresh {params.thresh} \
            --nmotifs {params.nmotifs} \
            {params.extra}        \
        > {log} 2>&1
        """


rule motif_parse:
    """Parse both STREME XML outputs into a single summary TSV."""
    input:
        xml = outdir + '/motif/streme/streme.xml',
        xml_bg = outdir + '/motif/streme_bg/streme.xml'
    output:
        tsv = outdir + '/motif/motif_summary.tsv'
    log:
        outdir + '/logs/motif_parse.log'
    conda:
        '../envs/motif.yaml'
    script:
        '../scripts/motif_parse.py'



rule extract_hotspot_regions:
    """
    Extract FASTA sequences for hotspot and coldspot regions.
    The regions input contains the coordinates of genic and intergenic regions in the reference genome.

    """
    input:
        hotspot_summary = outdir + '/detettore/hotspot_summary.tsv',
        reference = config['detettore']['reference'],
        regions = config['annotations']['mtbc0_regions'] 
    output:
        hotspot_seqs = outdir + '/motif/hotspot_sequences.fasta',
        coldspot_seqs = outdir + '/motif/coldspot_sequences.fasta'
    params:
        hotspot_threshold = config['motif']['hotspot_threshold'],
        outdir = outdir + '/motif'
    log:
        outdir + '/logs/extract_hotspot_regions.log'
    conda:
        '../envs/motif.yaml'          # biopython, pandas
    script:
        '../scripts/extract_hotspot_regions.py'







rule streme_to_meme:
    input:
        streme = outdir + '/motif/streme/streme.xml'
    output:
        meme = outdir + '/motif/streme/streme.meme'
    conda:
        '../envs/motif.yaml'
    shell:
        """
        streme2meme {input.streme} > {output.meme}
        """

rule run_fimo_hotspots:
    input: 
        motifs = outdir + '/motif/streme/streme.meme',
        hotspot_seqs = outdir + '/motif/hotspot_sequences.fasta'
    output:
        pass
    conda: '../envs/motif.yaml'
    shell:
        """
        fimo --oc hotspot_fimo {input.motifs} {input.hotspot_seqs}
        """

rule run_fimo_coldspots:
    input: 
        motifs = outdir + '/motif/streme/streme.meme',
        coldspot_seqs = outdir + '/motif/coldspot_sequences.fasta'
    output:
        pass
    conda: '../envs/motif.yaml'
    shell:
        """
        fimo --oc hotspot_fimo {input.motifs} {input.coldspot_seqs}
        """

