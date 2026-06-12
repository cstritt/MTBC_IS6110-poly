
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
        bg  = expand(outdir + '/motif/background_{i}.fasta', i=range(config.get('motif', {}).get('n_backgrounds', 10)))
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
        bg      = outdir + '/motif/background_{i}.fasta'
    output:
        xml = outdir + '/motif/streme_bg_{i}/streme.xml',
        txt = outdir + '/motif/streme_bg_{i}/streme.txt'
    params:
        outdir   = outdir + '/motif/streme_bg_{i}',
        minw     = config.get('motif', {}).get('minw',  4),
        maxw     = config.get('motif', {}).get('maxw',  12),
        thresh   = config.get('motif', {}).get('thresh', 0.05),
        nmotifs  = config.get('motif', {}).get('nmotifs', 10),  
        extra    = config.get('motif', {}).get('streme_extra', '')
    log:
        outdir + '/logs/motif_streme_{i}.log'
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
    """Parse all STREME XML outputs into a single summary TSV."""
    input:
        xml    = outdir + '/motif/streme/streme.xml',
        xml_bg = expand(outdir + '/motif/streme_bg_{i}/streme.xml',
                        i=range(config.get('motif', {}).get('n_backgrounds', 10)))
    output:
        tsv = outdir + '/motif/motif_summary.tsv'
    log:
        outdir + '/logs/motif_parse.log'
    conda:
        '../envs/motif.yaml'
    script:
        '../scripts/motif_parse.py'


rule fimo_top_motif:
    """Sensitive FIMO scan for top motif only."""
    input:
        motifs = outdir + '/motif/streme_bg_0/streme.xml'
    output:
        meme = outdir + '/motif/fimo_top_motif/top_motif.meme',
        tsv  = outdir + '/motif/fimo_top_motif/fimo.tsv'
    params:
        outdir     = outdir + '/motif/fimo_top_motif',
        motif_id   = config.get('motif', {}).get('top_motif_id', '1-TCTCAAAW'),
        thresh     = config.get('motif', {}).get('fimo_thresh_sensitive', 0.001),
        seqs       = config['annotations']['mtbc0_region_seqs']
    conda: '../envs/motif.yaml'
    shell:
        """
        meme-get-motif -id {params.motif_id} {input.motifs} \
            > {params.outdir}/top_motif.meme

        fimo \
            --thresh   {params.thresh} \
            --no-qvalue \
            --oc       {params.outdir} \
            {params.outdir}/top_motif.meme \
            {params.seqs}
        """

rule fimo:
    input: 
        motifs = outdir + '/motif/streme_bg_0/streme.xml'
    output:
        outdir + '/motif/fimo/fimo.tsv'
    params: 
        outdir = outdir + '/motif/fimo',
        seqs = config['annotations']['mtbc0_region_seqs']
    conda: '../envs/motif.yaml'
    shell:
        """
        fimo --oc {params.outdir} {input.motifs} {params.seqs}
        """
