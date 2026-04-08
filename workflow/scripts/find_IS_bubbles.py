#%%
bubbles = '../results/benchmarking/minigraph/bubbles.bed'
is_seq = '../software/detettore6110/resources/is_targets/IS6110.fasta'

import subprocess
import tempfile


def blastn(query, subject):
    cmd = ['blastn','-query', query,'-subject', subject,
           '-outfmt', '6 qseqid qstart qend qlen sseqid sstart send sstrand pident length']
    
    proc = subprocess.Popen(cmd, stdout=subprocess.PIPE)
    output = proc.stdout.read()
    
    hits = []
    for line in output.splitlines():
        line = line.decode('utf-8')
        fields = line.split('\t')
        hits.append(fields)

    hits = pandas.DataFrame(hits, columns=['qseqid', 'qstart', 'qend', 'qlen', 'sseqid', 'sstart', 'send', 'sstrand', 'pident', 'length'])
    for var in ['qstart', 'qend', 'qlen', 'sstart', 'send', 'length']:
        hits[var] = hits[var].astype(int)
    hits['pident'] = hits['pident'].astype(float)
    
    # Convert to 0-based index


with open(bubbles) as f:
    
    # Temporary fasta file
    tmp_fasta = tempfile.mkdtemp()
    fasta_handle = open(f'{tmp_fasta}/query.fa', 'w')
    
    for line in f:
        fields = line.strip().split('\t')
        
        len_shortest = int(fields[6]) 
        len_longest = int(fields[7]) 
        
        seq_shortest = fields[12]
        seq_longest = fields[13]
        
        if len_shortest >=1000:
            bl_res = blastn()
        
        
        

 
# %%
