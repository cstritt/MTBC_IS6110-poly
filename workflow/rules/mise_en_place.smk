configfile: 'config/config.yaml'

# Map assembly names to assembly paths
assemblies = {}
with open(config['benchmarking']['assemblies']) as f:
    next(f)
    for line in f:
        fields = line.strip().split("\t")
        assemblies[fields[0]] = fields[1]

# Create dictionary of G numbers with path to cram file as entries
gnrs = {}
with open(config['gnumbers']) as f:
    for line in f:
        gnr = line.strip().split('\t')[0]
        a = gnr[0:3]
        b = gnr[3:5]
        c = gnr[5:]
        gnrs[gnr] = f'{config['path_to_db']}/{a}/{b}/{c}/{gnr}.cram'


# Output directory
import os
outdir = config['outdir']['main']
if not os.path.exists(outdir):
    os.makedirs(outdir)

outdir_bm = config['outdir']['benchmarking']
if not os.path.exists(outdir_bm):
    os.makedirs(outdir_bm)

# Create output folder structure for G numbers
detettore_gnr_subdirs= []
for gnr in gnrs.keys():
    detettore_gnr_subdirs.append(os.path.join(gnr[:3], gnr[3:5], gnr[5:]))

# Create lineage dictionary
lineages = {}

with open(config["metadata"]) as f:
    header = next(f).strip().split("\t")
    gnumberi = header.index("GNUMBER")
    lineagei = header.index("LINEAGE_x")

    for line in f:
        fields = line.strip().split("\t")
        gnumber = fields[gnumberi]
        lineage = fields[lineagei]
        if lineage not in lineages:
            lineages[lineage] = []
        lineages[lineage].append(gnumber)
