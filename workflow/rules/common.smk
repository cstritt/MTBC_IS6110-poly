configfile: "config/config.yaml"

# Map assembly names to fastq paths
assemblies = {}

with open(config["assemblies"]) as f:
    next(f)
    for line in f:
        fields = line.strip().split("\t")
        assemblies[fields[0]] = fields[1]

# Create lineage dictionary
lineages = {}

with open(config["metadata"]) as f:
    header = next(f).strip().split("\t")
    gnumberi = header.index("GNUMBER")
    lineagei = header.index("LINEAGE")

    for line in f:
        fields = line.strip().split("\t")
        gnumber = fields[gnumberi]
        lineage = fields[lineagei]
        if lineage not in lineages:
            lineages[lineage] = []
        lineages[lineage].append(gnumber)
