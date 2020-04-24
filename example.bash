# generate pan-genome prediction using PEPPA
PEPPA -P examples/GCF_000010485.combined.gff.gz --min_cds 60 --incompleteCDS s -p examples/ST131 examples/*.gff.gz

# generate summaries for PEPPA predicted CDSs and pseudogenes
PEPPA_parser -g examples/ST131.PEPPA.gff -s examples/PEPPA_out -m -t -c -a 95 

