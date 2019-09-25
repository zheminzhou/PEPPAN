# generate pan-genome prediction using PEPPA
python PEPPA.py -P examples/GCF_000010485.combined.gff.gz --min_cds 60 --incompleteCDS s -p examples/ST131 examples/*.gff.gz
# generate summaries for PEPPA results
python PEPPA_parse.py -g examples/ST131.PEPPA.gff -s examples/PEPPA_gffs -m -t -c 3
