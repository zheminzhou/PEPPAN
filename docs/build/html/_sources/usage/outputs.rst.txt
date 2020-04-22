Outputs
*******
## Outputs for PEPPA.py
There are two final outputs for PEPPA.py:

1. &lt;prefix&gt;.PEPPA.gff
   
   This file includes all pan-genes predicted by PEPPA in GFF3 format. Intact CDSs are assigned as "CDS", disrupted genes (potential pseudogenes) are assigned as "pseudogene" and suspicious annotations that are removed are described as "misc_feature" entries. 

* If any of the predicted CDSs and psueogenes overlaps with old gene predictions in the original GFF files, the old gene is described in an attribute named "old_locus_tag" of the entry. 

* Each gene and pseudogene is assigned into one of the orthologous groups. This orthologous group is described in "inference" field in a format of: 
~~~~~~~~
inference=ortholog_group:<source_genome>:<exemplar_gene>:<allele_ID>:<start & end coordinates of alignment in the exemplar gene>:<start & end coordinates of alignmenet in the genome>
~~~~~~~~

2. &lt;prefix&gt;.alleles.fna
   
   This file contains all the unique alleles of all pan genes predicted by PEPPA. 

## Outputs for PEPPA_parse.py
PEPPA_parse.py generates: 

1. &lt;prefix&gt;.gene_content.matrix or &lt;prefix&gt;.CDS_content.matrix

A matrix of gene presence/absence in all genomes.

2. &lt;prefix&gt;.gene_content.nwk or &lt;prefix&gt;.CDS_content.nwk

A tree built based on gene presence/absence in all genomes.

3. &lt;prefix&gt;.gene_content.curve or &lt;prefix&gt;.CDS_content.curve

The rare-fraction curves for the pan-genome and core-genome

4. &lt;prefix&gt;.gene_CGAV.tree or &lt;prefix&gt;.CDS_CGAV.tree

Core Genome Allelic Variation trees based on the sequence differences of the core genes. This is similar to  but should not be treated as a cgMLST scheme, because the genes included in the analysis depend on the genomes. The result of CGAV analysis is not comparable across different analyses. 

