#################################################################################
#################################################################################
#################################################################################

pairs.Rdata contains the gene pairs for use in the AuPairWise method.
There are 6 variables in the binary file.
The housekeeping interactions or stoichiometric pairs are contained in the "stoich.pairs" variable.
All pairs are labelled by their human entrez gene ID.

##  List of variables
- stoich.pairs
  - dimensions: 2669 2
  - These gene pairs were those "highly" co-expressed in RNA-seq and microarray aggregate co-expression networks.
  - These were filtered on genes that also were annotated to the GO protein complex term (GO:0043234).

- ppin.pairs
  - dimensions: 7594    2
  - These are protein pairs from BIOGRID (version BIOGRID-ALL-3.2.106, June 2014) where the pairs were from physical interaction studies, affinity-MS experiments, and also intersected with the GO protein complex annotated genes.

- yeastmouse.bottom.pairs
  - dimensions: 247   2
  - These gene pairs were those "poorly" co-expressed in mouse and yeast.
  - The human homologs of these genes were used to get the final list. Homologs were found using Homologene.

- yeastmouse.top.pairs
  - dimensions: 352   2
  - These gene pairs were those "highly" co-expressed in mouse and yeast.
  - The human homologs of these genes were used to get the final list.

- multiple.complexes.pairs
  - dimensions: 11218     3
  - First column is the protein complex identifier
  - Columns 2 and 3 are the gene IDs
  - These protein complexes were taken from mips and parsed for human complexes. Only the protein complexes that were annotated to GO (similar to the previous lists) were kept.
  - http://mips.helmholtz-muenchen.de/genre/proj/corum

- single.complexes.pairs
  - dimensions: 11218     3
  - First column is the protein complex identifier
  - Columns 2 and 3 are the gene IDs
  - As with the multiple.complexes.pairs, but only the complexes with single pairs of genes.