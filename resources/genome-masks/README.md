# Genome masks for builds GRCh37 and GRCh38

The get the GRCh37 mask:

```
wget http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/supporting/accessible_genome_masks/20140520.strict_mask.autosomes.bed
```

Get the GRCh38 mask:

```
wget http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000_genomes_project/working/20160622_genome_mask_GRCh38/StrictMask/20160622.allChr.mask.bed
```

### NOTE:

- Mask 20140520.strict_mask.autosomes.bed is for GRCh37 and it is in build GRCh37. This is only for the autosomes.
- Mask 20160622.allChr.bed is for GRCh38 and it is in build GRCh38. This is for all chromosomes which includes the autosomes.
