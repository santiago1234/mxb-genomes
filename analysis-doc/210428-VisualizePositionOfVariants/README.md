
# Overview


What is the overlap of regions with no variant sites between the 1TGP (high coverage) data
and the 50 MXB genomes?


# Protocol


I computed the variant density in windows across the genome. See the script [variant_density.py](variant_density.py).


The variant density in the 1TGP data is higher, the reason for this is because
there are more samples. To make the visual comparison easier I scaled the
data in the range 0 to 1.



# Results

- Points in orange are regions with low variant density.
- The top plot shows the variant density before I ran the filter to remove the masked sites.
- The bottom plot shows the variant density after I ran the filter to remove the masked sites.

![image](plots/variant_density.png)


# Conclusion

After removing the masked sites, there is a nice overlap
in the regions with low density between the two data sets.
