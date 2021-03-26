# Downlod resources an accesory data to run the pipeline
# every output of the rules defined here goes into: mxb-genomes/resources/

# Genetic maps
# these maps are used in phasing, local ancestry, etc.

rule download_genetic_maps:
    # download genetic maps from github repo:
    # https://github.com/odelaneau/shapeit4
    output:
        expand("resources/genetic-maps/chr{chrn}.b38.gmap", chrn=CHROMS)
    shell:
        """
        wget https://github.com/odelaneau/shapeit4/raw/master/maps/genetic_maps.b38.tar.gz
        tar -xzf genetic_maps.b38.tar.gz
        gunzip chr*.gz
        mkdir -p resources/genetic-maps/
        mv *gmap resources/genetic-maps/
        rm genetic_maps.b38.tar.gz
        """

rule oneT_genome_population_info:
    #Download the population information
    #for the 1TGP data
    #NOTE: this population info table does not
    #include the additional samples what were sequenced.
    #This sample not included are related to the samples
    #here and we won't use them in the project
    output:
        "resources/1TGP-samples-meta-data/integrated_call_samples_v3.20130502.ALL.panel"
    shell:
        """
        wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/integrated_call_samples_v3.20130502.ALL.panel
        mv integrated_call_samples_v3.20130502.ALL.panel {output}
        """
