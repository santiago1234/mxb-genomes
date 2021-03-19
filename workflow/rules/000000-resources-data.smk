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
