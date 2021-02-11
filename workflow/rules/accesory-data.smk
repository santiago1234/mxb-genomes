# Download accessory data

rule download_genetic_maps:
    # download genetic maps from github repo:
    # https://github.com/joepickrell/1000-genomes-genetic-maps
    output:
        "resources/genetic-maps/chr22.b37.gmap.txt"
    shell:
        """
        wget https://github.com/odelaneau/shapeit4/raw/master/maps/genetic_maps.b37.tar.gz
        tar -xzf genetic_maps.b37.tar.gz
        gunzip chr*.gz
        for f in *.gmap; do
          echo $f
          cat $f | awk '{{print $2 "\t" $1 "\t" $3}}' > $f.txt
        done
        mkdir -p resources/genetic-maps/
        mv *gmap.txt resources/genetic-maps/
        rm genetic_maps.b37.tar.gz *gmap
        """
