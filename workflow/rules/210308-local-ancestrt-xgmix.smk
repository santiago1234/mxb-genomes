rule download_XGMix:
    # this rules dowloads the software (git repo)
    # to run XGMix
    output:
        directory("workflow/scripts/XGMix-master")
    shell:
        """
        wget https://github.com/AI-sandbox/XGMix/archive/master.zip
        unzip master.zip
        rm master.zip
        mv XGMix-master workflow/scripts/
        """


#rule make_sample_map_file:
    # sample map file matching reference samples
    # to their respective reference populations
#    pass



