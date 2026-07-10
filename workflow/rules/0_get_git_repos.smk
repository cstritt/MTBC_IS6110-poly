rule clone_detettore:
    output: 'software/detettore6110/detettore6110.py'
    shell: 
        """
        git clone https://github.com/cstritt/detettore6110.git software/detettore6110
        touch {output}
        """

rule clone_lva:
    output: 
        'software/large_variable_alignment/get_alignment.py'
    shell: 
        """
        git clone https://github.com/cstritt/large_variable_alignment.git software/large_variable_alignment
        touch {output}
        
        """

rule clone_IS_simulation:
    output:
        'software/IS_simulation/abc_sampler.py'
    shell:
        """
        git clone https://github.com/cstritt/IS_simulation.git software/IS_simulation
        """

rule get_iqtree:
    output:
        'software/iqtree-3.0.1-Linux/bin/iqtree3'
    shell:
        """
        wget https://github.com/iqtree/iqtree3/releases/download/v3.0.1/iqtree-3.0.1-Linux.tar.gz software/

        tar -xzf software/iqtree-3.0.1-Linux.tar.gz
        rm software/iqtree-3.0.1-Linux.tar.gz

        """