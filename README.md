# Prerequsites

1. Install perl Expect.pm
2. clean.ex and hbplus are provided in the bin/ dir if you find they will not run you can get the source files from the EBI. Clean is part of the procheck suite (https://www.ebi.ac.uk/thornton-srv/software/PROCHECK/) and hbplus can be downloaded via (https://www.ebi.ac.uk/thornton-srv/software/HBPLUS/install.html) 
3. svm_classify from svm_light may need recompiled (in data/)

# Example use

1)
Check that the chains you want to analyse are present in the pdb file

`> cp ./examples/1IAR.pdb ./`

`> ./bin/checkchains.pl 1IAR.pdb AB ./data/charmm19_ha.aaa`

If it is happy your chains are present proceed to
the next step

2)
Run cleanexpect, pass the location of the resdefs.dat

`> ./bin/clean_expect.pl 1IAR.pdb ./data/resdefs.dat`

This will produce a file 1IAR.new

3)
Run hbplus over the cleaned pdb file

`> ./bin/hbplus 1IAR.new 1IAR.pdb`

This will produce a file 1IAR.hb2

4)
Finally you can run the predictor, hspred

`> ./bin/hs-pred_v0_1.pl 1IAR A B ./data ./`

Args are the pdb code, the chains for the 2 domains to be analyses a path to the data files and a final temp path for the outputs. The pdb and hb2 files must be in the temp dir provided.

You should find a range of new ouptput files

5)
If you want to visualise the predictions annotated you can run splitpdb. This will output a pdb file for each chain with the ATOM temperature set to indicate the hotspot predictions.

`> ./bin/split_pdb.pl 1IAR A B ./`

ARGS: pdb code, chain IDs and the working/temp directory
