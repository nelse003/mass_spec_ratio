import os
import sys
sys.path.append("/dls/science/groups/i04-1/elliot-dev/parse_xchemdb")
from refinement.prepare_scripts import write_exhaustive_csh

exhaustive_multiple_sampling = "/dls/science/groups/i04-1/elliot-dev/" \
                               "Work/exhaustive_search/run_exhaustive_multiple_sampling.py"

ccp4_path = "/dls/science/groups/i04-1/software/pandda_0.2.12/ccp4/" \
            "ccp4-7.0/bin/ccp4.setup-sh"

parse_xchemdb_script_dir = "/dls/science/groups/i04-1/elliot-dev/parse_xchemdb"

#TODO Loop over output files in copy_atoms folder

for xtal in :

    write_exhaustive_csh(pdb=pdb,
                     mtz=mtz,
                     refinement_script_dir=refinement_script_dir,
                     out_dir=out_dir,
                     crystal=xtal,
                     script_dir=parse_xchemdb_script_dir,
                     exhaustive_multiple_sampling=exhaustive_multiple_sampling,
                     ccp4_path=ccp4_path)