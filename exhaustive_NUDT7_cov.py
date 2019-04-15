import os
import sys
sys.path.append("/dls/science/groups/i04-1/elliot-dev/parse_xchemdb")
from refinement.prepare_scripts import write_exhaustive_csh

exhaustive_multiple_sampling = "/dls/science/groups/i04-1/elliot-dev/" \
                               "Work/exhaustive_search/run_exhaustive_multiple_sampling.py"

ccp4_path = "/dls/science/groups/i04-1/software/pandda_0.2.12/ccp4/" \
            "ccp4-7.0/bin/ccp4.setup-sh"

parse_xchemdb_script_dir = "/dls/science/groups/i04-1/elliot-dev/parse_xchemdb"

copy_atoms_dir = "/dls/science/groups/i04-1/elliot-dev/Work/" \
                 "NUDT7A_mass_spec_refinements/copy_atoms"

out_script_dir = "/dls/science/groups/i04-1/elliot-dev/Work/" \
                 "NUDT7A_mass_spec_refinements/copy_atoms/exhaustive_scripts"

out_dir = "/dls/science/groups/i04-1/elliot-dev/Work/" \
                 "NUDT7A_mass_spec_refinements/copy_atoms/exhaustive"

xtal_prefix = "NUDT7A-x"

if not os.path.isdir(out_script_dir):
    os.makedirs(out_script_dir)

if not os.path.isdir(out_dir):
    os.makedirs(out_dir)

xtal_folders = [f.name
                for f in os.scandir(copy_atoms_dir)
                if f.is_dir()
                and xtal_prefix in f.name]

for xtal in xtal_folders:

    pdb = os.path.join(copy_atoms_dir, xtal, "refine.pdb")
    mtz = os.path.join(copy_atoms_dir, xtal, "refine.mtz")

    write_exhaustive_csh(pdb=pdb,
                     mtz=mtz,
                     refinement_script_dir=out_script_dir,
                     out_dir=out_dir,
                     crystal=xtal,
                     script_dir=parse_xchemdb_script_dir,
                     exhaustive_multiple_sampling=exhaustive_multiple_sampling,
                     ccp4_path=ccp4_path)