import os
import sys
import argparse
import subprocess
import datetime

sys.path.append("/dls/science/groups/i04-1/elliot-dev/parse_xchemdb")
from refinement.prepare_scripts import write_refmac_csh
from refinement.prepare_scripts import write_exhaustive_csh
from refinement.prepare_scripts import write_phenix_csh
from refinement.prepare_scripts import write_buster_csh

if __name__ == "__main__":

    """
    Refine crystals from visits:
    
    mx19301-17 x2074 - x2121
    mx19301-18 x2130 - x2154
    mx19301-20 x2160 - x2232
    mx19301-26 x2157 - x2319
    lb19758-25 x2297 - x2367
    
    using the same pdb and cif file, and param file from same folder. 
    Note that the compound numbering has been superseeded in scarab
    
    Notes
    -----
    
    Write qsub jobs. y
    
    Run Refinement
    
    Checking that these are empty:
    
    'NUDT7A-x2074', y
    'NUDT7A-x2075', y
    'NUDT7A-x2077', y
    'NUDT7A-x2078', y
    'NUDT7A-x2080', y
    'NUDT7A-x2086', y
    'NUDT7A-x2099', y
    'NUDT7A-x2107', y
    'NUDT7A-x2114', y
    'NUDT7A-x2115', y
    'NUDT7A-x2117', y
    'NUDT7A-x2118', y
    'NUDT7A-x2119'  y
    
    Quite a few due to dropped puck (10 missing)
    NUDT7A-x2146 fails to refine as low-res therefore incorrect spacegroup
    
    """

    # Add ability to parse command line arguments
    parser = argparse.ArgumentParser(description="Plot NUDT7 mass spectroscopy ratios")
    # Add path argument
    parser.add_argument("--program", action="store", dest="program", default="refmac")
    args = parser.parse_args()

    in_dir = "/dls/labxchem/data/2017/lb18145-49/processing/analysis/initial_model"

    out_dir = os.path.join(
        "/dls/science/groups/i04-1/elliot-dev/Work/NUDT7A_mass_spec_refinements/copy_atoms",
        args.program,
        str(datetime.date.today())
    )

    refinement_script_dir = (
        "/dls/science/groups/i04-1/elliot-dev/Work/NUDT7A_mass_spec_refinements/scripts"
    )

    cif = (
        "/dls/science/groups/i04-1/elliot-dev/Work/exhaustive_search_data/"
        "covalent_ratios_refine/NUDT7A-x1907/NUDT7A-x1907_LIG_CYS.cif"
    )
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

    os.chdir(out_dir)

    if args.program == "refmac":
        input_folder = "/dls/science/groups/i04-1/elliot-dev/Work/NUDT7A_mass_spec_refinements/copy_atoms_190525_refmac"
        params_fname = "multi-state-restraints.refmac.params"

    elif args.program == "buster":
        input_folder = "/dls/science/groups/i04-1/elliot-dev/Work/NUDT7A_mass_spec_refinements/copy_atoms_190525_buster"
        params_fname = "params.gelly"

    elif args.program == "phenix":
        input_folder = "/dls/science/groups/i04-1/elliot-dev/Work/NUDT7A_mass_spec_refinements/copy_atoms_190525_phenix"
        params_fname = "multi-state-restraints.phenix.params"

    elif args.program == "exhaustive":
        input_folder = "/dls/science/groups/i04-1/elliot-dev/Work/NUDT7A_mass_spec_refinements/copy_atoms_190525_buster"
        params_fname = ""

    prefix = "NUDT7A-x"

    xtals = []
    for num in range(2074, 2121 + 1):
        xtal_name = prefix + "{0:0>4}".format(num)
        xtals.append(xtal_name)

    for num in range(2130, 2154 + 1):
        xtal_name = prefix + "{0:0>4}".format(num)
        xtals.append(xtal_name)

    for num in range(2157, 2367 + 1):
        xtal_name = prefix + "{0:0>4}".format(num)
        xtals.append(xtal_name)

    mtz_to_check = []
    for xtal in xtals:

        mtz = os.path.join(in_dir, xtal, "dimple.mtz")
        params = os.path.join(input_folder, xtal, params_fname )
        refine_pdb = os.path.join(input_folder, xtal, "refine.pdb")
        pdb = os.path.join(input_folder, xtal, "refine.split.bound-state.pdb" )
        pdb_adj = os.path.join(input_folder, xtal, "refine.split-bound-state_new_link_record.pdb" )

        if not os.path.isfile(pdb):
            continue

        with open(pdb, 'r') as pdb_file:
            lines = pdb_file.readlines()
        with open(pdb_adj,'w') as pdb_adj_file:
            for line in lines:
                if line.strip("\n").startswith("LINK"):
                    pdb_adj_file.write("LINK         SG ACYS A  73                 C   LIG E   1     1555   1555  1.80\n")
                    pdb_adj_file.write("LINK         SG BCYS A  73                 C   LIG E   1     1555   1555  1.77\n")
                if not line.strip("\n").startswith("LINK"):
                    pdb_adj_file.write(line)


        if not os.path.isfile(mtz):
            continue

        if not args.program == "exhaustive":
            if not os.path.isfile(params):
                continue

        if args.program == "refmac":
            write_refmac_csh(
                pdb=pdb_adj,
                crystal=xtal,
                cif=cif,
                mtz=mtz,
                out_dir=os.path.join(out_dir,xtal),
                refinement_script_dir=refinement_script_dir,
                script_dir="/dls/science/groups/i04-1/elliot-dev/parse_xchemdb",
                ncyc=50,
                ccp4_path="/dls/science/groups/i04-1/elliot-dev/ccp4/ccp4-7.0/bin/ccp4.setup-sh",
            )
            csh_file = os.path.join(
                refinement_script_dir, "{}_{}.csh".format(xtal, "bound")
            )

        elif args.program == "phenix":
            write_phenix_csh(
                pdb=pdb_adj,
                mtz=mtz,
                cif=cif,
                script_dir="/dls/science/groups/i04-1/elliot-dev/parse_xchemdb",
                refinement_script_dir=refinement_script_dir,
                out_dir=os.path.join(out_dir,xtal),
                crystal=xtal,
                ncyc=20
            )
            csh_file = os.path.join(
                refinement_script_dir, "{}_{}.csh".format(xtal, "phenix")
            )

        elif args.program == "buster":
            write_buster_csh(
                pdb=pdb_adj,
                mtz=mtz,
                cif=cif,
                out_dir=os.path.join(out_dir,xtal),
                script_dir="/dls/science/groups/i04-1/elliot-dev/parse_xchemdb",
                refinement_script_dir=refinement_script_dir,
                crystal=xtal,
            )
            csh_file = os.path.join(
                refinement_script_dir, "{}_{}.csh".format(xtal, "buster")
            )

        elif args.program == "exhaustive":
             write_exhaustive_csh(
                 pdb=refine_pdb,
                 mtz=mtz,
                 script_dir="/dls/science/groups/i04-1/elliot-dev/parse_xchemdb",
                 refinement_script_dir=refinement_script_dir,
                 out_dir=os.path.join(out_dir,xtal),
                 crystal=xtal,
                 exhaustive_multiple_sampling="/dls/science/groups/i04-1/elliot-dev/Work/exhaustive_search/exhaustive/exhaustive_multiple_sampling.py",
                 ccp4_path="/dls/science/groups/i04-1/elliot-dev/ccp4/ccp4-7.0/bin/ccp4.setup-sh",
             )
             csh_file = os.path.join(
                 refinement_script_dir, "{}_{}.csh".format(xtal, "exhaustive")
             )

        os.system("qsub {}".format(csh_file))




