import os
import sys

sys.path.append("/dls/science/groups/i04-1/elliot-dev/parse_xchemdb")
from refinement.prepare_scripts import write_quick_refine_csh


if __name__ == "__main__":

    """
    Refine crystals from visits:
    
    mx19301-17 x2074 - x2121
    mx19301-18 x2130 -x2154
    mx19301-20 x2160 - x2232
    
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

    out_dir = "/dls/science/groups/i04-1/elliot-dev/Work/NUDT7A_mass_spec_refinements/copy_atoms"

    refinement_script_dir = (
        "/dls/science/groups/i04-1/elliot-dev/Work/NUDT7A_mass_spec_refinements/scripts"
    )

    pdb = (
        "/dls/science/groups/i04-1/elliot-dev/Work/exhaustive_search_data/"
        "covalent_ratios_refine/NUDT7A-x1907/refine.pdb"
    )

    cif = (
        "/dls/science/groups/i04-1/elliot-dev/Work/exhaustive_search_data/"
        "covalent_ratios_refine/NUDT7A-x1907/NUDT7A-x1907_LIG_CYS.cif"
    )

    if args.program == "refmac":

        params = (
            "/dls/science/groups/i04-1/elliot-dev/Work/exhaustive_search_data/"
            "covalent_ratios_refine/NUDT7A-x1907/refine_0004/input.params"
        )

    elif args.program == "buster":
        subprocess.call(
            [
                "/dls/science/groups/i04-1/elliot-dev/ccp4/ccp4-7.0/bin/giant.make_restraints",
                pdb,
            ]
        )

    elif args.program == "phenix":
        pass

    exit()

    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

    prefix = "NUDT7A-x"

    xtals = []
    for num in range(2074, 2121 + 1):
        xtal_name = prefix + "{0:0>4}".format(num)
        xtals.append(xtal_name)

    for num in range(2130, 2154 + 1):
        xtal_name = prefix + "{0:0>4}".format(num)
        xtals.append(xtal_name)

    for num in range(2160, 2232 + 1):
        xtal_name = prefix + "{0:0>4}".format(num)
        xtals.append(xtal_name)

    mtz_to_check = []
    for xtal in xtals:

        mtz = os.path.join(in_dir, xtal, "{}.mtz".format(xtal))

        if not os.path.exists(mtz):
            print("mtz does not exist: {}".format(mtz))
            mtz_to_check.append(xtal)
            continue

        crystal_dir = os.path.join(out_dir, xtal)
        if not os.path.exists(crystal_dir):
            os.makedirs(crystal_dir)

        write_quick_refine_csh(
            crystal=xtal,
            refine_pdb=pdb,
            cif=cif,
            out_dir=crystal_dir,
            refinement_params=params,
            refinement_script_dir=refinement_script_dir,
            free_mtz=mtz,
        )

    print(mtz_to_check)
