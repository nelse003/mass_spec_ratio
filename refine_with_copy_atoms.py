import sys
from copy_atoms import copy_atoms
from copy_phil import copy_phil

"""
Refine crystals from visits:

mx19301-17 x2074 - x2121
mx19301-18 x2130 - x2154
mx19301-20 x2160 - x2232
mx19301-26 x2157 - x2319
lb19758-25 x2297 - x2367

using the dimple files generated by running these crystals through XChemExplorer.
Copy atoms from known structure. 

Uses copy_atom function from exhaustive search
"""

# pdb = (
#     "/dls/science/groups/i04-1/elliot-dev/Work/exhaustive_search_data/"
#     "covalent_ratios_refine/NUDT7A-x1907/refine.pdb"
# )
pdb = "/dls/science/groups/i04-1/elliot-dev/Work/NUDT7A_mass_spec_refinements/dimple_with_lig.pdb"
cif = (
    "/dls/science/groups/i04-1/elliot-dev/Work/exhaustive_search_data/"
    "covalent_ratios_refine/NUDT7A-x1907/NUDT7A-x1907_LIG_CYS.cif"
)

copy_params = copy_phil.extract()

copy_params.input.path = (
    "/dls/labxchem/data/2017/lb18145-49/processing/analysis/initial_model"
)
copy_params.input.prefix = "NUDT7A-x"
copy_params.input.base_pdb = pdb
copy_params.input.atoms_new = [["E", "1"], ["A", "73"]]
copy_params.input.atoms_remove = [
    ["A", "73"],
    ["B", "208"],
    ["B", "151"],
    ["B", "60"],
    ["B", "189"],
    ["B", "33"],
    ["B", "40"],
    ["B", "11"],
    ["B", "196"],
]
copy_params.input.cif = cif
copy_params.input.link_record_list = [
    "LINKR        C  FLIG E   1                 SG  FCYS A  73                LIG-CYS\n",
    "LINKR        C  ELIG E   1                 SG  ECYS A  73                LIG-CYS\n"
    "LINK         SG FCYS A  73                 C  FLIG E   1     1555   1555  1.93\n",
]
copy_params.settings.ccp4_path = (
    "/dls/science/groups/i04-1/elliot-dev/ccp4/ccp4-7.0/bin/ccp4.setup-sh"
)

copy_params.output.out_dir = "/dls/science/groups/i04-1/elliot-dev/Work/NUDT7A_mass_spec_refinements/copy_atoms_190525_refmac"
copy_params.input.extra_params = "NCYC 5"
copy_params.settings.program = "refmac"
copy_params.settings.qsub = True

# Run multiple times for different ranges of crystals

copy_params.input.start_xtal_number = 2160
copy_params.input.end_xtal_number = 2160

copy_atoms(copy_params)

# copy_params.input.start_xtal_number = 2074
# copy_params.input.end_xtal_number = 2121
#
# copy_atoms(copy_params)
#
# copy_params.input.start_xtal_number = 2130
# copy_params.input.end_xtal_number = 2154
#
# copy_atoms(copy_params)
#
# copy_params.input.start_xtal_number = 2157
# copy_params.input.end_xtal_number = 2367
#
# copy_atoms(copy_params)
#
# copy_atoms(copy_params)

# Phenix
copy_params.output.out_dir = "/dls/science/groups/i04-1/elliot-dev/Work/NUDT7A_mass_spec_refinements/copy_atoms_190523_phenix"
copy_params.settings.program = "phenix"
copy_params.settings.param_file = "multi-state-restraints.phenix.params"
copy_params.input.extra_params = "refinement.main.number_of_macro_cycles=20"

# copy_params.input.start_xtal_number = 2074
# copy_params.input.end_xtal_number = 2121
#
# copy_atoms(copy_params)
#
# copy_params.input.start_xtal_number = 2130
# copy_params.input.end_xtal_number = 2154
#
# copy_atoms(copy_params)
#
# copy_params.input.start_xtal_number = 2157
# copy_params.input.end_xtal_number = 2367
#
# copy_atoms(copy_params)

# Buster
copy_params.output.out_dir = "/dls/science/groups/i04-1/elliot-dev/Work/NUDT7A_mass_spec_refinements/copy_atoms_190524_buster"
copy_params.settings.program = "buster"
copy_params.settings.param_file = "params.gelly"
copy_params.input.extra_params = ""

# copy_params.input.start_xtal_number = 2074
# copy_params.input.end_xtal_number = 2121
#
# copy_atoms(copy_params)
#
# copy_params.input.start_xtal_number = 2130
# copy_params.input.end_xtal_number = 2154
#
# copy_atoms(copy_params)
#
# copy_params.input.start_xtal_number = 2157
# copy_params.input.end_xtal_number = 2367
#
# copy_atoms(copy_params)
