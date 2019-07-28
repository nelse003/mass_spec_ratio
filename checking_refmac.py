import os

refmac_folder = "/dls/science/groups/i04-1/elliot-dev/Work/NUDT7A_mass_spec_refinements/copy_atoms/refmac/2019-05-29/"
refmac_superposed_folder = "/dls/science/groups/i04-1/elliot-dev/Work/NUDT7A_mass_spec_refinements/copy_atoms_190525_refmac/"
new_refmac_folder = "/dls/science/groups/i04-1/elliot-dev/Work/NUDT7A_mass_spec_refinements/copy_atoms/refmac/2019-07-08"

folder_interest = new_refmac_folder

for folder in os.listdir(folder_interest):
    if os.path.isdir(os.path.join(folder_interest, folder)):

        refine_pdb = os.path.join(folder_interest, folder, "refine.pdb")

        if os.path.isfile(refine_pdb):

            with open(refine_pdb) as f:
                contents = f.read()

        if "ALIG" in contents:
            print(refine_pdb)
