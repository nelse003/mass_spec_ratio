from exhaustive.utils.utils import get_minimum_fofc

exh_b_fix = "/dls/science/groups/i04-1/elliot-dev/Work/NUDT7A_mass_spec_refinements/copy_atoms/exhaustive_b_fix/2019-06-17/"

folders = [os.path.join(exh_b_fix, d) for d
           in os.path.listdir(exh_b_fix)
           if os.path.isdir(os.path.join(exh_b_fix,d))]

for folder in exh_b_fix:

    exh_csv = os.path.join(exh_b_fix, folder, "exhaustive_search.csv")

    if os.path.exists(exh_b_fix, folder) == exh_csv:

        print(get_minimum_fofc(csv_name = exh_csv))