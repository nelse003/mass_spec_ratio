import sys
import os
import pandas as pd

sys.path.append("/dls/science/groups/i04-1/elliot-dev/Work/exhaustive_search")
from exhaustive.utils.utils import u_iso_to_b_fac, get_minimum_fofc

out_dir = "/dls/science/groups/i04-1/elliot-dev/Work/" \
          "NUDT7A_mass_spec_refinements/copy_atoms/exhaustive/2019-05-29"

xtal_prefix = "NUDT7A-x"

xtal_folders = [f.name
                for f in os.scandir(out_dir)
                if f.is_dir()
                and xtal_prefix in f.name]

summary = {}

for xtal in xtal_folders:

    csv = os.path.join(out_dir, xtal, "exhaustive_search.csv")

    occ, u_iso, fo_fc = get_minimum_fofc(csv)
    summary[xtal] = (occ, u_iso_to_b_fac(u_iso), fo_fc)

df = pd.DataFrame.from_dict(summary,
                            orient='index',
                            columns=["occupancy",
                                     "b_factor",
                                     "fo_fc"])
df.index.name = "crystal"

out_csv = os.path.join(out_dir,"exhaustive_minima.csv")
df.to_csv(out_csv)
