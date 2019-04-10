sys.path.append("/dls/science/groups/i04-1/elliot-dev/parse_xchemdb")
from utils.filesystem import parse_refinement_folder

refinement_dir = "/dls/science/groups/i04-1/elliot-dev/Work/NUDT7A_mass_spec_refinements/
refinement_csv = os.path.join(refinement_dir, "NUDT7A_cov_log_pdb_mtz.csv")

parse_refinement_folder(refinement_csv=refinement_csv,
                        refinement_dir=refinement_dir,
                        refinement_type="superposed")
