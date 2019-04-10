import os
import sys
import luigi

sys.path.append("/dls/science/groups/i04-1/elliot-dev/parse_xchemdb")
from utils.filesystem import parse_refinement_folder
from tasks.plotting import PlotBoundOccHistogram
from path_config import Path

refinement_dir = "/dls/science/groups/i04-1/elliot-dev/Work/NUDT7A_mass_spec_refinements/"
refinement_csv = os.path.join(refinement_dir, "NUDT7A_cov_log_pdb_mtz.csv")

parse_refinement_folder(refinement_csv=refinement_csv,
                        refinement_dir=refinement_dir,
                        refinement_type="superposed")

out_dir = "/dls/science/groups/i04-1/elliot-dev/Work/NUDT7A_mass_spec_refinements"

paths = Path()
paths.log_pdb_mtz = os.path.join(out_dir, 'NUDT7A_cov_log_pdb_mtz.csv')
paths.superposed = os.path.join(out_dir, 'superposed.csv')
paths.refine = os.path.join(out_dir, 'refine.csv')
paths.log_occ_csv = os.path.join(out_dir, 'log_occ.csv')
paths.log_occ_resname = os.path.join(out_dir, 'log_occ_resname.csv')
paths.occ_state_comment_csv = os.path.join(out_dir, 'occ_state_comment.csv')
paths.refinement_summary = os.path.join(out_dir, 'refinement_summary.csv')
paths.occ_correct_csv = os.path.join(out_dir, 'occ_correct.csv')
#paths.refinement_summary_plot = os.path.join(out_dir, 'refinement_summary.png')
#paths.convergence_histogram = os.path.join(out_dir, 'convergence_hist.png')
#paths.occ_conv_scatter = os.path.join(out_dir, 'occ_conv_scatter.png')
#paths.ground_occ_histogram = os.path.join(out_dir, 'occ_ground_histogram.png')
paths.bound_occ_histogram = os.path.join(out_dir, 'occ_bound_histogram.png')

luigi.build([
    PlotBoundOccHistogram(occ_state_comment_csv=paths.occ_state_comment_csv,
           log_occ_resname=paths.log_occ_resname,
           log_occ_csv=paths.log_occ_csv,
           log_pdb_mtz_csv=paths.log_pdb_mtz,
           occ_correct_csv=paths.occ_correct_csv,
           plot_path=paths.bound_occ_histogram)
    ],
    local_scheduler=False, workers=10)
