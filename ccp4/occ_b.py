import pandas as pd
import numpy as np
from iotbx.pdb import hierarchy
import argparse

def get_occ_b(pdb, chain, resid, altloc=''):
    """
    Get occupancy and b factor of a single residue

    Parameters
    ----------
    pdb: str
        path to pdb file

    chain: str
        chain of interest

    resid: str
        residue of interest

    altloc: str
        altloc of interest

    Returns
    -------
    mean_occ: float
        mean occupancy of residue

    mean_b: float
        mean b factor of residue

    std_b: float
        standard deviation of b factor of refisude

    """

    # read into iotbx.hierarchy
    pdb_in = hierarchy.input(file_name=pdb)
    # read into iotbx.selection cache
    sel_cache = pdb_in.hierarchy.atom_selection_cache()

    # Get selection object which corresponds to supplied chain residue id and altloc
    if altloc == '':
        sel = sel_cache.selection(
            "chain {} resid {}".format(chain, resid)
        )

    else:
        sel = sel_cache.selection(
            "chain {} resid {} altloc {}".format(chain, resid, altloc)
        )

    # Select that residue from main hierarchy
    hier = pdb_in.hierarchy.select(sel)
    resnames = []
    for chain in hier.only_model().chains():
        for residue_group in chain.residue_groups():
            for atom_group in residue_group.atom_groups():
                resnames.append(atom_group.resname)

                # Get B factor and occ information on residue by looking a individual atoms
                b = []
                occ = []
                for atom in atom_group.atoms():
                    b.append(atom.b)
                    occ.append(atom.occ)

                mean_occ = np.mean(occ)
                mean_b = np.mean(b)
                std_b = np.std(b)

                return mean_occ, mean_b, std_b

if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument("pdb", help="path to pdb file", type=str)
    parser.add_argument("out", help="output file name", type=str)
    parser.add_argument("chain", help="chain of pdb file to get occ and b from", type=str)
    parser.add_argument("resid", help="resid of pdb file to get occ and b from", type=str)
    parser.add_argument("--altloc", help="resid of pdb file to get occ and b from", type=str, default='')

    args = parser.parse_args()

    with open(args.out,"w") as f:

        mean_occ, mean_b, std_b = get_occ_b(pdb=args.pdb, chain=args.chain, resid=args.resid, altloc=args.altloc)

        f.write("{},{},{}".format(mean_occ, mean_b, std_b))