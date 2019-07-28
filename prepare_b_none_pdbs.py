import os
import sys
import argparse
from iotbx.pdb import hierarchy


def copy_b(pdb, ref_pdb, out_pdb, chain, resid, altloc=""):

    """
    Copy b factor of a single residue to another pdb file

    Parameters
    ----------
    pdb: str
        path to pdb file

    ref_pdb: str
        path to reference pdb file

    chain: str
        chain of interest

    resid: str
        residue of interest

    altloc: str
        altloc of interest

    Returns
    -------


    """

    # read into iotbx.hierarchy
    ref_pdb_in = hierarchy.input(file_name=ref_pdb)
    # read into iotbx.selection cache
    sel_cache = ref_pdb_in.hierarchy.atom_selection_cache()

    # Get selection object which corresponds to supplied chain residue id and altloc
    if altloc == "":
        ref_sel = sel_cache.selection("chain {} resid {}".format(chain, resid))
    else:
        ref_sel = sel_cache.selection(
            "chain {} resid {} altloc {}".format(chain, resid, altloc)
        )
    # Select that residue from main hierarchy
    ref_hier = ref_pdb_in.hierarchy.select(ref_sel)
    ref_lig = {}
    for ref_chain in ref_hier.only_model().chains():
        for residue_group in ref_chain.residue_groups():
            for atom_group in residue_group.atom_groups():
                # Get B factor and occ information on residue by looking a individual atoms
                for atom in atom_group.atoms():
                    ref_lig[atom.name] = atom.b

    # read into iotbx.hierarchy
    pdb_in = hierarchy.input(file_name=pdb)

    # read into iotbx.selection cache
    sel_cache = pdb_in.hierarchy.atom_selection_cache()

    # Get selection object which corresponds to supplied chain residue id and altloc
    if altloc == "":
        sel = sel_cache.selection("chain {} resid {}".format(chain, resid))

    else:
        sel = sel_cache.selection(
            "chain {} resid {} altloc {}".format(chain, resid, altloc)
        )

    hier = pdb_in.hierarchy.select(sel)
    # Select that residue from main hierarchy
    for current_chain in hier.only_model().chains():
        for residue_group in current_chain.residue_groups():
            for atom_group in residue_group.atom_groups():
                for atom in atom_group.atoms():

                    atom.b = ref_lig[atom.name]

    if not os.path.isdir(os.path.dirname(out_pdb)):
        os.makedirs(os.path.dirname(out_pdb))

    with open(out_pdb, "w") as out:
        out.write(
            pdb_in.hierarchy.as_pdb_string(
                crystal_symmetry=pdb_in.input.crystal_symmetry()
            )
        )


if __name__ == "__main__":

    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--ref_pdb",
        action="store",
        help="reference pdb to copy B factors of selected residue from",
        dest="ref_pdb",
        type=str,
    )

    parser.add_argument(
        "--pdb", action="store", help="pdb to copy hierarchy from", dest="pdb", type=str
    )

    parser.add_argument(
        "--out_pdb",
        action="store",
        help="pdb to copy B factors to",
        dest="out_pdb",
        type=str,
    )

    parser.add_argument(
        "--chain",
        action="store",
        help="chain to select from",
        dest="chain",
        default="E",
        type=str,
    )

    parser.add_argument(
        "--resid",
        action="store",
        help="residue id to select",
        dest="resid",
        default="1",
        type=str,
    )

    parser.add_argument(
        "--altloc",
        action="store",
        help="altloc id to select",
        dest="altloc",
        default="",
        type=str,
    )

    args = parser.parse_args()

    copy_b(
        ref_pdb=args.ref_pdb,
        pdb=args.pdb,
        out_pdb=args.out_pdb,
        chain=args.chain,
        resid=args.resid,
        altloc=args.altloc,
    )

    # copy_b(ref_pdb="/dls/science/groups/i04-1/elliot-dev/Work/"\
    #                    "NUDT7A_mass_spec_refinements/copy_atoms_190525_buster/NUDT7A-x2160/refine.split.bound-state.pdb",
    #        pdb="/dls/science/groups/i04-1/elliot-dev/Work/"\
    #                    "NUDT7A_mass_spec_refinements/copy_atoms_190525_buster/NUDT7A-x2161/refine.split.bound-state.pdb",
    #        out_pdb="/dls/science/groups/i04-1/elliot-dev/Work/"\
    #                    "NUDT7A_mass_spec_refinements/copy_atoms/buster_b_none/NUDT7A-x2161/b_none.pdb",
    #         chain='E',
    #         resid='1',
    #         altloc="")
