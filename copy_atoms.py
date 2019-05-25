import os

from giant.structure.utils import transfer_residue_groups_from_other
from iotbx.pdb import hierarchy


def copy_atoms(copy_params):
    """ Copy atoms from one pdb file to many, then refine.

    Copy dimple pdb, mtz and cif with cys bond
    Copy ligand atoms from existing coordinates 
    Run giant.merge_conformations to generate a multi state model
    Copy link records suitable for both conformers of the ligand
    Run quick refine to generate refined ligand 
    """,

    # generate output directory if it doesn't exist
    if not os.path.exists(copy_params.output.out_dir):
        os.mkdir(copy_params.output.out_dir)

    # read in PDB file from which atoms are to be taken from (ground structure)
    pdb_in = hierarchy.input(file_name=copy_params.input.base_pdb)
    sel_cache = pdb_in.hierarchy.atom_selection_cache()

    # produce a hierarchy with atoms to copied
    selection_string_list = []
    chains_new = set()
    for atom_new in copy_params.input.atoms_new:
        selection_string = "(resid {} and chain {})".format(atom_new[1], atom_new[0])
        selection_string_list.append(selection_string)
        chains_new.add(atom_new[0])
    selection_string = "or".join(selection_string_list)
    new_atoms_sel = sel_cache.selection(selection_string)
    new_atoms_hier = pdb_in.hierarchy.select(new_atoms_sel)

    # Produce a selection string to determine which atoms are removed
    selection_string_list = []
    if copy_params.input.atoms_remove is not None:
        for atom_remove in copy_params.input.atoms_remove:
            selection_string = "(resid {} and chain {})".format(
                atom_remove[1], atom_remove[0]
            )
            selection_string_list.append(selection_string)

        selection_string = "or".join(selection_string_list)
        not_selection_string = "not ({})".format(selection_string)

    # Define xtals to loop over
    xtals = copy_params.input.xtal_list
    for num in range(
        copy_params.input.start_xtal_number, copy_params.input.end_xtal_number + 1
    ):
        xtal_name = copy_params.input.prefix + "{0:0>4}".format(num)
        xtals.append(xtal_name)

    # Loop over all xtals
    for xtal_name in xtals:

        # For quick rerun
        if (
            os.path.exists(
                os.path.join(
                    copy_params.output.out_dir, xtal_name, copy_params.output.refine_pdb
                )
            )
            and not copy_params.settings.overwrite
        ):
            continue

        # Run only if sufficent input data
        if not os.path.exists(
            os.path.join(copy_params.input.path, xtal_name, copy_params.input.pdb_style)
        ):
            print(
                "pdb does not exist: {}".format(
                    os.path.join(
                        copy_params.input.path, xtal_name, copy_params.input.pdb_style
                    )
                )
            )
            continue

        pdb_in_refine = hierarchy.input(
            file_name=os.path.join(
                copy_params.input.path, xtal_name, copy_params.input.pdb_style
            )
        )

        acceptor_hierarchy = pdb_in_refine.construct_hierarchy()

        # remove atoms from xtal
        if copy_params.input.atoms_remove is not None:
            refine_sel_cache = pdb_in_refine.hierarchy.atom_selection_cache()
            remove_atoms_sel = refine_sel_cache.selection(not_selection_string)
            removed_hier = acceptor_hierarchy.select(remove_atoms_sel)
            working_hier = removed_hier
        else:
            working_hier = acceptor_hierarchy

        # Add atoms from base_pdb
        donor_hierarchy = new_atoms_hier
        acceptor_hier = transfer_residue_groups_from_other(
            working_hier, donor_hierarchy, in_place=False, verbose=False
        )

        # Generate output xtal directories
        if not os.path.exists(os.path.join(copy_params.output.out_dir, xtal_name)):
            os.mkdir(os.path.join(copy_params.output.out_dir, xtal_name))

        # Write output pdb with changed atoms
        f = open(
            os.path.join(copy_params.output.out_dir, xtal_name, copy_params.output.pdb),
            "w+",
        )
        f.write(
            acceptor_hier.as_pdb_string(
                crystal_symmetry=pdb_in_refine.input.crystal_symmetry()
            )
        )
        f.close()

        # Copy the input pdb to output directory
        os.chdir(os.path.join(copy_params.output.out_dir, xtal_name))
        os.system(
            "cp {} {}".format(
                os.path.join(
                    copy_params.input.path, xtal_name, copy_params.input.pdb_style
                ),
                os.path.join(
                    copy_params.output.out_dir, xtal_name, copy_params.input.pdb_style
                ),
            )
        )

        # Copy the input cif to output_directory
        os.system(
            "cp {} {}".format(
                copy_params.input.cif,
                os.path.join(
                    copy_params.output.out_dir,
                    xtal_name,
                    os.path.basename(copy_params.input.cif),
                ),
            )
        )

        # Copy the input mtz to output directory
        os.system(
            "cp -rL {} {}".format(
                os.path.join(
                    copy_params.input.path, xtal_name, copy_params.input.mtz_style
                ),
                os.path.join(
                    copy_params.output.out_dir, xtal_name, copy_params.input.mtz_style
                ),
            )
        )
        # Run giant.merge_conforamtions
        os.system(
            "giant.merge_conformations major={} minor={}".format(
                os.path.join(
                    copy_params.output.out_dir, xtal_name, copy_params.input.pdb_style
                ),
                os.path.join(
                    copy_params.output.out_dir, xtal_name, copy_params.output.pdb
                ),
            )
        )

        # Add link record strings into multimodel pdb file, prior to refinement
        if copy_params.input.link_record_list is not None:

            with open(
                os.path.join(
                    copy_params.output.out_dir,
                    xtal_name,
                    copy_params.output.multi_state_model_pdb,
                ),
                "r",
            ) as original:

                multi_model = original.read()

            with open(
                os.path.join(
                    copy_params.output.out_dir,
                    xtal_name,
                    copy_params.output.multi_state_model_pdb,
                ),
                "w",
            ) as modified:

                for link_record in copy_params.input.link_record_list:
                    modified.write(link_record)

                modified.write(multi_model)

        # Add extra params
        if copy_params.input.extra_params is not None:
            with open(
                "multi-state-restraints.{}.params".format(copy_params.settings.program),
                "a+",
            ) as param_file:
                if copy_params.input.extra_params not in param_file.read():
                    param_file.write(copy_params.input.extra_params)

        # Run giant.quick_refine
        cmds = "source {}\n".format(copy_params.settings.ccp4_path)

        if copy_params.settings.program == "phenix":
            cmds += "module load phenix\n"
        elif copy_params.settings.program == "buster":
            cmds += "module load buster\n"

        cmds += "giant.quick_refine {} {} {} params={} program={}\n".format(
            os.path.join(
                copy_params.output.out_dir,
                xtal_name,
                copy_params.output.multi_state_model_pdb,
            ),
            os.path.join(
                copy_params.output.out_dir, xtal_name, copy_params.input.mtz_style
            ),
            os.path.join(copy_params.output.out_dir, xtal_name, copy_params.input.cif),
            os.path.join(
                copy_params.output.out_dir, xtal_name, copy_params.settings.param_file
            ),
            copy_params.settings.program,
        )
        cmds += "giant.split_conformations refine.pdb"

        if copy_params.settings.qsub:
            f = open(
                os.path.join(
                    copy_params.output.out_dir,
                    xtal_name,
                    "{}_quick_refine.sh".format(xtal_name),
                ),
                "w",
            )

            f.write(cmds)
            f.close()

            os.system(
                "qsub {}".format(
                    os.path.join(
                        copy_params.output.out_dir,
                        xtal_name,
                        "{}_quick_refine.sh".format(xtal_name),
                    )
                )
            )
        else:
            os.system(cmds)
