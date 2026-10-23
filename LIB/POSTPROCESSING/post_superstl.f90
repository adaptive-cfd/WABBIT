subroutine post_superstl()
    use module_globals
    use module_mesh
    use module_params
    use module_mpi
    use module_operators
    use module_physics_metamodule
    use module_time_step
    use module_stl_file_reader
    use module_helpers
    use module_ini_files_parser_mpi
    use module_t_files

    implicit none

    ! The wabbit implementation was not very elegant and used brute force. The script in python is much more efficient, and you can modify
    ! it easier for other purposes as well. The code can be found in the git repository history, should it be useful.
    call abort(2210261, "The computation of superSTL files is now a preprocessing script in https://github.com/adaptive-cfd/python-tools and no longer a part of wabbit")
end subroutine
