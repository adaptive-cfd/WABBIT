
program main_trimesh
    use module_trimesh
    implicit none
    type(surfaceTrimesh) :: Mesh

    print *, "TRIMESH module loaded successfully."

    call importGnuplotMesh("input.mesh", Mesh)
    call mergeNeighbours(Mesh, 1e-12_rk)


end program main_trimesh
