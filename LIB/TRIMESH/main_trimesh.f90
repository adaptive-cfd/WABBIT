
program main_trimesh
    use module_trimesh
    implicit none
    type(surfaceTrimesh) :: Mesh, Mesh3D
    integer :: n_displaced

    print *, "TRIMESH module loaded successfully."

    call import2DGnuplotMesh("input.mesh", Mesh)
    call mesh_2Dto3D(Mesh, Mesh3D)
    call group_identical_verts(Mesh3D, 1e-12_rk)
    call displace_vertex(Mesh3D, 3, (/0.0_rk, 0.0_rk, 1.0_rk/), n_displaced)
    print *, n_displaced, " points displaced"


end program main_trimesh
