! Module t3dmod defined in file t3dMod.F90

subroutine f90wrap_init(this, comm3d, nx, ny, nz, px, py, pz, periodic_, &
    reorder, fail, nghosts, createcrosscommunicators)
    use t3dmod, only: init, t3d
    implicit none
    
    type t3d_ptr_type
        type(t3d), pointer :: p => NULL()
    end type t3d_ptr_type
    type(t3d_ptr_type) :: this_ptr
    integer, intent(out), dimension(2) :: this
    integer, intent(in) :: comm3d
    integer, intent(in) :: nx
    integer, intent(in) :: ny
    integer, intent(in) :: nz
    integer, intent(in) :: px
    integer, intent(in) :: py
    integer, intent(in) :: pz
    logical, dimension(3), intent(in) :: periodic_
    logical, intent(in) :: reorder
    logical, intent(inout) :: fail
    integer, dimension(3), optional, intent(in) :: nghosts
    logical, optional, intent(in) :: createcrosscommunicators
    allocate(this_ptr%p)
    call init(this=this_ptr%p, comm3D=comm3d, nx=nx, ny=ny, nz=nz, px=px, py=py, &
        pz=pz, periodic_=periodic_, reorder=reorder, fail=fail, nghosts=nghosts, &
        createCrossCommunicators=createcrosscommunicators)
    this = transfer(this_ptr, this)
end subroutine f90wrap_init

subroutine f90wrap_destroy(this)
    use t3dmod, only: destroy, t3d
    implicit none
    
    type t3d_ptr_type
        type(t3d), pointer :: p => NULL()
    end type t3d_ptr_type
    type(t3d_ptr_type) :: this_ptr
    integer, intent(in), dimension(2) :: this
    this_ptr = transfer(this, this_ptr)
    call destroy(this=this_ptr%p)
    deallocate(this_ptr%p)
end subroutine f90wrap_destroy

subroutine f90wrap_transpose_3d_to_x(this, input, output, n0, n1, n2, n3, n4, &
    n5)
    use t3dmod, only: t3d, transpose_3d_to_x
    implicit none
    
    type t3d_ptr_type
        type(t3d), pointer :: p => NULL()
    end type t3d_ptr_type
    type(t3d_ptr_type) :: this_ptr
    integer, intent(in), dimension(2) :: this
    real(8), intent(in), dimension(n0,n1,n2) :: input
    real(8), intent(inout), dimension(n3,n4,n5) :: output
    integer :: n0
    !f2py intent(hide), depend(input) :: n0 = shape(input,0)
    integer :: n1
    !f2py intent(hide), depend(input) :: n1 = shape(input,1)
    integer :: n2
    !f2py intent(hide), depend(input) :: n2 = shape(input,2)
    integer :: n3
    !f2py intent(hide), depend(output) :: n3 = shape(output,0)
    integer :: n4
    !f2py intent(hide), depend(output) :: n4 = shape(output,1)
    integer :: n5
    !f2py intent(hide), depend(output) :: n5 = shape(output,2)
    this_ptr = transfer(this, this_ptr)
    call transpose_3d_to_x(this=this_ptr%p, input=input, output=output)
end subroutine f90wrap_transpose_3d_to_x

subroutine f90wrap_transpose_x_to_3d(this, input, output, n0, n1, n2, n3, n4, &
    n5)
    use t3dmod, only: t3d, transpose_x_to_3d
    implicit none
    
    type t3d_ptr_type
        type(t3d), pointer :: p => NULL()
    end type t3d_ptr_type
    type(t3d_ptr_type) :: this_ptr
    integer, intent(in), dimension(2) :: this
    real(8), intent(in), dimension(n0,n1,n2) :: input
    real(8), intent(inout), dimension(n3,n4,n5) :: output
    integer :: n0
    !f2py intent(hide), depend(input) :: n0 = shape(input,0)
    integer :: n1
    !f2py intent(hide), depend(input) :: n1 = shape(input,1)
    integer :: n2
    !f2py intent(hide), depend(input) :: n2 = shape(input,2)
    integer :: n3
    !f2py intent(hide), depend(output) :: n3 = shape(output,0)
    integer :: n4
    !f2py intent(hide), depend(output) :: n4 = shape(output,1)
    integer :: n5
    !f2py intent(hide), depend(output) :: n5 = shape(output,2)
    this_ptr = transfer(this, this_ptr)
    call transpose_x_to_3d(this=this_ptr%p, input=input, output=output)
end subroutine f90wrap_transpose_x_to_3d

subroutine f90wrap_transpose_3d_to_y(this, input, output, n0, n1, n2, n3, n4, &
    n5)
    use t3dmod, only: transpose_3d_to_y, t3d
    implicit none
    
    type t3d_ptr_type
        type(t3d), pointer :: p => NULL()
    end type t3d_ptr_type
    type(t3d_ptr_type) :: this_ptr
    integer, intent(in), dimension(2) :: this
    real(8), intent(in), dimension(n0,n1,n2) :: input
    real(8), intent(inout), dimension(n3,n4,n5) :: output
    integer :: n0
    !f2py intent(hide), depend(input) :: n0 = shape(input,0)
    integer :: n1
    !f2py intent(hide), depend(input) :: n1 = shape(input,1)
    integer :: n2
    !f2py intent(hide), depend(input) :: n2 = shape(input,2)
    integer :: n3
    !f2py intent(hide), depend(output) :: n3 = shape(output,0)
    integer :: n4
    !f2py intent(hide), depend(output) :: n4 = shape(output,1)
    integer :: n5
    !f2py intent(hide), depend(output) :: n5 = shape(output,2)
    this_ptr = transfer(this, this_ptr)
    call transpose_3d_to_y(this=this_ptr%p, input=input, output=output)
end subroutine f90wrap_transpose_3d_to_y

subroutine f90wrap_transpose_y_to_3d(this, input, output, n0, n1, n2, n3, n4, &
    n5)
    use t3dmod, only: transpose_y_to_3d, t3d
    implicit none
    
    type t3d_ptr_type
        type(t3d), pointer :: p => NULL()
    end type t3d_ptr_type
    type(t3d_ptr_type) :: this_ptr
    integer, intent(in), dimension(2) :: this
    real(8), intent(in), dimension(n0,n1,n2) :: input
    real(8), intent(inout), dimension(n3,n4,n5) :: output
    integer :: n0
    !f2py intent(hide), depend(input) :: n0 = shape(input,0)
    integer :: n1
    !f2py intent(hide), depend(input) :: n1 = shape(input,1)
    integer :: n2
    !f2py intent(hide), depend(input) :: n2 = shape(input,2)
    integer :: n3
    !f2py intent(hide), depend(output) :: n3 = shape(output,0)
    integer :: n4
    !f2py intent(hide), depend(output) :: n4 = shape(output,1)
    integer :: n5
    !f2py intent(hide), depend(output) :: n5 = shape(output,2)
    this_ptr = transfer(this, this_ptr)
    call transpose_y_to_3d(this=this_ptr%p, input=input, output=output)
end subroutine f90wrap_transpose_y_to_3d

subroutine f90wrap_transpose_3d_to_z(this, input, output, n0, n1, n2, n3, n4, &
    n5)
    use t3dmod, only: transpose_3d_to_z, t3d
    implicit none
    
    type t3d_ptr_type
        type(t3d), pointer :: p => NULL()
    end type t3d_ptr_type
    type(t3d_ptr_type) :: this_ptr
    integer, intent(in), dimension(2) :: this
    real(8), intent(in), dimension(n0,n1,n2) :: input
    real(8), intent(inout), dimension(n3,n4,n5) :: output
    integer :: n0
    !f2py intent(hide), depend(input) :: n0 = shape(input,0)
    integer :: n1
    !f2py intent(hide), depend(input) :: n1 = shape(input,1)
    integer :: n2
    !f2py intent(hide), depend(input) :: n2 = shape(input,2)
    integer :: n3
    !f2py intent(hide), depend(output) :: n3 = shape(output,0)
    integer :: n4
    !f2py intent(hide), depend(output) :: n4 = shape(output,1)
    integer :: n5
    !f2py intent(hide), depend(output) :: n5 = shape(output,2)
    this_ptr = transfer(this, this_ptr)
    call transpose_3d_to_z(this=this_ptr%p, input=input, output=output)
end subroutine f90wrap_transpose_3d_to_z

subroutine f90wrap_transpose_z_to_3d(this, input, output, n0, n1, n2, n3, n4, &
    n5)
    use t3dmod, only: t3d, transpose_z_to_3d
    implicit none
    
    type t3d_ptr_type
        type(t3d), pointer :: p => NULL()
    end type t3d_ptr_type
    type(t3d_ptr_type) :: this_ptr
    integer, intent(in), dimension(2) :: this
    real(8), intent(in), dimension(n0,n1,n2) :: input
    real(8), intent(inout), dimension(n3,n4,n5) :: output
    integer :: n0
    !f2py intent(hide), depend(input) :: n0 = shape(input,0)
    integer :: n1
    !f2py intent(hide), depend(input) :: n1 = shape(input,1)
    integer :: n2
    !f2py intent(hide), depend(input) :: n2 = shape(input,2)
    integer :: n3
    !f2py intent(hide), depend(output) :: n3 = shape(output,0)
    integer :: n4
    !f2py intent(hide), depend(output) :: n4 = shape(output,1)
    integer :: n5
    !f2py intent(hide), depend(output) :: n5 = shape(output,2)
    this_ptr = transfer(this, this_ptr)
    call transpose_z_to_3d(this=this_ptr%p, input=input, output=output)
end subroutine f90wrap_transpose_z_to_3d

subroutine f90wrap_fill_halo_x(this, array, n0, n1, n2)
    use t3dmod, only: fill_halo_x, t3d
    implicit none
    
    type t3d_ptr_type
        type(t3d), pointer :: p => NULL()
    end type t3d_ptr_type
    type(t3d_ptr_type) :: this_ptr
    integer, intent(in), dimension(2) :: this
    real(8), intent(inout), dimension(n0,n1,n2) :: array
    integer :: n0
    !f2py intent(hide), depend(array) :: n0 = shape(array,0)
    integer :: n1
    !f2py intent(hide), depend(array) :: n1 = shape(array,1)
    integer :: n2
    !f2py intent(hide), depend(array) :: n2 = shape(array,2)
    this_ptr = transfer(this, this_ptr)
    call fill_halo_x(this=this_ptr%p, array=array)
end subroutine f90wrap_fill_halo_x

subroutine f90wrap_fill_halo_y(this, array, n0, n1, n2)
    use t3dmod, only: t3d, fill_halo_y
    implicit none
    
    type t3d_ptr_type
        type(t3d), pointer :: p => NULL()
    end type t3d_ptr_type
    type(t3d_ptr_type) :: this_ptr
    integer, intent(in), dimension(2) :: this
    real(8), intent(inout), dimension(n0,n1,n2) :: array
    integer :: n0
    !f2py intent(hide), depend(array) :: n0 = shape(array,0)
    integer :: n1
    !f2py intent(hide), depend(array) :: n1 = shape(array,1)
    integer :: n2
    !f2py intent(hide), depend(array) :: n2 = shape(array,2)
    this_ptr = transfer(this, this_ptr)
    call fill_halo_y(this=this_ptr%p, array=array)
end subroutine f90wrap_fill_halo_y

subroutine f90wrap_fill_halo_z(this, array, n0, n1, n2)
    use t3dmod, only: fill_halo_z, t3d
    implicit none
    
    type t3d_ptr_type
        type(t3d), pointer :: p => NULL()
    end type t3d_ptr_type
    type(t3d_ptr_type) :: this_ptr
    integer, intent(in), dimension(2) :: this
    real(8), intent(inout), dimension(n0,n1,n2) :: array
    integer :: n0
    !f2py intent(hide), depend(array) :: n0 = shape(array,0)
    integer :: n1
    !f2py intent(hide), depend(array) :: n1 = shape(array,1)
    integer :: n2
    !f2py intent(hide), depend(array) :: n2 = shape(array,2)
    this_ptr = transfer(this, this_ptr)
    call fill_halo_z(this=this_ptr%p, array=array)
end subroutine f90wrap_fill_halo_z

subroutine f90wrap_get_sz3d(this, sz3d)
    use t3dmod, only: get_sz3d, t3d
    implicit none
    
    type t3d_ptr_type
        type(t3d), pointer :: p => NULL()
    end type t3d_ptr_type
    type(t3d_ptr_type) :: this_ptr
    integer, intent(in), dimension(2) :: this
    integer, dimension(3), intent(inout) :: sz3d
    this_ptr = transfer(this, this_ptr)
    call get_sz3d(this=this_ptr%p, sz3D=sz3d)
end subroutine f90wrap_get_sz3d

subroutine f90wrap_get_st3d(this, st3d)
    use t3dmod, only: t3d, get_st3d
    implicit none
    
    type t3d_ptr_type
        type(t3d), pointer :: p => NULL()
    end type t3d_ptr_type
    type(t3d_ptr_type) :: this_ptr
    integer, intent(in), dimension(2) :: this
    integer, dimension(3), intent(inout) :: st3d
    this_ptr = transfer(this, this_ptr)
    call get_st3d(this=this_ptr%p, st3D=st3d)
end subroutine f90wrap_get_st3d

subroutine f90wrap_get_en3d(this, en3d)
    use t3dmod, only: t3d, get_en3d
    implicit none
    
    type t3d_ptr_type
        type(t3d), pointer :: p => NULL()
    end type t3d_ptr_type
    type(t3d_ptr_type) :: this_ptr
    integer, intent(in), dimension(2) :: this
    integer, dimension(3), intent(inout) :: en3d
    this_ptr = transfer(this, this_ptr)
    call get_en3d(this=this_ptr%p, en3D=en3d)
end subroutine f90wrap_get_en3d

subroutine f90wrap_get_sz3dg(this, sz3dg)
    use t3dmod, only: t3d, get_sz3dg
    implicit none
    
    type t3d_ptr_type
        type(t3d), pointer :: p => NULL()
    end type t3d_ptr_type
    type(t3d_ptr_type) :: this_ptr
    integer, intent(in), dimension(2) :: this
    integer, dimension(3), intent(inout) :: sz3dg
    this_ptr = transfer(this, this_ptr)
    call get_sz3dg(this=this_ptr%p, sz3Dg=sz3dg)
end subroutine f90wrap_get_sz3dg

subroutine f90wrap_get_st3dg(this, st3dg)
    use t3dmod, only: t3d, get_st3dg
    implicit none
    
    type t3d_ptr_type
        type(t3d), pointer :: p => NULL()
    end type t3d_ptr_type
    type(t3d_ptr_type) :: this_ptr
    integer, intent(in), dimension(2) :: this
    integer, dimension(3), intent(inout) :: st3dg
    this_ptr = transfer(this, this_ptr)
    call get_st3dg(this=this_ptr%p, st3Dg=st3dg)
end subroutine f90wrap_get_st3dg

subroutine f90wrap_get_en3dg(this, en3dg)
    use t3dmod, only: t3d, get_en3dg
    implicit none
    
    type t3d_ptr_type
        type(t3d), pointer :: p => NULL()
    end type t3d_ptr_type
    type(t3d_ptr_type) :: this_ptr
    integer, intent(in), dimension(2) :: this
    integer, dimension(3), intent(inout) :: en3dg
    this_ptr = transfer(this, this_ptr)
    call get_en3dg(this=this_ptr%p, en3Dg=en3dg)
end subroutine f90wrap_get_en3dg

subroutine f90wrap_get_szx(this, szx)
    use t3dmod, only: t3d, get_szx
    implicit none
    
    type t3d_ptr_type
        type(t3d), pointer :: p => NULL()
    end type t3d_ptr_type
    type(t3d_ptr_type) :: this_ptr
    integer, intent(in), dimension(2) :: this
    integer, dimension(3), intent(inout) :: szx
    this_ptr = transfer(this, this_ptr)
    call get_szx(this=this_ptr%p, szX=szx)
end subroutine f90wrap_get_szx

subroutine f90wrap_get_szy(this, szy)
    use t3dmod, only: get_szy, t3d
    implicit none
    
    type t3d_ptr_type
        type(t3d), pointer :: p => NULL()
    end type t3d_ptr_type
    type(t3d_ptr_type) :: this_ptr
    integer, intent(in), dimension(2) :: this
    integer, dimension(3), intent(inout) :: szy
    this_ptr = transfer(this, this_ptr)
    call get_szy(this=this_ptr%p, szY=szy)
end subroutine f90wrap_get_szy

subroutine f90wrap_get_szz(this, szz)
    use t3dmod, only: get_szz, t3d
    implicit none
    
    type t3d_ptr_type
        type(t3d), pointer :: p => NULL()
    end type t3d_ptr_type
    type(t3d_ptr_type) :: this_ptr
    integer, intent(in), dimension(2) :: this
    integer, dimension(3), intent(inout) :: szz
    this_ptr = transfer(this, this_ptr)
    call get_szz(this=this_ptr%p, szZ=szz)
end subroutine f90wrap_get_szz

subroutine f90wrap_optimize_decomposition(comm3d, nx, ny, nz, periodic, &
    ret_this, nghosts)
    use t3dmod, only: t3d, optimize_decomposition
    implicit none
    
    type t3d_ptr_type
        type(t3d), pointer :: p => NULL()
    end type t3d_ptr_type
    integer, intent(in) :: comm3d
    integer, intent(in) :: nx
    integer, intent(in) :: ny
    integer, intent(in) :: nz
    logical, dimension(3), intent(in) :: periodic
    type(t3d_ptr_type) :: ret_this_ptr
    integer, intent(out), dimension(2) :: ret_this
    integer, dimension(3), optional, intent(in) :: nghosts
    allocate(ret_this_ptr%p)
    ret_this_ptr%p = optimize_decomposition(comm3D=comm3d, nx=nx, ny=ny, nz=nz, &
        periodic=periodic, nghosts=nghosts)
    ret_this = transfer(ret_this_ptr, ret_this)
end subroutine f90wrap_optimize_decomposition

! End of module t3dmod defined in file t3dMod.F90

