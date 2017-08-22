! Module gaussianstuff defined in file gaussian.F90

subroutine f90wrap_init(this, n_, periodic_)
    use gaussianstuff, only: init, gaussian
    implicit none
    
    type gaussian_ptr_type
        type(gaussian), pointer :: p => NULL()
    end type gaussian_ptr_type
    type(gaussian_ptr_type) :: this_ptr
    integer, intent(out), dimension(2) :: this
    integer, intent(in) :: n_
    logical, intent(in) :: periodic_
    allocate(this_ptr%p)
    call init(this=this_ptr%p, n_=n_, periodic_=periodic_)
    this = transfer(this_ptr, this)
end subroutine f90wrap_init

subroutine f90wrap_destroy(this)
    use gaussianstuff, only: destroy, gaussian
    implicit none
    
    type gaussian_ptr_type
        type(gaussian), pointer :: p => NULL()
    end type gaussian_ptr_type
    type(gaussian_ptr_type) :: this_ptr
    integer, intent(in), dimension(2) :: this
    this_ptr = transfer(this, this_ptr)
    call destroy(this=this_ptr%p)
    deallocate(this_ptr%p)
end subroutine f90wrap_destroy

subroutine f90wrap_filter1(this, f, fil, nb, nc, bc1_, bcn_, n0, n1, n2, n3, n4, &
    n5)
    use gaussianstuff, only: filter1, gaussian
    implicit none
    
    type gaussian_ptr_type
        type(gaussian), pointer :: p => NULL()
    end type gaussian_ptr_type
    type(gaussian_ptr_type) :: this_ptr
    integer, intent(in), dimension(2) :: this
    real(8), intent(in), dimension(n0,n1,n2) :: f
    real(8), intent(inout), dimension(n3,n4,n5) :: fil
    integer, intent(in) :: nb
    integer, intent(in) :: nc
    integer, optional, intent(in) :: bc1_
    integer, optional, intent(in) :: bcn_
    integer :: n0
    !f2py intent(hide), depend(f) :: n0 = shape(f,0)
    integer :: n1
    !f2py intent(hide), depend(f) :: n1 = shape(f,1)
    integer :: n2
    !f2py intent(hide), depend(f) :: n2 = shape(f,2)
    integer :: n3
    !f2py intent(hide), depend(fil) :: n3 = shape(fil,0)
    integer :: n4
    !f2py intent(hide), depend(fil) :: n4 = shape(fil,1)
    integer :: n5
    !f2py intent(hide), depend(fil) :: n5 = shape(fil,2)
    this_ptr = transfer(this, this_ptr)
    call filter1(this=this_ptr%p, f=f, fil=fil, nb=nb, nc=nc, bc1_=bc1_, bcn_=bcn_)
end subroutine f90wrap_filter1

subroutine f90wrap_filter2(this, f, fil, na, nc, bc1_, bcn_, n0, n1, n2, n3, n4, &
    n5)
    use gaussianstuff, only: gaussian, filter2
    implicit none
    
    type gaussian_ptr_type
        type(gaussian), pointer :: p => NULL()
    end type gaussian_ptr_type
    type(gaussian_ptr_type) :: this_ptr
    integer, intent(in), dimension(2) :: this
    real(8), intent(in), dimension(n0,n1,n2) :: f
    real(8), intent(inout), dimension(n3,n4,n5) :: fil
    integer, intent(in) :: na
    integer, intent(in) :: nc
    integer, optional, intent(in) :: bc1_
    integer, optional, intent(in) :: bcn_
    integer :: n0
    !f2py intent(hide), depend(f) :: n0 = shape(f,0)
    integer :: n1
    !f2py intent(hide), depend(f) :: n1 = shape(f,1)
    integer :: n2
    !f2py intent(hide), depend(f) :: n2 = shape(f,2)
    integer :: n3
    !f2py intent(hide), depend(fil) :: n3 = shape(fil,0)
    integer :: n4
    !f2py intent(hide), depend(fil) :: n4 = shape(fil,1)
    integer :: n5
    !f2py intent(hide), depend(fil) :: n5 = shape(fil,2)
    this_ptr = transfer(this, this_ptr)
    call filter2(this=this_ptr%p, f=f, fil=fil, na=na, nc=nc, bc1_=bc1_, bcn_=bcn_)
end subroutine f90wrap_filter2

subroutine f90wrap_filter3(this, f, fil, na, nb, bc1_, bcn_, n0, n1, n2, n3, n4, &
    n5)
    use gaussianstuff, only: gaussian, filter3
    implicit none
    
    type gaussian_ptr_type
        type(gaussian), pointer :: p => NULL()
    end type gaussian_ptr_type
    type(gaussian_ptr_type) :: this_ptr
    integer, intent(in), dimension(2) :: this
    real(8), intent(in), dimension(n0,n1,n2) :: f
    real(8), intent(inout), dimension(n3,n4,n5) :: fil
    integer, intent(in) :: na
    integer, intent(in) :: nb
    integer, optional, intent(in) :: bc1_
    integer, optional, intent(in) :: bcn_
    integer :: n0
    !f2py intent(hide), depend(f) :: n0 = shape(f,0)
    integer :: n1
    !f2py intent(hide), depend(f) :: n1 = shape(f,1)
    integer :: n2
    !f2py intent(hide), depend(f) :: n2 = shape(f,2)
    integer :: n3
    !f2py intent(hide), depend(fil) :: n3 = shape(fil,0)
    integer :: n4
    !f2py intent(hide), depend(fil) :: n4 = shape(fil,1)
    integer :: n5
    !f2py intent(hide), depend(fil) :: n5 = shape(fil,2)
    this_ptr = transfer(this, this_ptr)
    call filter3(this=this_ptr%p, f=f, fil=fil, na=na, nb=nb, bc1_=bc1_, bcn_=bcn_)
end subroutine f90wrap_filter3

! End of module gaussianstuff defined in file gaussian.F90

