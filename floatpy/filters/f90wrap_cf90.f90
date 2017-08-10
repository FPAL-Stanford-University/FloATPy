! Module cf90stuff defined in file cf90.F90

subroutine f90wrap_init(this, n_, periodic_)
    use cf90stuff, only: cf90, init
    implicit none
    
    type cf90_ptr_type
        type(cf90), pointer :: p => NULL()
    end type cf90_ptr_type
    type(cf90_ptr_type) :: this_ptr
    integer, intent(out), dimension(2) :: this
    integer, intent(in) :: n_
    logical, intent(in) :: periodic_
    allocate(this_ptr%p)
    call init(this=this_ptr%p, n_=n_, periodic_=periodic_)
    this = transfer(this_ptr, this)
end subroutine f90wrap_init

subroutine f90wrap_destroy(this)
    use cf90stuff, only: cf90, destroy
    implicit none
    
    type cf90_ptr_type
        type(cf90), pointer :: p => NULL()
    end type cf90_ptr_type
    type(cf90_ptr_type) :: this_ptr
    integer, intent(in), dimension(2) :: this
    this_ptr = transfer(this, this_ptr)
    call destroy(this=this_ptr%p)
    deallocate(this_ptr%p)
end subroutine f90wrap_destroy

subroutine f90wrap_filter1(this, f, df, na, nb, bc1_, bcn_, n0, n1, n2, n3, n4, &
    n5)
    use cf90stuff, only: cf90, filter1
    implicit none
    
    type cf90_ptr_type
        type(cf90), pointer :: p => NULL()
    end type cf90_ptr_type
    type(cf90_ptr_type) :: this_ptr
    integer, intent(in), dimension(2) :: this
    real(8), intent(in), dimension(n0,n1,n2) :: f
    real(8), intent(inout), dimension(n3,n4,n5) :: df
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
    !f2py intent(hide), depend(df) :: n3 = shape(df,0)
    integer :: n4
    !f2py intent(hide), depend(df) :: n4 = shape(df,1)
    integer :: n5
    !f2py intent(hide), depend(df) :: n5 = shape(df,2)
    this_ptr = transfer(this, this_ptr)
    call filter1(this=this_ptr%p, f=f, df=df, na=na, nb=nb, bc1_=bc1_, bcn_=bcn_)
end subroutine f90wrap_filter1

subroutine f90wrap_filter2(this, f, df, na, nb, bc1_, bcn_, n0, n1, n2, n3, n4, &
    n5)
    use cf90stuff, only: cf90, filter2
    implicit none
    
    type cf90_ptr_type
        type(cf90), pointer :: p => NULL()
    end type cf90_ptr_type
    type(cf90_ptr_type) :: this_ptr
    integer, intent(in), dimension(2) :: this
    real(8), intent(in), dimension(n0,n1,n2) :: f
    real(8), intent(inout), dimension(n3,n4,n5) :: df
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
    !f2py intent(hide), depend(df) :: n3 = shape(df,0)
    integer :: n4
    !f2py intent(hide), depend(df) :: n4 = shape(df,1)
    integer :: n5
    !f2py intent(hide), depend(df) :: n5 = shape(df,2)
    this_ptr = transfer(this, this_ptr)
    call filter2(this=this_ptr%p, f=f, df=df, na=na, nb=nb, bc1_=bc1_, bcn_=bcn_)
end subroutine f90wrap_filter2

subroutine f90wrap_filter3(this, f, df, na, nb, bc1_, bcn_, n0, n1, n2, n3, n4, &
    n5)
    use cf90stuff, only: cf90, filter3
    implicit none
    
    type cf90_ptr_type
        type(cf90), pointer :: p => NULL()
    end type cf90_ptr_type
    type(cf90_ptr_type) :: this_ptr
    integer, intent(in), dimension(2) :: this
    real(8), intent(in), dimension(n0,n1,n2) :: f
    real(8), intent(inout), dimension(n3,n4,n5) :: df
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
    !f2py intent(hide), depend(df) :: n3 = shape(df,0)
    integer :: n4
    !f2py intent(hide), depend(df) :: n4 = shape(df,1)
    integer :: n5
    !f2py intent(hide), depend(df) :: n5 = shape(df,2)
    this_ptr = transfer(this, this_ptr)
    call filter3(this=this_ptr%p, f=f, df=df, na=na, nb=nb, bc1_=bc1_, bcn_=bcn_)
end subroutine f90wrap_filter3

! End of module cf90stuff defined in file cf90.F90

