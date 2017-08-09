! Module cd10stuff defined in file cd10.F90

subroutine f90wrap_init(this, n_, dx_, periodic_, bc1_, bcn_)
    use cd10stuff, only: init, cd10
    implicit none
    
    type cd10_ptr_type
        type(cd10), pointer :: p => NULL()
    end type cd10_ptr_type
    type(cd10_ptr_type) :: this_ptr
    integer, intent(out), dimension(2) :: this
    integer, intent(in) :: n_
    real(8), intent(in) :: dx_
    logical, intent(in) :: periodic_
    integer, intent(in) :: bc1_
    integer, intent(in) :: bcn_
    allocate(this_ptr%p)
    call init(this=this_ptr%p, n_=n_, dx_=dx_, periodic_=periodic_, bc1_=bc1_, &
        bcn_=bcn_)
    this = transfer(this_ptr, this)
end subroutine f90wrap_init

subroutine f90wrap_destroy(this)
    use cd10stuff, only: destroy, cd10
    implicit none
    
    type cd10_ptr_type
        type(cd10), pointer :: p => NULL()
    end type cd10_ptr_type
    type(cd10_ptr_type) :: this_ptr
    integer, intent(in), dimension(2) :: this
    this_ptr = transfer(this, this_ptr)
    call destroy(this=this_ptr%p)
    deallocate(this_ptr%p)
end subroutine f90wrap_destroy

subroutine f90wrap_dd1(this, f, df, na, nb, bc1_, bcn_, n0, n1, n2, n3, n4, n5)
    use cd10stuff, only: dd1, cd10
    implicit none
    
    type cd10_ptr_type
        type(cd10), pointer :: p => NULL()
    end type cd10_ptr_type
    type(cd10_ptr_type) :: this_ptr
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
    call dd1(this=this_ptr%p, f=f, df=df, na=na, nb=nb, bc1_=bc1_, bcn_=bcn_)
end subroutine f90wrap_dd1

subroutine f90wrap_dd2(this, f, df, na, nb, bc1_, bcn_, n0, n1, n2, n3, n4, n5)
    use cd10stuff, only: dd2, cd10
    implicit none
    
    type cd10_ptr_type
        type(cd10), pointer :: p => NULL()
    end type cd10_ptr_type
    type(cd10_ptr_type) :: this_ptr
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
    call dd2(this=this_ptr%p, f=f, df=df, na=na, nb=nb, bc1_=bc1_, bcn_=bcn_)
end subroutine f90wrap_dd2

subroutine f90wrap_dd3(this, f, df, na, nb, bc1_, bcn_, n0, n1, n2, n3, n4, n5)
    use cd10stuff, only: dd3, cd10
    implicit none
    
    type cd10_ptr_type
        type(cd10), pointer :: p => NULL()
    end type cd10_ptr_type
    type(cd10_ptr_type) :: this_ptr
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
    call dd3(this=this_ptr%p, f=f, df=df, na=na, nb=nb, bc1_=bc1_, bcn_=bcn_)
end subroutine f90wrap_dd3

subroutine f90wrap_d2d1(this, f, df, na, nb, bc1_, bcn_, n0, n1, n2, n3, n4, n5)
    use cd10stuff, only: cd10, d2d1
    implicit none
    
    type cd10_ptr_type
        type(cd10), pointer :: p => NULL()
    end type cd10_ptr_type
    type(cd10_ptr_type) :: this_ptr
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
    call d2d1(this=this_ptr%p, f=f, df=df, na=na, nb=nb, bc1_=bc1_, bcn_=bcn_)
end subroutine f90wrap_d2d1

subroutine f90wrap_d2d2(this, f, df, na, nb, bc1_, bcn_, n0, n1, n2, n3, n4, n5)
    use cd10stuff, only: d2d2, cd10
    implicit none
    
    type cd10_ptr_type
        type(cd10), pointer :: p => NULL()
    end type cd10_ptr_type
    type(cd10_ptr_type) :: this_ptr
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
    call d2d2(this=this_ptr%p, f=f, df=df, na=na, nb=nb, bc1_=bc1_, bcn_=bcn_)
end subroutine f90wrap_d2d2

subroutine f90wrap_d2d3(this, f, df, na, nb, bc1_, bcn_, n0, n1, n2, n3, n4, n5)
    use cd10stuff, only: d2d3, cd10
    implicit none
    
    type cd10_ptr_type
        type(cd10), pointer :: p => NULL()
    end type cd10_ptr_type
    type(cd10_ptr_type) :: this_ptr
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
    call d2d3(this=this_ptr%p, f=f, df=df, na=na, nb=nb, bc1_=bc1_, bcn_=bcn_)
end subroutine f90wrap_d2d3

! End of module cd10stuff defined in file cd10.F90

