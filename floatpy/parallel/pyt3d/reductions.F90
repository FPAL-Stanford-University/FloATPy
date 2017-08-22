module reductions
    use mpi
    use kind_parameters, only: rkind,mpirkind
    use exits, only: GracefulExit
    
    implicit none

    private
    public :: P_MAXVAL, P_MINVAL, P_SUM, P_OR

    interface P_MAXVAL
        module procedure P_MAXVAL_arr4, P_MAXVAL_arr3, P_MAXVAL_arr2, P_MAXVAL_sca
    end interface

    interface P_MINVAL
        module procedure P_MINVAL_arr4, P_MINVAL_arr3, P_MINVAL_arr2, P_MINVAL_sca
    end interface

    interface P_SUM
        module procedure P_SUM_arr3, P_SUM_arr2, P_SUM_arr1, P_SUM_sca, P_SUM_sca_INT
    end interface

contains
    
    function P_MAXVAL_arr4(x) result(maximum)
        real(rkind), dimension(:,:,:,:), intent(in) :: x
        real(rkind) :: maximum
        real(rkind) :: mymax
        integer :: ierr

        mymax = MAXVAL(x)
        call MPI_Allreduce(mymax, maximum, 1, mpirkind, MPI_MAX, MPI_COMM_WORLD, ierr)

    end function

    function P_MAXVAL_arr3(x) result(maximum)
        real(rkind), dimension(:,:,:), intent(in) :: x
        real(rkind) :: maximum
        real(rkind) :: mymax
        integer :: ierr

        mymax = MAXVAL(x)
        call MPI_Allreduce(mymax, maximum, 1, mpirkind, MPI_MAX, MPI_COMM_WORLD, ierr)

    end function

    function P_OR(x, comm) result(maximum)
        logical, intent(in) :: x
        integer, intent(in), optional :: comm
        logical :: maximum 
        integer :: ierr

        if (present(comm)) then
            call MPI_Allreduce(x,maximum, 1, MPI_LOGICAL,MPI_LOR, comm, ierr)
        else
            call MPI_Allreduce(x,maximum, 1, MPI_LOGICAL,MPI_LOR, MPI_COMM_WORLD, ierr)
        end if 
    end function 

    function P_MAXVAL_arr2(x) result(maximum)
        real(rkind), dimension(:,:), intent(in) :: x
        real(rkind) :: maximum
        real(rkind) :: mymax
        integer :: ierr

        mymax = MAXVAL(x)
        call MPI_Allreduce(mymax, maximum, 1, mpirkind, MPI_MAX, MPI_COMM_WORLD, ierr)

    end function

    function P_MAXVAL_sca(x) result(maximum)
        real(rkind), intent(in) :: x
        real(rkind) :: maximum
        integer :: ierr

        call MPI_Allreduce(x, maximum, 1, mpirkind, MPI_MAX, MPI_COMM_WORLD, ierr)

    end function
    
    function P_MINVAL_arr4(x) result(minimum)
        real(rkind), dimension(:,:,:,:), intent(in) :: x
        real(rkind) :: minimum
        real(rkind) :: mymin
        integer :: ierr

        mymin = MINVAL(x)
        call MPI_Allreduce(mymin, minimum, 1, mpirkind, MPI_MIN, MPI_COMM_WORLD, ierr)

    end function

    function P_MINVAL_arr3(x) result(minimum)
        real(rkind), dimension(:,:,:), intent(in) :: x
        real(rkind) :: minimum
        real(rkind) :: mymin
        integer :: ierr

        mymin = MINVAL(x)
        call MPI_Allreduce(mymin, minimum, 1, mpirkind, MPI_MIN, MPI_COMM_WORLD, ierr)

    end function

    function P_MINVAL_arr2(x) result(minimum)
        real(rkind), dimension(:,:), intent(in) :: x
        real(rkind) :: minimum
        real(rkind) :: mymin
        integer :: ierr

        mymin = MINVAL(x)
        call MPI_Allreduce(mymin, minimum, 1, mpirkind, MPI_MIN, MPI_COMM_WORLD, ierr)

    end function

    function P_MINVAL_sca(x) result(minimum)
        real(rkind), intent(in) :: x
        real(rkind) :: minimum
        integer :: ierr

        call MPI_Allreduce(x, minimum, 1, mpirkind, MPI_MIN, MPI_COMM_WORLD, ierr)

    end function

    function P_SUM_arr3(x) result(summation)
        real(rkind), dimension(:,:,:), intent(in) :: x
        real(rkind) :: summation
        real(rkind) :: mysum
        integer :: ierr

        mysum = SUM(x)
        call MPI_Allreduce(mysum, summation, 1, mpirkind, MPI_SUM, MPI_COMM_WORLD, ierr)

    end function

    function P_SUM_arr2(x) result(summation)
        real(rkind), dimension(:,:), intent(in) :: x
        real(rkind) :: summation
        real(rkind) :: mysum
        integer :: ierr

        mysum = SUM(x)
        call MPI_Allreduce(mysum, summation, 1, mpirkind, MPI_SUM, MPI_COMM_WORLD, ierr)

    end function

    function P_SUM_arr1(x) result(summation)
        real(rkind), dimension(:), intent(in) :: x
        real(rkind) :: summation
        real(rkind) :: mysum
        integer :: ierr

        mysum = SUM(x)
        call MPI_Allreduce(mysum, summation, 1, mpirkind, MPI_SUM, MPI_COMM_WORLD, ierr)

    end function

    function P_SUM_sca(x) result(summation)
        real(rkind), intent(in) :: x
        real(rkind) :: summation
        integer :: ierr

        call MPI_Allreduce(x, summation, 1, mpirkind, MPI_SUM, MPI_COMM_WORLD, ierr)

    end function

    function P_SUM_sca_INT(x) result(summation)
        integer, intent(in) :: x
        integer :: summation
        integer :: ierr

        call MPI_Allreduce(x, summation, 1, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, ierr)

    end function
    

end module
