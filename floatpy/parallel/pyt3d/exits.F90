module exits

    use kind_parameters, only: stdout,stderr
    
    implicit none
    private
    public :: GracefulExit
        
contains

    subroutine GracefulExit(message, errcode)
        use mpi
        use kind_parameters, only: stderr, stdout
        integer, intent(in) :: errCode
        character(len=*), intent(in) :: message
        integer :: rank, ierr

        call mpi_comm_rank(mpi_comm_world, rank, ierr)
        if (rank == 0) then
            write(stderr,'(A)') '========== ERROR =========='
            write(stderr,'(A)') message
            write(stderr,'(A)') '==========================='
        end if 
        call mpi_abort(mpi_comm_world, errCode, ierr)
        if (ierr /= 0) then
            print*, "SHIT! It won't abort!"
        end if 

    end subroutine
        
end module 
