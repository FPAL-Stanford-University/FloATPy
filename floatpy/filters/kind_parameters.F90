! Set the floating point precision

module kind_parameters

    implicit none
    
    private
    public :: rkind, clen, stdin, stdout, stderr

    integer, parameter :: rkind=kind(0.d0)
   
    integer, parameter :: clen = 120

    integer, parameter :: stdin  = 5
    integer, parameter :: stdout = 6
    integer, parameter :: stderr = 0

end module
