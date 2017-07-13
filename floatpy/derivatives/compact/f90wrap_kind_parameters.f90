! Module kind_parameters defined in file kind_parameters.F90

subroutine f90wrap_kind_parameters__get__rkind(f90wrap_rkind)
    use kind_parameters, only: kind_parameters_rkind => rkind
    implicit none
    integer, intent(out) :: f90wrap_rkind
    
    f90wrap_rkind = kind_parameters_rkind
end subroutine f90wrap_kind_parameters__get__rkind

subroutine f90wrap_kind_parameters__get__clen(f90wrap_clen)
    use kind_parameters, only: kind_parameters_clen => clen
    implicit none
    integer, intent(out) :: f90wrap_clen
    
    f90wrap_clen = kind_parameters_clen
end subroutine f90wrap_kind_parameters__get__clen

subroutine f90wrap_kind_parameters__get__stdin(f90wrap_stdin)
    use kind_parameters, only: kind_parameters_stdin => stdin
    implicit none
    integer, intent(out) :: f90wrap_stdin
    
    f90wrap_stdin = kind_parameters_stdin
end subroutine f90wrap_kind_parameters__get__stdin

subroutine f90wrap_kind_parameters__get__stdout(f90wrap_stdout)
    use kind_parameters, only: kind_parameters_stdout => stdout
    implicit none
    integer, intent(out) :: f90wrap_stdout
    
    f90wrap_stdout = kind_parameters_stdout
end subroutine f90wrap_kind_parameters__get__stdout

subroutine f90wrap_kind_parameters__get__stderr(f90wrap_stderr)
    use kind_parameters, only: kind_parameters_stderr => stderr
    implicit none
    integer, intent(out) :: f90wrap_stderr
    
    f90wrap_stderr = kind_parameters_stderr
end subroutine f90wrap_kind_parameters__get__stderr

! End of module kind_parameters defined in file kind_parameters.F90

