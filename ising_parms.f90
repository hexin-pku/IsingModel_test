module ising_parms
implicit none
private
    public ::  betamin, betamax, Nsize, Nstep, Nsplit, samp_type, parmJ, parmh
    real(8) :: betamin, betamax, parmJ, parmh
    integer :: Nsize, Nsplit, Nstep, samp_type
    public :: read_parms
contains
subroutine read_parms(parmsfile)
    character(*) , intent(in) :: parmsfile
    logical :: my_exist
    
    inquire(file=parmsfile,exist=my_exist)
    if(.not. my_exist) stop 'parmsfile loss'
    open(unit=22, file=parmsfile, status='old')
    read(22,*)
    read(22,*) Nsize
    read(22,*)
    read(22,*) parmJ, parmh
    read(22,*)
    read(22,*) betamin, betamax, Nsplit
    read(22,*)
    read(22,*) Nstep
    read(22,*)
    read(22,*) samp_type
    close(unit=22)
end subroutine read_parms
end module ising_parms
