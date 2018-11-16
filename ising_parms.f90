module ising_parms
implicit none
private
    public ::  Nsize, Nstep, samp_type, work_type, parmb, parmJ, parmh
    real(8) :: parmb, parmJ, parmh
    integer :: Nsize, Nstep, samp_type, work_type
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
    read(22,*) Nstep
    read(22,*)
    read(22,*) parmb, parmJ, parmh
    read(22,*)
    read(22,*) samp_type
    read(22,*)
    read(22,*) work_type
    close(unit=22)
end subroutine read_parms
end module ising_parms
