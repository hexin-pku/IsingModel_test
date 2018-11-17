module ising_utils
implicit none
contains
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!-- init. random seed
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine init_seed()
  integer :: n, ival(8), v(3), i
  integer, allocatable :: seed(:)
  call date_and_time(values=ival)
  v(1) = ival(8) + 2048*ival(7)
  v(2) = ival(6) + 64*ival(5)     ! value(4) isn't really 'random'
  v(3) = ival(3) + 32*ival(2) + 32*8*ival(1)
  call random_seed(size=n)
  allocate(seed(n))
  call random_seed()   ! Give the seed an implementation-dependent kick
  call random_seed(get=seed)
  do i=1, n
     seed(i) = seed(i) + v(mod(i-1, 3) + 1)
  enddo
  call random_seed(put=seed)
  i = irand(seed(1)+1)
  deallocate(seed)
end subroutine init_seed

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!-- return random integer in [i1,i2]
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function rand_int(i1,i2) result(rn)
	integer, intent(in) :: i1,i2
	integer :: rn
	rn = irand()
	rn = mod(rn,i2-i1+1) + i1
end function rand_int

end module

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!-- test
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!program test
!use ising_utils
!implicit none
!    integer :: i, rn1
!    real(8) :: rn2
!    call init_seed()
!    do i=1, 10
!        call random_number(rn2)
!        rn1 = rand_int(0,10)
!        print *, rn1, rn2
!    enddo
!end program



