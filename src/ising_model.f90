!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!-- Ising Simulation procedure
!-- (Copyright) Xin He <1500011805@pku.edu.cn>)
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


module ising_model
use ising_utils
use ising_parms
implicit none
	logical :: PBC = .true.     !-- if boundary condition
	logical :: PRT = .false.    !-- if print info at each step
	
	integer :: out_unit = 20
	
	real(8), dimension(:), allocatable :: array_H, array_S
	real(8), dimension(:), allocatable :: array_G
	

type ising
	type(integer), dimension(:,:), allocatable :: grid
	integer :: sz
	real(8) :: avgH, avgS
end type ising

contains
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!-- constructor of the ising type
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine alloc_ising(myising,isize)
	type(ising), intent(inout) :: myising
	integer, intent(in) :: isize
	integer :: i,j
	Nsize = isize
	allocate(myising%grid(isize,isize))
	do j=1,isize
		do i=1,isize
			myising%grid(i,j)= 1-2*rand_int(0,1)
		enddo
	enddo
	!-- global allocating
	if(.not. allocated(array_H)) then
	    allocate(array_H(Nstep))
	    allocate(array_S(Nstep))
	    allocate(array_G(Nsize/2))
	endif
end subroutine alloc_ising


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!-- guess a distribution of ising state
!----- not a samll moving
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine guess_ising(myising)
	type(ising), intent(inout) :: myising
	integer :: i,j
	if(.not. allocated(myising%grid)) stop "overstack error at guess_ising" 
	do j=1,Nsize
		do i=1,Nsize
			myising%grid(i,j)= 1-2*rand_int(0,1)
		enddo
	enddo
end subroutine guess_ising


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!-- distructor of ising type
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine del_ising(myising)
	type(ising), intent(inout) :: myising
	if(allocated(myising%grid)) then
	    deallocate(myising%grid)
	endif
	!-- gloabal deallocating
	if(allocated(array_H)) then
	    deallocate(array_H)
	    deallocate(array_S)
	    deallocate(array_G)
	endif
end subroutine del_ising


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!-- return local Hamiltonian of the grid(i,j)
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function local_H(myising, i1, i2) result(H)
    real(8) :: H
	type(ising), intent(in) :: myising
	integer, intent(in) :: i1,i2
	integer :: i1m,i1p,i2m,i2p
	if(PBC) then
		i1m = mod(i1-2+Nsize, Nsize) + 1
		i2m = mod(i2-2+Nsize, Nsize) + 1
		i1p = mod(i1+Nsize, Nsize) + 1
		i2p = mod(i2+Nsize, Nsize) + 1
		
		H = - parmh * myising%grid(i1,i2)
		H = H - parmJ * myising%grid(i1,i2)* myising%grid(i1p,i2)
		H = H - parmJ * myising%grid(i1,i2)* myising%grid(i1m,i2)
		H = H - parmJ * myising%grid(i1,i2)* myising%grid(i1,i2p)
		H = H - parmJ * myising%grid(i1,i2)* myising%grid(i1,i2m)
	else
		H = - parmh * myising%grid(i1,i2)
		if(i1>1) then
			H = H - parmJ * myising%grid(i1,i2)* myising%grid(i1-1,i2)
		endif
		if(i2>1) then
			H = H - parmJ * myising%grid(i1,i2)* myising%grid(i1,i2-1)
		endif
		if(i1 < Nsize) then
			H = H - parmJ * myising%grid(i1,i2)* myising%grid(i1+1,i2)
		endif
		if(i2 < Nsize) then
			H = H - parmJ * myising%grid(i1,i2)* myising%grid(i1,i2+1)
		endif
	endif
end function local_H


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!-- return Sum of the Hamiltonian of the grid
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine avg_H(H, myising)
	real(8), intent(out):: H
	type(ising), intent(inout) :: myising
	integer :: i,j

    H = -parmh * sum(myising%grid)
    H = H - parmJ * sum( myising%grid(1:Nsize-1,:) * myising%grid(2:Nsize,:) )
    if(PBC) H = H - parmJ * sum( myising%grid(1,:) * myising%grid(Nsize,:) )
    H = H - parmJ * sum( myising%grid(:,1:Nsize-1) * myising%grid(:,2:Nsize) )
    if(PBC) H = H - parmJ * sum( myising%grid(1,:) * myising%grid(Nsize,:) )

	H = H/(Nsize**2)
	myising%avgH = H
	return
end subroutine avg_H


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!-- return local Spin of the grid
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function local_S(myising, i1, i2) result(S)
	integer :: S
	type(ising), intent(in) :: myising
	integer, intent(in) :: i1,i2
	integer :: i1m,i1p,i2m,i2p,Nsize
	S = 0
	if(PBC) then
		i1m = mod(i1-1+Nsize, Nsize)
		i2m = mod(i2-1+Nsize, Nsize)
		i1p = mod(i1+1+Nsize, Nsize)
		i1p = mod(i2+1+Nsize, Nsize)
		S = myising%grid(i1,i2)
		S = S + myising%grid(i1p,i2)
		S = S + myising%grid(i1m,i2)
		S = S + myising%grid(i1,i2p)
		S = S + myising%grid(i1,i2m)
	else
		S = myising%grid(i1,i2)
		if(i1>0) then
			S = S + myising%grid(i1-1,i2)
		endif
		if(i2>0) then
			S = S + myising%grid(i1,i2-1)
		endif
		if(i1 < Nsize) then
			S = S + myising%grid(i1+1,i2)
		endif
		if(i2 < Nsize) then
			S = S + myising%grid(i1,i2+1)
		endif
	endif
end function local_S


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!-- return Sum of the Spin of the grid
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine avg_S(S, myising)
	real(8), intent(out) :: S
	type(ising), intent(inout) :: myising
	integer :: i,j
	S = sum(myising%grid)
	S = S/(Nsize**2)
	myising%avgS = S
end subroutine avg_S


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!-- Metropolis algorithm process (Gibbs sampling)
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine ising_gibbs(myising)
	type(ising), intent(inout) :: myising
	integer :: rn1, rn2
	real(8) :: rn
	real(8) :: dH, tmp_S, tmp_H, old_H
	
	!-- generating preselection
	rn1 = rand_int(1,Nsize)
	rn2 = rand_int(1,Nsize)
	
	!-- metropolis algorithm
	dH = -2 * local_H(myising,rn1,rn2)
	
    if(dH<=0) then
		myising%avgS = myising%avgS - 2 * myising%grid(rn1,rn2) / real(Nsize**2)
		myising%avgH = myising%avgH + dH / real(Nsize**2)
		myising%grid(rn1,rn2) = - myising%grid(rn1,rn2)	
	else
		call random_number(rn)
		if(rn < exp(-parmb * dH)) then
		    myising%avgS = myising%avgS - 2 * myising%grid(rn1,rn2) / real(Nsize**2)
		    myising%avgH = myising%avgH + dH / real(Nsize**2)
		    myising%grid(rn1,rn2) = - myising%grid(rn1,rn2)
		endif
	endif
	
end subroutine ising_gibbs

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!-- Swendsen-Wang sampling
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
recursive subroutine adjoint_makesame(myising,i,j,bonds)
    type(ising), intent(inout) :: myising
    integer, intent(in) :: i,j
    integer, dimension(Nsize,Nsize,5), intent(inout) :: bonds
    integer :: newi, newj

    !-- check up
    if(bonds(i,j,1) .eq. 1) then
        newi = mod(i-2+Nsize,Nsize)+1
        if(bonds(newi,j,5).eq. 1) then
            bonds(newi,j,5) = 0
            myising%grid(newi,j) = myising%grid(i,j)
            call adjoint_makesame(myising,newi,j,bonds)
        endif
    endif
    !-- check right
    if(bonds(i,j,2) .eq. 1) then
        newj = mod(j,Nsize)+1
        if(bonds(i,newj,5).eq. 1) then
            bonds(i,newj,5) = 0
            myising%grid(i,newj) = myising%grid(i,j)
            call adjoint_makesame(myising,i,newj,bonds)
        endif
    endif
    !-- check down
    if(bonds(i,j,3) .eq. 1) then
        newi = mod(i,Nsize)+1
        if(bonds(newi,j,5).eq. 1) then
            bonds(newi,j,5) = 0
            myising%grid(newi,j) = myising%grid(i,j)
            call adjoint_makesame(myising,newi,j,bonds)
        endif
    endif
    !-- check left
    if(bonds(i,j,4) .eq. 1) then
        newj = mod(j-2+Nsize,Nsize)+1
        if(bonds(i,newj,5).eq. 1) then
            bonds(i,newj,5) = 0
            myising%grid(i,newj) = myising%grid(i,j)
            call adjoint_makesame(myising,i,newj,bonds)
        endif
    endif
end subroutine adjoint_makesame


subroutine randby_area(myising, bonds)
    type(ising), intent(inout) :: myising
    integer, dimension(Nsize,Nsize,5), intent(inout) :: bonds
    integer :: i,j,rn0

    do i=1,Nsize
        do j=1,Nsize
            if(bonds(i,j,5).eq.1) then
                bonds(i,j,5) = 0
                rn0 = 1 - 2*rand_int(0,1)
                myising%grid(i,j) = rn0
                call adjoint_makesame(myising,i,j,bonds)
            endif
        enddo
    enddo
end subroutine randby_area

subroutine ising_swendsenwang(myising)
    type(ising), intent(inout) :: myising
    integer, dimension(Nsize, Nsize,5) :: bonds
    real(8) :: irn, exp2bJ
    integer :: i,j,newi,newj

    exp2bJ = dexp(2*parmb*parmJ)
    bonds = 0
    if(PBC) then
        do i=1,Nsize
            do j=1,Nsize                                
                !-- bondij(2) --- (i,j) ~ (i,j+1)
                !-- bondij(4) --- (i,j) ~ (i,j-1)
                newj = mod(j,Nsize)+1
                if(myising%grid(i,j).eq.myising%grid(i,newj)) then
                    call random_number(irn)
                    if(irn*exp2bJ > 1) then
                        bonds(i,j,2) = 1
                        bonds(i,newj,4) = 1
                    endif
                endif
                
                !-- bondij(1) --- (i,j) ~ (i-1,j)
                !-- bondij(3) --- (i,j) ~ (i+1,j)
                newi = mod(i,Nsize)+1
                if(myising%grid(i,j).eq.myising%grid(newi,j)) then
                    call random_number(irn)
                    if(irn*exp2bJ > 1) then
                        bonds(i,j,3) = 1
                        bonds(newi,j,1) = 1
                    endif
                endif
                                               
                !-- bondij(5) --- record (i,j) flip status
                bonds(i,j,5) = 1
            enddo
        enddo
        call randby_area(myising,bonds)
    else
        stop "not support"
    endif
end subroutine ising_swendsenwang


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!-- Wolff sampling
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
recursive subroutine next_wolff(myising,i,j,bonds)
    type(ising), intent(inout) :: myising
    integer, intent(in) :: i,j
    integer, dimension(Nsize, Nsize,5), intent(inout) :: bonds
    integer :: newi,newj
    real(8) :: exp2bJ,irn
    
    exp2bJ = dexp(2*parmb*parmJ)
    bonds(i,j,5) = 1 !-- mark origin point
    
    !-- up
    newi = mod(i-2+Nsize,Nsize)+1
    if(myising%grid(i,j).ne.myising%grid(newi,j) .and. bonds(newi,j,5).eq. 0) then
        call random_number(irn)
        if(irn*exp2bJ > 1) then
            !print *, 'rand suc'
            bonds(i,j,1) = 1
            bonds(newi,j,3) = 1
            myising%grid(newi,j) = myising%grid(i,j)
            call next_wolff(myising,newi,j,bonds)
        endif
    endif
    
    !-- right
    newj = mod(j,Nsize)+1
    if(myising%grid(i,j).ne.myising%grid(i,newj) .and. bonds(i,newj,5).eq. 0) then
        !print *, 'right'
        call random_number(irn)
        if(irn*exp2bJ > 1) then
            !print *, 'rand suc'
            bonds(i,j,2) = 1
            bonds(i,newj,4) = 1
            myising%grid(i,newj) = myising%grid(i,j)
            call next_wolff(myising,i,newj,bonds)
        endif
    endif
    
    !-- down
    newi = mod(i,Nsize)+1
    if(myising%grid(i,j).ne.myising%grid(newi,j) .and. bonds(newi,j,5).eq. 0) then
        !print *, 'down'
        call random_number(irn)
        if(irn*exp2bJ > 1) then
            !print *, 'rand suc'
            bonds(i,j,3) = 1
            bonds(newi,j,1) = 1
            myising%grid(newi,j) = myising%grid(i,j)
            call next_wolff(myising,newi,j,bonds)
        endif
    endif
    
    !-- left
    newj = mod(j-2+Nsize,Nsize)+1
    if(myising%grid(i,j).ne.myising%grid(i,newj) .and. bonds(i,newj,5).eq. 0) then
        !print *, 'left'
        call random_number(irn)
        if(irn*exp2bJ > 1) then
            !print *, 'rand suc'
            bonds(i,j,4) = 1
            bonds(i,newj,2) = 1
            myising%grid(i,newj) = myising%grid(i,j)
            call next_wolff(myising,i,newj,bonds)
        endif
    endif 
end subroutine next_wolff

subroutine ising_wolff(myising)
    type(ising), intent(inout) :: myising
    integer, dimension(Nsize, Nsize,5) :: bonds
    integer :: rn1, rn2

    rn1 = rand_int(1,Nsize)
    rn2 = rand_int(1,Nsize)
    
    !-- flip a random site
    myising%grid(rn1,rn2) = - myising%grid(rn1,rn2)
    bonds = 0
    !-- (recursively) flip its neighborhood site according bonds
    call next_wolff(myising,rn1,rn2,bonds)

end subroutine ising_wolff


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!-- Kinetic Monte Carlo (BKL algorithm)
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine ising_kmc(myising)
    type(ising), intent(inout) :: myising
    integer :: rn1, rn2, rn
    stop 'not support now'
end subroutine ising_kmc


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!-- Simulation at a specific beta 
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine ising_simul(myising)
use ising_parms
    type(ising), intent(inout) :: myising
	real(8) :: rcdS, rcdH
	real(8) :: avgH, stdH, avgM, stdM
	integer :: i
	
	call guess_ising(myising)
	call avg_H(rcdH, myising)
	call avg_S(rcdS, myising)
	!print *, myising%grid
	!stop
	if(samp_type .eq. 0) then
	    do i=1,Nstep
		    call ising_gibbs(myising)
		    !-- needn't, for if we recise H,S at one step moving
            !call avg_H(rcdH, myising) !-- for metropolis, needn't
	        !call avg_S(rcdS, myising) !-- for metropolis, needn't
            array_H(i) = myising%avgH
            array_S(i) = myising%avgS
		    if(PRT) write(out_unit,*) i, myising%avgH , myising%avgS 
	    enddo
	elseif(samp_type .eq. 1) then
	    do i=1,Nstep
            !-- swendsen-wang
            call ising_swendsenwang(myising)
            call avg_H(rcdH, myising)
	        call avg_S(rcdS, myising)
	        array_H(i) = myising%avgH
            array_S(i) = myising%avgS
		    if(PRT) write(out_unit,*) i, myising%avgH , myising%avgS 
	    enddo
	elseif(samp_type .eq. 2) then
	    do i=1,Nstep
            !-- wolff
            call ising_wolff(myising)
            call avg_H(rcdH, myising)
	        call avg_S(rcdS, myising)
	        array_H(i) = myising%avgH
            array_S(i) = myising%avgS
		    if(PRT) write(out_unit,*) i, myising%avgH , myising%avgS 
	    enddo
	else
	    stop 'unknown args'
	endif
	
	avgH = sum(array_H) / Nstep
	stdH = sum(array_H**2/Nstep) - avgH**2
	avgM = sum(array_S) / Nstep
	stdM = sum(array_S**2/Nstep) - avgM**2
	
	print *, avgH, stdH, avgM, stdM
end subroutine ising_simul

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!-- Simulation at a specific magnetic field
!-- for problem 2
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

subroutine ising_simul2(myising)
use ising_parms
    type(ising), intent(inout) :: myising
	real(8) :: rcdS
	real(8) :: output
	integer :: i
	
	output = 0
	
	call guess_ising(myising)
	call avg_S(rcdS, myising)
	!print *, myising%grid
	!stop
	if(samp_type .eq. 0) then
	    do i=1,Nstep
		    call ising_gibbs(myising)
	        !call avg_S(rcdS, myising) !-- for metropolis, needn't
		    if(PRT) write(out_unit,*) i, myising%avgS
		    output = ( (i-1)*output + myising%avgS ) / real(i) 
	    enddo
	elseif(samp_type .eq. 1) then
	    do i=1,Nstep
            !-- swendsen-wang
            call ising_swendsenwang(myising)
	        call avg_S(rcdS, myising)
		    if(PRT) write(out_unit,*) i, myising%avgS
		    output = ( (i-1)*output + myising%avgS ) / real(i) 
	    enddo
	elseif(samp_type .eq. 2) then
	    do i=1,Nstep
            !-- wolff
            call ising_wolff(myising)
	        call avg_S(rcdS, myising)
		    if(PRT) write(out_unit,*) i, myising%avgS 
		    output = ( (i-1)*output + myising%avgS ) / real(i)
	    enddo
	else
	    stop 'unknown args'
	endif
	print *, output
end subroutine ising_simul2

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!-- Spatial Correlation
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine ising_spincorrel(gd,i)
    integer, dimension(Nsize,Nsize), intent(in) :: gd
    integer, intent(in) :: i
    integer, dimension(Nsize/2) :: tmp
    integer :: j,k,effk
    
    if(i.eq.1) array_G = 0
    
    tmp = 0
    do j=1,Nsize
        do k=1,Nsize/2
            effk = mod(j+k-1, Nsize) + 1
            tmp(k) = tmp(k) + sum(gd(:,j)*gd(:,effk))
            tmp(k) = tmp(k) + sum(gd(j,:)*gd(effk,:))
        enddo
    enddo
    
    array_G = ( array_G*(i-1) + tmp/real(2*Nsize**2) ) / real(i)
    
end subroutine ising_spincorrel

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!-- Configuration (Spin) correlation
!-- for problem 3 
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine ising_configcorrel(myising)
use ising_parms
    type(ising), intent(inout) :: myising
	real(8) :: rcdS, rcdH
	real(8) :: avgH, stdH, avgM, stdM
	integer :: i
	
	call guess_ising(myising)
	call avg_H(rcdH, myising)
	call avg_S(rcdS, myising)
	!print *, myising%grid
	!stop
	if(samp_type .eq. 0) then
	    do i=1,Nstep
		    call ising_gibbs(myising)
            call ising_spincorrel(myising%grid, i)
	    enddo
	elseif(samp_type .eq. 1) then
	    do i=1,Nstep
            !-- swendsen-wang
            call ising_swendsenwang(myising)
            call ising_spincorrel(myising%grid, i)
	    enddo
	elseif(samp_type .eq. 2) then
	    do i=1,Nstep
            !-- wolff
            call ising_wolff(myising)
            call ising_spincorrel(myising%grid, i)
	    enddo
	else
	    stop 'unknown args'
	endif

	print *, array_G
end subroutine ising_configcorrel

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!-- Simulation at a specific beta 
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine ising_simul4(myising)
use ising_parms
    type(ising), intent(inout) :: myising
	real(8) :: rcdS, rcdH
	real(8) :: avgH, stdH, avgM, stdM
	integer :: i
	
	call guess_ising(myising)
	call avg_H(rcdH, myising)
	call avg_S(rcdS, myising)
	!print *, myising%grid
	!stop
	if(samp_type .eq. 0) then
	    do i=1,Nstep
		    call ising_gibbs(myising)
		    !-- needn't, for if we recise H,S at one step moving
            !call avg_H(rcdH, myising) !-- for metropolis, needn't
	        !call avg_S(rcdS, myising) !-- for metropolis, needn't
	        call ising_spincorrel(myising%grid, i)
            array_H(i) = myising%avgH
            array_S(i) = myising%avgS
		    if(PRT) write(out_unit,*) i, myising%avgH , myising%avgS 
	    enddo
	elseif(samp_type .eq. 1) then
	    do i=1,Nstep
            !-- swendsen-wang
            call ising_swendsenwang(myising)
            call avg_H(rcdH, myising)
	        call avg_S(rcdS, myising)
	        call ising_spincorrel(myising%grid, i)
	        array_H(i) = myising%avgH
            array_S(i) = myising%avgS
		    if(PRT) write(out_unit,*) i, myising%avgH , myising%avgS 
	    enddo
	elseif(samp_type .eq. 2) then
	    do i=1,Nstep
            !-- wolff
            call ising_wolff(myising)
            call avg_H(rcdH, myising)
	        call avg_S(rcdS, myising)
	        call ising_spincorrel(myising%grid, i)
	        array_H(i) = myising%avgH
            array_S(i) = myising%avgS
		    if(PRT) write(out_unit,*) i, myising%avgH , myising%avgS 
	    enddo
	else
	    stop 'unknown args'
	endif
	
	avgH = sum(array_H) / Nstep
	stdH = sum(array_H**2/Nstep) - avgH**2
	avgM = sum(array_S) / Nstep
	stdM = sum(array_S**2/Nstep) - avgM**2
	
	print *, avgH, stdH, avgM, stdM, array_G
end subroutine ising_simul4

end module ising_model

program main
use ising_model
use ising_parms
implicit none
    type(ising) :: is1
	character(len=20) :: tmp
	integer :: n
	
	PBC = .true.
	call init_seed()
	
	call read_parms('is.parms')
	n = command_argument_count()
    if (n .ne. 0 .and. n .ne. 3) then
        stop 'args number mismatch, for: beta, J, h'
	end if
	if(n .eq. 3) then
	    call get_command_argument(1,tmp)
	    read(tmp,*) parmb
	    call get_command_argument(2,tmp)
	    read(tmp,*) parmJ
	    call get_command_argument(3,tmp)
	    read(tmp,*) parmh
	endif
	
	call alloc_ising(is1, Nsize)
	
	if(work_type .eq. 1) then
	    if(PRT) open(unit=out_unit, file='example.dat', status='replace')
	    call ising_simul(is1)
	    if(PRT) close(out_unit)
	elseif(work_type .eq. 2) then
	    if(PRT) open(unit=out_unit, file='example.dat', status='replace')
	    call ising_simul2(is1)
	    if(PRT) close(out_unit)
	elseif(work_type .eq. 3) then
	    call ising_configcorrel(is1)
	elseif(work_type .eq. 4) then
	    call ising_simul4(is1)
	else
	    stop 'wrong work type'
	endif
	
	call del_ising(is1)
end program main



