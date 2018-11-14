!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!-- Ising Simulation procedure
!-- (Copyright) Xin He <1500011805@pku.edu.cn>)
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


module ising_model
use ising_utils
implicit none
	logical :: PBC = .true.
	integer :: out_unit = 20
	integer, private :: set_Nstep
	real(8), private :: set_beta

type ising
	type(integer), dimension(:,:), allocatable :: grid
	integer :: sz
	real(8) :: J, h, b
	real(8) :: sumH, sumS
end type ising

contains
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!-- constructor of the ising type
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine alloc_ising(myising,isize)
	type(ising), intent(inout) :: myising
	integer, intent(in) :: isize
	integer :: i,j
	myising%sz = isize
	allocate(myising%grid(isize,isize))
	do j=1,isize
		do i=1,isize
			myising%grid(i,j)= 1-2*rand_int(0,1)
		enddo
	enddo
end subroutine alloc_ising


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!-- guess a distribution of ising state
!----- not a samll moving
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine guess_ising(myising)
	type(ising), intent(inout) :: myising
	integer :: i,j,sz
	if(.not. allocated(myising%grid)) stop "overstack error at guess_ising" 
	sz = myising%sz
	do j=1,sz
		do i=1,sz
			myising%grid(i,j)= 1-2*rand_int(0,1)
		enddo
	enddo
end subroutine guess_ising


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!-- set parameters of ising type
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine set_ising(myising,iJ,ih,ib)
	type(ising), intent(inout) :: myising
	real(8), intent(in) :: iJ, ih, ib
	myising%J = iJ
	myising%h = ih
	myising%b = ib
	myising%sumH = 0.0
	myising%sumS = 0.0
end subroutine set_ising


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!-- distructor of ising type
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine del_ising(myising)
	type(ising), intent(inout) :: myising
	if(allocated(myising%grid)) then
	    deallocate(myising%grid)
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
		i1m = mod(i1-2+myising%sz, myising%sz) + 1
		i2m = mod(i2-2+myising%sz, myising%sz) + 1
		i1p = mod(i1+myising%sz, myising%sz) + 1
		i2p = mod(i2+myising%sz, myising%sz) + 1
		
		H = - myising%h * myising%grid(i1,i2)
		H = H - myising%J * myising%grid(i1,i2)* myising%grid(i1p,i2)
		H = H - myising%J * myising%grid(i1,i2)* myising%grid(i1m,i2)
		H = H - myising%J * myising%grid(i1,i2)* myising%grid(i1,i2p)
		H = H - myising%J * myising%grid(i1,i2)* myising%grid(i1,i2m)
	else
		H = - myising%h * myising%grid(i1,i2)
		if(i1>1) then
			H = H - myising%J * myising%grid(i1,i2)* myising%grid(i1-1,i2)
		endif
		if(i2>1) then
			H = H - myising%J * myising%grid(i1,i2)* myising%grid(i1,i2-1)
		endif
		if(i1 < myising%sz) then
			H = H - myising%J * myising%grid(i1,i2)* myising%grid(i1+1,i2)
		endif
		if(i2 < myising%sz) then
			H = H - myising%J * myising%grid(i1,i2)* myising%grid(i1,i2+1)
		endif
	endif
end function local_H


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!-- return Sum of the Hamiltonian of the grid
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine sum_H(H, myising, iflag)
	real(8), intent(out):: H
	type(ising), intent(inout) :: myising
	integer, intent(in), optional :: iflag
	integer :: flag,i,j,sz
	if(present(iflag)) then
		flag = iflag
	else
		flag = 1 ! default, not using sum-local-H approach, recommend myising
	endif
	
	sz = myising%sz
	if(flag.eq.0) then
		H = 0.0
		do j=1,sz
			do i=1,sz
			    !-- sum of local hamiltonians counts for 2 times of total Hamiltonian 
				H = H + 0.5 * local_H(myising,i,j)
			enddo
		enddo
		if(PBC) then
		    myising%sumH = H
		    return
		endif
		!-- if not PBC, subtract the boundary interaction (subtract means + J*si*sj)
		do i=1,sz
			H = H + myising%J * myising%grid(1,i)* myising%grid(sz,i)
			H = H + myising%J * myising%grid(i,1)* myising%grid(i,sz)
		enddo
	elseif(flag.eq.1) then
		H = 0.0
		do j=1,sz-1
			do i=1,sz-1
				H = H - myising%h * myising%grid(i,j)
				H = H - myising%J * myising%grid(i,j)* myising%grid(i+1,j)
				H = H - myising%J * myising%grid(i,j)* myising%grid(i,j+1)
			enddo
			H = H - myising%h * myising%grid(sz,j)
			H = H - myising%J * myising%grid(sz,j)* myising%grid(sz,j+1)
			H = H - myising%h * myising%grid(j,sz)
			H = H - myising%J * myising%grid(j,sz)* myising%grid(j+1,sz)
		enddo
		H = H - myising%h * myising%grid(sz,sz)
		!-- if not PBC, end here
		if(.not. PBC) then
		    myising%sumH = H
		    return
		endif
		!-- otherwise, adding boundary interaction
		do i=1,sz
			H = H - myising%J * myising%grid(1,i)* myising%grid(sz,i)
			H = H - myising%J * myising%grid(i,1)* myising%grid(i,sz)
		enddo
	else
		stop "flag error at sum_H"
	endif
	myising%sumH = H
	return
end subroutine sum_H


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!-- return local Spin of the grid
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function local_S(myising, i1, i2) result(S)
	integer :: S
	type(ising), intent(in) :: myising
	integer, intent(in) :: i1,i2
	integer :: i1m,i1p,i2m,i2p,sz
	S = 0
	if(PBC) then
		i1m = mod(i1-1+myising%sz, myising%sz)
		i2m = mod(i2-1+myising%sz, myising%sz)
		i1p = mod(i1+1+myising%sz, myising%sz)
		i1p = mod(i2+1+myising%sz, myising%sz)
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
		if(i1 < myising%sz) then
			S = S + myising%grid(i1+1,i2)
		endif
		if(i2 < myising%sz) then
			S = S + myising%grid(i1,i2+1)
		endif
	endif
end function local_S


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!-- return Sum of the Spin of the grid
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine sum_S(S, myising)
	integer, intent(out) :: S
	type(ising), intent(inout) :: myising
	integer :: i,j,sz
	sz = myising%sz
	S = 0
	do j=1,sz
		do i=1,sz
			S = S + myising%grid(i,j)
		enddo
	enddo
	myising%sumS = S
end subroutine sum_S


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!-- Metropolis algorithm process (Gibbs sampling)
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine ising_gibbs(myising, iflag)
	type(ising), intent(inout) :: myising
	integer, optional :: iflag
	type(ising) :: tmp_ising
	integer :: sz, rn1, rn2, flag
	real(8) :: rn
	real(8) :: dH, tmp_H, old_H
	integer :: tmp_S
	
	if(present(iflag)) then
	    flag = iflag
	else
	    flag = 0  !-- move only one step
	endif
	
	sz = myising%sz
	if(flag .eq. 0) then
	    !-- generating preselection
	    rn1 = rand_int(1,sz)
	    rn2 = rand_int(1,sz)
	    
	    !-- metropolis algorithm
	    dH = -2 * local_H(myising,rn1,rn2)
	    
    	if(dH<=0) then
	    	myising%sumS = myising%sumS - 2 * myising%grid(rn1,rn2)
	    	!old_H = myising%sumH
	    	myising%sumH = myising%sumH + dH
	    	myising%grid(rn1,rn2) = - myising%grid(rn1,rn2)
	    	!call sum_H(tmp_H, myising)
	    	!if(tmp_H - old_H .ne. dH) then
	    	!    print *, tmp_H, old_H, dH
	    	!endif   	
	    else
	    	call random_number(rn)
	    	if(rn < exp(-myising%b * dH)) then
	    	    myising%sumS = myising%sumS - 2 * myising%grid(rn1,rn2)
	    	    !old_H = myising%sumH
	    	    myising%sumH = myising%sumH + dH
	    		myising%grid(rn1,rn2)= - myising%grid(rn1,rn2)
	    		!call sum_H(tmp_H, myising)
	        	!if(tmp_H - old_H .ne. dH) then
	        	!    print *, tmp_H, old_H, dH
	        	!endif
	    	endif
	    endif
	else if(flag .eq. 1) then
	    call alloc_ising(tmp_ising, sz)
	    call sum_H(tmp_H, tmp_ising)
	    call sum_S(tmp_S, tmp_ising)
	    dH = tmp_ising%sumH - myising%sumH
	    
	    if(dH<=0) then
	    	myising = tmp_ising
	    else
	    	call random_number(rn)
	    	if(rn < exp(-myising%b * dH)) then
	    	    myising = tmp_ising
	    	endif
	    endif
	    
	    call del_ising(tmp_ising)
	else
	    stop "flag error"
	endif
end subroutine ising_gibbs

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!-- Swendsen-Wang sampling
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
recursive subroutine adjoint_makesame(myising,i,j,bonds)
    type(ising), intent(inout) :: myising
    integer, intent(in) :: i,j
    integer, dimension(myising%sz,myising%sz,5), intent(inout) :: bonds
    integer :: sz, newi, newj
    sz = myising%sz
    !-- check up
    if(bonds(i,j,1) .eq. 1) then
        newi = mod(i-2+sz,sz)+1
        if(bonds(newi,j,5).eq. 1) then
            bonds(newi,j,5) = 0
            myising%grid(newi,j) = myising%grid(i,j)
            call adjoint_makesame(myising,newi,j,bonds)
        endif
    endif
    !-- check right
    if(bonds(i,j,2) .eq. 1) then
        newj = mod(j,sz)+1
        if(bonds(i,newj,5).eq. 1) then
            bonds(i,newj,5) = 0
            myising%grid(i,newj) = myising%grid(i,j)
            call adjoint_makesame(myising,i,newj,bonds)
        endif
    endif
    !-- check down
    if(bonds(i,j,3) .eq. 1) then
        newi = mod(i,sz)+1
        if(bonds(newi,j,5).eq. 1) then
            bonds(newi,j,5) = 0
            myising%grid(newi,j) = myising%grid(i,j)
            call adjoint_makesame(myising,newi,j,bonds)
        endif
    endif
    !-- check left
    if(bonds(i,j,4) .eq. 1) then
        newj = mod(j-2+sz,sz)+1
        if(bonds(i,newj,5).eq. 1) then
            bonds(i,newj,5) = 0
            myising%grid(i,newj) = myising%grid(i,j)
            call adjoint_makesame(myising,i,newj,bonds)
        endif
    endif
end subroutine adjoint_makesame

subroutine randby_area(myising, bonds)
    type(ising), intent(inout) :: myising
    integer, dimension(myising%sz,myising%sz,5), intent(inout) :: bonds
    integer :: i,j,sz,rn0

    sz=myising%sz
    do i=1,sz
        do j=1,sz
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
    integer, dimension(myising%sz, myising%sz,5) :: bonds
    real(8) :: irn, exp2bJ
    integer :: i,j,sz,newi,newj
    sz = myising%sz
    exp2bJ = dexp(2*myising%b*myising%J)
    bonds = 0
    if(PBC) then
        do i=1,sz
            do j=1,sz                                
                !-- bondij(2) --- (i,j) ~ (i,j+1)
                !-- bondij(4) --- (i,j) ~ (i,j-1)
                newj = mod(j,sz)+1
                if(myising%grid(i,j).eq.myising%grid(i,newj)) then
                    call random_number(irn)
                    if(irn*exp2bJ > 1) then
                        bonds(i,j,2) = 1
                        bonds(i,newj,4) = 1
                    endif
                endif
                
                !-- bondij(1) --- (i,j) ~ (i-1,j)
                !-- bondij(3) --- (i,j) ~ (i+1,j)
                newi = mod(i,sz)+1
                if(myising%grid(i,j).eq.myising%grid(newi,j)) then
                    call random_number(irn)
                    if(irn*exp2bJ > 1) then
                        bonds(i,j,3) = 1
                        bonds(newi,j,1) = 1
                    endif
                endif
                                               
                !-- bondij(5) --- (i,j) is flip or not
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
    integer, dimension(myising%sz, myising%sz,5), intent(inout) :: bonds
    integer :: newi,newj,sz
    real(8) :: exp2bJ,irn
    
    sz = myising%sz
    exp2bJ = dexp(2*myising%b*myising%J)
    bonds(i,j,5) = 1 !-- mark origin point
    
    !-- up
    newi = mod(i-2+sz,sz)+1
    if(myising%grid(i,j).ne.myising%grid(newi,j) .and. bonds(newi,j,5).eq. 0) then
        !print *, 'up'
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
    newj = mod(j,sz)+1
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
    newi = mod(i,sz)+1
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
    newj = mod(j-2+sz,sz)+1
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
    integer, dimension(myising%sz, myising%sz,5) :: bonds
    integer :: rn1, rn2, sz
    sz = myising%sz
    rn1 = rand_int(1,sz)
    rn2 = rand_int(1,sz)
    !print *, myising%grid+1
    !print *, 'choose ', rn1, rn2
    
    myising%grid(rn1,rn2) = - myising%grid(rn1,rn2)
    bonds = 0
    !print *, myising%grid+1
    call next_wolff(myising,rn1,rn2,bonds)
    !print *, myising%grid+1
    !stop 'debug'
end subroutine ising_wolff

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!-- KMC
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine ising_kmc(myising)
    type(ising), intent(inout) :: myising
    integer :: rn1, rn2, rn
end subroutine ising_kmc

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!-- Quench simulation
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!-- Tempering
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!-- Kinetic Monte Carlo (BKL algorithm)
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!-- Simulation at a specific beta 
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine ising_simul(myising, N, ibeta)
use ising_parms, only:samp_type
    type(ising), intent(inout) :: myising
	integer, intent(in) :: N
	real(8), intent(in) :: ibeta
	real(8) :: nowJ, nowh, rcdH, N2
	integer :: rcdS, i
	
	!-- set model paramaters
	set_Nstep = N
	set_beta = ibeta
	N2 = myising%sz ** 2
	
	!-- init. guess
	nowJ = myising%J
	nowh = myising%h
	call set_ising(myising, nowJ, nowH, set_beta)
	call guess_ising(myising)
	call sum_H(rcdH, myising)
	call sum_S(rcdS, myising)
	!print *, myising%grid
	!stop
	if(samp_type .eq. 0) then
	    do i=1,set_Nstep
		    call ising_gibbs(myising)
		    !-- needn't, for if we recise H,S at one step moving
            !call sum_H(rcdH, myising) !-- for metropolis, needn't
	        !call sum_S(rcdS, myising) !-- for metropolis, needn't

		    write(out_unit,*) i, myising%sumH / N2, myising%sumS / N2
	    enddo
	elseif(samp_type .eq. 1) then
	    do i=1,set_Nstep
            !-- swendsen-wang
            call ising_swendsenwang(myising)
            call sum_H(rcdH, myising)
	        call sum_S(rcdS, myising)
		    write(out_unit,*) i, myising%sumH / N2, myising%sumS / N2
	    enddo
	elseif(samp_type .eq. 2) then
	    do i=1,set_Nstep
            !-- wolff
            call ising_wolff(myising)
            call sum_H(rcdH, myising)
	        call sum_S(rcdS, myising)
		    write(out_unit,*) i, myising%sumH / N2, myising%sumS / N2
	    enddo
	else
	    stop 'unknown args'
	endif
end subroutine ising_simul

end module ising_model

program main
use ising_model
use ising_parms
implicit none
    type(ising) :: is1
	real(8) :: beta
	character(len=10) :: idxcs
	integer :: i
	
	PBC = .true.
	call init_seed()
	
	call read_parms('ising.parms')
	
	call alloc_ising(is1, Nsize)
	
	call set_ising(is1, parmJ, parmh, betamax)
	
	do i=1,Nsplit
		beta = betamin + betamax * (i/real(Nsplit))
		write(idxcs,'(i4)') i
		open(unit=out_unit, file='beta_i'//trim(adjustl(idxcs))//'.dat', status='replace')
		call ising_simul(is1, Nstep, beta)
		close(out_unit)
	enddo
	
	call del_ising(is1)
end program



