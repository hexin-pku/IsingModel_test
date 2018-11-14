FF=gfortran
src=ising_utils.f90 ising_parms.f90 ising_model.f90
exe=ising.run
default:
	${FF} ${src} -o ${exe}
clean:
	rm *.mod ${exe}	 
