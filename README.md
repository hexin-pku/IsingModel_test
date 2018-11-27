
<h4 style="text-align:right">CCME,&nbsp; Peking University</h4>
***
<p>
<h1 style="text-align:center"><font face="Times Roman" size=6> Ising Model Simulation </font></h1>
<br>
<h2 style="text-align:center"><font face="Times Roman">Xin He, 1500011805 </font> </h2>

# Source

the program contains gibbs sampling Metropolis, Swendsen-Wang, Wolff algorithms. (but the KMC is still need added.)   
For different problems, here we may use different algorithms.  
the source code can be found at the author's mainpage of github, for [IsingModel](https://github.com/XShinHe/IsingModel).

### the procedure implement

Here gives the `Makefile` by author, make sure `gfortran` is available. Then you can
compile the program by `make` to give a executable procedure `ising.run`.  

you can run the program in two ways: 
>* run without any arguments  
`./ising.run`  
the program will read $\beta, J, h$ from default parameter file --- `is.parms`.   
>* run with only 3 arguments   
`./ising.run  input_beta  input_J  input_h`  
so the relative arguments will be read from command line.

the ising simulation program has 4 parts.  
* ising_utils.f90:  
> gives some useful subroutines for generating rand number, contains **init_seed, rand_int**

* ising_parms.f90:  
> read parameters from the file (default **is.parms**), and sharing the public paramters.   
contains **Nsize, Nstep, parmb, parmJ, parmh, samp_type**  

the simulattion parameter:  

|name      | mean |
|----------|---|
|Nsize         | the size of 2-D grid of ising model |
|Nstep         | the total steps for simulation  |
|parmb         | beta ($\beta=1/k_B T$)  |
|parmJ         | coupling coefficient    |
|parmh         | the external magnetic field |


the sampling setting: 

|samp_type | method |
|----------|---|
|0         | Metropolis with gibbs sampling |
|1         | Swendsen-Wang sampling  |
|2         | Wolff sampling  |

the work type setting (suit for specific problem)

|work_type | method |
|----------|---|
|1         | calc both c and m, is for problem 1 |
|2         | calc only m, is for problem 2  |
|3         | calc only correlation, is for problem 3|
|4         | calc all, not support?|

* ising_model.f90:  
>* define **ising** type
>* support Metropolis with gibbs sampling
>* support Swendsen-Wang sampling
>* support Wolff sampling

* the control scripts/anaysis scripts  
>* here use shell scripts to give a series simulation, ref `run.sh.example`   
>* here use python3 scripts to analysis the result
  
