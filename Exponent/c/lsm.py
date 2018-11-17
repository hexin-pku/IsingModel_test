import numpy
import scipy.optimize as scimin
import matplotlib.pyplot as mpl

datax=numpy.array([1,2,3,4,5]) # data coordinates
datay=numpy.array([2.95,6.03,11.2,17.7,26.8])
constraintmaxx=numpy.array([0]) # list of maximum constraints
constraintmaxy=numpy.array([1.2])

# least square fit without constraints
def fitfunc(x,p): # model $f(x)=a x^2+c
    a,c=p
    return c+a*x**2
def residuals(p): # array of residuals
    return datay-fitfunc(datax,p)
p0=[1,2] # initial parameters guess
pwithout,cov,infodict,mesg,ier=scimin.leastsq(residuals, p0,full_output=True) #traditionnal least squares fit

# least square fir with constraints
def sum_residuals(p): # the function we want to minimize
    return sum(residuals(p)**2)
def constraints(p): # the constraints: all the values of the returned array will be >=0 at the end
    return constraintmaxy-fitfunc(constraintmaxx,p)
pwith=scimin.fmin_slsqp(sum_residuals,pwithout,f_ieqcons=constraints) # minimization with constraint

# plotting
ax=mpl.figure().add_subplot(1,1,1)
ax.plot(datax,datay,ls="",marker="x",color="blue",mew=2.0,label="Datas")
ax.plot(constraintmaxx,constraintmaxy,ls="",marker="x",color="red",mew=2.0,label="Max points")
morex=numpy.linspace(0,6,100)
ax.plot(morex,fitfunc(morex,pwithout),color="blue",label="Fit without constraints")
ax.plot(morex,fitfunc(morex,pwith),color="red",label="Fit with constraints")
ax.legend(loc=2)
mpl.show()
