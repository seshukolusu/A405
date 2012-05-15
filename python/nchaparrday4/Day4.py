import numpy as np
import matplotlib.pyplot as plt
import scipy.io

def probfv(C,a,v, T):    
    f_v = (C*v**2)*np.exp(1.0*(-a*v**2)/(2))
    return f_v 

def calc_testmean(v, fv):
    testmean = 0
    for i in range(len(v)):
        testmean = testmean + v[i]*fv[i]
    return testmean

#Constants and Factors
k_B = 1.38065*10**(-23) # m^2 kg s^-2 K^-1, Boltzman Constant
Na = 6.0221415*10**(23) #avogadro's number
m =(1.0*28/Na)*10**(-3) #kg, mass of Nitrogen molecule (N2)
T = np.array([270,370]) #K, the two given temperatures
a = 1.0*m/(k_B*T) #factor
C = np.sqrt((1.0*2/np.pi)*(a)**3) #factor 

#Calculations
v = np.arange(0, 2500, 1) #molecular velocity array
fv0 = probfv(C[0], a[0], v, T[0]) #probability distribution function of velocity at 270K
fv1 = probfv(C[1], a[1], v, T[1]) #probability distribution function of velocity at 370K
mean = np.array([np.dot(fv0, v), np.dot(fv1, v)]) #mean velocties at 270K and 370K, calculated by dot product of velocity and pdf matrices
testmean = np.array([calc_testmean(v, fv0), calc_testmean(v, fv1)])# as above except by elemental multiplication and summing
anmean = np.array([2*C[0]/a[0]**2, 2*C[1]/a[1]**2]) #analyitic mean from integration by parts of vf(v)


#Printouts
print 'mean velocities', mean
print 'test on mean velocities', testmean
print 'analytic mean velocities', anmean

#plots
fig = plt.figure(1)
ax1=fig.add_subplot(111)
ax1.set_title('Molecular Velocity vs its Probability Function')
ax1.set_xlabel('v (m/s)')
ax1.set_ylabel('f(v)')
ax1.plot(v, fv0, 'b', label='at 270K')
ax1.plot(v, fv1, 'r', label='at 370K')
ax1.legend(loc='upper right')
plt.show()






    


