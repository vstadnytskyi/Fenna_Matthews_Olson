import numpy as np
import sys
import os
import matplotlib.pyplot as plt
import math
import scipy.linalg
plt.close("all") #close all preexisting figures

DEBUG = False #this is a debug variable: True - will show debug values


def _debug(statement):
    if DEBUG:
        exec(statement)

def gaussian_func(x,amp,mean,sigma):
    y = x*0
    for i in range(len(x)):
        y[i] = amp*np.exp(-(x[i]-mean)**2/(2*(sigma/2.355)**2))
    return y
def exp_func(x,amp,zero,sigma):
    y = x*0
    for i in range(len(x)):
        y[i] = amp*np.exp(-(x[i]-zero)*sigma)
    return y
def linear_func(y_exp,x):
    y = x*0
    for i in range(len(x)):
        y[i] = -(2900.-y_exp[len(x)-1])*(x[i]-500)/(x[len(x)-1]-500) + 2900
    return y
data = np.loadtxt('BChl_a_extinction.txt')
dipole_moments_hat = np.loadtxt('dipole_moments.txt')

_debug('print dipole_moments_hat[1,:]')
muhat = dipole_moments_hat
_debug('print "sum muhat=,", sum(muhat)')
_debug('print np.linalg.norm(muhat[0,:])')
data.shape
data[:,0]



wavelength_lst = np.array([783., 770., 750., 730., 710., 690.,580.])
amp_lst = np.array([94000., 60000., 20000., 12000.000, 11000., 8000.,30380.])
x = data[:,0]
y_exp = data[:,1]
width = 25
y_fit_temp = x*0
for i in range(len(amp_lst)):
    y_fit_temp = y_fit_temp + gaussian_func(x,amp_lst[i],wavelength_lst[i],width)
    
y_fit_temp = y_fit_temp + linear_func(y_exp,x)
y_fit = (y_fit_temp/np.max(y_fit_temp))*np.max(y_exp)

plt.figure(1)
plt.ion()
plt.plot(x,y_exp)
plt.plot(x,y_fit)

plt.legend()
markerline, stemlines, baseline  = plt.stem(wavelength_lst,amp_lst, markerfmt=' ')
plt.setp(stemlines, 'color', 'r', 'linewidth', 5)
plt.setp(baseline, 'color', 'k', 'linewidth', 1)
plt.legend(['Measured absorption spectrum','fitted spectrum','BChl a extinction vibrational modes'])
plt.ylabel('BChl a extinction')
#plt.plot(x,y_fit-y_exp)
plt.show()



# plt.figure()
# plt.plot(wavelength_lst,amp_lst)
# plt.plot(wavelength_lst,25.*np.exp((wavelength_lst-750.)/20),'-o')
# plt.show()


wavenumber_lst = 10**7/wavelength_lst-10**7/wavelength_lst[0]
_debug('print wavenumber_lst')


_debug('print len(wavenumber_lst)')


Hamiltonian = np.asarray(np.loadtxt('FMO.txt'))

Hamiltonian_new = np.zeros((7*len(wavenumber_lst), 7*len(wavenumber_lst)))

for k_col in range(0,7):
    for k_row in range(k_col,7):
        for j in range(len(wavenumber_lst)):
            for i in range(len(wavenumber_lst)):
                Hamiltonian_new[k_col*len(wavenumber_lst)+i,k_row*len(wavenumber_lst)+j] = round(Hamiltonian[k_col,k_row],0)
for k_col in range(0,7):
    for k_row in range(k_col,7):
        for j in range(len(wavenumber_lst)):
            for i in range(len(wavenumber_lst)):
                Hamiltonian_new[k_row*len(wavenumber_lst)+i,k_col*len(wavenumber_lst)+j] = round(Hamiltonian[k_col,k_row],0) 
for k in range(7):
    for i in range(len(wavenumber_lst)):
        for j in range(len(wavenumber_lst)):
            Hamiltonian_new[len(wavenumber_lst)*k+i,len(wavenumber_lst)*k+j] = 0
        
for i in range(7):
    for j in range(len(wavenumber_lst)):
        Hamiltonian_new[len(wavenumber_lst)*i+j,len(wavenumber_lst)*i+j] = math.floor(Hamiltonian[i,i] + wavenumber_lst[j])
      
print 'Saving Hamiltonian in a file new_hamiltonian.csv'      
np.savetxt('new_hamiltonian.csv',Hamiltonian_new,delimiter=',')

               


dipole_moment_site  = np.asarray(amp_lst)/sum(amp_lst)
dipole_moment = np.zeros((len(wavenumber_lst)*7,3))
#print dipole_moment.shape
#print dipole_moment_site
#print 'sum of all amplidutes = ', sum(amp_lst)



#print sum(dipole_moment_site)
#print dipole_moment_site
for k in range(7):
    for i in range(len(wavenumber_lst)):
        #print ((np.asarray(amp_lst)[i]/sum(amp_lst)))
        dipole_moment[k*len(wavenumber_lst)+i,:] = (np.asarray(amp_lst)[i]/sum(amp_lst))*muhat[k,:]
#print muhat[0,:]
#print 'amplitudes', (np.asarray(amp_lst)[:]/sum(amp_lst))
#print len(wavenumber_lst)
#print 'sum dipole moments=', sum(dipole_moment)
#print 'sum mu hat=', sum(muhat)
#print 'mu hat [0]  -0.7410 -0.5606 -0.3696', muhat[0,:]
#print 'sum dipole moment 0-5 -0.7410 -0.5606 -0.3696' , sum(dipole_moment[0:6,:])
#print 'sum dipole moment 6-11 -0.8571 0.5038 -0.1073 ' , sum(dipole_moment[6:11,:])
#print 'sum dipole moment 6-11 -0.197121 0.957410 -0.210972 ' , sum(dipole_moment[12:18,:])
#print dipole_moment

temp_check = np.zeros(3)
for i in range(5):
    temp_check[:] = temp_check[:] + dipole_moment[i,:]
#print temp_check/np.linalg.norm(temp_check)
#print 'dipole_moment_sum (-0.7410, -0.5606, -0.3696)', sum(dipole_moment[0:5,:]) 



eigenenergy,eigenvector = scipy.linalg.eigh(Hamiltonian_new)
_debug('print eigenenergy')
_debug('print eigenvector')



dipole_strength = np.zeros(len(wavenumber_lst)*7)
temp = np.zeros((len(wavenumber_lst)*7,3))
for i in range(len(wavenumber_lst)*7):
    for j in range(len(wavenumber_lst)*7):
        temp[i,:] = eigenvector[i,j]*dipole_moment[j,:]
    dipole_strength[i] = temp[i,0]**2 + temp[i,1]**2 + temp[i,2]**2


_debug('print dipole_strength, dipole_strength.shape')
#dipole_strength[5] = 0
#dipole_strength[11] = 0
_debug('print temp')



y_theory = x*0
for i in range(len(wavenumber_lst)*7):
    y_theory = y_theory+ gaussian_func(x,dipole_strength[i],10**7/eigenenergy[i],10)

plt.figure()
plt.plot(10**7/eigenenergy,dipole_strength/np.max(dipole_strength),'o')
plt.plot(x,y_theory/np.max(y_theory))
plt.xlim((730, 850))
plt.ylim((-0.05, 1.05))
plt.show()


#print eigenvector[1,:]
#print 'sum', sum(eigenvector[1,:])
#does eigenvecotr needs to sum up to 1?



plt.figure()
plt.plot(eigenvector[:,5],'-o')
plt.show()





