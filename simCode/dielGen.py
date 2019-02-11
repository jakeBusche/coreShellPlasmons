##############################################################################
'''
This program creates a folder full of dielectric data files for use in DDA 
simulations. The dielectric data files describe a Drude-model metal with a
radially variable carrier density (plasma frequency) designed to mimic the 
carrier densities observed in ITO nanoparticles. 
'''
##############################################################################

##############################################################################
# importing packages
##############################################################################

import numpy as np
import matplotlib.pyplot as plt
import cmath as cm
import pprint as pp
import os
import shutil
import glob
import sys

##############################################################################
# set matplotlib.rc parameters to plot using LaTeX font
##############################################################################

plt.rc('text', usetex=True)
plt.rc('font', family = 'serif', size = 18,
	serif = ['computer modern roman'])

#############################################################################
#declaring universal constants
#############################################################################

c = 2.998*10**10. # speed of light (cm/s)
q = 4.8032042510*10**(-10.) # electron charge (statC)
nm = 10**(-7.) # (nm/cm)
hbar = 1.055*10**(-27) # reduced Planck constant (cm^2*g/s)
hbareV = 6.58211951*10**(-16) # reduced Planck constant (eV*s)
me = 9.10938*10**(-28) # mass of electron (g)
meV = 510998.946 # mass of electron (eV/c^2)
hceVnm = 1239.841984 # Planck constant times speed of light (eV*nm)

##############################################################################
# declaring system parameters
##############################################################################

hgamITO = 0.14 #Drude damping rate of ITO (eV)
epsbITO = 4 #Background dielectric constant of ITO (utls.)

##############################################################################
# defining useful functions
##############################################################################

'''
sig :   A sigmoid function built from an inverse tangent that has a defined
        minimum, maximum, and width. The distribution is minimum at x = -inf
        for sgn(dsig) > 0 and is minimum at x = inf when sgn(dsig) < 0.

        The arguments are:

        x       :    (scalar) The independent variable.
        x0      :    (scalar) The x-axis midpoint of the sigmoid distribution.
        dx      :    (scalar) Half of the distance between the x-axis points 
                        at which the distribution has values -pi/4 and pi/4. 
                        Defined such that sig(x = x0 + dx, args) = pi/4 * 
                        sgn(dsig), and sig(x = x0 - dx, args) = -pi/4 * 
                        sgn(dsig).
        dsig    :    (scalar) The y-axis distance between the minimum and 
                        maximum values of sig(dsig, args).
        sig0    :    (scalar) The minimum of the sigmoidal distribution. The
                        maximum is then sig0 + dsig.

'''
def sig(x,x0,dx,dsig,sig0):
	func = ( ( np.arctan( (x - x0)/dx * np.sign(dsig) ) + np.pi/2 ) * 
	    np.abs(dsig) / np.pi + sig0)
	return func


'''
X   :   An inverse sigmoid function built from a tangent that returns the x-
        axis value associated with a y-axis value of the distribution. 

        The arguments are:

        sig     :    (scalar) The independent variable.
        x0      :    (scalar) The x-axis midpoint of the sigmoid distributiom.
        dx      :    (scalar) Half of the distance between the x-axis points 
                        at which the distribution has values -pi/4 and pi/4. 
                        Defined such that sig(x = x0 + dx, args) = pi/4 * 
                        sgn(dsig), and sig(x = x0 - dx, args) = -pi/4 * 
                        sgn(dsig).
        dsig    :    (scalar) The y-axis distance between the minimum and 
                        maximum values of sig(dsig, args).
        sig0    :    (scalar) The minimum of the sigmoidal distribution. The
                        maximum is then sig0 + dsig.

'''
def X(sig,x0,dx,dsig,sig0):
	func = ( ( np.tan( np.pi/np.abs(dsig)*(sig - sig0) - np.pi/2 ) * 
		dx/np.sign(dsig) ) + x0 )
	if np.abs(func) > x0 + 3*dx:
		print('X : Desired value of x may be far from distribution midpoint.')
	else:
		pass
	return func


'''
lor :   An amplitude-normalized Lorentzian function with a FWHM of 2*gam, a
        minimum at lor0, and a maximum at lor0 + dlor.

        The arguments are:

        x       :   (scalar) The independent variable.
        x0      :   (scalar) The center of the distribution. The x-value of 
                        the maximum point.
        gam     :   (scalar) The half-width-at-half-maximum of the
                        distribution.
        dlor    :   (scalar) The height (peak value minus minimum value) of
                        the distribution.
        lor0    :    (scalar) The minimum value of the distribution.

'''
def lor(x,x0,gam,dlor,lor0):
	func = dlor * gam**2/( (x - x0)**2 + gam**2 ) + lor0
	return func


'''
sqr :   A quadratic function with a defined curvature and no linear 
        component. Defined such that the roots are x0 +/- sqrt(-y0/k). 

        The arguments are:

        x   :   (scalar) The independent variable.
        x0  :   (scalar) The (double) root when y0 = 0. The minimum of the 
                    function when k > 0, and the maximum when k < 0.

        k   :   (scalar) The curvature.
        y0  :   (scalar) The y-intercept when x0 = 0.

'''
def sqr(x,x0,k,y0):
	func = y0 + k*(x - x0)**2
	return func


'''
lin :   A linear function with a defined slope. Defined such that the root is
        x0 - y0/m. 

        The arguments are:

        x   :   (scalar) The independent variable.
        x0  :   (scalar) The root when y0 = 0. 
        m   :   (scalar) The slope.
        y0  :   (scalar) The y-intercept when x0 = 0.

'''
def lin(x,x0,m,y0):
	func = y0 + m*(x - x0)
	return func


'''
epsDrude    :   A (complex) dielectric function that models the response of a 
                noninteracting electron gas to electric excitations.

                The arguments are:

                hw          :   (scalar) The independent variable.
                hwp         :   (scalar) The plasma energy frequency of the 
		                            gas.
                hgam        :   (scalar) The plasma damping energy of the gas.
                epsInf      :   (scalar) The constant high-energy 
	                                contribution to the dielectric function.
'''
def epsDrude(hw,hwp,hgam,epsInf):
	func = epsInf  - hwp**2/(hw**2 + 1j * hgam * hw)
	return func

##############################################################################
# Initialize shape parameters and relevant arrays
##############################################################################

# The minimum and maximum plasma energies to be used in the calculation
hwpMin = 0.5 
hwpMax = 2.0

# The radius of the sphere
R = 80

# The number of shells of distinct material in the sphere
nMat = 14

# The radii of the boundary surfaces between shells
aRad = np.linspace(0,R,nMat + 1)
aRad = aRad[1::]

# The radii of the midpoints (along the radial direction) of each shell
aRadMid = np.zeros( len(aRad) )
for i in range(len(aRadMid)):
	if i == 0:
		aRadMid[i] = aRad[i] / 2
	else:
		aRadMid[i] = (aRad[i] + aRad[i-1])/2

# The energy values at which the dielectric function will be defined
aEn = np.arange(0.1,6.0 + 0.1,0.1)

''' 
Arrays containing DISCRETE plasma energy values that fall along a sigmoidal
    distribution, a Lorentzian distribution, and a linear distribution. Each
    distribution begins at hwpMin at r = 0 and ends at hwpMax at r = R. The 
    sigmoid is centered about r = R/2.
'''
aHwpSig = np.zeros( len(aRad) )
aHwpLor = np.zeros( len(aRad) )
aHwpLin = np.zeros( len(aRad) )

for i in range(len(aRadMid)):
	aHwpSig[i] = sig(aRadMid[i],aRad[len(aRad)//2 - 1],2,
		hwpMax - hwpMin,hwpMin)
	aHwpLor[i] = lor(aRadMid[i],aRad[len(aRad)//2 - 1],R/8,(hwpMax - hwpMin),
		hwpMin)
	aHwpLin[i] = lin(aRadMid[i],0,(hwpMax - hwpMin)/R,hwpMin)

'''
Arrays containing QUASI-CONTINOUS energy values that all along a sigmoidal
    distribution, a Lorentzian distribution, and a linear distribution. Each
    distribution begins at hwpMin at r = 0 and ends at hwpMax at r = R. The 
    sigmoid is centered about r = R/2.
''' 
aR = np.arange(0,R,0.01)
aSig = np.zeros( len(aR) )
aLor = np.zeros( len(aR) )
aLin = np.zeros( len(aR) )
for i in range( len(aR) ):
	aSig[i] = sig(aR[i],aRad[len(aRad)//2 - 1],2,hwpMax - hwpMin,hwpMin)
	aLor[i] = lor(aR[i],aRad[len(aRad)//2 - 1],R/8,(hwpMax - hwpMin),
		hwpMin)
	aLin[i] = lin(aR[i],0,(hwpMax - hwpMin)/R,hwpMin)

##############################################################################
# plotting 
##############################################################################

# Enable interactive plotting to avoid plots blocking code interpretation
plt.ion()

# Plot the discrete values along the continuous "ideal" distributions
fig1 = plt.figure(1)
ax1 = fig1.add_subplot(1,1,1)

# Demarcate the boundary radii between shells
for i in range( len(aRad) ):
	ax1.axvline(x = aRad[i],lw = 1,ls = '-',color = 'LightGrey')

ax1.plot(aR,aSig,lw = 2,color = 'k',label = 'Sigmoid')
ax1.plot(aR,aLor,lw = 2,color = 'b',label = 'Lorentzian')
ax1.plot(aR,aLin,lw = 2,color = 'r',label = 'Linear')

ax1.scatter(aRadMid,aHwpSig,c = 'k',edgecolor = 'k',s = 20)
ax1.scatter(aRadMid,aHwpLor,c = 'b',edgecolor = 'b',s = 20)
ax1.scatter(aRadMid,aHwpLin,c = 'r',edgecolor = 'r',s = 20)

ax1.legend(loc = 'best',prop = {'size':12})

ax1.set_xlabel('Radius [nm]')
ax1.set_ylabel('Plasma Energy [eV]')

ax1.set_xlim([0.0,80.0])
ax1.set_ylim([0.49,2.01])

ax1.set_title('Ideal Carrier Distribution')

ax1.minorticks_on()
fig1.tight_layout()


# Plot the discrete values along the step-wise "discretized" distributions
fig2 = plt.figure(2)
ax2 = fig2.add_subplot(1,1,1)

# Demarcate the boundary radii between shells
for i in range( len(aRad) ):
	ax2.axvline(x = aRad[i],lw = 1,ls = '-',color = 'LightGrey')

ax2.step(aRadMid,aHwpSig,lw = 2,color = 'k',where = 'mid',
	label = 'Sigmoid')
ax2.step(aRadMid,aHwpLor,lw = 2,color = 'b',where = 'mid',
	label = 'Lorentzian')
ax2.step(aRadMid,aHwpLin,lw = 2,color = 'r',where = 'mid',
	label = 'Linear')

ax2.scatter(aRadMid,aHwpSig,c = 'k',edgecolor = 'k',s = 20)
ax2.scatter(aRadMid,aHwpLor,c = 'b',edgecolor = 'b',s = 20)
ax2.scatter(aRadMid,aHwpLin,c = 'r',edgecolor = 'r',s = 20)

ax2.legend(loc = 'best',prop = {'size' : 12})

ax2.set_xlabel('Radius [nm]')
ax2.set_ylabel('Plasma Energy [eV]')

ax2.set_xlim([0.0,80.0])
ax2.set_ylim([0.49,2.01])

ax2.set_title('Discretized Carrier Distribution')

ax2.minorticks_on()
fig2.tight_layout()


# Save figures for later presentations
pwd = os.getcwd()
pArg = input('Save plots to in directory "' + pwd + '"? (y/n) ')
if pArg is 'y':
	fig1.savefig(pwd + '/idealDist.png',dpi = 300)
	fig2.savefig(pwd + '/discDist.png',dpi = 300)
	print('Figures successfully saved!')
else:
	dArg = input('Enter name of desired directory' + 
		' (enter "n" to forgo saving): ')
	if dArg is 'n':
		print('Figures were not saved.')
	else:
		fig1.savefig(dArg + '/idealDist.png',dpi = 300)
		fig2.savefig(dArg + '/discDist.png',dpi = 300)
		print('Figures saved successfully!')


plt.show()

##############################################################################
# printing data to files
##############################################################################

# Define paths to directories in which dielectric data will be saved
pwd = os.getcwd()
dataPathSig = '/mnt/c/DDA/pwExperiments/coreShellPlasmons/coreShellDielsSig'
dataPathLor = '/mnt/c/DDA/pwExperiments/coreShellPlasmons/coreShellDielsLor'
dataPathLin = '/mnt/c/DDA/pwExperiments/coreShellPlasmons/coreShellDielsLin'

# Define lists of directory paths and plasma energy arrays
pathList = [dataPathSig,dataPathLor,dataPathLin]
aHwpList = [aHwpSig,aHwpLor,aHwpLin]

# Give the user a heads up
print('Generating dielectric data in directories',pathList)

'''
For each directory in pathList:

1.  Create directory if it does not exist. Delete it and recreate it if it
        does exist.
2.  Create arrays of dielectric function values. Each array will describe a
        dielectric function characterized by a plasma energy contained in 
        aHwp"Sig,Lor,Lin,..."
3.  Save each array in the proper directory as a text file with proper DDA 
        dielectric file formatting and a name given with formatting required
        by ddscatInputGen.py
'''
for dataPath in pathList:
	print('---')
	print('Current directory is : ',dataPath)
	print('---')

	if os.path.isdir(dataPath):
		pArg = input('Dielectric folder ' + dataPath + ' already exists.'+ 
			' Delete and recreate? (y/n) ')
		if pArg is 'y':
			shutil.rmtree(dataPath)
			os.mkdir(dataPath)
			print('Folder deleted and recreated.')
		else:
			print('Folder was not deleted.')
			print('Exiting program.')
			sys.exit()
	else:
		print('Creating dielectric folder.')
		os.mkdir(dataPath)

	pArg2 = input('Write dielectric data? (y/n) ')
	if pArg2 is 'y':
		print('Writing data...')
		
		hwpArray = aHwpList[ pathList.index(dataPath) ]

		for i in range( len(hwpArray) ):
			aEps = np.zeros( (len(aEn) ),dtype = complex)
			for j in range(len(aEn)):
				aEps[j] = epsDrude(aEn[j],hwpArray[i],hgamITO,epsbITO)

			epsRealArray = np.real(aEps)
			epsImagArray = np.imag(aEps)

			wlArray = 1239.841984 / aEn / 1000

			outMat = np.vstack( (wlArray,epsRealArray) )
			outMat = np.vstack( (outMat,epsImagArray) )
			outMat = np.transpose(outMat)
			outMat = outMat[::-1]

			print(dataPath + '/diel' + str(i+1) + '.tab' + ' written.')

			hdr1 = ('Real and imaginary epsilon of ITO with hwp = ' +
				str(hwpArray[i]))
			hdr2 = (str(1) + ' ' + str(0) + ' ' + str(0) + ' ' + 
				str(2) + ' ' + str(3) + ' = specifies whether n,k or ' + 
				'epsilon' + ' are read in (see manual)')
			hdr3 = (3*'%12.12s') % ('lambda (um)','real(eps)','imag(eps)')
			hdr = '\n'.join([hdr1,hdr2,hdr3])

			with open(dataPath + '/diel' + str(i+1) + '.tab','w+') as file:
				np.savetxt(dataPath + '/diel' + str(i+1) + '.tab',outMat,
					comments = '',header = hdr,fmt = '%12.12s')
				file.close()
	else:
		print('Data was not written.')
		print('Exiting program.')
		sys.exit()


##############################################################################
# testing
##############################################################################

