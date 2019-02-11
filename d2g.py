##############################################################################
# importing packages
##############################################################################

import numpy as np
import numpy.linalg as linalg
import matplotlib.pyplot as plt
import math as m
import cmath as cm
import pprint as pp
import matplotlib.cm as colormap
from scipy import interpolate

##############################################################################
# set matplotlib.rc parameters to plot using LaTeX font
##############################################################################

plt.rc('text', usetex=True)
plt.rc('font', family = 'serif', size = 18,
	serif = ['computer modern roman'])

##############################################################################
# defining universal constants
##############################################################################

q = 1.6021766208*10**(-19.) #electron charge (C)
c = 2.998*10**8. #speed of light (m/s)
hbareV = 1.054571628251774 *10**(-34.)/q #reduced Planck constant (eV*s)

##############################################################################
# importing data
##############################################################################

dataPath = '/mnt/c/DDA/pwExperiments/coreShellPlasmons/Calculations/job5/'
dataName = 'specData.txt'

aData = np.genfromtxt(dataPath + dataName,delimiter = '',comments = '#')

aEn = aData[:,0]
aExt = aData[:,1]
aAbs = aData[:,2]
aSca = aData[:,3]

##############################################################################
# plotting data
##############################################################################

plt.ion()

fig1 = plt.figure(1)
ax1 = fig1.add_subplot(1,1,1)

ax1.plot(aEn,aAbs,color = 'r',lw = 2,label = 'Absorption')
ax1.plot(aEn,aSca,color = 'b',lw = 2,label = 'Scattering')
ax1.plot(aEn,aExt,color = 'k',lw = 2,ls = '--',label = 'Extinction')

ax1.legend(loc = 'best',prop = {'size' : 12})

ax1.set_xlabel('Energy (eV)')
ax1.set_ylabel(r'Cross-Section (nm$^2$)')
ax1.set_title('Optical cross-sections of a\n' + 
	'10 nm radius ITO of sigmoidally\n' + r' radially varying $\omega_p$')

ax1.minorticks_on()
fig1.tight_layout()

pArg = input('Save figure to file "fig.png"? (y/n) ')

if pArg is 'y':
	fig1.savefig('fig.png',dpi = 300)
	print('Figure was saved successfully!')
else:
	fArg = input('Please enter desired filename (enter "n" to cancel): ')
	if fArg is 'n':
		print('Figure was not saved.')
	else:
		fig1.savefig(fArg,dpi = 300)
		print('Figure was saved successfully!')

plt.show()
