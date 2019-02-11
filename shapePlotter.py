##############################################################################
# importing packages
##############################################################################

import numpy as np
import os as os
import sys as sys
import shutil as shutil
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import matplotlib.cm as cmap

##############################################################################
# import shape file
##############################################################################

pwd = os.getcwd()
fSh = 'radInhomSph80_1.dat'

print('Plotting file : ' + fSh)
pArg = input('Continue? (y/n) ')
if pArg is 'y':
	pass
else:
	print('Exiting program.')
	sys.exit()

aData = []
with open(pwd + '/' + fSh,'r') as shFile:
	shLines = shFile.readlines()[7:]
	for shLine in shLines:
		shLine = shLine.split()
		shLine = [int(i) for i in shLine]
		aData.append(shLine)
aData = np.array(aData,dtype = int)

##############################################################################
# plot shape file data
##############################################################################

coords = aData[:,1:4]
# print(np.amax(coords[:,0]),np.amax(coords[:,1]),np.amax(coords[:,2]))
mats = aData[:,5]
colors = cmap.rainbow(np.linspace(0,1,len(np.unique(mats))))

print('Number of points in file is :',str(len(coords[:,0])))
pArg1 = input('Thin data? (y/n) ')
if pArg1 is 'y':
	pArg2 = input('By what factor? ')
	print('Displaying thinned data set')
	coords = coords[::int(pArg2),:]
	mats = aData[::int(pArg2),5]
	colors = cmap.rainbow(np.linspace(0,1,len(np.unique(mats))))
else:
	print('Displaying complete data set.')

colorList = []
for i in range( len(mats) ):
	colorList.append(colors[mats[i]-1])
# print(colorList)

aXidx = coords[:,0]
aYidx = coords[:,1]
aZidx = coords[:,2]



pArg2 = input('Print 3d shape? (y/n) ')
if pArg2 is 'y':
	fig1 = plt.figure(1)
	ax1 = fig1.add_subplot(111, projection = '3d')
	ax1.set_aspect('equal','box')

	ax1.scatter(aXidx,aYidx,aZidx,s=20,c = colorList,edgecolors = colorList)

	ax1.set_xlabel('x index')
	ax1.set_ylabel('y index')
	ax1.set_zlabel('z index')
else:
	pass




aYidxSl1 = np.array([])
aZidxSl1 = np.array([])
matsSl1 = np.array([],dtype = int)
for i in range(len(aXidx)):
	if aXidx[i] == 0.:
		aYidxSl1 = np.append(aYidxSl1,aYidx[i])
		aZidxSl1 = np.append(aZidxSl1,aZidx[i])
		matsSl1 = np.append(matsSl1,mats[i])
	else:
		pass

colorListSl1 = []
for i in range( len(matsSl1) ):
	colorListSl1.append(colors[matsSl1[i]-1])




fig2 = plt.figure(2)
ax2 = fig2.add_subplot(1,1,1)

cax2 = ax2.imshow( matsSl1, interpolation = None, cmap = cmap.rainbow )
cbar2 = fig2.colorbar( cax2, ticks = matsSl1.tolist() )

ax2.scatter(aYidxSl1,aZidxSl1,s = 20,c = colorListSl1,
	edgecolors = colorListSl1)

ax2.set_xlabel('y (nm)')
ax2.set_ylabel('z (nm)')

ax2.minorticks_on()
ax2.axis('equal')
fig2.tight_layout()





plt.show()