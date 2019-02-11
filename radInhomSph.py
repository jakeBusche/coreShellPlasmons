##############################################################################
'''
This program creates a sphere enveloped by a number of shells, each with a
different isotropic dielectric function.
'''
##############################################################################

##############################################################################
# import packages
##############################################################################

import numpy as np
import matplotlib.pyplot as plt
import time
import os

##############################################################################
# define functions
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

def X(sig,x0,dx,dsig,sig0):
	func = ( ( np.tan( np.pi/np.abs(dsig)*(sig - sig0) - np.pi/2 ) * 
		dx/np.sign(dsig) ) + x0 )
	if np.abs(func) > x0 + 3*dx:
		print('X : Desired value of x may be far from distribution midpoint.')
	else:
		pass
	return func

##############################################################################
# define shape parameters
##############################################################################

# Time checkpoint 1
ckpt1 = time.time()

'''
Define sphere parameters:

	R 		:	radius of sphere in index units
	Didx 	:	integer index spacing of shape
	aXidx 	:	list of x indices to be considered
	aYidx 	: 	list of y indices to be considered
	aZidx 	: 	list of z indices to be considered
'''
R = 80
Didx = 1
aXidx = np.arange(-R, R+Didx, Didx,dtype = int)
aYidx = np.arange(-R, R+Didx, Didx,dtype = int)
aZidx = np.arange(-R, R+Didx, Didx,dtype = int)

##############################################################################
# build shape
##############################################################################


'''
Create three Nx by Ny by Nz mesh matrices such that:

	aXXidx	:	x index values vary by sheet
	aYYidx	:	y index values vary by row
	aZZidx	:	z index values vary by column
'''
aXXidx,aYYidx,aZZidx = np.meshgrid(aXidx,aYidx,aZidx,indexing = 'ij')


'''
Flatten the 3d mesh matrices and concatenate them into an array in which each
row is a location vector (Xidx,Yidx,Zidx). The Zidx cycles each Nz rows, the
Yidx cylces every Ny*Nz rows, and the Xidx cycles every Nx*Ny*Nz rows (once).
'''
coords = np.vstack((aXXidx.flatten(),aYYidx.flatten()))
coords = np.vstack((coords,aZZidx.flatten()))
coords = np.transpose(coords)


'''
Iterate through list of coordinate vectors and select those that are outside 
of the sphere's volume for deletion.

	killList 	:	List of row indices of coords to pass to np.delete() and
					delete extraneous points from coords.
'''
killList = []
for i in range(len(coords[:,0])):
	if coords[i,0]**2 + coords[i,1]**2 + coords[i,2]**2 > R**2:
		killList.append(i)
	else:
		pass

 #    # Uncomment to monitor shape creation progress
	if i % 1000 == 0:
		print('Shape construction',"{0:4.1f}".format(i/len(coords)*100),
			'percent complete.\r',end = '')

coords = np.delete(coords,killList,axis = 0)
nPts = len(coords[:,0])

# Time checkpoint 2
ckpt2 = time.time()


'''
Print out total number of points in shape and total time taken to demarcate 
the indices inside the shape.
'''
print('Total number of points in shape :',nPts)
print('Total time (in seconds) to demarcate shape :',ckpt2 - ckpt1)


'''
Concatenate array of dummy indices to coords for each point in shape. Also, 
concatenate array of material indices to coords for each point in shape.
'''
coords = np.transpose(
	np.vstack((
	np.arange(1,nPts + 1,1),
	np.transpose(coords)
	))
	)

# Build an array of ones to concatenate with aData and edit in place 
aMat = np.ones( (nPts,3),dtype = int)

aData = np.hstack((
	coords,
	aMat
	))

##############################################################################
# define material values
##############################################################################

# number of regions of distinct material in shape
nMat = 14

# make list of material boundary radii
aRad = np.linspace(0,R,nMat + 1)
aRad = aRad[1::]

# Uncomment to debug material assignment
print('List of boundary radii:',aRad)

for i in range(len(aData[:,0])):
	for j in range(len(aRad)):
		if j == 0:
			if aData[i,1]**2 + aData[i,2]**2 + aData[i,3]**2 <= aRad[j]**2:
			    aData[i,4],aData[i,5],aData[i,6] = j+1,j+1,j+1
			else:
				pass
		elif j > 0:
			if ( aData[i,1]**2 + aData[i,2]**2 + aData[i,3]**2 > 
	    		aRad[j - 1]**2 and 
	    		aData[i,1]**2 + aData[i,2]**2 + aData[i,3]**2 <= 
	    		aRad[j]**2 ):

			    aData[i,4],aData[i,5],aData[i,6] = j+1,j+1,j+1
			else:
				pass

		# Uncomment to debug material assignment
		# print('Latest assignment: material',str(j+1),'at radius',
		# 	"{0:4.1f}".format(np.sqrt(aData[i,1]**2 + aData[i,2]**2 +
		# 	 aData[i,3]**2)),'\r',end = '')

	# Uncomment to monitor shape creation progress
	if i % 1000 == 0:
		print('Material assignment',"{0:4.1f}".format(i/len(coords)*100),
			'percent complete.\r',end = '')

ckpt3 = time.time()
print('Total time (in seconds) to assign material indices :',ckpt3 - ckpt2)

dataList = aData.tolist()

##############################################################################
# write shape.dat file 
##############################################################################

'''
!!! All data printed to file must be printed as a string to be read with a 
    FORTRAN 90 program !!!

Printing shape file header:

	hdr1 	:	[text read as comments]
	hdr2 	:	nPts [text read as comments]
	hdr3 	:	a1_x a1_y a1_z [text read as comments]
	hdr4	:	a2_x a2_y a2_z [text read as comments]
	hdr5 	: 	d_x/d d_y/d d_z/d [text read as comments]
	hdr6 	:	x0 y0 z0 [text read as comments]
	hdr7	:	[text read as comments]
'''
hdr1 = ' '.join(['---Radially inhomogeneous sphere with index radius',
	str(R), 'and index spacing',str(Didx),'---'])
hdr2 = ' '.join([str(int(nPts)),' = number of dipoles in shape'])
hdr3 = ' '.join(['{:.3f}'.format(1.000),' ','{:.3f}'.format(0.000),' ',
	'{:.3f}'.format(0.000),' ','= target vector a1 (in TF)'])
hdr4 = ' '.join(['{:.3f}'.format(0.000),' ','{:.3f}'.format(1.000),' ',
	'{:.3f}'.format(0.000),' ','= target vector a2 (in TF)'])
hdr5 = ' '.join([str(1.000),' ',str(1.000),' ',str(1.000),' ', 
	'= d_x/d d_y/d d_z/d (normally 1. 1. 1.)'])
hdr6 = ' '.join([str(0.000),' ',str(0.000),' ',str(0.000),' ', 
	'= X0(1-3) = location in lattice of target origin'])
hdr7 = (7*'%-8.8s') % ('I', 'IX', 'IY', 'IZ', 'ICOMPX', 'ICOMPY', 'ICOMPZ')
# hdr = '\n'.join([hdr1,hdr2,hdr3,hdr4,hdr5,hdr6,hdr7])
hdr = '\n'.join([hdr1,hdr2,hdr3,hdr4,hdr5,hdr6,hdr7])

print('File header is :\n\n',hdr,'\n')
print('Shape data is :\n\n',aData,'\n')

dataName = 'radInhomSph' + str(int(R)) + '_' + str(int(Didx)) + '.dat'
pArg = input('Save data to file: ' + dataName + '? (y/n) ')

if pArg is 'y':
	np.savetxt(dataName,
		aData,delimiter = ' ',header = hdr,
		comments = '',fmt = '%-8.8s')
	print('Data saved succesfully!')
elif pArg is not 'y':
	dataName = input('Please enter desired filename for data' + 
		' (including extension): ')
	np.savetxt(dataName,
		aData,delimiter = ' ',header = hdr,
		comments = '',fmt = '%-8.8s')
	print('Data saved succesfully!')
else:
	print('Data was not saved.')

##############################################################################
# testing
##############################################################################
