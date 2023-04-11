import numpy as np
import matplotlib.pyplot as plt
from scipy import interpolate
import os
import shutil

#############################################################################################################################
############################################# start input data  #############################################################
#############################################################################################################################

t=5      # you just need to specify in what time you want your results to be plotted.

############################################################################################################################
############################################# end input data ###############################################################
############################################################################################################################


t_name = str(t)

directory = 't=' + t_name	
  
# Parent Directory path
parent_dir = "."

path = directory
isExist = os.path.exists(path)

if isExist :
	shutil.rmtree(path, ignore_errors=False, onerror=None)
  
path = os.path.join(parent_dir, directory)
  
# Create the directory
os.mkdir(path)


contourNumber_of_levels = 31
colorbarNumber_of_levels = 11
        
x, y, u, v = np.loadtxt('../solver/results/t=' + t_name + '/Streamlines.txt', delimiter=', ', unpack=True)

xi = np.linspace(x.min(), x.max(), 100)                                
yi = np.linspace(y.min(), y.max(), 100) 

Xi, Yi = np.meshgrid(xi,yi)
uC = interpolate.Rbf(x, y, u, function='cubic')
vC = interpolate.Rbf(x, y, v, function='cubic')

uCi = uC(Xi, Yi)
vCi = vC(Xi, Yi)

fig,ax=plt.subplots(1,1)          

ax.streamplot(xi, yi, uCi, vCi, density=4, linewidth=1, color='k', arrowsize=0.5)            

ax.set_title('Streamline')
ax.set_xlabel('x (m)')
ax.set_ylabel('y (m)')
ax.axis(xmin=x.min(),xmax=x.max())
ax.axis(ymin=y.min(),ymax=y.max())

fig.savefig(directory + '/Streamlines.png')
	
#####################################################################################################
#####################################################################################################

x, y, z = np.loadtxt('../solver/results/t=' + t_name + '/P.txt', delimiter=', ', unpack=True)
	
fig,ax=plt.subplots(1,1)

ax.tricontour(x, y, z)

levels1 = np.linspace(z.min(), z.max(), contourNumber_of_levels)
levels2 = np.linspace(z.min(), z.max(), colorbarNumber_of_levels)
cp = ax.tricontourf(x, y, z, levels=levels1)

fig.colorbar(cp, ticks=levels2) # Add a colorbar to a plot

ax.set_title('Pressure')
ax.set_xlabel('x (m)')
ax.set_ylabel('y (m)')

plt.savefig(directory + '/Pressure.png')

#####################################################################################################
#####################################################################################################

x, y, z = np.loadtxt('../solver/results/t=' + t_name + '/U.txt', delimiter=', ', unpack=True)
	
fig,ax=plt.subplots(1,1)

ax.tricontour(x, y, z)

levels1 = np.linspace(z.min(), z.max(), contourNumber_of_levels)
levels2 = np.linspace(z.min(), z.max(), colorbarNumber_of_levels)
cp = ax.tricontourf(x, y, z, levels=levels1)

fig.colorbar(cp, ticks=levels2) # Add a colorbar to a plot

ax.set_title('u_Velocity')
ax.set_xlabel('x (m)')
ax.set_ylabel('y (m)')

plt.savefig(directory + '/u_Velocity.png')

#####################################################################################################
#####################################################################################################

x, y, z = np.loadtxt('../solver/results/t=' + t_name + '/V.txt', delimiter=', ', unpack=True)
	
fig,ax=plt.subplots(1,1)

ax.tricontour(x, y, z)

levels1 = np.linspace(z.min(), z.max(), contourNumber_of_levels)
levels2 = np.linspace(z.min(), z.max(), colorbarNumber_of_levels)
cp = ax.tricontourf(x, y, z, levels=levels1)

fig.colorbar(cp, ticks=levels2) # Add a colorbar to a plot

ax.set_title('v_Velocity')
ax.set_xlabel('x (m)')
ax.set_ylabel('y (m)')

plt.savefig(directory + '/v_Velocity.png')


#####################################################################################################
#####################################################################################################

path = '../solver/results/t=' + t_name + '/Mu_turbulence.txt'
isExist = os.path.exists(path)

if isExist :

	x, y, z = np.loadtxt('../solver/results/t=' + t_name + '/Mu_turbulence.txt', delimiter=', ', unpack=True)
		
	fig,ax=plt.subplots(1,1)

	ax.tricontour(x, y, z)

	levels1 = np.linspace(z.min(), z.max(), contourNumber_of_levels)
	levels2 = np.linspace(z.min(), z.max(), colorbarNumber_of_levels)
	cp = ax.tricontourf(x, y, z, levels=levels1)

	fig.colorbar(cp, ticks=levels2) # Add a colorbar to a plot

	ax.set_title('TurbulentViscosity (mu)')
	ax.set_xlabel('x (m)')
	ax.set_ylabel('y (m)')

	plt.savefig(directory + '/Mu_turbulence.png')

	
#####################################################################################################
#####################################################################################################

x, y = np.loadtxt('../solver/results/t=' + t_name + '/ShearStress_bottomWall.txt', delimiter=', ', unpack=True)

fig, ax = plt.subplots(1,1,figsize=(10,10))

ax.plot(x, y)
ax.set_xlabel('Length (m)')
ax.set_ylabel('Wall Shear Stress (Kg/m.s^2)')
ax.set_title('Bottom Wall')

plt.savefig(directory + '/ShearStress_bottomWall.png')

#####################################################################################################
#####################################################################################################

x, y = np.loadtxt('../solver/results/t=' + t_name + '/ShearStress_leftWall.txt', delimiter=', ', unpack=True)

fig, ax = plt.subplots(1,1,figsize=(10,10))

ax.plot(x, y)
ax.set_xlabel('Length (m)')
ax.set_ylabel('Wall Shear Stress (Kg/m.s^2)')
ax.set_title('Left Wall')

plt.savefig(directory + '/ShearStress_leftWall.png')

#####################################################################################################
#####################################################################################################

x, y = np.loadtxt('../solver/results/t=' + t_name + '/ShearStress_rightWall.txt', delimiter=', ', unpack=True)

fig, ax = plt.subplots(1,1,figsize=(10,10))

ax.plot(x, y)
ax.set_xlabel('Length (m)')
ax.set_ylabel('Wall Shear Stress (Kg/m.s^2)')
ax.set_title('Right Wall')

plt.savefig(directory + '/ShearStress_rightWall.png')

#####################################################################################################
#####################################################################################################

path = '../solver/results/t=' + t_name + '/yPlus_bottomWall.txt'
isExist = os.path.exists(path)

if isExist :

	x, y = np.loadtxt('../solver/results/t=' + t_name + '/yPlus_bottomWall.txt', delimiter=', ', unpack=True)

	fig, ax = plt.subplots(1,1,figsize=(10,10))

	ax.plot(x, y)
	ax.set_xlabel('Length (m)')
	ax.set_ylabel('Y Plus')
	ax.set_title('Bottom Wall')

	plt.savefig(directory + '/yPlus_bottomWall.png')
	
#####################################################################################################
#####################################################################################################

path = '../solver/results/t=' + t_name + '/yPlus_leftWall.txt'
isExist = os.path.exists(path)

if isExist :

	x, y = np.loadtxt('../solver/results/t=' + t_name + '/yPlus_leftWall.txt', delimiter=', ', unpack=True)

	fig, ax = plt.subplots(1,1,figsize=(10,10))

	ax.plot(x, y)
	ax.set_xlabel('Length (m)')
	ax.set_ylabel('Y Plus')
	ax.set_title('Left Wall')

	plt.savefig(directory + '/yPlus_leftWall.png')
	
#####################################################################################################
#####################################################################################################

path = '../solver/results/t=' + t_name + '/yPlus_rightWall.txt'
isExist = os.path.exists(path)

if isExist :

	x, y = np.loadtxt('../solver/results/t=' + t_name + '/yPlus_rightWall.txt', delimiter=', ', unpack=True)

	fig, ax = plt.subplots(1,1,figsize=(10,10))

	ax.plot(x, y)
	ax.set_xlabel('Length (m)')
	ax.set_ylabel('Y Plus')
	ax.set_title('Right Wall')

	plt.savefig(directory + '/yPlus_rightWall.png')
	
	

