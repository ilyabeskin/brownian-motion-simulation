#-------------------------------------------------------------------------------
# This code is based on the method described in
# 	"Simulation of a Brownian particle in an optical trap" 
# 	by Giorgio Volpe and Giovanni Volpe
#	
#	American Journal of Physics 81, 224 (2013); doi: 10.1119/1.4772632
#
#-------------------------------------------------------------------------------
#
# Code written by Ilya M. Beskin
# Updated: Jun 2, 2024
#
#-------------------------------------------------------------------------------

import matplotlib.pyplot as plt
import numpy as np
import random

# WhiteNoise returns a single random number with a Gaussian distribution.
# 	delT -> time step given in seconds (s)
def WhiteNoise(delT):
	return random.gauss(0,1)/(np.sqrt(delT))
	# Above line insures that the variance will be 1/(delta T)

# generateATraj returns a 1D numpy array of positions representing a spherical Brownian
# 		particle in a harmonic potential.
#	totTime -> total time of trajectory in seconds.
#	delT 	-> time interval between points on trajectory.
#	gamma	-> Stokes drag on particle.
#	m 		-> Mass of particle.
#	spK		-> Spring coefficient of particle.
#	showPlots -> Boolean variable. If set to true, MatPlotLib plots of generated trajectory will appear.
def generateATraj(totTime,delT, gamma, m, spK, showPlots):
	totSteps = int(np.ceil(totTime/delT))

	randNums = np.zeros(totSteps) # randNums will record the sequence of generated white noise
	brownPos = np.zeros(totSteps) # brownPos will record the position of your particle

									# Particle always starts for the first two time steps at 0.
	randNums[1] = WhiteNoise(delT)	# This line generates the first push on the second time step.

	kB = 1.38e-23	#Boltzmann constant (m^2 kg / s^2 K)
	T = 300			#Temp in Kelvin

	tauInv = gamma/m 					# tauInv is 1/tau where tau is the characteristic time constant
	#print(tau)							#		t>>tau -> diffusive motion
										#		t<<tau -> ballistic motion
	denom = 1+tauInv*delT+spK/m*(delT**2)	# This simplifies the calculation in the loop.

	for i in range(2,totSteps):
		brownPos[i] = ((2+tauInv*delT)/(denom))*brownPos[i-1]
		brownPos[i] -= (1/(denom))*brownPos[i-2]
		brownPos[i] +=(np.sqrt(2*kB*T*gamma)/(m*(denom)))*(delT**(2))*randNums[i-1]	# Update position

		randNums[i] = WhiteNoise(delT)		#Generate a new random push.

	if showPlots is True:
		fig, (ax1, ax2) = plt.subplots(2)
		ax1.plot(range(totSteps), randNums)
		ax2.plot(range(totSteps), brownPos)
		plt.show()

	return brownPos


# findMSD generates multiple trajectories and graphs/outputs the mean squared displacement at each timestep.
#	numRuns	-> Integer giving number of trajectories to generate.
#	runTime	-> Time of each trajectory in [s]
#	delT 	-> Time step between points on trajectory
#	gamma	-> Stokes drag on spherical particle.
#	m 		-> Mass of spherical particle.
#	spK		-> Spring coefficient.
def findMSD(numRuns, runTime, delT, gamma, m, spK):
	totSteps = int(runTime/delT)
	allTheTraj = np.zeros((numRuns,totSteps))

	for i in range(numRuns):
		allTheTraj[i] = generateATraj(runTime, delT, gamma, m, spK, showPlots=False)
		#print("Your particle has a standard deviation of" + str(np.std(allTheTraj[i])) + " m")
		plt.plot(range(totSteps),allTheTraj[i])
	plt.show()

	msd = np.zeros(totSteps)

	for j in range(totSteps):
		msd[j] = np.mean(allTheTraj[:,j]**2)
	plt.plot(range(totSteps),msd)
	plt.show()

	return msd

	#print(allTheTraj)

runTime = 1					# Seconds
delT = 1e-5					# Inverse of sampling rate [s]
eta = 0.001 				# viscocity (H2O) [Ns/m^2]
r = 100e-9					# radius [m]
gamma = 6*np.pi*eta*r 		# Stokes drag
m = 4/3*np.pi*(r**3)*1000	# mass in Kg (polystyrene microparticle)
print("Mass is "+ str(m))
spK = 15e-7					# Spring coefficient of trap in N/m
print("Your gamma is " + str(gamma))

#findMSD(100, runTime, delT, gamma, m, spK)
traj = generateATraj(runTime, delT, gamma, m, spK, showPlots=False)
plt.plot(np.arange(0,runTime,delT), traj*1e9, color='black')
plt.xlabel("Time (s)")
plt.ylabel("X Displacement (nm)")
plt.show()

#np.savetxt('brownSimData.csv',traj,delimiter=',')