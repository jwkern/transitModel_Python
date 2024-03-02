import matplotlib.pyplot as plt
import math
import numpy as np
import pylab

Ls = 3.846E26  #luminosity of the star
rp = 1.35*69E6     #radius of the planet in Jupiter radii
rs = 1.14*696E6    #radius of the star in solar radii
days = 3.52474541 # 5.579  #period of orbit in days
P = days*86400.0    #period of orbit in seconds
a = 0.045*1.496E11	#semi-major axis
w = (2*math.pi)/P     #orbital frequency
t_steps = 1146		#number of time steps
s_steps = 100
L = np.zeros(t_steps)    #Luminosity array
x = np.zeros(t_steps)     #x array

time_transit_days = 2451370.048
time_transit = time_transit_days*86400.0

n = (86.1/180)*math.pi   #incliniation
b = (a*math.cos(n))/rs     #impact parameter
D = 3.784E16     	#distance to the star
time, flux = np.loadtxt('100binned.txt', unpack=True)
# t = np.arange (0, P, P/t_steps) #time array
# convert times to seconds, apparently.
time=time*86400.0
t=time
flux=1.0-flux
chi2 = np.zeros(t_steps)

length = 2*rs + 4*rp     #length of the box used to check whether xp is inside or outside of the star/planet
area = np.zeros(t_steps) #array of area

for i in range (0, t_steps):
	yp = b*rs
	xp = a*math.cos(w*(t[i]-time_transit))   #position of the planet at time t[i]
	zp = a*math.sin(w*(t[i]-time_transit))   #position of the planet at time t[i]
	x[i]=xp
	dist = (abs((xp**2)+(yp**2)))**(0.5)
	if (dist >= (rp + rs)) or (zp< 0):
		L[i] = 1.0 #Ls
	else:
		theta = math.asin((((xp**2)+(yp**2))**(0.5))/(rs+(2*rp)))
		countstar = 0        #counts inside the star
		countplanet = 0       #counts inside the planet
		for j in range (0, (s_steps)-1):
			for k in range (0, (s_steps)-1):
				xtestj = (j*length) / s_steps - length/2.0
				ytestk = (k*length) / s_steps - length/2.0
				ls = ((xtestj**2) + (ytestk**2))**(0.5)
				if ls < rs:
					countstar = countstar + 1
					lp = (((xtestj - xp)**2) + ((ytestk - yp)**2))**(0.5)
					if lp < rp:
						countplanet = countplanet + 1
		area[i] = 1.0*countplanet/ countstar
		L[i] = 1-area[i] # (math.cos(theta)*area[i]) #*Ls
chi2 =np.sum(((flux-L)**2)/L)
print(chi2)

plt.scatter(time,flux, linewidth=2.0)
plt.plot (time, L, label='HD 209458')
plt.legend (loc='upper right')
plt.xlabel('t')
plt.ylabel('L')
plt.title('Exoplanet Transit')
#plt.xlim (120000, 220000)
#plt.ylim (-0.01, 0.021)
plt.show()
