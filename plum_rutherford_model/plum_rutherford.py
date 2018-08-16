"""
Variables that can be changed
plum_or_rutherford
n_alpha,x0, y0, vx0, vy0, theta0, n_alpha_screen
animation_speed, tmax, n_atom_high, n_atom_wide, n_atom_high_m, t_frames
"""
import numpy as np
import scipy.integrate as spi
import matplotlib.pyplot as plt
import matplotlib.lines as matlines
import matplotlib.animation as an
import numpy.random as rn

def electric(r,t):	
	"""
	Solving differential equations for electric fields for Rutherford's atoms
	The sum of every electric field created by each gold nucleus is the resultant electric field.
	When alpha is inside the nucleus, then use equation F=kqQr/Rn^3,. Elsewhere, F = kqQ/r^2 where k=1/4(pi)Eo. F/m = a
	Every atom away from the edge will be adjacent to 6 other atoms (if the size of the foil allows it).
	"""
	x=r[0]
	vx=r[1]
	y=r[2]
	vy=r[3] 
	axx=np.array([])
	ayy=np.array([])
	if  n_atom_high_parity==0:	#If there are even number of rows
		for i in range(n_atom_wide):
			for j in range(int(n_atom_high/2)):
				k = 2*j	#Used for every even row
				l = 2*j+1	#Used for every odd row
				if np.sqrt((x-i*gold_space)**2+(y-k*gold_space)**2)<Rn:	#If atom in nucleus
					ax_acc=const_in_nucleus*(x-i*gold_space)
					ay_acc=const_in_nucleus*(y-k*gold_space)
				else:
					ax_acc=(const*(x-i*gold_space))/((x-i*gold_space)**2.+(y-k*gold_space)**2.)**(3./2.)
					ay_acc=(const*(y-k*gold_space))/((x-i*gold_space)**2.+(y-k*gold_space)**2.)**(3./2.)	
				axx=np.append(axx,ax_acc)
				ayy=np.append(ayy,ay_acc)
				if i != range(n_atom_wide)[-1]:		#Every odd row has one less atom than the even rows. This is to allow the even rows to act as the edges of the foil.
					if np.sqrt((x-i*gold_space-0.5*gold_space)**2+(y-l*gold_space)**2)<Rn:	#If atom in nucleus
						ax_acc=const_in_nucleus*(x-i*gold_space-0.5*gold_space)
						ay_acc=const_in_nucleus*(y-l*gold_space)
					else:
						ax_acc=(const*(x-i*gold_space-0.5*gold_space))/((x-i*gold_space-0.5*gold_space)**2.+(y-l*gold_space)**2.)**(3./2.)
						ay_acc=(const*(y-l*gold_space))/((x-i*gold_space-0.5*gold_space)**2.+(y-l*gold_space)**2.)**(3./2.)
				axx=np.append(axx,ax_acc)
				ayy=np.append(ayy,ay_acc)
	elif n_atom_high == 1:	#If there is one row.
		for i in range(n_atom_wide):
			if np.sqrt((x-i*gold_space)**2+(y)**2)<Rn:	#If atom in nucleus
				ax_acc=const_in_nucleus*(x-i*gold_space)	
				ay_acc=const_in_nucleus*y
			else:
				ax_acc=(const*(x-i*gold_space))/((x-i*gold_space)**2.+(y)**2.)**(3./2.)
				ay_acc=(const*(y))/((x-gold_space)**2.+(y)**2.)**(3./2.)
			axx=np.append(axx,ax_acc)
			ayy=np.append(ayy,ay_acc)

	else:		#If there are odd number of rows thats not equal to 1.
		for i in range(n_atom_wide):
			for j in range(int(n_atom_high/2)):
				k = 2*j
				l = 2*j+1
				if np.sqrt((x-i*gold_space)**2+(y-k*gold_space)**2)<Rn:	#If atom in nucleus
					ax_acc=const_in_nucleus*(x-i*gold_space)
					ay_acc=const_in_nucleus*(y-k*gold_space)
				else:
					ax_acc=(const*(x-i*gold_space))/((x-i*gold_space)**2.+(y-k*gold_space)**2.)**(3./2.)
					ay_acc=(const*(y-k*gold_space))/((x-i*gold_space)**2.+(y-k*gold_space)**2.)**(3./2.)		
				axx=np.append(axx,ax_acc)
				ayy=np.append(ayy,ay_acc)
				if i != range(n_atom_wide)[-1]:		#Every odd row has one less atom than the even rows. This is to allow the even rows to act as the edges of the foil.
					if np.sqrt((x-i*gold_space-0.5*gold_space)**2+(y-l*gold_space)**2)<Rn:	#If atom in nucleus
						ax_acc=const_in_nucleus*(x-i*gold_space-0.5*gold_space)
						ay_acc=const_in_nucleus*(y-l*gold_space)
					else:
						ax_acc=(const*(x-i*gold_space-0.5*gold_space))/((x-i*gold_space-0.5*gold_space)**2.+(y-l*gold_space)**2.)**(3./2.)
						ay_acc=(const*(y-l*gold_space))/((x-i*gold_space-0.5*gold_space)**2.+(y-l*gold_space)**2.)**(3./2.)
					axx=np.append(axx,ax_acc)
					ayy=np.append(ayy,ay_acc)
			if np.sqrt((x-i*gold_space)**2.+(y-(k+2)*gold_space)**2.)<Rn:	#If atom in nucleus
				ax_acc=const_in_nucleus*(x-i*gold_space)
				ay_acc=const_in_nucleus*(y-(k+2)*gold_space)
			else:
				ax_acc=(const*(x-i*gold_space))/((x-i*gold_space)**2.+(y-(k+2)*gold_space)**2.)**(3./2.)
				ay_acc=(const*(y-(k+2)*gold_space))/((x-i*gold_space)**2.+(y-(k+2)*gold_space)**2.)**(3./2.)	
			axx=np.append(axx,ax_acc)
			ayy=np.append(ayy,ay_acc)
	ax=sum(axx)
	ay=sum(ayy)	
	return (vx,ax,vy,ay)
	
def electricplum(r,t):	
	"""
	Solving differential equations for electric fields for plum pudding atoms
	The sum of every electric field created by each gold nucleus is the resultant electric field.
	When alpha is inside the nucleus, then use equation F=kqQr/R^3,. Elsewhere, F = kqQ/r^2 where k=1/4(pi)Eo. F/m = a
	Every atom away from the edge will be adjacent to 6 other atoms (if the size of the foil allows it).
	"""
	x=r[0]
	vx=r[1]
	y=r[2]
	vy=r[3] 
	axx=np.array([])
	ayy=np.array([])
	if  n_atom_high_parity==0:
		for i in range(n_atom_wide):
			for j in range(int(n_atom_high/2)):
				k = 2*j
				l = 2*j+1
				if np.sqrt((x-i*gold_space)**2+(y-k*gold_space)**2)<R:
					ax_acc=const_in_atom*(x-i*gold_space)
					ay_acc=const_in_atom*(y-k*gold_space)
				else:
					ax_acc=(q*Q*(x-i*gold_space))/(m*4*np.pi*Eo*((x-i*gold_space)**2.+(y-k*gold_space)**2.)**(3./2.))
					ay_acc=(q*Q*(y-k*gold_space))/(m*4*np.pi*Eo*((x-i*gold_space)**2.+(y-k*gold_space)**2.)**(3./2.))			
				axx=np.append(axx,ax_acc)
				ayy=np.append(ayy,ay_acc)
				if i != range(n_atom_wide)[-1]:
					if np.sqrt((x-i*gold_space-0.5*gold_space)**2+(y-l*gold_space)**2)<R:
						ax_acc=const_in_atom*(x-i*gold_space-0.5*gold_space)
						ay_acc=const_in_atom*(y-l*gold_space)
					else:
						ax_acc=(q*Q*(x-i*gold_space-0.5*gold_space))/(m*4*np.pi*Eo*((x-i*gold_space-0.5*gold_space)**2.+(y-l*gold_space)**2.)**(3./2.))
						ay_acc=(q*Q*(y-l*gold_space))/(m*4*np.pi*Eo*((x-i*gold_space-0.5*gold_space)**2.+(y-l*gold_space)**2.)**(3./2.))
				axx=np.append(axx,ax_acc)
				ayy=np.append(ayy,ay_acc)
	elif n_atom_high == 1:
		for i in range(n_atom_wide):
			if np.sqrt((x-i*gold_space)**2+(y)**2)<R:
				ax_acc=const_in_atom*(x-i*gold_space)
				ay_acc=const_in_atom*y
			else:
				ax_acc=(q*Q*(x-i*gold_space))/(m*4*np.pi*Eo*((x-i*gold_space)**2.+(y)**2.)**(3./2.))
				ay_acc=(q*Q*(y))/(m*4*np.pi*Eo*((x-gold_space)**2.+(y)**2.)**(3./2.))			
			axx=np.append(axx,ax_acc)
			ayy=np.append(ayy,ay_acc)

	else:
		for i in range(n_atom_wide):
			for j in range(int(n_atom_high/2)):
				k = 2*j
				l = 2*j+1
				if np.sqrt((x-i*gold_space)**2+(y-k*gold_space)**2)<R:
					ax_acc=const_in_atom*(x-i*gold_space)
					ay_acc=const_in_atom*(y-k*gold_space)
				else:
					ax_acc=(q*Q*(x-i*gold_space))/(m*4*np.pi*Eo*((x-i*gold_space)**2.+(y-k*gold_space)**2.)**(3./2.))
					ay_acc=(q*Q*(y-k*gold_space))/(m*4*np.pi*Eo*((x-i*gold_space)**2.+(y-k*gold_space)**2.)**(3./2.))			
				axx=np.append(axx,ax_acc)
				ayy=np.append(ayy,ay_acc)
				if i != range(n_atom_wide)[-1]:
					if np.sqrt((x-i*gold_space-0.5*gold_space)**2+(y-l*gold_space)**2)<R:
						ax_acc=const_in_atom*(x-i*gold_space-0.5*gold_space)
						ay_acc=const_in_atom*(y-l*gold_space)
					else:
						ax_acc=(q*Q*(x-i*gold_space-0.5*gold_space))/(m*4*np.pi*Eo*((x-i*gold_space-0.5*gold_space)**2.+(y-l*gold_space)**2.)**(3./2.))
						ay_acc=(q*Q*(y-l*gold_space))/(m*4*np.pi*Eo*((x-i*gold_space-0.5*gold_space)**2.+(y-l*gold_space)**2.)**(3./2.))
					axx=np.append(axx,ax_acc)
					ayy=np.append(ayy,ay_acc)
			if np.sqrt((x-i*gold_space)**2.+(y-(k+2)*gold_space)**2.)<R:
				ax_acc=const_in_atom*(x-i*gold_space)
				ay_acc=const_in_atom*(y-(k+2)*gold_space)
			else:
				ax_acc=(q*Q*(x-i*gold_space))/(m*4*np.pi*Eo*((x-i*gold_space)**2.+(y-(k+2)*gold_space)**2.)**(3./2.))
				ay_acc=(q*Q*(y-(k+2)*gold_space))/(m*4*np.pi*Eo*((x-i*gold_space)**2.+(y-(k+2)*gold_space)**2.)**(3./2.))			
			axx=np.append(axx,ax_acc)
			ayy=np.append(ayy,ay_acc)
	ax=sum(axx)
	ay=sum(ayy)	
	return (vx,ax,vy,ay)

#The function that resets animations
def init():
	"""
	Resets animation after each frame.
	"""
	b=()
	for i in range(n_alpha_screen):
			alpha_screen[i][0].set_data([],[])
			b+=alpha_screen[i][0],
	return b

#Used to animate by varying location of each plot with frame. 
#This script animates several alpha particles at once by chosing the number of alpha particles visible at any time (n_alpha_screen).
def animate(i):
	"""
	Creates the animation by going through every array of motion for each alpha.
	Each alpha particle on the screen is equal spaced out.
	"""
	c=()
	z=int(n_alpha_screen*i/t_frames-n_alpha_screen+1)
	zz=n_alpha-n_alpha_screen+1
	zzz=n_alpha%n_alpha_screen
	if n_alpha == 1:
		u=x[i%t_frames]
		v=y[i%t_frames]
		alpha_screen[0][0].set_data(u,v)
		c+=alpha_screen[0][0],
	else:					
		for j in range(n_alpha_screen):
			new_i = i-(j)*fire_t
			if new_i>=0:
				u=x[(n_alpha_screen*(new_i/t_frames)+j)%n_alpha,(new_i)%t_frames]
				v=y[(n_alpha_screen*(new_i/t_frames)+j)%n_alpha,(new_i)%t_frames]
				alpha_screen[j][0].set_data(u,v)
				if z>=zz:
					for jj in range((z-zzz)%n_alpha_screen):
						alpha_screen[(jj+zzz)%n_alpha_screen][0].set_data([],[])
			elif new_i<0:
				alpha_screen[j][0].set_data([],[])
			c+=alpha_screen[j][0],
	return c

#Constants (including the properties of the gold leaf)
eV = 1.6e-19    #Charge of proton/electron
Z = au_atomic_number = 79   #Atomic number of gold
R = au_atom_radius = 1.35e-10   #Radius of gold atom
Rn = au_nucleus_radius = 5.e-15 #Radius of nucleus
gold_space = 2.88e-10   #Separation between gold centres
leaf_thickness = 8.e-8  #Thickness of gold foil
space_between_atom=gold_space-2*au_atom_radius
n_atom_wide =  int(leaf_thickness/gold_space)   #Thickness in terms of number of atoms
#n_atom_wide = 1		#Width of gold leaf in atoms.
actual_thickness = gold_space*n_atom_wide-space_between_atom
#n_atom_high_m = 0.011					#height of the gold leaf in metres.
#n_atom_high = int(n_atom_high_m/gold_space)	#Works out the number of atoms when the height is chosen in metres.
n_atom_high = 3	#The height of the gold leaf in atoms.
n_atom_high_parity =  n_atom_high%2	#Odd or even number of rows?

Q=au_atomic_number*eV		#Charge of gold nucleus
q = 2*eV		#Charge of alpha particle
Eo=8.85418782e-12	#Permittivity of free space
KE0 = 5.e6*1.6e-19	#Initial kinetic energy of alpha particles
m = 6.644e-27		#Mass of alpha particle
const = (q*Q)/(m*4*np.pi*Eo)		#Work out constants to reduce calculation times. (For outside charge region)
const_in_nucleus = (q*Q)/(m*4*np.pi*Eo*Rn**3)	#Work out constants to reduce calculation times. (For rutherford model)
const_in_atom = (q*Q)/(m*4*np.pi*Eo*R**3) 	#Work out constants to reduce calculation times. (For plum pudding model)

plum = 0
rutherford = 1

#choose to have the plum pudding model or the rutherford model. 
#0 or plum for plum pudding. 1 or rutherford for rutherford model.
plum_or_rutherford = 1 

#initial conditions of each alpha. 
#Used random uniform distribution on the height of each alpha particle due to thickness of alpha ray.
gold_centre = ((n_atom_high-1)*gold_space)/2		#The height at the centre of the gold leaf.
theta0 = 0		#Angle of initial velocity of alpha particles anti-clockwise from east.
n_alpha = 3	#Number of alpha particles
alpha_thickness = n_atom_high*gold_space	#Thickness of alpha particles. Programmed so the thickness is the height of the gold leaf.
x0=-10*R		#Starting horizontal distance from the edge of the gold leaf.

#All must be arrays
x0=np.ones(n_alpha)*x0	#Starting horizontal distance in array
y0=np.array(map(lambda x: rn.uniform(-alpha_thickness/2,alpha_thickness/2), range(n_alpha)))+gold_centre	#Random initial heights of alpha particles
#y0=np.array([0]*n_alpha)+gold_centre 	#Used to test if model worked by firing at centre of nuclei.
vx0 = np.ones(n_alpha)*np.sqrt(2*KE0/m)*np.cos(theta0)	#Horizontal initial velocity
vy0=np.ones(n_alpha)*np.sqrt(2*KE0/m)*np.sin(theta0)		#Vertical initial velocity

#Time for for each alpha particle's movement. Limited for clarity when plotting and animating.
tmax = (actual_thickness+4*np.absolute(x0[0]))/np.sqrt(2*KE0/m)	#Works out distance/max velocity. 

#Change these values to modify the speed of the animation.
animation_speed=1  
t_frames=30
t=np.linspace(0.,tmax,t_frames)    #time elapsed in seconds


if plum_or_rutherford == 0:
	#Creates an array with values for variables of motion for the first alpha particle.
	soln=spi.odeint(electricplum,[x0[0],vx0[0],y0[0],vy0[0]],t,mxstep=5000) #mxstep=5000 to increase number of steps in intergration

	#Stacks the array with values for variables of motion for the alpha particles afer the first. This makes a 2D array.
	for i in range(n_alpha-1):
		soln=np.vstack((soln,spi.odeint(electricplum,[x0[i+1],vx0[i+1],y0[i+1],vy0[i+1]],t,mxstep=5000))) 
		
else:
	#Creates an array with values for variables of motion for the first alpha particle.
	soln=spi.odeint(electric,[x0[0],vx0[0],y0[0],vy0[0]],t,mxstep=5000) #mxstep=5000 to increase number of steps in intergration

	#Stacks the array with values for variables of motion for the alpha particles afer the first. This makes a 2D array.
	for i in range(n_alpha-1):
		soln=np.vstack((soln,spi.odeint(electric,[x0[i+1],vx0[i+1],y0[i+1],vy0[i+1]],t,mxstep=5000))) 
		
#Assigns each variable of motion from the 2D array. This array has only the data for the first alpha particle.
x = soln[:t_frames,0]
y = soln[:t_frames,2]
vx = soln[:t_frames,1]
vy = soln[:t_frames,3]

#Stacks each array with values for other particles. This makes a 2D array and has a row for each particle.
for i in range(n_alpha-1):	
	x = np.vstack((x,soln[t_frames*(i+1):t_frames*(i+2),0]))	
	y = np.vstack((y,soln[t_frames*(i+1):t_frames*(i+2),2]))
	vx = np.vstack((vx,soln[t_frames*(i+1):t_frames*(i+2),1]))
	vy = np.vstack((vy,soln[t_frames*(i+1):t_frames*(i+2),3]))

#Plotting
fig = plt.figure(num=1,figsize=(15,9))	#Create figure for alpha particles firing at gold.

#Modifing figure and axes for clarity.
plot1 = fig.add_subplot(1,1,1)
#plot1.set_xlim(100*x0[0],actual_thickness-x0[0])
plot1.set_ylim(-2.*gold_space,(n_atom_high+1)*gold_space)
plot1.set_xlim(-4e-8,actual_thickness-x0[0])

#n_alpha_screen is the number of alpha particles visible on the screen at any one time.
n_alpha_screen = 5
#Number of alphas on the screen at one time cannot be more than the number of frames for each alpha particle and not more than the number of alpha particles in total.
#This is because each alpha particle are equally seperated improve presentation.
if n_alpha_screen > t_frames:
	n_alpha_screen = t_frames
if n_alpha_screen > n_alpha:
	n_alpha_screen = 1
alpha_screen = []
for i in range(n_alpha_screen):
	alpha_screen.append(plt.plot([],[],'ro'))

#Used for animation. Creating these variables will reduce the calculation times while the function 'animate' is active.
fire_t = t_frames/n_alpha_screen
zz=n_alpha-n_alpha_screen+1
zzz=n_alpha%n_alpha_screen

#Generating the circles for the atoms. Tried generating Nuclei circles, but they were too small to be visible.
#Every atom away from the edge will be adjacent to 6 other atoms (if the size of the foil allows it).
if  n_atom_high_parity==0:
	for i in range(n_atom_wide):
		for j in range(n_atom_high/2):
			k = 2*j
			l = 2*j+1
			circles = plt.Circle((i*gold_space,k*gold_space),au_atom_radius,color='y')
			plot1.add_patch(circles)
			if i != range(n_atom_wide)[-1]:
				circles = plt.Circle((0.5*gold_space+i*gold_space,l*gold_space),au_atom_radius,color='y')
				plot1.add_patch(circles)
elif n_atom_high ==1:
	for i in range(n_atom_wide):
		circles = plt.Circle((i*gold_space,0),au_atom_radius,color='y')
		plot1.add_patch(circles)
else:
	for i in range(n_atom_wide):
		for j in range(int(n_atom_high/2)):
			k = 2*j
			l = 2*j+1
			circles = plt.Circle((i*gold_space,k*gold_space),au_atom_radius,color='y')
			plot1.add_patch(circles)
			if i != range(n_atom_wide)[-1]:
				circles = plt.Circle((0.5*gold_space+i*gold_space,l*gold_space),au_atom_radius,color='y')
				plot1.add_patch(circles)
		circles = plt.Circle((i*gold_space,(k+2)*gold_space),au_atom_radius,color='y')
		plot1.add_patch(circles)
		
#Starts the animation.
anim=an.FuncAnimation(fig,animate,init_func=init,frames=int(t_frames*(1.+(n_alpha-1.)/n_alpha_screen)),interval=animation_speed,blit=True)

#Labels.
plt.title("Path of Alpha Particle due to Gold Nucleus")
plt.xlabel("Horizontal Displacement (m)")
plt.ylabel("Vertical Displacement (m)")

#Used to save the animation.
#writer = an.writers['ffmpeg'](fps=30)
#dpi=150
#anim.save('real atoms - 2.mp4',writer=writer,dpi=dpi)

all_deflected_angle = np.array([])
n_reflected=0

#Counts number of reflected alphas. Due to Python being unable to calculate the correct angle
#when reflected, adding or subtracting pi is necessary.
if n_alpha == 1:
	deflected_angle = np.arctan(vy[-1]/vx[-1])
	if vx[-1]<0:
		n_reflected+=1
		if vy[-1]>=0:
			deflected_angle+=np.pi
		elif vy[-1]<0:
			deflected_angle-=np.pi
	all_deflected_angle = np.append(all_deflected_angle,deflected_angle)

else:
	for i in range(n_alpha):
		deflected_angle = np.arctan(vy[i][-1]/vx[i][-1])
		if vx[i][-1]<0:
			n_reflected+=1
			if vy[i][-1]>=0:
				deflected_angle+=np.pi
			elif vy[i][-1]<0:
				deflected_angle-=np.pi
		all_deflected_angle = np.append(all_deflected_angle,deflected_angle)

#Changes to degrees.
all_deflected_angle = all_deflected_angle*180/np.pi

#Velocity graph of the first alpha.
fig3 = plt.figure(num=2)
plt.plot(t,vx[0])

#Plotting histograms to show the distribution of angle deflections.
fig2 = plt.figure(num=3,figsize=(10,6))
ax = fig2.add_subplot(1,2,1)
plt.hist(all_deflected_angle,n_alpha)
plt.xlabel("Angle of Deflection")
plt.ylabel("Number of Alpha Particles")
plt.title("Gold Leaf Experiment Model")

#Shows the number of alphas released initially and the number reflected.
txt = fig2.add_subplot(1,2,2)
txt0 = plt.text(0,0.78, "Number of Alpha Particles Released", fontsize="large")
txt01 = plt.text(0,0.70, n_alpha,weight="bold", fontsize="large")
txt1 = plt.text(0,0.58,"Number of Alpha Particles Reflected",fontsize="large")
txt2 = plt.text(0,0.5,n_reflected,weight="bold",fontsize="large")
plt.axis("off")

plt.show()