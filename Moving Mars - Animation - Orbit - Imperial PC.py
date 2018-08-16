"""
Variables that can be changed:
escape
r0, x0, y0, m
vx0, vy
theta0, initangle
tmax, tframes, animation_speed

The code near the bottom allows for animation to be recorded.
"""

import numpy as np
import scipy.integrate as spi
import matplotlib.pyplot as plt
import matplotlib.lines as matlines
import matplotlib.animation as an

def gravity(r,t):	
	"""
	Solves differential equations for gravity in cartesian coordinaes
	Stops moving when crash on Mars.
	Includes the speed of Mars for inertial frame.
	"""
	x=r[0]
	vx=r[1]
	y=r[2]
	vy=r[3]     
	x_mars=r[4]
	y_mars=r[5]
	if np.sqrt((x-x_mars)**2.+(y-y_mars)**2.)<=R:
		vx=vxmars
		vy=vymars
		ax=0
		ay=0
	else:
		ax=-G*M*(x-x_mars)/((x-x_mars)**2.+(y-y_mars)**2.)**(3./2.)
		ay=-G*M*(y-y_mars)/((x-x_mars)**2.+(y-y_mars)**2.)**(3./2.)
	vx_mars=vxmars
	vy_mars=vymars
	return (vx,ax,vy,ay,vx_mars,vy_mars)

def gravitymove(r,t):	
	"""
	Solves differential equations for gravity in cartesian coordinaes
	Stops moving when crash on Mars.
	Does not include the speed of Mars for non-inertial frame
	"""
	xr=r[0]
	vxr=r[1]
	yr=r[2]
	vyr=r[3] 
	if np.sqrt(xr**2.+yr**2.)<=R:
		vxr=0
		vyr=0
		axr=0
		ayr=0
	else:
		axr=-G*M*xr/(xr**2.+yr**2.)**(3./2.)
		ayr=-G*M*yr/(xr**2.+yr**2.)**(3./2.)
	return (vxr,axr,vyr,ayr)

def init():
	"""
	Used to clear each frame.
	After each frame, Mars resets to near infinity.
	"""
	orbit.set_data([],[])
	mars.center = (1e308,1e308)
	ax.add_patch(mars)	
	GPEanim.set_data([],[])
	KEanim.set_data([],[])
	GPE_KEanim.set_data([],[])
	orbitmove.set_data([],[])
	return orbit, mars, orbitmove,


def animate(i):				#Animating the orbit. The animation is at frame i.
	"""
	Animates orbits in both frames, Mars and also the energy and allows for repeated orbits.
	The animation is at frame i.
	Function goes through all elements of arrays of x and y coordinates of the satellite, Mars and the energies in moving frame.
	"""
	uorbit=x[i%t_frames]		#Goes through x-coordinates of satellite as animation continues.
	vorbit=y[i%t_frames]		#Goes through y-coordinates of satellite as animation continues.
	orbit.set_data(uorbit,vorbit)	#Sets coordinates of satellite at each frame.
	umars=xmars[i%t_frames]	#Goes through x-coordinates of Mars as animation continues.
	vmars=ymars[i%t_frames]	#Goes through y-coordinates of Mars as animation continues.
	mars.center = (umars, vmars)	#Sets coordinates of Mars's centre.
	uGPE = t[:i%t_frames]
	vGPE = GPE[:i%t_frames]
	GPEanim.set_data(uGPE,vGPE)
	uKE = t[:i%t_frames]
	vKE = KE[:i%t_frames]
	KEanim.set_data(uKE,vKE)
	uGPE_KE = t[:i%t_frames]
	vGPE_KE = GPE_KE[:i%t_frames]	
	GPE_KEanim.set_data(uGPE_KE,vGPE_KE)
	uorbitmove=xmove[i%t_frames]
	vorbitmove=ymove[i%t_frames]
	orbitmove.set_data(uorbitmove,vorbitmove)
	return orbit, mars, GPEanim, KEanim, GPE_KEanim, orbitmove

#Constants
G=6.67e-11	#gravitational constant
M=6.4e23		#mass of Mars
m=260.		#mass of satellite
R=3.4e6		#Radius of Mars
mars_speed=24100.	#Speed of Mars

#initial conditions
r0=5*R      #The satellite's initial distance from Mars's initial position (0,0)
theta0=0		#The satellite's position's angle anti-clockwise from east of Mars's initial position (0,0)
if theta0>=np.pi or theta0<=-np.pi:		#Makes sure the angle is between pi and -pi.
	n_pi=round(theta0/(2*np.pi),0)
	theta0=theta0-n_pi*2*np.pi
	
initangle=3*np.pi/4	#The direction of the satellite's velocity with angle anti-clockwise from east of the satellite's position
if initangle>=np.pi or initangle<=-np.pi:	#Makes sure the angle is between pi and -pi.
	n_pi=round(initangle/(2*np.pi),0)	
	initangle=initangle-n_pi*2*np.pi
	
escape = 0	#Type 0 for the satellite to have a stable orbit, or type 1 to escape

vmarsangle=np.pi/2					#Mars's velocity's angle anti-clockwise from east of Mars's initial position (0,0).
vxmars=mars_speed*np.cos(vmarsangle)	#Mars's velocty in x-direction.
vymars=mars_speed*np.sin(vmarsangle)	#Mars's velocity in y-direction.

psi0=initangle-np.pi/2				#The angle anti-clockwise from north to the direction of velocity. Used to calculated x,y vectors of escape velocity (Only possible if Mars is stationary).
x0=r0*np.cos(theta0)    				#initial displacement in x-direction from Mars's initial position (0,0)
y0=r0*np.sin(theta0)    				#initial displacement in y-direction from Mars's initial position (0,0)
v_escape=np.sqrt((2*G*M)/r0)			#escape velocity for the satellite in Mars's moving frame.
vx0_escape=-v_escape*np.sin(psi0)    	#escape velocity in x-direction in Mars's moving frame.
vy0_escape=v_escape*np.cos(psi0)    	#escape velocity in y-direction in Mars's moving frame.

#Using this ensures that there will be a stable orbit around mars. (Only works when the velocity of Mars is taken into account, as shown below).
if vx0_escape<-1 and vy0_escape<-1:		#Direction away from Mars is Bottom Left.
	vxsign=1							
	vysign=1
elif vx0_escape<-1 and vy0_escape>1:		#Direction away from Mars is Top Left
	vxsign=1
	vysign=-1
elif vx0_escape>1 and vy0_escape>1:		#Direction away from Mars is Top Right
	vxsign=-1
	vysign=-1
elif vx0_escape>1 and vy0_escape<-1:		#Direction away from Mars is Bottom Right
	vxsign=-1
	vysign=1
elif vx0_escape<1 and vx0_escape>-1 and vy0_escape>1:	#Direction away from Mars is Up
	vxsign=0
	vysign=-1
elif vx0_escape<1 and vx0_escape>-1 and vy0_escape<-1:	#Direction away from Mars is Down
	vxsign=0
	vysign=1
elif vx0_escape>1 and vy0_escape<1 and vy0_escape>-1:	#Direction away from Mars is Right
	vxsign=-1
	vysign=0
elif vx0_escape<-1 and vy0_escape<1 and vy0_escape>-1:	#Direction away from Mars is Left
	vxsign=1
	vysign=0

#From previous if statements, vxsign and vysign will have a sign that will allow for a stable orbit.
#If an unstable orbit is chosen (i.e. escape == 1), the signs will be reversed.
if escape==1:
	vxsign=vxsign*-1	
	vysign=vysign*-1
	
#Initial velocity of the satellite. Takes into account the velocity of mars to allow vxsign/vysign to work.
vx0=vx0_escape+vxsign*200+vxmars	
vy0=vy0_escape+vysign*200+vymars	

animation_speed=5	#The 'interval' argument for matplotlib.animation.FuncAnimation. Number of milliseconds between each frame
t_frames=400		#Number of frames.
tmax = 500000.
t=np.linspace(0.,tmax,t_frames)    #time elapsed in seconds

soln=spi.odeint(gravity,[x0,vx0,y0,vy0,0,0],t)	#Integration to calculate distance and speed in inertial frame

solnmove = spi.odeint(gravitymove,[x0,vx0-vxmars,y0,vy0-vymars],t)	#Integration to calculate distance and speed of satellite in Mars's moving frame.

#Inertial Frame
x=soln[:,0]				#x-coordinates of satellite
y=soln[:,2]				#y-coordinates of satelite
vx=soln[:,1]				#x-velocity coordinates of satellite
vy=soln[:,3]				#y-velocity coordinates of satellite
xmars=soln[:,4]			#x-coordinates of Mars
ymars=soln[:,5]			#y-coordinates of Mars
vxmars=np.array([vxmars]*t_frames)	#Array of mars's horizontal velocity
vymars=np.array([vymars]*t_frames)	#Array of mars's vertical velocity

#Mars's Moving Frame
xmove =solnmove[:,0]				#x-coordinates of satellite
ymove =solnmove[:,2]				#y-coordinates of satelite
vxmove =solnmove[:,1]				#x-velocity coordinates of satellite
vymove =solnmove[:,3]				#y-velocity coordinates of satellite

crash="The satellite does not crashed on Mars"	#Shows a message to say the satellite does not crash, if the satellite doesn't crash.
L=0		#Changes to 1 when satellite crashes.

#Detects when satellite crashes on Mars. If the satellite crashes on Mars, it remains stationary on Mars.
for i in range(x.size):	
	if x[i]<=xmars[i]+R and x[i]>=xmars[i]-R:
		if y[i]<=ymars[i]+R and y[i]>ymars[i]-R:
			crash="The satellite crashes on Mars"
			L=1
			break

r=np.sqrt((x-xmars)**2+(y-ymars)**2)		#radius from Mars's centre
GPE=-G*M*m/r			#Gravitational Potential Energy
KE=0.5*m*((vx-vxmars)**2+(vy-vymars)**2)	#Kinetic Energy
GPE_KE=GPE+KE			#GPE + KE
KEinert=0.5*m*(vx**2+vy**2) #Kinetic Energy in Inertial Frame
GPE_KEinert=GPE+KEinert		#Total Energy in Inertial Frame

#Lines for apoareion and periareion for the satellite
peri_x=np.linspace(xmove[np.where(r==min(r))[0][0]],0,10000)	#x-coordinates for periareion
peri_y=np.linspace(ymove[np.where(r==min(r))[0][0]],0,10000)	#y-coordinates for periareion

if escape == 0 and L==0:
	apo_x=np.linspace(0,xmove[np.where(r==max(r))[0][0]],10000)	#x-coordinates for apoareion
	apo_y=np.linspace(0,ymove[np.where(r==max(r))[0][0]],10000)	#y-coordinates for apoareion
	a_line_x = np.append(peri_x,apo_x)					#x-coordinates for the major axis line
	a_line_y = np.append(peri_y,apo_y)					#y-coordinates for the major axis line
	grad_a = (apo_y[-1]-peri_y[0])/(apo_x[-1]-peri_x[0])	#gradient of the major axis
	grad_b = -1/grad_a								#gradient of the minor axis
	angle_b = np.arctan(grad_b)						#anti-clockwise angle between the minor axis line and the horizontal.

#Using this allows the centre of the orbit to be correctly recorded in Mars's moving frame.
if peri_x[0]<0 and peri_y[0]<0:	#Periareion is at Bottom Left
	centrex_sign=1
	centrey_sign=1
	orbit_legend="upper left"		#Location of Legend is Upper Left
elif peri_x[0]<0 and peri_y[0]>0:	#Periareion is at Top Left
	centrex_sign=1
	centrey_sign=-1
	orbit_legend="upper right"	#Location of Legend is Upper Right
elif peri_x[0]>0 and peri_y[0]>0:	#Periareion is at Top Right
	centrex_sign=-1
	centrey_sign=-1
	orbit_legend="upper left"		#Location of Legend is Upper Left
elif peri_x[0]>0 and peri_y[0]<0:	#Periareion is at Bottom Right
	centrex_sign=-1
	centrey_sign=1
	orbit_legend="upper right"	#Location of Legend is Upper Right
elif peri_x[0]<1 and peri_x[0]>-1 and peri_y[0]>1:		#Periareion is at Up
	centrex_sign=0
	centrey_sign=-1
	orbit_legend="upper right"			#Location of Legend is Upper Right
elif peri_x[0]<1 and peri_x[0]>-1 and peri_y[0]<-1:	#Periareion is at Down
	centrex_sign=0
	centrey_sign=1
	orbit_legend="upper right"			#Location of Legend is Upper Right
elif peri_x[0]>1 and peri_y[0]<1 and peri_y[0]>-1:		#Periareion is at Right
	centrex_sign=-1
	centrey_sign=0
	orbit_legend="upper right"			#Location of Legend is Upper Right
elif peri_x[0]<-1 and vy0_escape<1 and vy0_escape>-1:	#Periareion is at Left
	centrex_sign=1
	centrey_sign=0
	orbit_legend="upper right"			#Location of Legend is Upper Right

#Calculations to plot the major and the minor axis lines
if escape == 0 and L==0:
	centre_x = centrex_sign*(np.absolute(apo_x[-1])+np.absolute(peri_x[0]))/2+peri_x[0]	#x-coordinates of the centre of the orbit
	centre_y = centrey_sign*(np.absolute(apo_y[-1])+np.absolute(peri_y[0]))/2+peri_y[0]	#y-coordinates of the centre of the orbit

	apo_dist = max(r)			#Distance of Apoareion
	a_dist = (max(r)+min(r))/2	#semi-major axis
	b_dist = np.sqrt(max(r)*min(r)) #semi-minor axis	
	bx1 = centre_x+b_dist*np.cos(angle_b)	#x-coordinate of one end of the minor axis line
	by1 = centre_y+b_dist*np.sin(angle_b)	#y-coordinate of one end of the minor axis line
	bx2 = centre_x-b_dist*np.cos(angle_b)	#x-coordinate of the other end of the minor axis line
	by2 = centre_y-b_dist*np.sin(angle_b)	#y-coordinate of the other end of the minor axis line
	bx = np.linspace(bx2,bx1,100)		#x-coordinates for the minor axis line
	by = np.linspace(by2,by1,100)		#y-coordinates for the minor axis line

	T=np.sqrt(4*np.pi**2*a_dist**3/(G*M)) #Time Period
	
peri_dist = min(r) #Distance of periareion	

#Plot seperate energy graphs to see how each type of energy changes in inertial and in Mars's moving frame.
fig2 = plt.figure("Energy",figsize=(15,9))
zerolength=np.linspace(t[-1]*-0.1,t[-1]*1.1,100)

#KE
energy1 = plt.subplot(4,2,1)
plt.plot(zerolength,np.zeros([zerolength.size]),'k--')
KEplot1, = plt.plot(t,KE,'r-', label="KE",lw=2)
energy1handles, energy1labels = energy1.get_legend_handles_labels()
plt.legend(energy1handles, energy1labels,loc="upper right")
plt.axis([t[-1]*-0.1,t[-1]*1.1,min(KE)*0.9,max(KE)*1.1])
plt.ylabel("Energy (J)", weight="bold")
plt.title("Energy of the Satellite in Mars's Moving Frame",weight="bold",fontsize="small")

#KE in Inertial Frame
energy2 = plt.subplot(4,2,2)
plt.plot(zerolength,np.zeros([zerolength.size]),'k--')
KEplot2, = plt.plot(t,KEinert,'r-',label="KE",lw=2)
energy2handles, energy2labels = energy2.get_legend_handles_labels()
plt.legend(energy2handles, energy2labels,loc="upper right")
plt.title("Energy of the Satellite in Inertial Frame",weight="bold",fontsize="small")
plt.ylabel("Energy (J)", weight="bold")
plt.axis([t[-1]*-0.1,t[-1]*1.1,min(KEinert)*0.9,max(KEinert)*1.1])

#Total Energy
energy3 = plt.subplot(4,2,3)
plt.plot(zerolength,np.zeros([zerolength.size]),'k--')
GPE_KEplot1, = plt.plot(t,GPE_KE,'g-', label="Total Energy",lw=2)
energy3handles, energy3labels = energy3.get_legend_handles_labels()
plt.legend(energy3handles, energy3labels,loc="upper right")
plt.ylabel("Energy (J)", weight="bold")
plt.axis([t[-1]*-0.1,t[-1]*1.1,min(GPE)*1.1,max(KE)*1.1])

#Total Energy in Inertial Frame
energy4 = plt.subplot(4,2,4)
plt.plot(zerolength,np.zeros([zerolength.size]),'k--')
GPE_KEplot2, = plt.plot(t,GPE_KEinert,'g-',label="Total Energy",lw=2)
energy4handles, energy4labels = energy4.get_legend_handles_labels()
plt.legend(energy4handles, energy4labels,loc="upper right")
plt.axis([t[-1]*-0.1,t[-1]*1.1,min(GPE_KEinert)*0.9,max(GPE_KEinert)*1.1])
plt.xlabel("Time (s)", weight="bold")
plt.ylabel("Energy (J)", weight="bold")

#GPE
energy5 = plt.subplot(4,2,5)
plt.plot(zerolength,np.zeros([zerolength.size]),'k--')
GPEplot1, = plt.plot(t,GPE,'b-',label="GPE",lw=2)
energy5handles, energy5labels = energy5.get_legend_handles_labels()
plt.legend(energy5handles, energy5labels,loc="upper right")
plt.axis([t[-1]*-0.1,t[-1]*1.1,min(GPE)*1.1,max(GPE)*0.9])
plt.ylabel("Energy (J)", weight="bold")
plt.xlabel("Time (s)", weight="bold")

#Text, showing information of the model.
plot1 = fig2.add_subplot(4,2,6)
titletxt1 = plot1.text(-0.10,-1.00,"Mars in Inertial Frame")
txt11 = plt.text(-0.10,-1.15,'Mass = '+"%0.4g"%(M)+"kg")
txt12 = plt.text(-0.10,-1.30,'Radius = '+"%0.4g"%(R)+"m")
txt13 = plt.text(-0.10,-1.45,'Speed of Mars = '+"%0.4g"%(mars_speed)+"m/s")

titletxt2 = plt.text(0.40,0.90,"Satellite in Inertial Frame")
txt21 = plt.text(0.4,0.75,'Mass = '+"%0.4g"%(m)+"kg")
txt22 = plt.text(0.4,0.60,"Initial Distance from Mars's Centre = "+"%0.4g"%(r0)+"m")
txt23 = plt.text(0.4,0.45,"Initial Speed = "+"%0.4g"%(round(np.sqrt(vx0**2+vy0**2),3))+"m/s")
txt24 = plt.text(0.4,0.30,"Slowest Speed = "+"%0.4g"%(min(np.sqrt(vx**2+vy**2)))+"m/s")
txt25 = plt.text(0.4,0.15,"Fastest Speed = "+"%0.4g"%(max(np.sqrt(vx**2+vy**2)))+"m/s")
txt26 = plt.text(0.4,0.00,"Initial KE = "+"%0.4g"%(KEinert[0])+"J")
txt27 = plt.text(0.4,-0.15,"Minimum KE = "+"%0.4g"%(min(KEinert))+"J")
txt28 = plt.text(0.4,-0.30,"Maximum KE = "+"%0.4g"%(max(KEinert))+"J")

titletxt3 = plt.text(-0.10,0.90,"Satellite in Mars's Moving Frame")
txt31 = plt.text(-0.10,0.75,"Initial Speed = "+"%0.4g"%(np.sqrt((vx0-vxmars[0])**2+(vy0-vymars[0])**2))+"m/s")
txt32 = plt.text(-0.10,0.60,"Slowest Speed = "+"%0.4g"%(min(np.sqrt((vx-vxmars[0])**2+(vy-vymars[0])**2)))+"m/s")
txt33 = plt.text(-0.10,0.45,"Fastest Speed = "+"%0.4g"%(max(np.sqrt((vx-vxmars[0])**2+(vy-vymars[0])**2)))+"m/s")
txt34 = plt.text(-0.10,0.30,"Initial KE = "+"%0.4g"%(KE[0])+"J")
txt35 = plt.text(-0.10,0.15,"Minimum KE = "+"%0.4g"%(min(KE))+"J")
txt36 = plt.text(-0.10,0.00,"Maximum KE = "+"%0.4g"%max(KE)+"J")
txt37 = plt.text(-0.10,-0.15,"Initial GPE = "+"%0.4g"%(GPE[0])+"J")
txt38 = plt.text(-0.10,-0.30,"Minimum GPE = "+"%0.4g"%(min(GPE))+"J")
txt39 = plt.text(-0.10,-0.45,"Maximum GPE = "+"%0.4g"%(max(GPE))+"J")
txt310 = plt.text(-0.10,-0.60,"Total Energy = "+"%0.4g"%(np.mean(GPE_KE))+"J")
txt311 = plt.text(-0.10,-0.75,"Escape Velocity = "+"%0.4g"%(np.sqrt(vx0_escape**2+vy0_escape**2))+"m/s")

#The if statement allows more information to be added if possible
if escape == 0 and L==0:
	txt401 = plt.text(0.40,-0.55, 'Semi-Major Axis = '+"%0.4g"%(a_dist)+"m")
	txt402 = plt.text(0.40,-0.70,'Semi-Minor Axis = '+"%0.4g"%(b_dist)+"m")
	txt403 = plt.text(0.40,-0.85,'Apoareion = '+"%0.4g"%(apo_dist)+"m")
	txt404 = plt.text(0.40,-1.00,'Time Period = '+"%0.4g"%(T)+"s")
	
txt42 = plt.text(0.40,-1.15,'Periareion = '+"%0.4g"%(peri_dist)+"m")
txt43 = plt.text(0.40,-1.30,'Speed of Mars = '+"%0.4g"%(mars_speed)+"m/s")
txt44 = plt.text(0.40,-1.45,"Time Elapsed = "+"%0.4g"%(t[-1])+"s")
txt45 = plt.text(0.40,-1.60,crash)
plt.axis("off")

#This makes the font of the text easier to modify.
txt_set=[txt11,txt12,txt13,txt21,txt22,txt23,txt24,txt25,txt26,txt27,txt28,txt31,txt32,txt33,txt34,txt35,txt36,txt37,txt38,txt39,txt310,txt311,txt42,txt43,txt44,txt45]
if escape == 0 and L==0: #If satellite makes orbit and doesn't crash.
	txt_set+=[txt401,txt402,txt403,txt404]
title_set=[titletxt1,titletxt2,titletxt3]
for text in txt_set:
	text.set_fontsize("medium")
for title in title_set:
	title.set_fontsize("small")
	title.set_weight("bold")
plt.tight_layout(w_pad=3, h_pad=0.2, pad=5)

#Plotting
fig=plt.figure("Main",figsize=(16,9))	#Figure size

#Trajectory (Inertial Frame)
ax = fig.add_subplot(1,3,3)	#Right
ax.set_xlim(min(x),max(x))	#Change x axis for clarity
ax.set_ylim(min(y),max(y))	#Change y axis for clarity
orbit, = plt.plot([],[],'bo',label="Satellite",linestyle='--')	#Allows for animation of satellite
orbitpath, = plt.plot(x,y,'b--')				#Plot Path of Satellite
mars = plt.Circle((1e308,1e308),R,color="r",)	#Add Mars to Diagram		
marspath, = plt.plot(xmars,ymars,'r--')		#Plot Path of Mars
mars_circle=matlines.Line2D(range(1),range(1),color='red',marker='o',markersize=15,label='Mars',linestyle="None")		#Mars's Legend
plt.legend([mars_circle,orbit],("Mars", "Satellite"),numpoints=1,handlelength=2,loc="upper right")		#Add legend, and modified for clarity.
if L==1:
	ax.set_xlim(min(xmars)-2*R,max(x)+2*R)	#Change x axis for clarity
	ax.set_ylim(min(y)-2*R,max(y)+2*R)	#Change y axis for clarity
	
plt.title("Inertial Frame")
plt.xlabel("Horizontal Displacement (m)")
plt.ylabel("Vertical Displacement (m)")

#Trajectory (Moving Frame)
axmove = fig.add_subplot(1,3,1)	#Left
axmove.set_xlim(min(xmove),max(xmove))	#Change x axis for clarity
axmove.set_ylim(min(ymove),max(ymove))	#Change y axis for clarity
orbitmove, = plt.plot([],[],'bo',label="Satellite",linestyle='--')	#Allows for animation of satellite
orbitpathmove, = plt.plot(xmove,ymove,'b--')				#Plot Path of Satellite
marsmove = plt.Circle((0,0),R,color="r",)	#Add Mars to Diagram		
mars_circle=matlines.Line2D(range(1),range(1),color='red',marker='o',markersize=15,label='Mars',linestyle="None")		#Mars's Legend
axmove.add_patch(marsmove)
if L==0 and escape == 0:	#If satellite makes orbit and doesn't crash.
	Major_Axis, = plt.plot(a_line_x,a_line_y,'k--', label='Major Axis')	#Plot Major Axis
	Minor_Axis, = plt.plot(bx,by,'g--', label='Minor Axis')			#Plot Minor Axis
	plt.legend([mars_circle,orbitmove,Major_Axis,Minor_Axis],("Mars", "Satellite", "Major Axis", "Minor Axis"),numpoints=1,handlelength=2,loc=orbit_legend)	#Add Legend
elif L==1:
	axmove.set_xlim(-2*R,2*R)	#Change x axis for clarity
	axmove.set_ylim(-2*R,2*R)	#Change y axis for clarity
	plt.legend([mars_circle,orbit],("Mars", "Satellite"),numpoints=1,handlelength=2,loc=orbit_legend)		#Add Legend	
else:
    plt.legend([mars_circle,orbit],("Mars", "Satellite"),numpoints=1,handlelength=2,loc=orbit_legend)		#Add Legend	
plt.title("Mars's Moving Frame")
plt.xlabel("Horizontal Displacement (m)")
plt.ylabel("Vertical Displacement (m)")

#Graph of Energy against Time
plot2 = fig.add_subplot(1,3,2)				#Middle	
GPEplot, = plt.plot(t,GPE,'b--',label="GPE")				#Plot GPE
GPEanim, = plt.plot([],[],'b-',lw=2)			#Animation of GPE graph
KEplot, = plt.plot(t,KE,'r--',label="KE")				#Plot KE
KEanim, = plt.plot([],[],'r-',lw=2)			#Animation of KE graph
GPE_KEplot, = plt.plot(t,GPE_KE,'g--',label="Total Energy")	#Plot Total Energy
GPE_KEanim, = plt.plot([],[],'g-',lw=2)		#Animation of Total Energy Graph
plot2.set_xlim(t[-1]*-0.1,t[-1]*1.1)			#Change time-axis so all of each energy type will be clearly visible
plt.legend([GPEplot,KEplot,GPE_KEplot],("GPE", "KE", "Total Energy"),loc="upper right")		#Legend
plt.title("Satellite's Energy in Moving Frame")					
plt.xlabel("Time (s)")								
plt.ylabel("Energy (J)")

plt.tight_layout() #Ensures all plots and animations is as big as possible for clarity
	
anim=an.FuncAnimation(fig,animate,init_func=init,frames=t_frames,interval=animation_speed,blit=True)		#Animate Orbit

#Can save animation if wanted. Uses ffmpeg for the library of codecs.
#writer = an.writers['ffmpeg'](fps=30)
#dpi=300
#anim.save('orbit.mp4',writer=writer,dpi=dpi)

plt.show()