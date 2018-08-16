import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as an

"""
method_list (only single pendulum has sin(theta)=theta approximation and Implicit Euler method):
'0. Explicit Euler'
'1. Leapfrog'
'2. RK4'
'3. Implicit Euler'

type_list:
'0. Single Pendulum'
'1. Double Pendulum'
"""

############ CONFIGURATION ###############

method_list = ['Explicit Euler', 'Leapfrog', 'RK4', 'Implicit Euler']
type_list = ['Single Pendulum', 'Double Pendulum']
approx = 'n'             #Only used for single pendulum. Assumes sin(theta) = theta. 'Y' or 'y' for yes, anything else for no.
method = method_list[1]  #Choose your FDM from the list
Type = type_list[0]      #Choose your problem type in the list. Single or Double pendulum?

if Type == type_list[1] and method == method_list[3]:
    raise TypeError("Double Pendulum does not have Implicit Euler method")

"""        
These lines are used to find the critical step length. The algorithm finds the largest step length which satisfies the conditions below.
If the system is undamped: (Only available for single pendulum)
-If the system is approximated with sin(theta) = theta, then 'theta_difference' is the criteria for maximum theta difference between FDM and analytical solution.
-If the system is not approximated, the 'energy_percentage' value is the maximum/minimum mean average energy allowed as a percentage of initial energy.
If the system is damped, the 'energy_percentage' value is the percentage of initial energy the final energy value must fall below.
"""
find_critical_step_length = 'n'    #'Y' or 'y' for yes, anything else for no.
critical_step_start = 10       
critical_step_end = 9
number_of_steps = 10
energy_percentage = 1 #Units of %. Usually 1 is a good value. Only used for undamped non-appoximated single-pendulum and undamped double pendulum.
theta_difference = 0.001 #Only used for undamped approximated single pendulum system.

#Precautions because algorithm does not find critical step size for double pendulum. 
if Type == type_list[1]:
    if find_critical_step_length.upper() == 'Y':
        raise TypeError("Critical step size algorithm is only available for single pendulum system. You have chosen double pendulum.")
        
#Needed for both single and double pendulum.
h = 1             #Ignore this if you are finding critical step-size.
tmax = 40          #Time when model ends. The time is in terms of t hat, not t.
theta0 = 1
omega0 = 0

m = 1.
l = 1.
g = 1.
tau = np.sqrt(float(l)/g) #t = tau * t_hat (where t is the usual time value)
D = 1

#D_hat = D/(m*np.sqrt(l*g))
D_hat = D*tau/(m*l)

#Only for double pendulum.
phi0 = 3
nu0 = 0
M = 1
R = float(M)/m
G = D_hat

############ FUNCTIONS ##################

def initiate(tmax, h):
    """
    Creates the array for t, theta, omega, phi and nu.
    The array already includes the initial conditions.
    """
    t = np.arange(0,tmax,h)    #t is t_hat (the dimensionless time)
    tsize = len(t)
    theta = np.zeros(tsize)
    omega = np.zeros(tsize)
    phi = np.zeros(tsize) 
    nu = np.zeros(tsize)
    theta[0], omega[0], phi[0], nu[0] = theta0, omega0, phi0, nu0
    return t, tsize, theta, omega, phi, nu

def solveODE_true(D_hat, theta0, omega0, t):
    """
    Solves the ODE (specifically for the single pendulum).
    All workings out done on paper.
    This function finds the values of all constants as well as theta and omega.
    """
    modD = abs(D_hat)
    if modD < 2:
        x = -D_hat/2.
        y = np.sqrt(4.-D_hat**2)/2.
        C = theta0
        D_hat = (omega0-theta0*x)/y
        true_theta = np.exp(x*t)*(C*np.cos(y*t)+D_hat*np.sin(y*t))
        true_omega = np.exp(x*t)*((x*C+y*D_hat)*np.cos(y*t) + (x*D_hat-y*C)*np.sin(y*t))
    elif modD == 2:
        A = theta0
        B = omega0 + (D_hat/2.)*theta0
        true_theta = A*np.exp(-D_hat*t/2.) + B*t*np.exp(-D_hat*t/2.)
        true_omega = (-D_hat/2.)*(A*np.exp(-D_hat*t/2.) + B*t*np.exp(-D_hat*t/2.)) + B*np.exp(-D_hat*t/2.)
    elif modD > 2:
        m1 = (-D_hat+np.sqrt(D_hat**2-4.))/2.
        m2 = (-D_hat-np.sqrt(D_hat**2-4.))/2.
        B = (theta0 - omega0/m1)/(1 - m2/m1)
        A = theta0 - B
        true_theta = A*np.exp(m1*t) + B*np.exp(m2*t)
        true_omega = A*m1*np.exp(m1*t) + B*m2*np.exp(m2*t)
    return true_theta, true_omega

def dwdt(n, fargs, add_omega = 0):
    """
    d(omega)/dt ODE. f = f(t, omega)
    fargs = omega, theta, D_hat
    fargs is a tuple. fargs is alphabetical, with variables first.
    """
    omega, theta, D_hat = fargs
    omega = omega[n] + add_omega  
    theta = theta[n]
    return -D_hat*omega-theta
    
def dwdt_noapprox(n, fargs, add_omega = 0):
    """
    d(omega)/dt ODE without approximating sin(theta) = theta. . f = f(t, omega)
    fargs = omega, theta, D_hat
    fargs is a tuple. fargs is alphabetical, with variables first.
    """
    omega, theta, D_hat = fargs
    omega = omega[n] + add_omega  
    theta = theta[n]
    return -D_hat*omega-np.sin(theta)
    
def dwdt_double(n, fargs, add_omega = 0):
    """
    d(omega)/dt ODE for double pendulum. f = f(t, omega)
    fargs = omega, phi, theta,  R, G
    fargs is a tuple. fargs is alphabetical, with variables first.
    """
    omega, phi, theta, G, R = fargs
    omega = omega[n] + add_omega 
    theta, phi = theta[n], phi[n]
    return -(R+1.)*theta + R*phi - G*omega

def dnudt(n, fargs, add_nu = 0):
    """
    d(nu)/dt ODE for double pendulum. f = f(t, nu)
    fargs = omega, theta, phi, R, G
    fargs is a tuple.
    """
    omega, phi, theta, nu, G, R = fargs
    nu = nu[n] + add_nu
    omega, phi, theta = omega[n], phi[n], theta[n]
    return (R+1.)*theta - (R+1.)*phi + G*(1.-(1./R))*omega - (G/R)*nu
    
def dthetadt(n, fargs, add_theta = 0):
    """
    d(theta)/dt ODE. f = f(t, theta)
    fargs = omega
    fargs is a tuple.
    """
    omega = fargs
    omega = omega[n]
    return omega    

def dphidt(n, fargs, add_phi = 0):
    """
    d(theta)/dt ODE. f = f(t, phi)
    fargs = nu
    nu is a tuple.
    """
    nu = fargs
    nu = nu[n]
    return nu  

def expliciteuler(y, f, fargs, h, n):
    """
    Explicit Euler method
    """
    return y[n] + h*f(n,fargs)
    
def leapfrog(y, f, fargs, h, n):
    """
    Leapfrog method
    """
    if n == 0:
        return expliciteuler(y, f, fargs, h, n)
    return y[n-1] + 2*h*f(n, fargs)

def RK4(y, f, fargs, h, n):
    """
    Runge-Kutta 4 method
    """
    k1 = h*f(n, fargs)
    k2 = h*f(n, fargs, 0.5*k1)
    k3 = h*f(n, fargs, 0.5*k2)
    k4 = h*f(n, fargs, k3)
    y[n+1] = y[n] + (1./6.)*(k1 + 2*k2 + 2*k3 + k4)
    return y[n+1]

def impliciteulertheta(y, f, fargs, h, n):
    """
    Implicit Euler for Theta
    """
    omega = fargs
    y[n+1] = (y[n] + h*omega[n] + h*D_hat*y[n])/(1 + h*D_hat + h**2)
    return y[n+1]

def impliciteuleromega(y, f, fargs, h, n):
    """
    Implicit Euler for Omega
    """
    omega, theta, D_hat = fargs
    y[n+1] = (y[n] - h*theta[n])/(1+h*D_hat+h**2)
    return y[n+1]

def solveODE_single(method, tsize, omega, theta, D_hat, h, approx):
    """
    Uses FDMs to solve the single pendulum ODE.
    If approx = 'Y' or 'y', this will use sin(theta) = theta approximation.
    Type of FDM used will depends on the 'method' variable.
    method can equal to:
    'Explicit Euler'
    'Leapfrog'
    'RK4'
    'Implicit Euler'
    """
    if approx.upper() == 'Y':
        dwdt_func = dwdt
    else:
        dwdt_func = dwdt_noapprox

    for n in range(tsize-1):
        if method == 'Explicit Euler':
            omega[n+1] = expliciteuler(omega, dwdt_func, (omega, theta, D_hat), h, n)
            theta[n+1] = expliciteuler(theta, dthetadt, (omega), h, n)
        elif method == 'Leapfrog':
            omega[n+1] = leapfrog(omega, dwdt_func, (omega, theta, D_hat), h, n)
            theta[n+1] = leapfrog(theta, dthetadt, (omega), h, n)
        elif method == 'RK4':
            omega[n+1] = RK4(omega, dwdt_func, (omega, theta, D_hat), h, n)
            theta[n+1] = RK4(theta, dthetadt, (omega), h, n)
        elif method == 'Implicit Euler':
            omega[n+1] = impliciteuleromega(omega, dwdt_func, (omega, theta, D_hat), h, n)
            theta[n+1] = impliciteulertheta(theta, dthetadt, (omega), h, n)
    return theta, omega

def solveODE_double(method, tsize, omega, phi, theta, nu, G, R, h):
    """
    Uses FDMs to solve the double pendulum ODE.
    Type of FDM used will depends on the 'method' variable.
    method can equal to:
    'Explicit Euler'
    'Leapfrog'
    'RK4'
    """
    for n in range(tsize-1):
        if method == 'Explicit Euler':
            omega[n+1] = expliciteuler(omega, dwdt_double, (omega, phi, theta, G, R), h, n)
            theta[n+1] = expliciteuler(theta, dthetadt, (omega), h, n)
            nu[n+1]    = expliciteuler(nu, dnudt, (omega, phi, theta, nu, G, R), h, n)
            phi[n+1]   = expliciteuler(phi, dphidt, (nu), h, n)
        elif method == 'Leapfrog':
            omega[n+1] = leapfrog(omega, dwdt_double, (omega, phi, theta, G, R), h, n)
            theta[n+1] = leapfrog(theta, dthetadt, (omega), h, n)
            nu[n+1]    = leapfrog(nu, dnudt, (omega, phi, theta, nu, G, R), h, n)
            phi[n+1]   = leapfrog(phi, dphidt, (nu), h, n)
        elif method == 'RK4':
            omega[n+1] = RK4(omega, dwdt_double, (omega, phi, theta, G, R), h, n)
            theta[n+1] = RK4(theta, dthetadt, (omega), h, n)
            nu[n+1]    = RK4(nu, dnudt, (omega, phi, theta, nu, G, R), h, n)
            phi[n+1]   = RK4(phi, dphidt, (nu), h, n)
    return theta, omega, phi, nu
    
def findstepsize(method, criteria, critical_step_start, critical_step_end, number_of_steps, approx):
    """
    It iterates through all of 'h' between the given intervals. Number of interations is defined by 'number_of_steps'.
    It will see if the energy will diverge to infinity. This is done by seeing if the final energy value at tmax is larger than a percentage of the initial value.
    For no damping, there are oscillations in the energy values due to errors in numerical calculations. 
    The criteria is a value which depends on if the system is approximation, damped or undamped. More information at top of the configurations.
    If a critical step size is not found, then it will return False.
    If a critical step size is found, it will return True.
    """
    h = np.linspace(critical_step_start, critical_step_end, number_of_steps)
    finalenergy = []
    for i in range(len(h)):
        t, tsize, theta, omega, phi, nu = initiate(tmax, h[i])
        theta, omega = solveODE_single(method, tsize, omega, theta, D_hat, h[i], approx)
        KE = 0.5*m*(omega*l)**2
        PE = m*g*l*(1-np.cos(theta))
        Etotal = KE + PE
        if D_hat == 0:
            if approx.upper() == 'Y':
                true_theta, true_omega = solveODE_true(D_hat, theta0, omega0, t)
                theta_diff = abs(true_theta - theta)
                if max(theta_diff) <= criteria:
                    return findstepsize_return(h, i, t)
            else:
                oscillation_criteria = Etotal[0]*criteria/100.
                average_energy = np.mean(Etotal)
                if abs(average_energy - Etotal[0]) <= oscillation_criteria:
                    print 'Average energy = ', average_energy
                    return findstepsize_return(h, i, t)
        else:
            energy_criteria = Etotal[0]*criteria/100.
            if Etotal[-1] <= energy_criteria:
                return findstepsize_return(h, i, t)
        finalenergy += Etotal[-1]
    return h[-1], False, t

def findstepsize_return(h, i, t):
    """
    Returns the critical step-size and the final energy value and whether the critical step size is the initial step size.
    True means the critical step size is not the intial step size.
    """
    if i == 0:
        return h[i], False, t
    return h[i], True, t

########### SCRIPT #####################
"""
This will find the critical step size by using the 'findstepsize' function defined above.
If findstepsize function returns False, then max/min of 'h' must be changed to find the critical step size.
"""
if find_critical_step_length.upper() == 'Y':
    if approx.upper() == 'Y' and D_hat == 0:
        criteria = theta_difference
    else:
        criteria = energy_percentage
    h, true_or_false, t = findstepsize(method, criteria, critical_step_start, critical_step_end, number_of_steps, approx)
    if true_or_false == True:
        print 'Critical step size = ', h
    elif true_or_false == False:
        print 'Critical step size = ', h
        print 'Initial step-size = ', critical_step_start
        print 'Final step-size = ', critical_step_end
        print 'The initial/final step-size is not high/low enough.'
    
if Type == 'Single Pendulum':
    t, tsize, theta, omega, phi, nu = initiate(tmax, h)
    #Calculates the actual solution to ODE (with sin(theta) = theta approximation)    
    true_theta, true_omega = solveODE_true(D_hat, theta0, omega0, t)
    
    #Calculates the FDM solution to ODE    
    theta, omega = solveODE_single(method, tsize, omega, theta, D_hat, h, approx)
   
    thetalist = [[],[],[],[]]
    omegalist = [[],[],[],[]]
    KElist = [[],[],[],[]]
    PElist = [[],[],[],[]]
    Etotallist = [[],[],[],[]]
    for i in range(4):
        t, tsize, thetalist[i], omegalist[i], phi, nu = initiate(tmax, h)
        thetalist[i], omegalist[i] = solveODE_single(method_list[i], tsize, omegalist[i], thetalist[i], D_hat, h, approx)
        KElist[i] = 0.5*m*(omegalist[i]*l)**2
        PElist[i] = m*g*l*(1-np.cos(thetalist[i]))
        Etotallist[i] = KElist[i]+ PElist[i]
    
    theta_diff = abs(true_theta - theta)
    
    #Energy
    KE = 0.5*m*(omega*l)**2
    PE = m*g*l*(1-np.cos(theta))
    Etotal = KE + PE
    
    KE_true = 0.5*m*(true_omega*l)**2
    PE_true = m*g*l*(1-np.cos(true_theta))
    Etotal_true = KE_true + PE_true
    
elif Type == 'Double Pendulum':
    #Calculates the FDM solution to ODE
    t, tsize, theta, omega, phi, nu = initiate(tmax, h)
    theta, omega, phi, nu = solveODE_double(method, tsize, omega, phi, theta, nu, G, R, h)
    
    KE = 0.5*m*(omega*l)**2 + 0.5*M*(nu*l+omega*l)**2
    PE = m*g*l*(1-np.cos(theta)) + M*g*l*(1 - np.cos(theta) + 1 - np.cos(phi))
    Etotal = KE + PE
    
#Plots the difference between theta and the analytically solved theta.
if Type == 'Single Pendulum' and approx.upper() == 'Y':
    fig2 = plt.figure('Secondary', figsize=(8,6))
    plt.plot(t, theta_diff, 'g-')
    plt.title('Theta difference between solution to approx ODE and using FDM')
    plt.ylabel('Theta Difference')
    plt.xlabel('Time')
    
fig = plt.figure("Main", figsize=(13,10))
"""
#Plots theta and the analytical solution (and phi if double pendulum) aginst time
subplot1 = fig.add_subplot(2,2,1)
thetaplot, = plt.plot(t, theta, 'g-', label="Theta")
if Type == 'Single Pendulum':
    plot2, = plt.plot(t, true_theta, 'k-', label="Solution to approx ODE")
elif Type == 'Double Pendulum':
    plot2, = plt.plot(t, phi, 'r-', label="Phi")
plt.title(method + ': h = %.4g' %h)
plt.xlabel('Time')
plt.ylabel('Angle')
plt.legend(handles=[thetaplot, plot2])

#Plots omega (and nu if double pendulum)  against time.
subplot2 = fig.add_subplot(2,2,2)
omegaplot, = plt.plot(t, omega, 'g-', label="Omega")
if Type == 'Single Pendulum':
    plot3, = plt.plot(t, true_omega, 'k-', label="Solution to approx ODE")
elif Type == 'Double Pendulum':
    plot3, = plt.plot(t, nu, 'r-', label="Nu")
plt.title("Angular Velocities")
plt.xlabel('Time')
plt.ylabel('Angular Velocity')
plt.legend(handles=[omegaplot, plot3])
"""
#Plots total energy against time.
colours = ['g-', 'r-', 'b-', 'k-']
energyplots = []
for i in range(4):
    energyplots += plt.plot(t, Etotallist[i], colours[i], linewidth = 2, label=method_list[i])
plt.legend(handles=energyplots, loc='Upper Left')

#plt.plot(t, Etotal)
plt.title('Total energy against time')
plt.xlabel('Time')
plt.ylabel('Energy (J)')

plt.show()