# Code to solve the diffusion equation to represent the variation of wind speed
# with height in the atmospheric boundary layer

import numpy as np
import matplotlib.pyplot as plt

def FTCS_diffuse(phi, d, nt, dtForcing):
    '''Solve the one dimensional diffusion coeffiecient with the FTCS finite
    difference scheme starting from initial conditions in array phi.
    d is the non-dimensional diffusion coefficient: d = dt K/dx^2
    dtForcing is the time step multiplied by an additional term on the right
    hand side of the equation.
    The start boundary condition is zero and the end boundary condition is
    zero gradient.
    phi after nt time steps is returned.'''

    # Create an array for phi at the next time step
    phiNew = phi.copy()
    
    # Short cut to the length of the array phi
    nx = len(phi)
    
    # Loop through all time steps. phiNew is always phi at time step n+1
    # and phi is at time step n. 
    for it in range(nt):
        # Loop over space away from the boundaries
        for i in range(1,nx-1):
            phiNew[i] = phi[i]

        # Update the boundary conditions
        ...
        
        # Update phi for the next time step
        ...
    
    return phi

def FTCS_turbulentDiffuse(u, z, Kmin, nt, dt, forcing):
    '''Solve the one dimensional diffusion coeffiecient with the FTCS finite
    difference scheme starting from initial conditions in array u (wind speed).
    z is height above the ground.
    Kmin is the minimum diffusion coefficient.
    dt is the time step.
    forcing is the additional term on the right hand side of the equation.
    The start boundary condition is zero and the end boundary condition is
    zero gradient.
    u after nt time steps is returned.
    The diffusion coefficient is calculated as
    K = L^2 |dudz|
    where the length scale L = 0.4z'''

    # von-Karmen's constant
    kappa = 0.4

    # Calculate dz and declare an array for the height of the K locations
    dz = z[1] - z[0]
    zK = np.arange(0.5*dz, z[-1], dz)
    
    # Declare array for the diffusion coefficients, K and the wind shear
    K = np.zeros_like(zK)
    dudz = np.zeros_like(zK)
    
    # Loop over all time steps
    for it in range(nt):
        # Calculate dudz and K for each level
        for k in range(len(zK)):
            dudz[k] = (u[k+1] - u[k])/dz
            length = kappa*zK[k]
            K[k] = length**2*abs(dudz[k])
    
        # Update u for each internal level based on dudz and K
        for k in range(1,len(u)-1):
            u[k] = u[k] + ..
    
        # Update the boundary conditions
    
    return u

def diffuse():
    # Problem setup, solution and plotting
    
    # Grid and time
    dz = 5
    zTop = 100.
    z = np.arange(0,zTop+1, dz)
    dt = 0.5
    Tend = 3600
    nTimes = int(np.round(Tend/dt))

    # Constant pressure gradient
    dpdx = -5e-3

    # Diffusion coefficient
    K = 2.5
    # Minimum diffusion coefficient for turbulent simulation 
    Kmin = 1e-6
    
    # Non-dimensional diffusion coefficient
    d = K*dt/dz**2
    
    print('Solving diffusion equation for ', nTimes, \
          ' time steps with non-dimensional diffusion coefficient ', d)

    # Initial wind
    uInit = 10.
    u = uInit*np.ones_like(z)

    # Solution
    u_FTCS = FTCS_diffuse(u.copy(), d, nTimes, -dt*dpdx)
    u_turb = FTCS_turbulentDiffuse(u, z, Kmin, nTimes, dt, -dpdx)

    # Plot with height
    plt.plot(u_FTCS, z, 'b-', label='u fixed K')
    plt.plot(u_turb, z, 'r-', label='u turbulent')
    plt.xlabel('Wind speed [m/s]')
    plt.ylabel('Height [m]')
    plt.xlim([0,10])
    plt.ylim([0,100])
    plt.legend(loc='best')
    plt.show()

diffuse()
