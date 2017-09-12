"""
transport class
---
Stefan Ringe (stefan.ringe.tum@gmail.com)
"""
import scipy.integrate as integrate
import numpy as np
import matplotlib.pyplot as plt
#from ase import units
from scipy.sparse import diags
from units import *
from itertools import cycle


class Transport:
#just a few variables which can be hopefully later taken from catmap


#self.diffusion_dict{['

    def __init__(self):
        
        self.dt=2.e-2
        self.dx=2.e-2
        self.tmax=10
        self.xmax=0.5
        self.tmesh=np.arange(0,self.tmax+self.dt,self.dt)
        self.xmesh=np.arange(0,self.xmax+self.dx,self.dx)


        self.bulk_concentrations=np.array([0.1,0.1])
        self.dOHP = 1.0
        self.eps = 80.0
        self.charges=np.array([1,-1])
        self.T = 300
        self.u = [0.0,0.0] #001]
        self.D = np.array([0.276,0.5])/100**2 #in cm2/s
        self.nspecies=len(self.D)
        
        c0 = np.zeros((len(self.tmesh),)) # initial condition
        c0 = 10 * 10**3 #in mol/L

        colorlist=cycle(['b','k'])

        #cout,s=self.diffusion_Crank_Nicolson(self.dx,len(self.xmesh),self.dt,\
        #        len(self.tmesh),self.D,c0,10)
        cout=self.diffusion_FTCS_odeint(self.dx,len(self.xmesh),self.dt,\
                len(self.tmesh),self.D,c0,10)
        color_offset=0.1
        for k in range(0,self.nspecies):
            color=next(colorlist)
            i=-1
            for c in cout:
                i+=1
                brightness=1.-(color_offset+(i*1.)/len(cout)*(1.-color_offset))
                if k==0:
                    color=str(brightness)
                elif k==1:
                    color=(brightness,1.,1.)
                plt.plot(self.xmesh[:-1],c[k*len(self.xmesh):k*len(self.xmesh)+len(self.xmesh)-1] /10**3,'-',color=color,linewidth=1.5)
        #c=self.integrate_FTCS(self.dt,self.dx)
        #for t in np.arange(0.0,1.,0.1):
        #    plt.plot(self.xmesh,c[int(t/self.dt),:],'-o',label=str(t))
        plt.legend()
        plt.xlabel('x')
        plt.ylabel('c')
        plt.show()
    #    self.integrate_pb()


    def diffusion_FTCS_odeint(self,dx,nx,dt,nt,D,c0,ntout):

        self.bc_outputted=False

        def ode_func(c,t):
            dc_dt = np.zeros([nx*self.nspecies])

            #boundary condition for dc/dt
            bc_kind='neumann'
            bc_pos=self.xmesh[-2] #x-position of boundary condition
            bc_val=[0.0,0.0] #value in mol/L /(s or m)

            z=1*unit_e
            eps=80*unit_eps0
            beta = self.T * unit_R
            T=300
            j=-2
            for k in range(0,self.nspecies):
                j+=1
                for i in range(0,nx-1): # space
                    j+=1
                    dc_dx=(c[j+1]-c[j-1])/(2.*dx)
                    dc_dx_2=(c[j-1]-2*c[j]+c[j+1])/dx**2
                    if i==int(((bc_pos-min(self.xmesh))/dx)):
                        if not self.bc_outputted:
                            print 'species ',k,': Applying bc dc/dt(x=',round(bc_pos,3),',t) = ', bc_val, ' (mol/L)/s'
                            if k==self.nspecies-1:
                                self.bc_outputted=True
                        if bc_kind=='neumann':
                            #von-Neumann boundary condition in t-dimension
                            dc_dt[j] =bc_val[k] *10**3
                            continue
                    dc_dt[j] =\
                            (\
                            -self.u[k]*dc_dx \
                            +self.D[k]*(\
                                dc_dx_2# \
                        #        -dc_dx**2/c[j]\
                        #        +z*beta*c[j]/eps*z*c[j]\
                            )\
                            )
            return dc_dt

        #boundary condition for c
        cout=[]
        bc_kind='dirichlet'
        bc_pos=self.xmesh[-2]
        bc_val=[10.0,10.0]
        c0=np.array([0.0]*nx*self.nspecies) #value at t=0 for all x-values
        for k in range(0,self.nspecies):
            c0[k*nx+int(((bc_pos-min(self.xmesh))/dx))]=bc_val[k]*10**3 #value at t=0 and x=0
        print 'Applying bc c(x=',round(bc_pos,3),',t=0) = ',bc_val, 'mol/L (zero for all other x)'
        #solve time problem
        sol = integrate.odeint(ode_func, c0, range(0,nt-1)) #, args=(b, c))

        #return the results
        for n in range(0,nt-1):
            for j in range(1,nx-1):
                if n % int(nt/float(ntout)) == 0 or n==nt-1:
                    cout.append(sol[n,:].copy()) # numpy arrays are mutable, 
        return cout

    def diffusion_FTCS(self,dx,nx,dt,nt,D,c0,ntout):
        # diffusion number (has to be less than 0.5 for the
        # solution to be stable):
        cout = []
        s = D*dt/dx**2
        c = np.zeros([nt,nx])
        c[:,0] = c0
        u=0.0 #1e-4
        z=1*unit_e
        eps=80*unit_eps0
        beta = self.T * unit_R
        T=300
        for n in range(0,nt-1): # time
            for j in range(1,nx-1): # space
                dc_dx=(c[n,j+1]-c[n,j-1])/(2.*dx)
                dc_dx_2=(c[n,j-1]-2*c[n,j]+c[n,j+1])/dx**2
                c[n+1,j] = c[n,j] + dt*\
                        (\
                        -u*dc_dx \
                        +D*(\
                            dc_dx_2 \
                            -dc_dx**2/c[n,j]\
                            +z*beta*c[n,j]/eps*z*c[n,j]\
                        )\
                        )
                if n % int(nt/float(ntout)) == 0 or n==nt-1:
                    cout.append(c[n,:].copy()) # numpy arrays are mutable, 
        return cout,s
    # note that this can be written without the for-loop
    # in space, but it is easier to read it this way

    def diffusion_Crank_Nicolson(self,dx,nx,dt,nt,D,c,ntout):
        cout = [] # list for storing c arrays at certain time steps
        c0 = c[0] # boundary condition on left side
        c1 = c[-1] # boundary condition on right side
        s = D*dt/dx**2  # diffusion number
        # create coefficient matrix:
        A = diags([-0.5*s, 1+s, -0.5*s], [-1, 0, 1], 
              shape=(nx-2, nx-2)).toarray() 
        B1 = diags([0.5*s, 1-s, 0.5*s],[-1, 0, 1], shape=(nx-2, nx-2)).toarray()
        for n in range(1,nt): # time is going from second time step to last
            cn = c
            B = np.dot(cn[1:-1],B1) 
            B[0] = B[0]+0.5*s*(c0+c0)
            B[-1] = B[-1]+0.5*s*(c1+c1)
            c[1:-1] = np.linalg.solve(A,B)
            if n % int(nt/float(ntout)) == 0 or n==nt-1:
                cout.append(c.copy()) # numpy arrays are mutable, 
                #so we need to write out a copy of c, not c itself
        return cout,s


    def integrate_FTCS(self,dt,dx):
        c=np.zeros([len(self.tmesh),len(self.xmesh)])
        #c[0,:] = 1.0 #initial concentrations
        c[:,0] = 10.0 #fix concentrations at x=0
#        dc_dv=self.bulk_concentrations * (-self.charges/units.kB/T) * np.exp(-self.charges * potential/units.kB/T)
        for it in range(0,len(self.tmesh)-1):
            for ix in range(1,len(self.xmesh)-1):
                c[it+1,ix] = c[it,ix] \
                    + self.D * (c[it,ix+1]-2*c[it,ix]+c[it,ix-1]) * dt/dx**2 #\
                # - 0.5 * self.u * (c[it,ix+1]-c[it,ix-1]) * dt/dx \
                  #      + self.D*self.charges[k]/units.kB/self.T * (c[t,x+1]-c[t,x-1])**2 / dc_dv * dt/dx**2
#c[t+dt,x] = c[t,x] + 1/2 * u * (c[t,x+dx]-c[t,x-dx]) * dt/dx - D * (c[t,x+dx]-2c[t,x]+c[t,x-dx]) * dt/dx^2 
#	+ const. * (c[t,x+dx]-c[t,x-dx])^2/dcv * dt/(dx)^2  + const * c[t,x]^2 = R
        return c
#    def update():
#
#        self.rates = []
#        self.rates += self.migration_rates
#        self.rates += self.diffusion_rates
#        self.rates += self.convection_rates
#        self.rates += self.reaction_rates

    #def get_rates(self):
    #        for i in self.xmesh:
    #            rates[i] = self.D*self.z/unit.kB/self.T*concentraitions*
    #    return rates

    def ionic_charge(self,pot):
        charge=0.0
        for j,z in enumerate(self.charges):
            charge+=z*self.bulk_concentrations[j]*\
            np.exp(-self.charges[j]*pot/(units.kB * self.T))
        return charge

    def integrate_pbe(self):
        """This function integrates the 1d PBE and solves for the electric field:
            int (sum_i z_i q c_i) dz
            solution is performed with finite differences"""

        def d(z, t):
            #z = [v',v]
            return np.array((
                -1./self.eps*self.ionic_charge(z[1]),       #this is z[0]'
                z[0]                                        #this is z[1]'
                ))

        sol = integrate.odeint(d,[0.,0.],self.xmesh)
        
        plt.plot(self.xmesh, sol[:, 0], 'b', label='v(x)')
        plt.plot(self.xmesh, sol[:, 1], 'g', label='E(x)')
        plt.legend(loc='best')
        plt.xlabel('t')
        plt.ylim([-2,2])

        plt.grid()
        plt.show()



#    def read_catmap():
#        """Purpose: read the catmap input file"""
         
