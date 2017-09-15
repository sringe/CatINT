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
import sys

class Transport:
#just a few variables which can be hopefully later taken from catmap


#self.diffusion_dict{['

    def __init__(self):

        #THE MESH
        self.dt=5.e-2
        self.dx=5.e-2 #*1e-10
        self.tmax=10
        self.xmax=10 #*1e-10
        self.tmesh=np.arange(0,self.tmax+self.dt,self.dt)
        self.xmesh=np.arange(0,self.xmax+self.dx,self.dx)

        #A FEW VARIABLES
        self.bulk_concentrations=np.array([0.1,0.1])
        self.dOHP = 1.0
        self.eps = 80.0*unit_eps0
        #eps=80*unit_eps0
        self.charges=np.array([1,-1])*unit_F
        self.T = 300
        self.u = np.array([0.0,0.0])
        self.D = np.array([0.276,0.5])/100**2 #in cm2/s
        self.nspecies=len(self.D)
        self.beta = 1./(self.T * unit_R)
        self.external_charge=np.zeros([len(self.xmesh)])
        self.external_charge=self.gaussian(sigma=0.3,z=1.0,mu=3.0,cmax=0.1) 


        #BOUNDARY AND INITIAL CONDITIONS
        self.set_initial_conditions(
                c_initial_general={'all':0.1},
                c_initial_specific={'0':{'0':0.5},'1':{'0':0.3}})

        self.set_boundary_conditions(\
                dc_dt_boundary={'all':{'all':0.0}},     #in mol/l/s
                efield_boundary={'r':0.0})    #in V/Ang

    def run(self):

        self.cout=self.integrate_pnp(self.dx,len(self.xmesh),self.dt,\
                len(self.tmesh),self.D,10)
       
    def plot(self):

        colorlist=cycle(['b','k'])

        ax1=plt.subplot('221')
        ax2=plt.subplot('222')
        ax3=plt.subplot('223')
        ax4=plt.subplot('224')
    
        color_offset=0.3
        for k in range(0,self.nspecies):
            color=next(colorlist)
            i=-1
            for c in self.cout:
                i+=1
                if i!=len(self.cout)-1:
                    lw=0.5
                    zorder=100
                else:
                    lw=2.0
                    zorder=0
                brightness=1.-(color_offset+(i*1.)/len(self.cout)*(1.-color_offset))
                if k==0:
                    color=str(brightness)
                elif k==1:
                    color=(brightness,1.,1.)
                ax1.plot(self.xmesh[:-1],c[k*len(self.xmesh):k*len(self.xmesh)+len(self.xmesh)-1] /10**3,'-',color=color,linewidth=lw,zorder=zorder)
        ax1.legend()
        ax1.set_xlabel('x (m)')
        ax1.set_ylabel('c (mol/L)')
        #c=self.integrate_FTCS(self.dt,self.dx)
        #for t in np.arange(0.0,1.,0.1):
        #    plt.plot(self.xmesh,c[int(t/self.dt),:],'-o',label=str(t))
        ax2.plot(self.xmesh[:-1],self.efield[:-1]/1e10,'-')
        ax2.set_title('Electric field')
        ax2.set_xlabel('x (m)')
        ax2.set_ylabel('E (V/Ang)')

        ax3.set_title('Potential')
        self.potential=np.zeros([len(self.xmesh)])
        integral=0.0
        for i in range(len(self.xmesh)):
            integral-=self.efield[i]*self.dx
            self.potential[i]+=integral
        ax3.plot(self.xmesh[:-1],self.potential[:-1],'-')
        ax3.set_ylabel('v (V)')
        ax3.set_xlabel('x (m)')
        ax4.plot(self.xmesh[:-1],self.external_charge[:-1]/10**3/unit_F, '-',label='n_ext')
        ax4.plot(self.xmesh[:-1],self.total_charge[:-1]/10**3/unit_F,'-',label='n_tot')
        ax4.set_ylabel('n_ion (e*mol/L)')
        ax4.legend()
        ax4.set_title('Charge Density')
        plt.tight_layout()
        plt.show()
    #    self.integrate_pb()

    def get_initial_conditions(self):
        return self.c0 

    def set_initial_conditions(self,\
        c_initial_general={},\
        c_initial_specific={}):
        """accepts the initial concentrations of all species in the c_general list (all in mol/L), 
        the c_specific dictionary allows to parse specific values of the concentrations at 
        specific grid points (1st key is species index, 2nd is xmesh index, 
        value is concentration)."""
        c0=np.zeros([self.nspecies*len(self.xmesh)])
        j=-1

        for k in range(self.nspecies):
            for i in range(len(self.xmesh)):
                j+=1
                if str(k) in c_initial_general:
                    c0[j]=c_initial_general[k]
                elif 'all' in c_initial_general:
                    c0[j]=c_initial_general['all']
                if (str(k) in c_initial_specific and str(i) in c_initial_specific[str(k)]):
                    c0[j]=c_initial_specific[str(k)][str(i)]
                if ('all' in c_initial_specific and str(i) in c_initial_specific['all']):
                    c0[j]=c_initial_specific['all'][str(i)]

        c0=np.array(c0)
#        c0=[0.1+0.1*self.xmesh]+[0.1+0.1*self.xmesh]
#        c0=np.array([item for sublist in c0 for item in sublist])
        c0*=10**3 #multiply by 1000 to get m^-3
        self.c0=c0


    def get_boundary_conditions(self):
        return self.c_bound,self.dc_dt_bound,self.efield_bound

    def set_boundary_conditions(self,\
        dc_dt_boundary={},\
        efield_boundary={}):

        """allows to set boundary conditions in the concentrations themselves (c_boundary) 
        and the derivatives dc_dx_boundary and the efield. 
        species, direction, value"""

        dc_dt_bound=np.zeros([self.nspecies,2])
        efield_bound=np.array([None,None])

        def key_to_index(key):
            if key=='l':
                index=0
            elif key=='r':
                index=1
            return index

        def species_to_index(key):
            if key.isdigit():
                index=int(key)
            else:
                index=self.nspecies+1

        for key in efield_boundary:
            efield_bound[key_to_index(key)]=efield_boundary[key]

        for key1_raw in dc_dt_boundary:
            if key1_raw=='all':
                keys1=map(str,range(self.nspecies))
            else:
                keys1=[str(key1_raw)]
            for key1 in keys1:
                for key2_raw in dc_dt_boundary[key1_raw]:
                    if key2_raw=='all':
                        keys2=map(str,range(2))
                    else:
                        keys2=[str(key2_raw)]
                    for key2 in keys2:
                        dc_dt_bound[int(key1),int(key2)]=dc_dt_boundary[key1_raw][key2_raw]

        for l in [dc_dt_bound,efield_bound]:
            l=np.array(l)

        self.dc_dt_bound=dc_dt_bound*10**3
        self.efield_bound=[]
        for e in efield_bound:
            if e==None:
                self.efield_bound+=[None]
            else:
                self.efield_bound+=[e*1e10]
        self.efield_bound=np.array(self.efield_bound)

    def gaussian(self,sigma=0.01,z=1,mu=0.2,cmax=1):
        """define gaussian charge density. c is in molar"""
#        sigma=1./(np.sqrt(2.*np.pi))*1/cmax
        gaussian=np.zeros([len(self.xmesh)])
        for i in range(0,len(self.xmesh)):
            gaussian[i]=\
                    1./(sigma*np.sqrt(2.*np.pi))\
                    *np.exp(-0.5*((self.xmesh[i]-mu)/sigma)**2)
        gaussian=gaussian/max(gaussian)*cmax*10**3*z*unit_F
        return gaussian

    def integrate_pnp(self,dx,nx,dt,nt,D,ntout):

        def ode_func(c,t):

            dc_dt = np.zeros([nx*self.nspecies])

            #first determine E = int E' = 1/eps int sum c_i z_i by numerical integration over the grid
            self.efield=np.zeros([nx])
            self.defield_dx=np.zeros([nx])

            self.total_charge=np.zeros([nx])
            #INTEGRATION OF PBE: get electric field
            if self.efield_bound[1] !=None:
                #integration from the right
                ll=reversed(range(0,nx))
            elif self.efield_bound[0] !=None:
                #integration from the left
                ll=range(0,nx)

            self.total_concentrations=np.zeros([self.nspecies+1,nx])

            for k in range(0,self.nspecies):
                total_int=0.0
                m=-1
                for i in ll:
                    m+=1
                    j=k*nx+i
                    current_charge=self.charges[k]*c[j]
                    if k==0:
                        current_charge+=self.external_charge[i]
                    total_int+=current_charge*dx/self.eps
                    self.efield[i]+=total_int
                    self.defield_dx[i]+=current_charge/self.eps
                    #others:
                    self.total_charge[i]+=current_charge
                    self.total_concentrations[k,i]+=c[j]

            ##integrate charge with Gaussian quadrature, default order=5
            #for k in range(0,self.nspecies):
            #    func=self.charges[k]*c[k*nx:(k+1)*nx]
            #    #integrate.fixed_quad(func,min(func),max(func))
            #    integral+=integrate.simps(func,self.xmesh)

            self.total_concentrations[-1,:]=self.external_charge/unit_F

            #plt.plot(self.xmesh,self.total_charge/10**3/unit_F)
            #for t in self.total_concentrations:
            #    plt.plot(self.xmesh,t/10**3)
            #plt.show()
            #sys.exit()

            if self.efield_bound[1] != None:
                #we integrated from right to left, so we have to switch the sign
                #of the efield here
                self.efield = -self.efield
            #constant shift of electric field due to boundary conditions
            self.efield += [e for e in self.efield_bound if e!=None][0]


            #IMPLEMENTATION OF FLOWS
            #go over all point in space and all species

            for k in range(0,self.nspecies):
                m=-1
                for i in range(0,nx):
                    m+=1
                    j=k*nx+i

                    if m==0:
                        dc_dt[j] = self.dc_dt_bound[k,0]
                        continue
                    elif m==nx-1:
                        dc_dt[j] = self.dc_dt_bound[k,1]
                        continue

                    dc_dx=(c[j+1]-c[j-1])/(2.*dx)
                    dc_dx_2=(c[j-1]-2*c[j]+c[j+1])/dx**2

                    dc_dt[j] =\
                        (\
                        -self.u[k]*dc_dx \
                        +self.D[k]*(\
                            dc_dx_2 \
                            +self.charges[k]*self.beta*dc_dx*self.efield[i]\
                            +self.charges[k]*self.beta*c[j]*self.defield_dx[i]\
                        )\
                        )
            return dc_dt

        sol = integrate.odeint(ode_func, self.c0, range(0,nt-1)) #, args=(b, c))

        #return the results
        cout = []
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
         
