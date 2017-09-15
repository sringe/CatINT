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

    def __init__(self,integrator='FTCS-odeint'):

        #THE MESH
        self.dt=5.e-3
        self.dx=5.e-2 #*1e-10
        self.tmax=10
        self.xmax=10 #*1e-10
        self.tmesh=np.arange(0,self.tmax+self.dt,self.dt)
        self.xmesh=np.arange(0,self.xmax+self.dx,self.dx)

        #A FEW VARIABLES
        self.integrator=integrator
        self.bulk_concentrations=np.array([0.1,0.1])
        self.dOHP = 1.0
        self.eps = 80.0*unit_eps0
        #eps=80*unit_eps0
        self.charges=np.array([1,-1])*unit_F
        self.T = 300
        self.u = np.array([0.0,0.0])
        self.D = np.array([0.276])/100**2 #in cm2/s
        self.nspecies=len(self.D)
        self.beta = 1./(self.T * unit_R)
        self.external_charge=np.zeros([len(self.xmesh)])
        self.external_charge=self.gaussian(sigma=0.3,z=1.0,mu=3.0,cmax=0.1) 
        self.count=1
        self.ax1=plt.subplot('211')
        self.ax2=plt.subplot('212')

        #BOUNDARY AND INITIAL CONDITIONS
        self.set_initial_conditions(
                c_initial_general={'all':0.1},
                c_initial_specific={'0':{'0':0.5}})#,'1':{'0':0.3}})

        self.set_boundary_conditions(\
                dc_dt_boundary={'all':{'all':0.0}},     #in mol/l/s
                efield_boundary={'r':0.0})    #in V/Ang

    def run(self):

        self.cout=self.integrate_pnp(self.dx,len(self.xmesh),self.dt,\
                len(self.tmesh),10,method=self.integrator)
       
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

    def integrate_pnp(self,dx,nx,dt,nt,ntout,method='FTCS'):

        def ode_func(c,t):

            dc_dt = np.zeros([nx*self.nspecies])

            #first determine E = int E' = 1/eps int sum c_i z_i by numerical integration over the grid
            self.efield=np.zeros([nx])
            self.defield_dx=np.zeros([nx])

            self.total_charge=np.zeros([nx])

            self.total_concentrations=np.zeros([self.nspecies+1,nx])


#            #position of test charge
#            for i in range(0,nx):
#                #integrate over grid
#                for ii in range(0,nx):
#                    for k in range(0,self.nspecies):
#                        if ii==i:
#                            continue
#                        j=k*nx+ii
#                        dist_sq=(self.xmesh[ii]-self.xmesh[i])**2
#                        current_charge=self.charges[k]*c[j]
#                        if k==0:
#                            current_charge+=self.external_charge[j]
#                        self.efield[i]+=current_charge*dx/self.eps/dist_sq
#            
#            for k in range(0,self.nspecies):
#                for i in range(0,nx):
#                    j=k*nx+i
#                    current_charge=self.charges[k]*c[j]
#                    if k==0:
#                        current_charge+=self.external_charge[j]
#                    self.defield_dx[i]+=current_charge/self.eps

            for k in range(self.nspecies):
                #INTEGRATION OF PBE: get electric field
                if self.efield_bound[1] !=None:
                    #integration from the right
                    ll=reversed(range(0,nx))
                elif self.efield_bound[0] !=None:
                    #integration from the left
                    ll=range(0,nx)
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

#            plt.plot(self.xmesh,(self.charges[0]*c[:len(self.xmesh)]+self.charges[1]*c[len(self.xmesh):])/self.eps,'o')
#            plt.plot(self.xmesh,self.external_charge/self.eps,'o')
##            plt.plot(self.xmesh,self.efield,label='e')
#            plt.plot(self.xmesh,self.defield_dx,label='de/dx')
#            plt.legend()

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
        #                    +self.charges[k]*self.beta*dc_dx*self.efield[i]\
        #                    +self.charges[k]*self.beta*c[j]*self.defield_dx[i]\
                        )\
                        )

         #   if self.count in range(50,100):
         #       self.ax1.plot(self.xmesh,dc_dt[:len(self.xmesh)]/10**3*dt,'-',linewidth=0.3,color='r')
         #       self.ax1.plot(self.xmesh,dc_dt[len(self.xmesh):]/10**3*dt,'-',linewidth=0.3,color='k')
         #       self.ax2.plot(self.xmesh,c[:len(self.xmesh)]/10**3,'-',linewidth=0.3,color='r')
         #       self.ax2.plot(self.xmesh,c[len(self.xmesh):]/10**3,'-',linewidth=0.3,color='k')
         #   self.count+=1
         #   if self.count==100:
         #       plt.show()
         #       sys.exit()
            return dc_dt

        def integrate_FTCS_odeint(dx,nx,dt,nt,c0,ntout):
            print np.shape(c0)
            sol = integrate.odeint(ode_func, c0, range(0,nt)) #, args=(b, c))
            print np.shape(sol)
            #return the results
            cout = []
            for n in range(0,nt):
                if n % int(nt/float(ntout)) == 0 or n==nt-1:
                    cout.append(sol[n,:].copy()) # numpy arrays are mutable, 
            return cout

        def integrate_FTCS(dt,dx,nt,nx,V0):
            # diffusion number (has to be less than 0.5 for the 
            # solution to be stable):
            s = self.D[0]*dt/dx**2
            V = np.zeros([nt,nx])
            V[:,0] = [V0[0]]*(nt)       #boundary left  (dc_dt=0)
            V[:,-1] = [V0[-1]]*(nt)     #boundary right (dc_dt=0)
            V[0,:]=V0                   #initial
            for n in range(0,nt-1): # time
                for j in range(1,nx-1): # space
                    V[n+1,j] = V[n,j] + s*(V[n,j-1] -
                        2*V[n,j] + V[n,j+1]) 
            cout=V[-1,:]
            return cout

        def integrate_Crank_Nicolson(dx,nx,dt,nt,c,ntout):
            """only for dc_dt=0 at both boundary sites implemented yet"""
            cout = [] # list for storing c arrays at certain time steps
            c0 = c[0] # boundary condition on left side
            c1 = c[-1] # boundary condition on right side
            s = self.D[0]*dt/dx**2  # diffusion number
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

        if method=='FTCS-odeint':
            cout=integrate_FTCS_odeint(dx,nx,dt,nt,self.c0,ntout)
            print np.shape(cout)
            dataplot=cout[-1][:nx]
        elif method=='FTCS':
            cout=integrate_FTCS(dt,dx,nt,nx,self.c0[:nx])
            dataplot=cout
        elif method=='Crank-Nicolson':
            cout,s=integrate_Crank_Nicolson(dx,nx,dt,nt,self.c0[:nx],ntout)
            dataplot=cout[-1]
        plt.plot(self.xmesh,dataplot/10**3,'-')
        plt.show()
        sys.exit()
        return cout

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



#    def read_catmap():
#        """Purpose: read the catmap input file"""
         
