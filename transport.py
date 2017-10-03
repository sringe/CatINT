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
from scipy.linalg import block_diag
from units import *
from itertools import cycle
import sys
from copy import deepcopy
from scipy.optimize import minimize
from scipy.integrate import ode

class Transport:
#just a few variables which can be hopefully later taken from catmap


#self.diffusion_dict{['

    def __init__(self,integrator='FTCS-odeint'):

        self.eps = 80.0*unit_eps0 #1.1e11 #*unit_eps0
        self.T = 300
        self.beta = 1./(self.T * unit_R)
        self.vzeta=-0.025 #Volt
        self.charges=np.array([-1,1])*unit_F
        self.bulk_concentrations=np.array([0.001,0.001])*10**3 #in mol/m^3
        self.debye_length = np.sqrt( self.eps/self.beta / sum(self.charges**2*self.bulk_concentrations)) #in m

        self.u = np.array([0.0,0.0]) #
        self.D = np.array([1.96e-9,1.2e-9])  #diffusion coefficient in m^2/s
        self.mu = self.D * self.charges *self.beta  #ion mobilities according to Einstein relation

        #THE MESH
        self.dt=1e-12 #1e-11 #5e-6 #1.0 #5.e-3
        self.tmax=1e-10#1e-8 #100e-3 #2024.2369851 #10
        self.tmesh=np.arange(0,self.tmax+self.dt,self.dt)
        self.dx=self.debye_length/10.
        #self.xmax=80.e-6 #*1e-10
        self.xmax=10*self.debye_length
        self.xmesh=np.arange(0,self.xmax+self.dx,self.dx)

        #A FEW VARIABLES
        self.integrator=integrator
        #eps=80*unit_eps0
        self.nspecies=len(self.D)
        self.external_charge=np.zeros([len(self.xmesh)])
        #self.external_charge=self.gaussian(sigma=5e-6,z=1.0,mu=4.e-5,cmax=5e-5)+self.gaussian(sigma=5e-6,z=-1.0,mu=5.e-5,cmax=5e-5)
        self.count=1
        self.ax1=plt.subplot('311')
        self.ax2=plt.subplot('312')
        self.ax3=plt.subplot('313')

        #BOUNDARY AND INITIAL CONDITIONS

        self.set_initial_conditions(
                c_initial_general={'all':self.bulk_concentrations[0]})
                #c_initial_specific={'0':{'0':0.1},'1':{'0':0.1}})

        self.set_boundary_conditions(\
                flux_boundary={'0':{'r':5e-8},'1':{'r':0.0}},         #in mol/s/cm^2
                dc_dt_boundary={'all':{'l':0.0}},     #in mol/l/s #give either left OR right boundary condition here
                #integration will start at the site where dc_dt is defined
                efield_boundary={'l':0.0})    #in V/Ang

        self.initialize='bla' #diffusion' #initialize with with diffusion or nothing
        #c_inf=self.get_static_concentrations()


    def get_static_concentrations(self):
        """solves the PBE in order to get the static limit for the concentrations"""
        def d(z, x):
            zold=z
            z=[]
            #z[0] = phi'
            #z[1] = phi
            #first block defines z[0]' which is d^2/dx^2 phi as function of
            #second block defines z[1]' which is phi, so z[0]
            z.append(-1./self.eps*self.external_charge_func(x)) #(sum([charge*unit_F*0.1*10**3*np.exp(-self.beta*unit_F*charge*zold[1]) for charge in self.charges])+self.external_charge_func(x)))
            z.append(zold[0])
            return z

        v0=np.zeros([2]) #this is phi' and phi at x=0
        sol,output = integrate.odeint(d, v0, self.xmesh,full_output=True) #, args=(b, c))
        print np.shape(sol)
        print sol[:,0]
        #b_min = minimize(func, b0, args=(z0, m, k, g, a), constraints=cons)
        self.ax1.plot(self.xmesh,sol[:,0]/1e10,'-',label='E (V/Ang)')
        self.ax2.plot(self.xmesh,sol[:,1],'-',label='phi (V)')
        self.ax1.legend()
        self.ax2.legend()
        plt.show()
        sys.exit()
        return sol

    def run(self):

        self.cout=self.integrate_pnp(self.dx,len(self.xmesh),self.dt,\
                len(self.tmesh),5,method=self.integrator)
       
    def gouy_chapman(self,x):
        term1 = 1.+np.tanh(self.vzeta*self.beta*unit_F/4.)*\
            np.exp(-1./self.debye_length*x)
        term2 = 1.-np.tanh(self.vzeta*self.beta*unit_F/4.)*\
            np.exp(-1./self.debye_length*x)
        return 2./(self.beta*unit_F)*np.log(term1/term2)

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
        #self.potential=np.zeros([len(self.xmesh)])
        #integral=0.0
        #for i in range(len(self.xmesh)):
        #    integral-=self.efield[i]*self.dx
        #    self.potential[i]+=integral
        ax3.plot(self.xmesh[:-1],self.potential[:-1],'-')
        ax3.plot(self.xmesh[:-1],[self.gouy_chapman(x) for x in self.xmesh[:-1]],'-',color='k',linewidth=lw)
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
        #c0*=10**3 #multiply by 1000 to get m^-3
        self.c0=c0


    def get_boundary_conditions(self):
        return self.c_bound,self.dc_dt_bound,self.efield_bound

    def set_boundary_conditions(self,\
        flux_boundary={},\
        dc_dt_boundary={},\
        efield_boundary={}):

        """allows to set boundary conditions in the concentrations themselves (c_boundary) 
        and the derivatives dc_dx_boundary and the efield. 
        species, direction, value"""

        if len(flux_boundary)>0:
            self.boundary_type='flux'
        elif len(dc_dt_boundary)>0:
            self.boundary_type='dc_dt'
        else:
            print('No boundary conditions defined, stopping here.')
            sys.exit()

        flux_bound=np.zeros([self.nspecies,2])
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

        def allocate_fluxes(array_in,array_out):
            for key1_raw in array_in:
                if key1_raw=='all':
                    keys1=map(str,range(self.nspecies))
                else:
                    keys1=[str(key1_raw)]
                for key1 in keys1:
                    for key2_raw in array_in[key1_raw]:
                        if key2_raw=='all':
                            keys2=map(str,range(2))
                        else:
                            keys2=[key_to_index(key2_raw)]
                        for key2 in keys2:
                            array_out[int(key1),int(key2)]=array_in[key1_raw][key2_raw]
            return array_out

        dc_dt_bound=allocate_fluxes(dc_dt_boundary,dc_dt_bound)
        flux_bound=allocate_fluxes(flux_boundary,flux_bound)

        for l in [dc_dt_bound,efield_bound,flux_bound]:
            l=np.array(l)

        self.flux_bound=flux_bound*100**2
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

    def external_charge_func(self,x):
        return self.gaussian_func(x,sigma=5e-6,z=1.0,mu=4.e-5,cmax=5e-5)+self.gaussian_func(x,sigma=5e-6,z=-1.0,mu=5.e-5,cmax=5e-5)

    def gaussian_func(self,x,sigma=0.01,z=1,mu=0.2,cmax=1):
        """define gaussian charge density. c is in molar"""
#        sigma=1./(np.sqrt(2.*np.pi))*1/cmax
        gaussian=\
             1./(sigma*np.sqrt(2.*np.pi))\
             *np.exp(-0.5*((x-mu)/sigma)**2)
        gaussian=gaussian*cmax*10**3*z*unit_F
        return gaussian

    def integrate_pnp(self,dx,nx,dt,nt,ntout,method='FTCS'):

        def calculate_efield_FD(nx,c):
            #solve A * u = f for u
            efield=np.zeros([nx])
            defield_dx=np.zeros([nx])
            #sum up all charges in order to get derivative of e-field
            for k in range(self.nspecies):
                defield_dx+=self.charges[k]*c[k*nx:(k+1)*nx]
            #add external charge
            defield_dx+=self.external_charge
            #devide by eps
            defield_dx/=self.eps
            #now perform FD scheme in order to solve d/dx efield = rho for efield:
            A = diags([-1, 0, 1], [-1, 0, 1], 
                  shape=(nx-2, nx-2)).toarray()
            A/=(2*dx)
            print('\n'.join([''.join(['{:10}'.format(item) for item in row])
      for row in A]))
            #RHS vector:
            f=[]
            for k in range(self.nspecies):
                f.extend(defield_dx[k*nx+1:(k+1)*nx-1])
            print f
            print np.linalg.det(A)
            tmp = np.linalg.solve(A,f) #this gives vector without initial and final elements
            for k in range(self.nspecies):
                efield[k*nx]=0.0
                efield[(k+1)*nx-1]=0.0
                efield[k*nx+1,(k+1)*nx-1]=tmp
            return efield, defield_dx

        def calculate_efield_both_sites(nx,c):
            #first determine E = int E' = 1/eps int sum c_i z_i by numerical integration over the grid
            efield=np.zeros([nx])
            defield_dx=np.zeros([nx])

            self.total_charge=np.zeros([nx])

            self.total_concentrations=np.zeros([self.nspecies+1,nx])
            for k in range(self.nspecies):
                defield_dx+=self.charges[k]*c[k*nx:(k+1)*nx]
            #add external charge
            defield_dx+=self.external_charge
            self.total_charge=defield_dx
            #devide by eps
            defield_dx/=self.eps

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
                    total_int+=current_charge*dx
                    efield[i]+=total_int

            for k in range(self.nspecies):
                #INTEGRATION OF PBE: get electric field
                if self.efield_bound[1] !=None:
                    #integration from the right
                    ll=range(0,nx)
                elif self.efield_bound[0] !=None:
                    #integration from the left
                    ll=reversed(range(0,nx))
                total_int=0.0
                m=-1
                for i in ll:
                    m+=1
                    j=k*nx+i
                    current_charge=self.charges[k]*c[j]
                    if k==0:
                        current_charge+=self.external_charge[i]
                    total_int+=current_charge*dx
                    efield[i]+=total_int
            efield/=self.eps
#            self.total_charge=(self.charges[0]*c[:nx]+self.charges[1]*c[nx:])/self.eps
#            self.defield_dx=self.total_charge
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
                efield = -efield
            #constant shift of electric field due to boundary conditions
            efield += [e for e in self.efield_bound if e!=None][0]
            #if self.count<200:
            #    self.ax1.plot(self.xmesh,efield,'-')
            #    self.ax2.plot(self.xmesh,self.charges[0]*c[:nx]+self.charges[1]*c[nx:],'-')
            #else:
            #    plt.show()
            #    sys.exit()
            #self.count+=1

            return efield,defield_dx

        def calculate_efield(nx,c):
            #first determine E = int E' = 1/eps int sum c_i z_i by numerical integration over the grid
            efield=np.zeros([nx])
            defield_dx=np.zeros([nx])

            self.total_charge=np.zeros([nx])

            self.total_concentrations=np.zeros([self.nspecies+1,nx])

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
                    total_int+=current_charge*dx
                    efield[i]+=total_int
                    defield_dx[i]+=current_charge
                    #others:
                    self.total_charge[i]+=current_charge
                    self.total_concentrations[k,i]+=c[j]
            efield/=self.eps
            defield_dx/=self.eps
#            self.total_charge=(self.charges[0]*c[:nx]+self.charges[1]*c[nx:])/self.eps
#            self.defield_dx=self.total_charge
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
                efield = -efield
            #constant shift of electric field due to boundary conditions
            efield += [e for e in self.efield_bound if e!=None][0]
            #if self.count<200:
            #    self.ax1.plot(self.xmesh,efield,'-')
            #    self.ax2.plot(self.xmesh,self.charges[0]*c[:nx]+self.charges[1]*c[nx:],'-')
            #else:
            #    plt.show()
            #    sys.exit()
            #self.count+=1

            return efield,defield_dx

        def ode_func_old(c,t):

            self.efield,self.defield_dx=calculate_efield(nx,c)

            dc_dt = np.zeros([nx*self.nspecies])


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
                            +self.charges[k]*self.beta*dc_dx*self.efield[i]\
                            +self.charges[k]*self.beta*c[j]*self.defield_dx[i]\
                        )\
                        )
#            self.ax1.plot(self.xmesh,dc_dt/10**3,'-')
##            plt.show()
#            sys.exit()
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

        def integrate_FTCS_odeint_old(dx,nx,dt,nt,c0,ntout):
            sol,output = integrate.odeint(ode_func, c0, range(0,nt),full_output=True) #, args=(b, c))
            print 'Used minimal time = ',min(output['tcur']),' and maximal time = ',max(output['tcur'])
            #return the results
            cout = []
            for n in range(0,nt):
                if n % int(nt/float(ntout)) == 0 or n==nt-1:
                    cout.append(sol[n,:].copy()) # numpy arrays are mutable, 
            return cout

        def integrate_Crank_Nicolson_pnp(dx,nx,dt,nt,C0,ntout):
            """Integrates PNP equations, BCs:
                        concentrations      potential
                left    ROBIN j=0           DIRICHLET vzeta
                right   DIRICHLET cbulk     DIRICHLET 0.0
            """
            #if true add an artifical diffusion term which enhances the stability of the solution
            self.use_lax_friedrichs=True

            # create coefficient matrix:
            def a_matrix(s):
                return diags([-0.5*s, 1+s, -0.5*s], [-1, 0, 1],\
                    shape=(nx-2, nx-2)).toarray()

            def b1_matrix(s):
                return diags([0.5*s, 1-s, 0.5*s],[-1, 0, 1],\
                    shape=(nx-2,nx-2)).toarray()

            def add_field(B1,A,grad_v,lapl_v,ee):
                for i in range(nx-2):
                    if i==0:
                        jvalues=[i,i+1]
                    elif i==nx-3:
                        jvalues=[i-1,i]
                    else:
                        jvalues=[i-1,i,i+1]
                    for j in jvalues:
                        if i==j:
                            B1[i,j]+=ee*lapl_v[i]
                        if abs(i-j)==1:
                            if j<i:
                                B1[i,j]-=ee*grad_v[i]/4./dx
                                A[i,j]+=ee*grad_v[i]/4./dx
                            elif i<j:
                                B1[i,j]+=ee*grad_v[i]/4./dx
                                A[i,j]-=ee*grad_v[i]/4./dx
                return B1,A

            def add_boundary_values(B,grad_v,s,ee,C0,C1,C0_old,C1_old):
                #C0 and C1 are the current boundary values
                #C0_old and C1_old were the last iteration's boundary values
                B[0] += (0.5*s+\
                    ee*grad_v[0]/4./dx)*(C0+C0_old)
                B[-1] += (0.5*s-\
                    ee*grad_v[-1]/4./dx)*(C1+C1_old)
                return B

            C = np.zeros([self.nspecies,nx])
            COLD = np.zeros([self.nspecies,nx])
            COUT=[]

            #initial conditions for concentrations
            for k in range(self.nspecies):
                C[k,:] = C0[k*nx:(k+1)*nx]
    

            #time iteration
            for n in range(1,nt):
                print 'time step = ',n
                v,grad_v,lapl_v=get_potential_and_gradient(C,dx,nx)
                for k in range(self.nspecies):
                    if n==1:
                        COLD[k,:]=deepcopy(C[k,:])

                    #Robin BC for concentrations on left side (wall)
                    C[k,0] =\
                        (-4*self.D[k]-self.mu[k]*(v[1]-self.vzeta))/\
                        (-4*self.D[k]+self.mu[k]*(v[1]-self.vzeta))*C[k,1]
                    #Dirichlet BC for concentrations on right side (bulk)
                    C[k,-1]=C0[(k+1)*nx-1]


                    s = self.D[k]*dt/dx**2  # diffusion number
                    if self.use_lax_friedrichs:
                        #add artificial diffusion term for better stability
                        s+=0.5
                    ee = self.charges[k]*self.beta*dt*self.D[k]

                    A=a_matrix(s)
                    B1=b1_matrix(s)
                    B1,A=add_field(B1,A,grad_v,lapl_v,ee)
                    B = np.dot(C[k,1:-1],B1)
                    B=add_boundary_values(B,grad_v,s,ee,C[k,0],C[k,-1],COLD[k,0],COLD[k,-1])

                    CTMP = np.linalg.solve(A,B) #this gives vector without initial and final elements
                    C[k,1:-1] = CTMP
                    COLD[k,:]=C[k,:]

                if n % int(nt/float(ntout)) == 0 or n==nt-1: # or True:
                    COUT.append(np.ndarray.flatten(C)) # numpy arrays are mutable, 
                    #so we need to write out a copy of c, not c itself
            return COUT,s

        def integrate_Crank_Nicolson_pnp_old(dx,nx,dt,nt,c,ntout):
            """only for dc_dt=0 at both boundary sites implemented yet"""
            cout = [] # list for storing c arrays at certain time steps
            c0 = np.zeros([self.nspecies])
            c1 = np.zeros([self.nspecies])
            s = np.zeros([self.nspecies])
            ee = np.zeros([self.nspecies])
            for k in range(self.nspecies):
                c0[k] = c[k*nx] # boundary condition on left side
                c1[k] = c[(k+1)*nx-1] # boundary condition on right side
                s[k] = self.D[k]*dt/dx**2  # diffusion number
                ee[k] = self.charges[k]*self.beta*dt*self.D[k]
            # create coefficient matrix:
            def a_matrix(s):
                return diags([-0.5*s, 1+s, -0.5*s], [-1, 0, 1], 
                  shape=(nx-2, nx-2)).toarray() 

            A=a_matrix(s[0])
            Atot=deepcopy(A)
            for k in range(self.nspecies-1):
                Atot=block_diag(Atot,a_matrix(s[k+1]))
            A=Atot
            def b1_matrix(s):
                return diags([0.5*s, 1-s, 0.5*s],[-1, 0, 1], shape=(nx-2,nx-2)).toarray()

            B1=b1_matrix(s[0])
            B1tot=deepcopy(B1)
            for k in range(self.nspecies-1):
                B1tot=block_diag(B1tot,b1_matrix(s[k+1]))
            B1=B1tot

            def unpack(ctmp,c0,c1,nx):
                cout=[]
                j=-1
                for k in range(self.nspecies):
                    for i in range(nx):
                        if i!=0 and i!=nx-1:
                            j+=1
                            cout.append(ctmp[j])
                        elif i==0:
                            cout.append(c0[k])
                        elif i==nx-1:
                            cout.append(c1[k])
                return np.array(cout)
            
            for n in range(1,nt): # time is going from second time step to last
                #self.ax1.plot(self.xmesh,c[:nx],'-')
                #if n<50 and abs(self.charges[0])>0.0:
                #    self.ax1.plot(self.xmesh,c[:nx],label='c1')
                #    self.ax1.plot(self.xmesh,c[nx:],label='c2')
                #    self.ax1.legend()
                #if n==50 and abs(self.charges[0])>0.0:
                #    plt.show()
                #    sys.exit()
                cn = c
                #PNP electric field modifications
                self.efield,self.defield_dx=calculate_efield(nx,cn)
                #if n<50 and abs(self.charges[0])>0.0:
                #    self.ax2.plot(self.xmesh,self.efield[:nx],'-')
                #    self.ax3.plot(self.xmesh,self.defield_dx[:nx],'-')
                cn_slice=[cc for ic,cc in enumerate(cn) if ((ic+1)%nx!=0 and (ic)%nx!=0)]

                for i in range((nx-2)*self.nspecies):
                    for j in [i-1,i,i+1]: #range((nx-2)*self.nspecies):
                        if (i%(nx-2)==0 and j==i-1):
                            #beginning of new species. should not calculate derivative
                            #using last position
                            continue
                        if (i+1)%(nx-2)==0 and j==i+1:
                            #end of species. should not use next position
                            continue
                        k=i//(nx-2)     #current species
                        ii=i-k*(nx-2)   #current i index
                        if i==j:
                            B1[i,j]+=ee[k]*self.defield_dx[ii+1]
                        if abs(i-j)==1:
                            if j<i:
                                B1[i,j]-=ee[k]*self.efield[ii+1]/4./dx
                                A[i,j]+=ee[k]*self.efield[ii+1]/4./dx
                            elif i<j:
                                B1[i,j]+=ee[k]*self.efield[ii+1]/4./dx
                                A[i,j]-=ee[k]*self.efield[ii+1]/4./dx
                #print '-'*50
                #print('\n'.join([''.join(['{:4}'.format(item) for item in row])
                #    for row in A]))
                #print('\n'.join([''.join(['{:4}'.format(item) for item in row])
                #    for row in B1]))
                #print '-'*50
                #end electric field modifications

                B = np.dot(cn_slice,B1)
                for k in range(self.nspecies):
                    B[k*(nx-2)] += (0.5*s[k]+\
                        ee[k]*self.efield[0]/4./dx)*(c0[k]+c0[k])
                    B[(k+1)*(nx-2)-1] += (0.5*s[k]-\
                        ee[k]*self.efield[-1]/4./dx)*(c1[k]+c1[k])
                    #B[k*(nx-2)] += 0.5*s[k]*(c0[k]+c0[k])
                    #B[(k+1)*(nx-2)-1] += 0.5*s[k]*(c1[k]+c1[k])
                #if abs(self.charges[0])>0.0 and n==1:
                #    for i in range(len(B)):
                #        print i, B[i]
                #    for i in range(len(A[:,0])):
                #        for j in range(len(A[0,:])):
                #            if A[i,j]!=0.0:
                #                print i,j,A[i,j]
                ctmp = np.linalg.solve(A,B) #this gives vector without initial and final elements
                c = unpack(ctmp,c0,c1,nx) #add left and right boundary values back
                if n % int(nt/float(ntout)) == 0 or n==nt-1: # or True:
                    cout.append(c.copy()) # numpy arrays are mutable, 
                    #so we need to write out a copy of c, not c itself
            return cout,s

        def get_potential_and_gradient(C,dx,nx):
            """Calculates the potential and its gradient (and laplacian=charge density)
                from the PBE by Jacobi relaxation (FD)"""
            # calculate RHS of PBE
            rhs=np.zeros([nx])
            for k in range(self.nspecies):
                rhs+=self.charges[k]*C[k,:]/self.eps
            lapl_v=-rhs
            v=np.zeros([nx])
            v[0]=self.vzeta
            v_old=deepcopy(v)
            tau_jacobi=1e-5
            error=np.inf
            i_step=0
            while error>tau_jacobi**2:# or i_step<3:
                i_step+=1
                for i in range(1,nx-1):
                    v[i] = 1/2.*(dx**2*rhs[i]+v_old[i+1]+v_old[i-1])
                error=sum((v_old-v)**2/len(v))
                v_old=deepcopy(v)
            grad_v=np.zeros([nx])
            for i in range(1,nx-1):
                grad_v[i] = 1./(2*dx)*(v[i+1]-v[i-1])
            #linearly extrapolate to get derivative at boundaries:
            grad_v[0] = grad_v[1]-(grad_v[2]-grad_v[1])
            grad_v[-1] = grad_v[-2]-(grad_v[-2]-grad_v[-3])
            self.efield=grad_v
            self.potential=v
            self.total_charge=lapl_v*self.eps
            return v,grad_v,lapl_v

        def integrate_FTCS_odeint(dx,nx,dt,nt,c0,ntout):
            """Integrates PNP equations, BCs:
                        concentrations      potential
                left    ROBIN j=0           DIRICHLET vzeta
                right   DIRICHLET cbulk     DIRICHLET 0.0
            """
            def ode_func(c,t,dx):

                #map concentrations onto 2D array:
                C = np.zeros([self.nspecies,nx])
                for k in range(self.nspecies):
                    C[k,:]=c[k*nx:(k+1)*nx]

                #get potential and field
                v,grad_v,lapl_v=get_potential_and_gradient(C,dx,nx)

                DC_DT = np.zeros([self.nspecies,nx])
#                flow = np.zeros([self.nspecies,nx])
#
#
#                #first tabulate D*(dc/dx + mu * c dv/dx) = flow
#                for k in range(0,self.nspecies):
#                    #no flow boundary condition on the left
#                    flow[k,0]=0.0
#                    for i in range(1,nx-1):
#                        dc_dx=(C[k,i+1]-C[k,i-1])/(2.*dx)
##                    if self.lax_friedrich:
##                        corr=dc_dx/2.*dx**2/dt
##                    else:
##                        corr=0.0
#                    corr=0.0
#                    flow[k,i] =\
#                        self.D[k]*\
#                            (\
#                            dc_dx\
#                            +self.beta*self.charges[k]*C[k,i]*grad_v[i]\
#                            )\
#                    #extrapolate flow on the right:
#                    flow[k,-1]=flow[k,-2]+(flow[k,-2]-flow[k,-3])
                self.lax_friedrich=False #True

                #go over the flow and calculate derivative
                for k in range(0,self.nspecies):
                    #set concentration to be constant on the right side
                    DC_DT[k,-1] = 0.0
                    for i in range(0,nx-1):
                        #dflow_dx=(flow[k,i+1]-flow[k,i-1])/(2.*dx)
                        dc_dx=(C[k,i+1]-C[k,i-1])/(2.*dx)
                        dc_dx_2=(C[k,i+1]-2*C[k,i]+C[k,i-1])/(dx**2)
                        dcgradv_dx=(C[k,i+1]*grad_v[i+1]-C[k,i-1]*grad_v[i-1])/(2.*dx)
                        corr=0.0
                        if i==0:
                            if self.lax_friedrich:
                                corr=2*(C[k,1]-C[k,0])/dt/2.
                            #we have to consider no flow boundary condition here
                            #setting C[k,-1]=C[k,1] considers this:
                            DC_DT[k,i]=\
                                corr+\
                                self.D[k]*\
                                    (\
                                    2*(C[k,1]-C[k,0])/dx**2\
                                    +self.beta*self.charges[k]*dcgradv_dx\
                                    )
                        else:
                            if self.lax_friedrich:
                                corr=dc_dx_2*dx**2/dt/2.
                            #we can directly implement the 2nd derivatives here
                            DC_DT[k,i]=\
                                corr+\
                                self.D[k]*\
                                    (\
                                    dc_dx_2\
                                    +self.beta*self.charges[k]*dcgradv_dx\
                                    )

                    #extrapolate derivative on the left side
#                    DC_DT[k,0]=DC_DT[k,1]+(DC_DT[k,1]-DC_DT[k,2])

                #map dc_dt's onto output format
                dc_dt = np.zeros([self.nspecies*nx])
                for k in range(self.nspecies):
                    dc_dt[k*nx:(k+1)*nx]=DC_DT[k,:]
                return dc_dt

            sol,output = integrate.odeint(ode_func, c0, self.tmesh, args=(dx,),full_output=True, ml=self.nspecies, mu=self.nspecies)
#            sol,output = integrate.odeint(ode_func, c0, range(0,nt),full_output=True) #, args=(b, c))
            print 'Used minimal time = ',min(output['tcur']),' and maximal time = ',max(output['tcur'])
            #return the results
            print np.shape(sol)
            cout = []
            for n in range(0,nt):
                if n % int(nt/float(ntout)) == 0 or n==nt-1:
                    cout.append(sol[n,:].copy()) # numpy arrays are mutable, 
            return cout


        def integrate_FTCS(dt,dx,nt,nx,C0,ntout):
            """Integrates PNP equations, BCs:
                        concentrations      potential
                left    ROBIN j=0           DIRICHLET vzeta
                right   DIRICHLET cbulk     DIRICHLET 0.0
            """

            C = np.zeros([self.nspecies,nx])
            COUT=[]

            #initial conditions for concentrations
            for k in range(self.nspecies):
                C[k,:] = C0[k*nx:(k+1)*nx]

            for n in range(0,nt):
                print 'time step = ',n
                v,grad_v,lapl_v=get_potential_and_gradient(C,dx,nx)
                for k in range(self.nspecies):
                    #Robin BC for concentrations on left side (wall)
                    C[k,0] =\
                        (-4*self.D[k]-self.mu[k]*(v[1]-self.vzeta))/\
                        (-4*self.D[k]+self.mu[k]*(v[1]-self.vzeta))*C[k,1]
                    #Dirichlet BC for concentrations on right side (bulk)
                    C[k,-1]=C0[(k+1)*nx-1]
                    temp = np.zeros([nx])
                    temp[0]=C[k,0]
                    temp[-1]=C[k,-1]
                    for i in range(1,nx-1):
                        W = self.D[k]*dt/dx**2-\
                            dt/(2.*dx)*self.mu[k]*grad_v[i+1]+0.5
                        M = -2.*self.D[k]*dt/dx**2
                        E = self.D[k]*dt/dx**2+\
                            dt/(2.*dx)*self.mu[k]*grad_v[i-1]+0.5
                        temp[i] = E * C[k,i-1] + M * C[k,i] + W * C[k,i+1]
                    C[k,:]=temp
                if n % int(nt/float(ntout)) == 0 or n==nt-1:
                    COUT.append(np.ndarray.flatten(C))
                #if n==100:
                #    return cout
            return COUT

        def integrate_FTCS_old(dt,dx,nt,nx,V0,ntout):
            # diffusion number (has to be less than 0.5 for the 
            # solution to be stable):
            V = np.zeros([nt,self.nspecies*nx])
            ee=np.zeros([self.nspecies])
            s=np.zeros([self.nspecies])
            for k in range(self.nspecies):
                s[k] = self.D[k]*dt/dx**2
                ee[k] = self.charges[k]*self.beta*dt*self.D[k]
                #print 'before',self.charges[k]*self.beta*dt*self.D[k]
                #ee[k] = 7.62e-8*dt
                #print 'after',ee[k]
                V[:,k*nx] = [V0[k*nx]]*(nt)       #boundary left (
                V[:,(k+1)*nx-1] = [V0[(k+1)*nx-1]]*(nt)     #boundary right 

            #initial conditions
            V[0,:]=V0                   #initia
            cout=[]

            #boundary conditions for change in concentrations
            dV=np.zeros([self.nspecies*nx])
            for k in range(self.nspecies):
                dV[k*nx]=self.dc_dt_bound[k,0]*dt
                dV[(k+1)*nx-1]=self.dc_dt_bound[k,1]*dt

            for n in range(0,nt-1): # time
                self.efield,self.defield_dx=calculate_efield(nx,V[n,:])
                if self.boundary_type=='flux':
                    if sum([fb**2 for fb in self.flux_bound[:,0]])!=0.0:
                        #left bound is defined by flux. this means that right bound will be defined by dc_dt_bound
                        #we have to start integrating from the right
                        #(we cannot use different flux integration for both species)
                        xiter=reversed(range(1,nx-1))
                    else:
                        xiter=range(1,nx-1)
                else:
                    xiter=range(1,nx-1)
                for k in range(self.nspecies):
                    #update concentrations at the boundaries (from boundary condition)
                    V[n+1,k*nx]=V[n,k*nx]+dV[k*nx]
                    V[n+1,(k+1)*nx-1]=V[n,(k+1)*nx-1]+dV[(k+1)*nx-1]
                for k in range(self.nspecies):
                    if self.boundary_type=='flux':
                        flux_bound=np.sqrt(max(self.flux_bound[k,:]**2))
                   # V[n+1,(k+1)*nx-1]=V[n,(k+1)*nx-1]*dV[(k+1)*nx-1]
                    ii=0
                    for i in xiter: #k*nx+1,(k+1)*nx-1): # space
                        ii+=1
                        j=k*nx+i
                        if self.boundary_type=='flux' and ii==nx-2:
                            #define V update including flux boundary
                            V[n+1,j] = V[n,j] + self.D[k]*dt/dx*(\
                                flux_bound-(V[n,j]-V[n,j-1])/dx) +\
                                ee[k]*\
                                (self.efield[i]*(V[n,j+1]-V[n,j-1])/(2.*dx)+\
                                self.defield_dx[i]*V[n,j])
                        else:
                            V[n+1,j] = V[n,j] + s[k]*(V[n,j-1]-\
                                2*V[n,j] + V[n,j+1]) + ee[k]*\
                                (self.efield[i]*(V[n,j+1]-V[n,j-1])/(2.*dx)+\
                                self.defield_dx[i]*V[n,j])
                if n % int(nt/float(ntout)) == 0 or n==nt-1:
                    cout.append(V[n,:].copy())
                #if n==100:
                #    return cout
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
            dataplot=cout
        elif method=='FTCS':
            cout=integrate_FTCS(dt,dx,nt,nx,self.c0,ntout)
            dataplot=cout
        elif method=='Crank-Nicolson':
            if self.initialize=='diffusion':
                charges_tmp=self.charges
                self.charges=np.array([0.0]*len(self.charges))
                cout,s=integrate_Crank_Nicolson_pnp(dx,nx,dt,nt,self.c0,ntout)
                self.charges=charges_tmp
                print np.shape(cout)
                self.ax1.plot(self.xmesh,cout[-1][:nx],'-r')
                self.ax1.set_ylim([min(cout[-1][:nx]),max(cout[-1][:nx])])
                cout,s=integrate_Crank_Nicolson_pnp(dx,nx,dt,nt,cout[-1],ntout)
            else:
                print 'not diff'
                cout,s=integrate_Crank_Nicolson_pnp(dx,nx,dt,nt,self.c0,ntout)
            dataplot=cout[-1]
        #print dataplot
        #plt.plot(self.xmesh,dataplot[:nx]/10**3,'-')
        #plt.plot(self.xmesh,dataplot[nx:]/10**3,'-')
        #plt.show()
        #sys.exit()
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
         
