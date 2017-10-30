"""
calculator class
uses finite difference techniques to solve the PNP equations
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
import sys
from copy import deepcopy
from scipy.optimize import minimize
from scipy.integrate import ode
from scipy import interpolate
import os
import odespy, numpy

class Calculator():

    def __init__(self,transport=None,dt=None,tmax=None,ntout=5,calc=None,
            scale_pb_grid=None,tau_jacobi=1e-7):

        if transport is None:
            self.tp.logger.error('No transport object provided for calculator. Stopping here.')
            sys.exit()
        else:
            self.tp=transport #transport object
        
        if calc is None:
            calc=self.tp.calc
        
#        if os.path.exists('results.txt'):
#            os.remove('results.txt')

        self.use_lax_friedrich=False
        string=calc.split('--')
        self.calc_method=None
        if len(string)>1:
            if string[-1]=='LF':
                self.use_lax_friedrich=True
            else:
                self.calc_method=string[-1]
            self.calc=string[0]
        else:
            self.calc=calc

        self.calc_list=['FTCS','Crank-Nicolson','odeint','vode','lsoda','dopri5','dop853','odeint','odespy']

        if self.calc not in self.calc_list:
            self.tp.logger.error('No calculator found with this name. Aborting.')
            sys.exit()
        
        #use a finer grid for PBE integration (factor of scale_pb_grid smaller than
        #normal xmesh
        self.scale_pb_grid=scale_pb_grid
        #accuracy of Jacoby iteration if potential needs to be outputted:
        self.tau_jacobi=tau_jacobi

        if dt is not None: # and hasattr(self.tp,'dt'):
            self.tp.dt=dt
        if tmax is not None: # and hasattr(self.tp,'tmax'):
            self.tp.tmax=tmax
       # if hasattr(self.tp,'dt') and hasattr(self.tp,'tmax'):
        if tmax is not None or dt is not None:
            self.tp.tmesh=np.arange(0,self.tp.tmax+self.tp.dt,self.tp.dt)
        self.oldtime=np.inf
        self.tp.ntout=ntout
        self.initialize='bla'
        return

    def get_rates(self,C):
        #get list of all species names:
        species_names=[sp for sp in self.tp.species]
        rates=np.zeros([self.tp.nspecies,self.tp.nx])
        for i in range(self.tp.nx):
            for r in self.tp.reactions:
                reaction=self.tp.reactions[r]
                if 'rates' not in reaction:
                    continue
                #rates of educts:
                for reactant in reaction['reactants'][0]:
                    if reactant not in self.tp.species:
                        continue
                    k=species_names.index(reactant)
                    rates[k,i]=0.0
                    prod=1.
                    for reactant2 in reaction['reactants'][0]:
                        if reactant2 not in self.tp.species:
                            continue
                        k2=species_names.index(reactant2)
                        prod*=C[k2,i]
                    rates[k,i]-=prod*reaction['rates'][0]
                    prod=1.
                    for reactant2 in reaction['reactants'][1]:
                        if reactant2 not in self.tp.species:
                            continue
                        k2=species_names.index(reactant2)
                        prod*=C[k2,i]
                    rates[k,i]+=prod*reaction['rates'][1]
                #rates of products
                for reactant in reaction['reactants'][1]:
                    if reactant not in self.tp.species:
                        continue
                    k=species_names.index(reactant)
                    rates[k,i]=0.0
                    prod=1.
                    for reactant2 in reaction['reactants'][0]:
                        if reactant2 not in self.tp.species:
                            continue
                        k2=species_names.index(reactant2)
                        prod*=C[k2,i]
                    rates[k,i]+=prod*reaction['rates'][0]
                    prod=1.
                    for reactant2 in reaction['reactants'][1]:
                        if reactant2 not in self.tp.species:
                            continue
                        k2=species_names.index(reactant2)
                        prod*=C[k2,i]
                    rates[k,i]-=prod*reaction['rates'][1]
        return rates

    def integrate_pnp(self,dx,nx,dt,nt,ntout,method):

        def calculate_efield_FD(nx,c):
            #solve A * u = f for u
            efield=np.zeros([nx])
            defield_dx=np.zeros([nx])
            #sum up all charges in order to get derivative of e-field
            for k in range(self.tp.nspecies):
                defield_dx+=self.tp.charges[k]*c[k*nx:(k+1)*nx]
            #add external charge
            defield_dx+=self.tp.external_charge
            #devide by eps
            defield_dx/=self.tp.eps
            #now perform FD scheme in order to solve d/dx efield = rho for efield:
            A = diags([-1, 0, 1], [-1, 0, 1], 
                  shape=(nx-2, nx-2)).toarray()
            A/=(2*dx)
            self.tp.logger_db.debug('\n'.join([''.join(['{:10}'.format(item) for item in row])
      for row in A]))
            #RHS vector:
            f=[]
            for k in range(self.tp.nspecies):
                f.extend(defield_dx[k*nx+1:(k+1)*nx-1])
            tmp = np.linalg.solve(A,f) #this gives vector without initial and final elements
            for k in range(self.tp.nspecies):
                efield[k*nx]=0.0
                efield[(k+1)*nx-1]=0.0
                efield[k*nx+1,(k+1)*nx-1]=tmp
            return efield, defield_dx

        def calculate_efield_both_sites(nx,c):
            #first determine E = int E' = 1/eps int sum c_i z_i by numerical integration over the grid
            efield=np.zeros([nx])
            defield_dx=np.zeros([nx])

            self.tp.total_charge=np.zeros([nx])

            self.tp.total_concentrations=np.zeros([self.tp.nspecies+1,nx])
            for k in range(self.tp.nspecies):
                defield_dx+=self.tp.charges[k]*c[k*nx:(k+1)*nx]
            #add external charge
            defield_dx+=self.tp.external_charge
            self.tp.total_charge=defield_dx
            #devide by eps
            defield_dx/=self.tp.eps

            for k in range(self.tp.nspecies):
                #INTEGRATION OF PBE: get electric field
                if self.tp.efield_bound[1] !=None:
                    #integration from the right
                    ll=reversed(range(0,nx))
                elif self.tp.efield_bound[0] !=None:
                    #integration from the left
                    ll=range(0,nx)
                total_int=0.0
                m=-1
                for i in ll:
                    m+=1
                    j=k*nx+i
                    current_charge=self.tp.charges[k]*c[j]
                    if k==0:
                        current_charge+=self.tp.external_charge[i]
                    total_int+=current_charge*dx
                    efield[i]+=total_int

            for k in range(self.tp.nspecies):
                #INTEGRATION OF PBE: get electric field
                if self.tp.efield_bound[1] !=None:
                    #integration from the right
                    ll=range(0,nx)
                elif self.tp.efield_bound[0] !=None:
                    #integration from the left
                    ll=reversed(range(0,nx))
                total_int=0.0
                m=-1
                for i in ll:
                    m+=1
                    j=k*nx+i
                    current_charge=self.tp.charges[k]*c[j]
                    if k==0:
                        current_charge+=self.tp.external_charge[i]
                    total_int+=current_charge*dx
                    efield[i]+=total_int
            efield/=self.tp.eps
#            self.tp.total_charge=(self.tp.charges[0]*c[:nx]+self.tp.charges[1]*c[nx:])/self.tp.eps
#            self.tp.defield_dx=self.tp.total_charge
            ##integrate charge with Gaussian quadrature, default order=5
            #for k in range(0,self.tp.nspecies):
            #    func=self.tp.charges[k]*c[k*nx:(k+1)*nx]
            #    #integrate.fixed_quad(func,min(func),max(func))
            #    integral+=integrate.simps(func,self.tp.xmesh)

            self.tp.total_concentrations[-1,:]=self.tp.external_charge/unit_F

            #plt.plot(self.tp.xmesh,self.tp.total_charge/10**3/unit_F)
            #for t in self.tp.total_concentrations:
            #    plt.plot(self.tp.xmesh,t/10**3)
            #plt.show()
            #sys.exit()

            if self.tp.efield_bound[1] != None:
                #we integrated from right to left, so we have to switch the sign
                #of the efield here
                efield = -efield
            #constant shift of electric field due to boundary conditions
            efield += [e for e in self.tp.efield_bound if e!=None][0]
            #if self.tp.count<200:
            #    self.tp.ax1.plot(self.tp.xmesh,efield,'-')
            #    self.tp.ax2.plot(self.tp.xmesh,self.tp.charges[0]*c[:nx]+self.tp.charges[1]*c[nx:],'-')
            #else:
            #    plt.show()
            #    sys.exit()
            #self.tp.count+=1

            return efield,defield_dx

        def calculate_efield(nx,c):
            #first determine E = int E' = 1/eps int sum c_i z_i by numerical integration over the grid
            efield=np.zeros([nx])
            defield_dx=np.zeros([nx])

            self.tp.total_charge=np.zeros([nx])

            self.tp.total_concentrations=np.zeros([self.tp.nspecies+1,nx])

            for k in range(self.tp.nspecies):
                #INTEGRATION OF PBE: get electric field
                if self.tp.efield_bound[1] !=None:
                    #integration from the right
                    ll=reversed(range(0,nx))
                elif self.tp.efield_bound[0] !=None:
                    #integration from the left
                    ll=range(0,nx)
                total_int=0.0
                m=-1
                for i in ll:
                    m+=1
                    j=k*nx+i
                    current_charge=self.tp.charges[k]*c[j]
                    if k==0:
                        current_charge+=self.tp.external_charge[i]
                    total_int+=current_charge*dx
                    efield[i]+=total_int
                    defield_dx[i]+=current_charge
                    #others:
                    self.tp.total_charge[i]+=current_charge
                    self.tp.total_concentrations[k,i]+=c[j]
            efield/=self.tp.eps
            defield_dx/=self.tp.eps
#            self.tp.total_charge=(self.tp.charges[0]*c[:nx]+self.tp.charges[1]*c[nx:])/self.tp.eps
#            self.tp.defield_dx=self.tp.total_charge
            ##integrate charge with Gaussian quadrature, default order=5
            #for k in range(0,self.tp.nspecies):
            #    func=self.tp.charges[k]*c[k*nx:(k+1)*nx]
            #    #integrate.fixed_quad(func,min(func),max(func))
            #    integral+=integrate.simps(func,self.tp.xmesh)

            self.tp.total_concentrations[-1,:]=self.tp.external_charge/unit_F

            #plt.plot(self.tp.xmesh,self.tp.total_charge/10**3/unit_F)
            #for t in self.tp.total_concentrations:
            #    plt.plot(self.tp.xmesh,t/10**3)
            #plt.show()
            #sys.exit()

            if self.tp.efield_bound[1] != None:
                #we integrated from right to left, so we have to switch the sign
                #of the efield here
                efield = -efield
            #constant shift of electric field due to boundary conditions
            efield += [e for e in self.tp.efield_bound if e!=None][0]
            #if self.tp.count<200:
            #    self.tp.ax1.plot(self.tp.xmesh,efield,'-')
            #    self.tp.ax2.plot(self.tp.xmesh,self.tp.charges[0]*c[:nx]+self.tp.charges[1]*c[nx:],'-')
            #else:
            #    plt.show()
            #    sys.exit()
            #self.tp.count+=1

            return efield,defield_dx

        def ode_func_old(c,t):

            self.tp.efield,self.tp.defield_dx=calculate_efield(nx,c)

            dc_dt = np.zeros([nx*self.tp.nspecies])


#            plt.plot(self.tp.xmesh,(self.tp.charges[0]*c[:self.tp.nx]+self.tp.charges[1]*c[self.tp.nx:])/self.tp.eps,'o')
#            plt.plot(self.tp.xmesh,self.tp.external_charge/self.tp.eps,'o')
##            plt.plot(self.tp.xmesh,self.tp.efield,label='e')
#            plt.plot(self.tp.xmesh,self.tp.defield_dx,label='de/dx')
#            plt.legend()

            #IMPLEMENTATION OF FLOWS
            #go over all point in space and all species

            for k in range(0,self.tp.nspecies):
                m=-1
                for i in range(0,nx):
                    m+=1
                    j=k*nx+i

                    if m==0:
                        dc_dt[j] = self.tp.dc_dt_bound[k,0]
                        continue
                    elif m==nx-1:
                        dc_dt[j] = self.tp.dc_dt_bound[k,1]
                        continue

                    dc_dx=(c[j+1]-c[j-1])/(2.*dx)
                    dc_dx_2=(c[j-1]-2*c[j]+c[j+1])/dx**2

                    dc_dt[j] =\
                        (\
                        -self.tp.u[k]*dc_dx \
                        +self.tp.D[k]*(\
                            dc_dx_2 \
                            +self.tp.charges[k]*self.tp.beta*dc_dx*self.tp.efield[i]\
                            +self.tp.charges[k]*self.tp.beta*c[j]*self.tp.defield_dx[i]\
                        )\
                        )
#            self.tp.ax1.plot(self.tp.xmesh,dc_dt/10**3,'-')
##            plt.show()
#            sys.exit()
         #   if self.tp.count in range(50,100):
         #       self.tp.ax1.plot(self.tp.xmesh,dc_dt[:self.tp.nx]/10**3*dt,'-',linewidth=0.3,color='r')
         #       self.tp.ax1.plot(self.tp.xmesh,dc_dt[self.tp.nx:]/10**3*dt,'-',linewidth=0.3,color='k')
         #       self.tp.ax2.plot(self.tp.xmesh,c[:self.tp.nx]/10**3,'-',linewidth=0.3,color='r')
         #       self.tp.ax2.plot(self.tp.xmesh,c[self.tp.nx:]/10**3,'-',linewidth=0.3,color='k')
         #   self.tp.count+=1
         #   if self.tp.count==100:
         #       plt.show()
         #       sys.exit()
            return dc_dt

        def integrate_FTCS_odeint_old(dx,nx,dt,nt,c0,ntout):
            sol,output = integrate.odeint(ode_func, c0, range(0,nt),full_output=True) #, args=(b, c))
            self.tp.logger.info('Used minimal time = ',min(output['tcur']),' and maximal time = ',max(output['tcur']))
            #return the results
            cout = []
            for n in range(0,nt):
                if n % int(nt/float(ntout)) == 0 or n==nt-1:
                    cout.append(sol[n,:].copy()) # numpy arrays are mutable, 
            return cout

        def integrate_Crank_Nicolson(dx,nx,dt,nt,C0,ntout):
            """Integrates PNP equations, BCs:
                        concentrations      potential
                left    ROBIN j=0           DIRICHLET vzeta
                right   DIRICHLET cbulk     DIRICHLET 0.0
            """

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

            C = np.zeros([self.tp.nspecies,nx])
            COLD = np.zeros([self.tp.nspecies,nx])
            COUT=[]

            #initial conditions for concentrations
            for k in range(self.tp.nspecies):
                C[k,:] = C0[k*nx:(k+1)*nx]
    

            #time iteration
            for n in range(1,nt):
                self.tp.logger.info('time step = {}'.format(n))
                if self.tp.use_migration:
                    v,grad_v,lapl_v=get_potential_and_gradient(C,dx,nx)
#                plt.plot(self.tp.xmesh,v,'o')
#                plt.show()
#                sys.exit()
#                with open('results.txt','a') as outfile:
#                    rmsd1=np.sqrt(sum((v-np.array([self.tp.gouy_chapman(x)[0] for x in self.tp.xmesh]))**2/len(v)))
#                    rmsd2=np.sqrt(sum((v-np.array([self.tp.gouy_chapman(x,vzeta=self.tp.vzeta_init)[0] for x in self.tp.xmesh]))**2/len(v)))
#                    outfile.write('{} {} {}\n'.format(n*dt,rmsd1,rmsd2))
                for k in range(self.tp.nspecies):
                    if n==1:
                        COLD[k,:]=deepcopy(C[k,:])

                    #Robin BC for concentrations on left side (wall)
                    C[k,0] =\
                        (-2*self.tp.D[k]-self.tp.mu[k]*(v[1]-self.tp.system['vzeta']))/\
                        (-2*self.tp.D[k]+self.tp.mu[k]*(v[1]-self.tp.system['vzeta']))*C[k,1]\
                        -2*self.tp.flux_bound[k,0]*dx/\
                        (-2*self.tp.D[k]+self.tp.mu[k]*(v[1]-self.tp.system['vzeta']))
                    #C[k,0] =\
                    #    (-4*self.tp.D[k]-self.tp.mu[k]*grad_v[0])/\
                    #    (-4*self.tp.D[k]+self.tp.mu[k]*grad_v[0])*C[k,1]\
                    #    -4*self.tp.flux_bound[k,0]*dx/\
                    #    (-4*self.tp.D[k]+self.tp.mu[k]*grad_v[0])

                    #Dirichlet BC for concentrations on right side (bulk)
                    C[k,-1]=C0[(k+1)*nx-1]


                    s = self.tp.D[k]*dt/dx**2  # diffusion number
                    if self.use_lax_friedrich:
                        #add artificial diffusion term for better stability
                        s+=0.5
                    ee = self.tp.charges[k]*self.tp.beta*dt*self.tp.D[k]

                    A=a_matrix(s)
                    B1=b1_matrix(s)
                    if self.tp.use_migration:
                        B1,A=add_field(B1,A,grad_v,lapl_v,ee)
                    B = np.dot(C[k,1:-1],B1)
                    B=add_boundary_values(B,grad_v,s,ee,C[k,0],C[k,-1],COLD[k,0],COLD[k,-1])

                    CTMP = np.linalg.solve(A,B) #this gives vector without initial and final elements
                    C[k,1:-1] = CTMP
                    COLD[k,:]=C[k,:]

                if n % int(nt/float(ntout)) == 0 or n==nt-1: # or True:
                    COUT.append(np.ndarray.flatten(C)) # numpy arrays are mutable, 
                    #so we need to write out a copy of c, not c itself.tp
            return COUT,s

        def integrate_Crank_Nicolson_pnp_old(dx,nx,dt,nt,c,ntout):
            """only for dc_dt=0 at both boundary sites implemented yet"""
            cout = [] # list for storing c arrays at certain time steps
            c0 = np.zeros([self.tp.nspecies])
            c1 = np.zeros([self.tp.nspecies])
            s = np.zeros([self.tp.nspecies])
            ee = np.zeros([self.tp.nspecies])
            for k in range(self.tp.nspecies):
                c0[k] = c[k*nx] # boundary condition on left side
                c1[k] = c[(k+1)*nx-1] # boundary condition on right side
                s[k] = self.tp.D[k]*dt/dx**2  # diffusion number
                ee[k] = self.tp.charges[k]*self.tp.beta*dt*self.tp.D[k]
            # create coefficient matrix:
            def a_matrix(s):
                return diags([-0.5*s, 1+s, -0.5*s], [-1, 0, 1], 
                  shape=(nx-2, nx-2)).toarray() 

            A=a_matrix(s[0])
            Atot=deepcopy(A)
            for k in range(self.tp.nspecies-1):
                Atot=block_diag(Atot,a_matrix(s[k+1]))
            A=Atot
            def b1_matrix(s):
                return diags([0.5*s, 1-s, 0.5*s],[-1, 0, 1], shape=(nx-2,nx-2)).toarray()

            B1=b1_matrix(s[0])
            B1tot=deepcopy(B1)
            for k in range(self.tp.nspecies-1):
                B1tot=block_diag(B1tot,b1_matrix(s[k+1]))
            B1=B1tot

            def unpack(ctmp,c0,c1,nx):
                cout=[]
                j=-1
                for k in range(self.tp.nspecies):
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
                #self.tp.ax1.plot(self.tp.xmesh,c[:nx],'-')
                #if n<50 and abs(self.tp.charges[0])>0.0:
                #    self.tp.ax1.plot(self.tp.xmesh,c[:nx],label='c1')
                #    self.tp.ax1.plot(self.tp.xmesh,c[nx:],label='c2')
                #    self.tp.ax1.legend()
                #if n==50 and abs(self.tp.charges[0])>0.0:
                #    plt.show()
                #    sys.exit()
                cn = c
                #PNP electric field modifications
                self.tp.efield,self.tp.defield_dx=calculate_efield(nx,cn)
                #if n<50 and abs(self.tp.charges[0])>0.0:
                #    self.tp.ax2.plot(self.tp.xmesh,self.tp.efield[:nx],'-')
                #    self.tp.ax3.plot(self.tp.xmesh,self.tp.defield_dx[:nx],'-')
                cn_slice=[cc for ic,cc in enumerate(cn) if ((ic+1)%nx!=0 and (ic)%nx!=0)]

                for i in range((nx-2)*self.tp.nspecies):
                    for j in [i-1,i,i+1]: #range((nx-2)*self.tp.nspecies):
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
                            B1[i,j]+=ee[k]*self.tp.defield_dx[ii+1]
                        if abs(i-j)==1:
                            if j<i:
                                B1[i,j]-=ee[k]*self.tp.efield[ii+1]/4./dx
                                A[i,j]+=ee[k]*self.tp.efield[ii+1]/4./dx
                            elif i<j:
                                B1[i,j]+=ee[k]*self.tp.efield[ii+1]/4./dx
                                A[i,j]-=ee[k]*self.tp.efield[ii+1]/4./dx
                #print '-'*50
                #self.tp.logger.info('\n'.join([''.join(['{:4}'.format(item) for item in row])
                #    for row in A]))
                #self.tp.logger.info('\n'.join([''.join(['{:4}'.format(item) for item in row])
                #    for row in B1]))
                #print '-'*50
                #end electric field modifications

                B = np.dot(cn_slice,B1)
                for k in range(self.tp.nspecies):
                    B[k*(nx-2)] += (0.5*s[k]+\
                        ee[k]*self.tp.efield[0]/4./dx)*(c0[k]+c0[k])
                    B[(k+1)*(nx-2)-1] += (0.5*s[k]-\
                        ee[k]*self.tp.efield[-1]/4./dx)*(c1[k]+c1[k])
                    #B[k*(nx-2)] += 0.5*s[k]*(c0[k]+c0[k])
                    #B[(k+1)*(nx-2)-1] += 0.5*s[k]*(c1[k]+c1[k])
                #if abs(self.tp.charges[0])>0.0 and n==1:
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
                    #so we need to write out a copy of c, not c itself.tp
            return cout,s


        def get_potential_and_gradient(C,dx,nx):
            """Calculates the potential and its gradient (and laplacian=charge density)
                from the PBE by Jacobi relaxation (FD)""" 

            if self.scale_pb_grid is not None:
                #define a finer x grid for PB integration
                if self.scale_pb_grid == 'linear':
                    dx/=self.scale_pb_grid
                    xmesh=np.arange(0,self.tp.xmax+dx,dx)
                elif self.scale_pb_grid == 'log':
                    if min(self.tp.xmesh)<=0:
                        min_x=1e-15
                    else:
                        min_x=min(self.tp.xmesh)
                    xmesh=np.logspace(np.log10(min_x),np.log10(self.tp.xmax-dx),self.tp.nx)
                nx=len(xmesh)
                CNEW=np.zeros([self.tp.nspecies,nx])
                #interpolate concentrations on this grid:
                for k in range(self.tp.nspecies):
                    f = interpolate.interp1d(self.tp.xmesh,C[k,:])
                    CNEW[k,:]=[f(x) for x in xmesh]
                #plt.plot(xmesh,CNEW[0,:],'x',label='interpolated 1')
                #plt.plot(xmesh,CNEW[1,:],'x',label='interpolated 2')
                #plt.plot(self.tp.xmesh,C[0,:],'-',label='before 1')
                #plt.plot(self.tp.xmesh,C[1,:],'-',label='before 2')
                #plt.plot(self.tp.xmesh,np.sum([self.tp.charges[k]*C[k,:] for k in range(self.tp.nspecies)],axis=0)/unit_F,'-',label='charge')
                #plt.legend()
                C=CNEW
                #plt.show()
                #sys.exit()
            bounds=self.tp.pb_bound

            if bounds['gradient']['wall'] is not None and bounds['gradient']['bulk'] is not None:
                self.tp.logger.error('Cannot use two boundary conditions for gradient, quitting.')
                sys.exit()

            def solve_poisson(q,nx,dx,sol0):
                #solve poisson equation with rhs = q
                #rescale xmesh to be linear:
                A=diags([1, -2, 1], [-1, 0, 1],\
                      shape=(nx-2, nx-2)).toarray()
                b=np.array(q[1:nx-1])*dx**2
                b[0]-=sol0[0]
                b[-1]-=sol0[-1]
                x=np.linalg.solve(A,b)
                print 'Check', np.allclose(np.dot(A, x), b)
                sol_final=[sol0[0]]
                for s in x:
                    sol_final.append(s)
                sol_final.append(sol0[-1])
                return sol_final

            def integrate_1d_func(var0,integrand,n=1,inv=False):
                #integrates any function "rhs" n times, var0 is the initial guess for the result
                var_old=deepcopy(var0)
                var=deepcopy(var0)
                if n==2:
                    var=solve_poisson(integrand,nx,self.tp.dx,var0)
                    #interpolate var and 
                    i_step=0
                    #error=np.inf
                    #i_step=0
                    #while error>self.tau_jacobi**2:# or i_step<3:
                    #    i_step+=1
                    #    if inv:
                    #        iterator=reversed(range(1,nx-1))
                    #    else:
                    #        iterator=range(1,nx-1)
                    #    for i in iterator:
                    #        var[i] = 1/2.*(-dx**2*integrand[i]+var_old[i+1]+var_old[i-1])
                    #    error=sum((var_old-var)**2/len(var))
                    #    var_old=deepcopy(var)
                elif n==1:
                    if inv:
                        iterator=reversed(range(1,nx-1))
                    else:
                        iterator=range(1,nx-1)
                    for i in iterator:
                        if inv:
                            var[i] = var[i+1] - integrand[i]*dx
                        else:
                            var[i] = var[i-1] + integrand[i]*dx
                if (n==2):
                    self.tp.logger.info('Converged in {} steps.'.format(i_step))
                return var

            # calculate RHS of PBE
            rhs=np.zeros([nx])
            for i in range(nx):
                for k in range(self.tp.nspecies):
                    rhs[i]-=self.tp.charges[k]*C[k,i]/self.tp.eps
            lapl_v=rhs

            v=np.zeros([nx])
            grad_v=np.zeros([nx])

            if bounds['potential']['wall'] is not None:
                v[0]=bounds['potential']['wall']
            if bounds['potential']['bulk'] is not None:
                v[-1]=bounds['potential']['bulk']
            if bounds['potential']['wall'] is not None and bounds['potential']['bulk'] is not None:
                v=integrate_1d_func(v,rhs,n=2)
                v=np.array(v)
                for i in range(1,nx-1):
                    grad_v[i] = 1./(2*dx)*(v[i+1]-v[i-1])
                grad_v[0]=grad_v[1]+(grad_v[1]-grad_v[2])
                grad_v[-1]=grad_v[-2]+(grad_v[-2]-grad_v[-3])
            else:
                #integrate once to get gradient first
                if bounds['gradient']['wall'] is not None:
                    grad_v[0]=bounds['gradient']['wall']
                    grad_v=integrate_1d_func(grad_v,rhs,n=1)
                    grad_v[-1]=grad_v[-2]+(grad_v[-2]-grad_v[-3])
                if bounds['gradient']['bulk'] is not None:
                    grad_v[-1]=bounds['gradient']['bulk']
                    grad_v=integrate_1d_func(grad_v,rhs,n=1,inv=True)
                    grad_v[0]=grad_v[1]+(grad_v[1]-grad_v[2])
                #extrapolate boundaries and calculate potential
                if bounds['potential']['wall'] is not None:
                    v=integrate_1d_func(v,grad_v,n=1)
                    v[-1]=v[-2]+(v[-2]-v[-3])
                if bounds['potential']['bulk'] is not None:
                    v=integrate_1d_func(v,grad_v,n=1,inv=True)
                    v[0]=v[1]+(v[1]-v[2])

            if self.scale_pb_grid is not None:
                #interpolate the functions back on the normal grid
                def reinterpolate(func):
                    f = interpolate.interp1d(xmesh,func,fill_value='extrapolate')
                    return np.array([f(x) for x in self.tp.xmesh])
                v=reinterpolate(v)
                grad_v=reinterpolate(grad_v)
                lapl_v=reinterpolate(lapl_v)


            #save results
            self.tp.efield=-grad_v
            self.tp.potential=v
            self.tp.total_charge=-lapl_v*self.tp.eps
            return v,grad_v,lapl_v

        def integrate_odeint(dx,nx,dt,nt,c0,ntout):
            """Integrates PNP equations, BCs:
                        concentrations      potential
                left    ROBIN j=0           DIRICHLET vzeta
                right   DIRICHLET cbulk     DIRICHLET 0.0
            """
            def ode_func(t,c,dx,dt):

                if self.oldtime!=t:
                    self.tp.logger.info('Current time = {}'.format(t,max(self.tp.tmesh)))
                self.oldtime=t

                #map concentrations onto 2D array:
                C = np.zeros([self.tp.nspecies,nx])
                for k in range(self.tp.nspecies):
                    C[k,:]=c[k*nx:(k+1)*nx]

                if self.tp.use_migration:
                    #get potential and field
                    v,grad_v,lapl_v=get_potential_and_gradient(C,dx,nx)

                #with open('results.txt','a') as outfile:
                #    rmsd1=np.sqrt(sum((v-np.array([self.tp.gouy_chapman(x)[0] for x in self.tp.xmesh]))**2/len(v)))
                #    rmsd2=np.sqrt(sum((v-np.array([self.tp.gouy_chapman(x,vzeta=self.tp.vzeta_init)[0] for x in self.tp.xmesh]))**2/len(v)))
                #    outfile.write('{} {} {}\n'.format(t,rmsd1,rmsd2))

                DC_DT = np.zeros([self.tp.nspecies,nx])
#                flux = np.zeros([self.tp.nspecies,nx])
#
#
#                #first tabulate D*(dc/dx + mu * c dv/dx) = flux
#                for k in range(0,self.tp.nspecies):
#                    #no flux boundary condition on the left
#                    flux[k,0]=0.0
#                    for i in range(1,nx-1):
#                        dc_dx=(C[k,i+1]-C[k,i-1])/(2.*dx)
##                    if self.tp.lax_friedrich:
##                        corr=dc_dx/2.*dx**2/dt
##                    else:
##                        corr=0.0
#                    corr=0.0
#                    flux[k,i] =\
#                        self.tp.D[k]*\
#                            (\
#                            dc_dx\
#                            +self.tp.beta*self.tp.charges[k]*C[k,i]*grad_v[i]\
#                            )\
#                    #extrapolate flux on the right:
#                    flux[k,-1]=flux[k,-2]+(flux[k,-2]-flux[k,-3])
                corr=0.0

                if self.tp.use_reactions:
                    rates=self.get_rates(C)
                    self.tp.logger_db.debug('x = {}'.format(self.tp.xmesh[0]))

                    for sp in ['CO2','HCO3-','CO32-','OH-','H2','CO','CH4','C2H4','HCOO-','etol','propol','allyl','metol','acet','etgly','K']:
                        for isp2,sp2 in enumerate(self.tp.species):
                            if sp==sp2:
                                break
                        self.tp.logger_db.debug('{} {} {}'.format(sp, rates[isp2,0], C[isp2,0]))
                else:
                    rates=np.zeros([self.tp.nspecies,nx])
                
                for k in range(0,self.tp.nspecies):
                    #set concentration to be constant on the right side
                    DC_DT[k,-1] = 0.0
                    for i in range(0,nx-1):
                        #dflux_dx=(flux[k,i+1]-flux[k,i-1])/(2.*dx)
                        dc_dx=(C[k,i+1]-C[k,i-1])/(2.*dx)
                        dc_dx_2=(C[k,i+1]-2*C[k,i]+C[k,i-1])/(dx**2)
                        if self.tp.use_migration:
                            dcgradv_dx=(C[k,i+1]*grad_v[i+1]-C[k,i-1]*grad_v[i-1])/(2.*dx)
                        else:
                            dcgradv_dx=0.0
                            lapl_v=[0.0]*nx
                            grad_v=[0.0]*nx
                        if i==0:
                            if self.use_lax_friedrich:
                                corr=(C[k,1]-C[k,0])/dt
                            #we have to consider no flux boundary condition here
                            #setting C[k,-1]=C[k,1] considers this:
                            flux=self.tp.flux_bound[k,0] #this is dc/dx at i=0
                            DC_DT[k,i]=\
                                corr+\
                                    (\
                                    #1st part: [dc/dx(i=1)-dc/dx(i=0)] / dx
                                    #2*(C[k,1]-C[k,0]-self.tp.flux_bound[k,0]*2.*dx)/dx**2\
                                    (self.tp.D[k]*((C[k,2]-C[k,0])/(2.*dx)+self.tp.beta*self.tp.charges[k]*\
                                    C[k,1]*grad_v[1])-flux)/dx\
                                    #2nd part
                                    #+self.tp.beta*self.tp.charges[k]*dcgradv_dx\
            #                        +self.tp.beta*self.tp.charges[k]*(\
            #                            C[k,0]*lapl_v[0]+flux*grad_v[0])
                                    )#+\
                                #rates[k,i]
                        else:
                            if self.use_lax_friedrich:
                                corr=dc_dx_2*dx**2/dt/2.
                            #we can directly implement the 2nd derivatives here
                            DC_DT[k,i]=\
                                corr+\
                                self.tp.D[k]*\
                                    (\
                                    dc_dx_2\
                                    +self.tp.beta*self.tp.charges[k]*dcgradv_dx\
                                    )+\
                                rates[k,i]

                    #extrapolate derivative on the left side
#                    DC_DT[k,0]=DC_DT[k,1]+(DC_DT[k,1]-DC_DT[k,2])
                #map dc_dt's onto ouself.tput format
                dc_dt = np.zeros([self.tp.nspecies*nx])
                for k in range(self.tp.nspecies):
                    dc_dt[k*nx:(k+1)*nx]=DC_DT[k,:]
                return dc_dt

            def ode_func_inv(t,c,dx,dt):
                return ode_func(c,t,dx,dt)

            def ode_func_odespy(c,t):
                print 'going', np.shape(c),np.shape(t)
                dx=self.tp.dx
                dt=self.tp.dt
                return ode_func(c,t,dx,dt)

            if self.tp.calc in ['lsoda','odeint']:
                sol,output = integrate.odeint(ode_func_inv, c0, self.tp.tmesh, args=(dx,dt),full_output=True, ml=self.tp.nspecies, mu=self.tp.nspecies)
                self.tp.logger.info('Used minimal time = ',min(output['tcur']),' and maximal time = ',max(output['tcur']))
            elif self.tp.calc == 'odespy':
                solver = odespy.Vode(ode_func_odespy, rtol=0.0, atol=1e-6,
                     adams_or_bdf='adams', order=10)
                solver.set_initial_condition(c0)
                sol, t = solver.solve(self.tp.tmesh)
            else:
                if self.calc_method is not None:
                    r = ode(ode_func).set_integrator(self.tp.calc) #('dopri5') #, method='bdf')
                else:
                    r = ode(ode_func).set_integrator(self.tp.calc,method=self.calc_method,nsteps=10000) #('dopri5') #, method='bdf')
                r.set_initial_value(c0).set_f_params(dx,dt) #c, t, dx, dt
                sol=[]
                while r.successful() and r.t < nt*dt:
                    sol.append(r.integrate(r.t+dt))
                sol=np.array(sol)


            #return the results
            self.tp.logger_db.debug(np.shape(sol))
            cout = []
            for n in range(0,nt):
                if n % int(nt/float(ntout)) == 0 or n==nt-1:
                    cout.append(sol[n,:].copy()) # numpy arrays are mutable, 
            return cout


        def integrate_FTCS(dt,dx,nt,nx,C0,ntout):
            """Integrates PNP equations, BCs:
                        concentrations      
                left    ROBIN j=j0          
                right   DIRICHLET cbulk     
            """

            C = np.zeros([self.tp.nspecies,nx])
            COUT=[]

            #initial conditions for concentrations
            for k in range(self.tp.nspecies):
                C[k,:] = C0[k*nx:(k+1)*nx]

            for n in range(0,nt):
                self.tp.logger.info('time step = {}'.format(n))
                if self.tp.use_migration:
                    v,grad_v,lapl_v=get_potential_and_gradient(C,dx,nx)
                else:
                    v=np.zeros([nx])
                    grad_v=np.zeros([nx])
                    lapl_v=np.zeros([nx])
                rates=self.get_rates(C)
                for k in range(self.tp.nspecies):
                    #get the flux
                    flux=self.tp.flux_bound[k,0] #this is dc/dx at i=0
                    #Robin BC for concentrations on left side (wall)
                    divisor=2*self.tp.D[k]-self.tp.mu[k]*(v[1]-self.tp.system['vzeta'])
                    C[k,0] =(\
                        (2*self.tp.D[k]+self.tp.mu[k]*(v[1]-self.tp.system['vzeta']))*C[k,1]+\
                        flux*2.*dx)/divisor
                    #Dirichlet BC for concentrations on right side (bulk)
                    C[k,-1]=C0[(k+1)*nx-1]
                    temp = np.zeros([nx])
                    temp[0]=C[k,0]
                    temp[-1]=C[k,-1]
                    for i in range(1,nx-1):
                        W = self.tp.D[k]*dt/dx**2-\
                            dt/(2.*dx)*self.tp.mu[k]*grad_v[i+1]+0.5
                        M = -2.*self.tp.D[k]*dt/dx**2
                        E = self.tp.D[k]*dt/dx**2+\
                            dt/(2.*dx)*self.tp.mu[k]*grad_v[i-1]+0.5
                        if not self.use_lax_friedrich:
                            W-=0.5
                            E-=0.5
                            M+=1
                        temp[i] = E * C[k,i-1] + M * C[k,i] + W * C[k,i+1] + rates[k,i]*dt
                    C[k,:]=temp
                if n % int(nt/float(ntout)) == 0 or n==nt-1:
                    COUT.append(np.ndarray.flatten(C))
                #if n==100:
                #    return cout
            return COUT

        def integrate_FTCS_old(dt,dx,nt,nx,V0,ntout):
            # diffusion number (has to be less than 0.5 for the 
            # solution to be stable):
            V = np.zeros([nt,self.tp.nspecies*nx])
            ee=np.zeros([self.tp.nspecies])
            s=np.zeros([self.tp.nspecies])
            for k in range(self.tp.nspecies):
                s[k] = self.tp.D[k]*dt/dx**2
                ee[k] = self.tp.charges[k]*self.tp.beta*dt*self.tp.D[k]
                #print 'before',self.tp.charges[k]*self.tp.beta*dt*self.tp.D[k]
                #ee[k] = 7.62e-8*dt
                #print 'after',ee[k]
                V[:,k*nx] = [V0[k*nx]]*(nt)       #boundary left (
                V[:,(k+1)*nx-1] = [V0[(k+1)*nx-1]]*(nt)     #boundary right 

            #initial conditions
            V[0,:]=V0                   #initia
            cout=[]

            #boundary conditions for change in concentrations
            dV=np.zeros([self.tp.nspecies*nx])
            for k in range(self.tp.nspecies):
                dV[k*nx]=self.tp.dc_dt_bound[k,0]*dt
                dV[(k+1)*nx-1]=self.tp.dc_dt_bound[k,1]*dt

            for n in range(0,nt-1): # time
                self.tp.efield,self.tp.defield_dx=calculate_efield(nx,V[n,:])
                if self.tp.boundary_type=='flux':
                    if sum([fb**2 for fb in self.tp.flux_bound[:,0]])!=0.0:
                        #left bound is defined by flux. this means that right bound will be defined by dc_dt_bound
                        #we have to start integrating from the right
                        #(we cannot use different flux integration for both species)
                        xiter=reversed(range(1,nx-1))
                    else:
                        xiter=range(1,nx-1)
                else:
                    xiter=range(1,nx-1)
                for k in range(self.tp.nspecies):
                    #update concentrations at the boundaries (from boundary condition)
                    V[n+1,k*nx]=V[n,k*nx]+dV[k*nx]
                    V[n+1,(k+1)*nx-1]=V[n,(k+1)*nx-1]+dV[(k+1)*nx-1]
                for k in range(self.tp.nspecies):
                    if self.tp.boundary_type=='flux':
                        flux_bound=np.sqrt(max(self.tp.flux_bound[k,:]**2))
                   # V[n+1,(k+1)*nx-1]=V[n,(k+1)*nx-1]*dV[(k+1)*nx-1]
                    ii=0
                    for i in xiter: #k*nx+1,(k+1)*nx-1): # space
                        ii+=1
                        j=k*nx+i
                        if self.tp.boundary_type=='flux' and ii==nx-2:
                            #define V update including flux boundary
                            V[n+1,j] = V[n,j] + self.tp.D[k]*dt/dx*(\
                                flux_bound-(V[n,j]-V[n,j-1])/dx) +\
                                ee[k]*\
                                (self.tp.efield[i]*(V[n,j+1]-V[n,j-1])/(2.*dx)+\
                                self.tp.defield_dx[i]*V[n,j])
                        else:
                            V[n+1,j] = V[n,j] + s[k]*(V[n,j-1]-\
                                2*V[n,j] + V[n,j+1]) + ee[k]*\
                                (self.tp.efield[i]*(V[n,j+1]-V[n,j-1])/(2.*dx)+\
                                self.tp.defield_dx[i]*V[n,j])
                if n % int(nt/float(ntout)) == 0 or n==nt-1:
                    cout.append(V[n,:].copy())
                #if n==100:
                #    return cout
            return cout

        def integrate_Crank_Nicolson_old(dx,nx,dt,nt,c,ntout):
            """only for dc_dt=0 at both boundary sites implemented yet"""
            cout = [] # list for storing c arrays at certain time steps
            c0 = c[0] # boundary condition on left side
            c1 = c[-1] # boundary condition on right side
            s = self.tp.D[0]*dt/dx**2  # diffusion number
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
                    #so we need to write out a copy of c, not c itself.tp
            return cout,s

        if method not in ['FTCS','Crank-Nicolson']: #in ['vode','lsoda','dopri5','dop853','odeint']:
            cout=integrate_odeint(dx,nx,dt,nt,self.tp.c0,ntout)
            dataplot=cout
        elif method=='FTCS':
            cout=integrate_FTCS(dt,dx,nt,nx,self.tp.c0,ntout)
            dataplot=cout
        elif method=='Crank-Nicolson':
            if self.initialize=='diffusion':
                charges_tmp=self.tp.charges
                self.tp.charges=np.array([0.0]*len(self.tp.charges))
                cout,s=integrate_Crank_Nicolson_pnp(dx,nx,dt,nt,self.tp.c0,ntout)
                self.tp.charges=charges_tmp
                self.tp.logger_db.debug(np.shape(cout))
                self.tp.ax1.plot(self.tp.xmesh,cout[-1][:nx],'-r')
                self.tp.ax1.set_ylim([min(cout[-1][:nx]),max(cout[-1][:nx])])
                cout,s=integrate_Crank_Nicolson_pnp(dx,nx,dt,nt,cout[-1],ntout)
            else:
                cout,s=integrate_Crank_Nicolson(dx,nx,dt,nt,self.tp.c0,ntout)
            dataplot=cout[-1]
        return cout

    def run(self,dt=None,tmax=None):
        if dt!=None:
            self.tp.dt=dt
        if tmax!=None:
            self.tp.tmax=tmax
        self.tp.tmesh=np.arange(0,self.tp.tmax+self.tp.dt,self.tp.dt)

        cout=self.integrate_pnp(self.tp.dx,self.tp.nx,self.tp.dt,\
                len(self.tp.tmesh),self.tp.ntout,method=self.calc)
        return cout
