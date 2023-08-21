#!/usr/bin/env python
# coding: utf-8

# In[1]:


import readfof
import numpy as np

from scipy.interpolate import interp1d
from scipy import integrate
import matplotlib.pyplot as plt#约定俗成的写法plt

from scipy.optimize import curve_fit
from scipy import special

from scipy.integrate import odeint,quad,dblquad,tplquad,simps
import scipy.special as special
from numpy.polynomial import polynomial, legendre
import sympy as sym
from colossus.cosmology import cosmology


from mcfit import xi2P
from scipy import linalg


def LegendrePolynomials(N,x):
    """
    N 为勒让德多项式的阶数
    x 为自变量 sym.symbol对象
    """
    if N == 0:
        return 1
    if N == 1:
        return x
    p0 = LegendrePolynomials(0,x)
    p1 = LegendrePolynomials(1,x)
    assert N>=2
    for i in range(1,N):
        p = (2*i+1)/(i+1)*x*p1 - i/(i+1)*p0
        p0 = p1
        p1 = p
    return sym.simplify(p1)


# In[2]:


##读取暗物质功率谱 
params = {'flat': True, 'H0': 67.72, 'Om0': 0.315, 'Ob0': 0.049, 'sigma8': 0.82, 'ns': 0.95}
cosmology.addCosmology('myCosmo', **params)
cosmo = cosmology.setCosmology('myCosmo')
k=10**np.linspace(-4,0.5,10000)
Pnl0 = cosmo.matterPowerSpectrum(k,z=0)
Pdm=interp1d(k,Pnl0,kind='cubic')


# In[3]:


rlist=np.linspace(25,120,20)
thetalist=np.linspace(0.1+19*0.15,0.1,20)
coslist=np.cos(thetalist)

legendre0list=np.zeros(20)
legendre2list=np.zeros(20)
legendre4list=np.zeros(20)
for i in range(20):
    legendre0list[i]=LegendrePolynomials(0,coslist[i])
    legendre2list[i]=LegendrePolynomials(2,coslist[i])
    legendre4list[i]=LegendrePolynomials(4,coslist[i])





dataori=np.fromfile('/data1/quijote/lzy_pos/result_rsd2pcf/qjthalorsd2pcf0.txt',dtype=np.float32)
dataori=dataori.reshape(20,20)

##calculate multipoles for samples
nhs=350
sample_multipole0=np.zeros((nhs,20))
sample_multipole2=np.zeros((nhs,20))
sample_multipole4=np.zeros((nhs,20))

for n in range(nhs):

    dataori=np.fromfile('/data1/quijote/lzy_pos/result_rsd2pcf/qjthalorsd2pcf'+str(n)+'.txt',dtype=np.float32)
    dataori=dataori.reshape(20,20)

    multipole0=np.zeros(20)
    multipole2=np.zeros(20)
    multipole4=np.zeros(20)

    for i in range(20):
        ##calculate multipole through a simpson intergrate
        ##simpsy
        simpsy=dataori[i]
        simpsx=coslist

        ##intergrate
        multipole0[i]=simps(simpsy*legendre0list,simpsx)*(2*0+1)/2
        multipole2[i]=simps(simpsy*legendre2list,simpsx)*(2*2+1)/2
        multipole4[i]=simps(simpsy*legendre4list,simpsx)*(2*4+1)/2

    sample_multipole0[n]=multipole0
    sample_multipole2[n]=multipole2
    sample_multipole4[n]=multipole4


# In[4]:


sample_total=np.zeros((nhs,40))
for i in range(nhs):
    sample_total[i]=np.append(sample_multipole0[i],sample_multipole2[i])

sample_total_avrg=np.zeros(40)
for i in range(nhs):
    sample_total_avrg=sample_total_avrg+sample_total[i]
sample_total_avrg=sample_total_avrg/nhs

##Rijtest
cov_matrix=np.zeros((40,40))
##再读取300个halo样本里的2pcf数据
for i in range(40):
    for j in range(40):
        sum_covij=0
        for k in range(nhs):
            total_data=sample_total[k]
            sum_covij=sum_covij+(total_data[i]-sample_total_avrg[i])*(total_data[j]-sample_total_avrg[j])

        cov_matrix[i][j]=(sum_covij/(nhs-1))

ncov=cov_matrix
invcov=linalg.inv(ncov)



# In[101]:


##th rsd2pcf

from nbodykit.lab import cosmology
import matplotlib.pyplot as plt
import numpy as np
from mcfit import P2xi
from scipy.integrate import quad 
from scipy.integrate import simps 




c = cosmology.Planck15
c=c.clone(P_k_max=10)
c=c.clone(Omega0_b=0.048)
c=c.clone(Omega0_cdm=(0.31-0.048))
c=c.clone(h=0.6777)
c=c.clone(n_s=0.965)
c=c.match(sigma8=0.834)




Plin = cosmology.LinearPower(c, redshift=0.)
Pnl = cosmology.HalofitPower(c, redshift=0)
PZel = cosmology.ZeldovichPower(c, redshift=0)
cf_nl = cosmology.CorrelationFunction(Pnl)
cf_zel = cosmology.CorrelationFunction(PZel)
cf_lin=cosmology.CorrelationFunction(Plin)


rl = np.logspace(-1, np.log10(200), 10000)
y=cf_zel(rl)

from scipy.interpolate import interp1d
cfzel1 = interp1d(rl, y, kind="cubic")

def twopcf1d(r):
    return cfzel1(r)


def r2_multiple_twopcf1d(r):
    return (r**2)*twopcf1d(r)

def r4_multiple_twopcf1d(r):
    return (r**4)*twopcf1d(r)


def barred_1(r):
    simp_x=np.linspace(0.1,r,100)
    simp_y=r2_multiple_twopcf1d(simp_x)
    return (3/(r**3))*simps(simp_y,simp_x)
    #return (3/(r**3))*quad(r2_multiple_twopcf1d,0.01,r)[0]

def barred_2(r):
    simp_x=np.linspace(0.1,r,100)
    simp_y=r4_multiple_twopcf1d(simp_x)
    return (5/(r**5))*simps(simp_y,simp_x)
    #return (5/(r**5))*quad(r4_multiple_twopcf1d,0.01,r)[0]









def multipole_0(r,f,b):
    bias=b
    sigma_8=0.834
    beta=f/b

    #return (1+(2/3)*0.519623*sigma_8*1+(1/5)*((0.519623*sigma_8)**2))*twopcf1d(r)/(sigma_8**2)
    return ((sigma_8*bias)**2+(2/3)*f*bias*sigma_8**2+(1/5)*((f*sigma_8)**2))*twopcf1d(r)/(sigma_8**2)
    ##return twopcf1d(r)


def multipole_2(r,f,b):
    bias=b
    sigma_8=0.834
    beta=f/b
    return ((4/3)*f*bias*sigma_8**2+(4/7)*(f*sigma_8)**2)*(twopcf1d(r)/(sigma_8**2)-barred_1(r)/(sigma_8**2))
    ##return (twopcf1d(r)-barred_1(r))

def multipole_4(r,f,b):
    bias=b
    sigma_8=0.834
    beta=f/b
    return (8/35)*((f*sigma_8)**2)*(twopcf1d(r)/(sigma_8**2)+2.5*barred_1(r)/(sigma_8**2)-3.5*barred_2(r)/(sigma_8**2))
    ##return (twopcf1d(r)+2.5*barred_1(r)-3.5*barred_2(r))

def thrsd2pcf(r,mnu,f,b):
    multipole0=multipole_0(r,f,b)
    multipole2=multipole_2(r,f,b)
    multipole4=multipole_4(r,f,b)
    return multipole0*LegendrePolynomials(0,mnu)+multipole2*LegendrePolynomials(2,mnu)+multipole4*LegendrePolynomials(4,mnu)

def thrsd2pcf_s(s_par,s_ver,f,b):
    r=np.sqrt(s_par**2+s_ver**2)
    mu=s_par/r
    return thrsd2pcf(r,mu,f,b)


# In[125]:


##add non-linear effect
##distribution function of random pairwise velocities
def random_pv(v,sigmav):
    return (1/(np.sqrt(2)*sigmav))*np.exp(-np.sqrt(2)*np.abs(v)/sigmav)

def thrsd2pcf_s_nl(s_par,s_ver,f,b,sigmav):
    def inner_f(v):
        return random_pv(v,sigmav)*thrsd2pcf_s(s_par-v/100,s_ver,f,b)
    

    return quad(inner_f,-2200,2200,epsabs=1e-1)[0]

def thrsd2pcf_nl(r,mnu,f,b,sigmav):
    is_par=r*mnu
    is_ver=r*np.sqrt(1-mnu**2)
    return thrsd2pcf_s_nl(is_par,is_ver,f,b,sigmav)


# In[126]:


thrsd2pcf_nl(30,1,0.5,2,300)


# In[127]:


def LegendrePolynomials(N,x):

    if N == 0:
        return 1
    if N == 1:
        return x
    p0 = LegendrePolynomials(0,x)
    p1 = LegendrePolynomials(1,x)
    assert N>=2
    for i in range(1,N):
        p = (2*i+1)/(i+1)*x*p1 - i/(i+1)*p0
        p0 = p1
        p1 = p
    return sym.simplify(p1)

def multipole_0_nl(r,f,b,sigmav):
    def inner_0(mnu):
        return thrsd2pcf_nl(r,mnu,f,b,sigmav)*LegendrePolynomials(0,mnu)
    simpsx=np.linspace(-1,1,20)
    simpsy=np.zeros(20)
    for i in range(20):
        simpsy[i]=inner_0(simpsx[i])
    return simps(simpsy,simpsx)*(2*0+1)/2

def multipole_2_nl(r,f,b,sigmav):
    def inner_2(mnu):
        return thrsd2pcf_nl(r,mnu,f,b,sigmav)*LegendrePolynomials(2,mnu)
    return quad(inner_2,-1,1,epsabs=1e-3)[0]*(2*2+1)/2

def multipole_4_nl(r,f,b,sigmav):
    def inner_4(mnu):
        return thrsd2pcf_nl(r,mnu,f,b,sigmav)*LegendrePolynomials(4,mnu)
    return quad(inner_4,-1,1,epsabs=1e-3)[0]*(2*4+1)/2
    


# In[ ]:





# In[128]:


multipole_2_nl(30,0.5,2,300)


# In[69]:





# In[129]:


r_log=np.linspace(25,25+19*5,20)

rsd2pcf_mock=np.zeros(40)
rsd2pcf_mock=sample_total_avrg

def loglike(theta):
    f,b,sigmav=theta
    flagup = theta>np.array([1, 5,1000])
    flagdown = theta<np.array([0, 1,0])
    if(sum(flagup) + sum(flagdown) != 0): return -np.inf
    th_vector=np.zeros(40)
    for i in range(20):
        th_vector[i]=multipole_0(r_log[i],f,b,sigmav)
        th_vector[i+20]=multipole_2(r_log[i],f,b,sigmav)
        
    errorvec=(th_vector-rsd2pcf_mock)



    chi2=np.matmul(errorvec,np.matmul((1-41/(nhs-1))*invcov,errorvec))

    return -0.5*chi2


# In[ ]:





# In[130]:





# In[115]:


import sys

import numpy as np
import pandas as pd
import os
import emcee
import time 
from multiprocessing import Pool

os.environ["OMP_NUM_THREADS"] = "1"
nsteps = 8000

savefile = 'rsd2pcf2par.h5'
initial = np.random.random([4, 2]) * 0.01 + np.array([0.6,2.4])
nwalkers, ndim = initial.shape
backend = emcee.backends.HDFBackend(savefile)
backend.reset(nwalkers, ndim)
sampler = emcee.EnsembleSampler(nwalkers, ndim, loglike, pool=Pool(processes=96),backend=backend,moves=emcee.moves.StretchMove()) 

sampler.run_mcmc(initial, nsteps)

