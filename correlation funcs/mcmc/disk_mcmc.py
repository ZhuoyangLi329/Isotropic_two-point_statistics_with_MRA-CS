#!/usr/bin/env python
# coding: utf-8

# In[2]:


from nbodykit.lab import cosmology
import matplotlib.pyplot as plt
import numpy as np
from mcfit import P2xi
from scipy.integrate import quad 
from scipy.integrate import simps 
from scipy import special

##parameter set
sigma_8=0.834
bias=2.2
beta=0.5196231/bias


c = cosmology.Planck15
c=c.clone(P_k_max=10)
c=c.clone(Omega0_b=0.049)
c=c.clone(Omega0_cdm=(0.3175-0.049))
c=c.clone(h=0.6711)
c=c.clone(n_s=0.9624)
c=c.match(sigma8=sigma_8)




Plin = cosmology.LinearPower(c, redshift=0.)
Pnl = cosmology.HalofitPower(c, redshift=0)
PZel = cosmology.ZeldovichPower(c, redshift=0)
cf_nl = cosmology.CorrelationFunction(Pnl)
cf_zel = cosmology.CorrelationFunction(PZel)
cf_lin=cosmology.CorrelationFunction(Plin)


# In[3]:


rl = np.logspace(-1, np.log10(200), 1000000)
y=cf_zel(rl)

from scipy.interpolate import interp1d
cfzel1 = interp1d(rl, y, kind="cubic")



# In[5]:


twopcf1d=cfzel1


# In[6]:


def r2_multiple_twopcf1d(r):
    return (r**2)*twopcf1d(r)

def r4_multiple_twopcf1d(r):
    return (r**4)*twopcf1d(r)


def barred_1(r):
    simp_x=np.linspace(0.1,r,20)
    simp_y=r2_multiple_twopcf1d(simp_x)
    return (3/(r**3))*simps(simp_y,simp_x)
    #return (3/(r**3))*quad(r2_multiple_twopcf1d,0.01,r)[0]

def barred_2(r):
    simp_x=np.linspace(0.1,r,20)
    simp_y=r4_multiple_twopcf1d(simp_x)
    return (5/(r**5))*simps(simp_y,simp_x)
    #return (5/(r**5))*quad(r4_multiple_twopcf1d,0.01,r)[0]


# In[7]:


##parameter f,b
parameter_fiducial=[0.51,2.32]


# In[8]:


def fofz(z):
    return 0.519623


def multipole_0(r,parameter):
    bias=parameter[1]
    f=parameter[0]
    
    #return (1+(2/3)*0.519623*sigma_8*1+(1/5)*((0.519623*sigma_8)**2))*twopcf1d(r)/(sigma_8**2)
    return ((sigma_8*bias)**2+(2/3)*f*bias*sigma_8**2+(1/5)*((f*sigma_8)**2))*twopcf1d(r)/(sigma_8**2)
    ##return twopcf1d(r)


def multipole_2(r,parameter):
    bias=parameter[1]
    f=parameter[0]
    return ((4/3)*f*bias*sigma_8**2+(4/7)*(f*sigma_8)**2)*(twopcf1d(r)/(sigma_8**2)-barred_1(r)/(sigma_8**2))
    ##return (twopcf1d(r)-barred_1(r))

def multipole_4(r,parameter):
    bias=parameter[1]
    f=parameter[0]
    return (8/35)*((f*sigma_8)**2)*(twopcf1d(r)/(sigma_8**2)+2.5*barred_1(r)/(sigma_8**2)-3.5*barred_2(r)/(sigma_8**2))
    ##return (twopcf1d(r)+2.5*barred_1(r)-3.5*barred_2(r))


# In[9]:


from numpy.polynomial import polynomial, legendre

import sympy as sym
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



def thrsd2pcf(r,mnu,parameter):
    multipole0=multipole_0(r,parameter)
    multipole2=multipole_2(r,parameter)
    multipole4=multipole_4(r,parameter)
    return float(multipole0*LegendrePolynomials(0,mnu)+multipole2*LegendrePolynomials(2,mnu)+multipole4*LegendrePolynomials(4,mnu))


# In[10]:


def pcf_cyls(rho,mnu,parameter):
    H=rho*mnu
    R=rho*np.sqrt(1-mnu**2)


    def innerfunc(z):
        return thrsd2pcf(np.sqrt(R**2+z**2),z/(np.sqrt(R**2+z**2)),parameter)

    nsimps=15
    simpsx=np.linspace(-H,H,nsimps)
    simpsy=np.zeros(nsimps)
    for i in range(nsimps):
        simpsy[i]=innerfunc(simpsx[i])
    return simps(simpsy,simpsx)/(2*H)



def pcf_disk(rho,mnu,parameter):
    H=rho*mnu
    R=rho*np.sqrt(1-mnu**2)
    def innerfunc(r):
        return thrsd2pcf(np.sqrt(H**2+r**2),H/(np.sqrt(H**2+r**2)),parameter)*r

    nsimps=15
    simpsx=np.linspace(0,R,nsimps)
    simpsy=np.zeros(nsimps)
    for i in range(nsimps):
        simpsy[i]=innerfunc(simpsx[i])
    return 2*simps(simpsy,simpsx)/(R**2)


# In[11]:


##multipole analysis
rlist=np.linspace(50,50+12*9,10)
coslist=(np.linspace(-0.98,-0.98+0.103*19,20))
thetalist=np.arccos(coslist)


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


# In[12]:


def pcf_cyls_multipole(rho,l,parameter):
    def innercyls(mnu):
        return pcf_cyls(rho,mnu,parameter)*LegendrePolynomials(l,mnu)
    
    nsimps=20
    simpsx=np.linspace(-0.98,0.98,nsimps)
    simpsy=np.zeros(nsimps)
    for i in range(nsimps):
        simpsy[i]=innercyls(simpsx[i])
    return simps(simpsy,simpsx)*(2*l+1)/2



def pcf_disk_multipole(rho,l,parameter):
    def innerdisk(mnu):
        return pcf_disk(rho,mnu,parameter)*LegendrePolynomials(l,mnu)
    
    nsimps=20
    simpsx=np.linspace(-0.98,0.98,nsimps)
    simpsy=np.zeros(nsimps)
    for i in range(nsimps):
        simpsy[i]=innerdisk(simpsx[i])
    return simps(simpsy,simpsx)*(2*l+1)/2

    


# In[38]:


legendre0list=np.zeros(20)
legendre2list=np.zeros(20)
legendre4list=np.zeros(20)
for i in range(20):
    legendre0list[i]=LegendrePolynomials(0,coslist[i])
    legendre2list[i]=LegendrePolynomials(2,coslist[i])
    legendre4list[i]=LegendrePolynomials(4,coslist[i])



nhs=120
sample_multipole0=np.zeros((nhs,10))
sample_multipole2=np.zeros((nhs,10))
sample_multipole4=np.zeros((nhs,10))

for n in range(nhs):

    dataori=np.fromfile('/data1/quijote/lzy_pos/res_disk/4/qjthalodisk'+str(n)+'.txt',dtype=np.float64)
    dataori=dataori.reshape(10,20)


    multipole0=np.zeros(10)
    multipole2=np.zeros(10)
    multipole4=np.zeros(10)

    for i in range(10):
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


##multipole avrg
sample_multipole0_avrg=np.zeros(10)
for n in range(nhs):
    sample_multipole0_avrg+=sample_multipole0[n]
sample_multipole0_avrg/=nhs

sample_multipole2_avrg=np.zeros(10)
for n in range(nhs):
    sample_multipole2_avrg+=sample_multipole2[n]
sample_multipole2_avrg/=nhs

sample_multipole4_avrg=np.zeros(10)
for n in range(nhs):
    sample_multipole4_avrg+=sample_multipole4[n]
sample_multipole4_avrg/=nhs











sample_total=np.zeros((nhs,20))
for i in range(nhs):
    sample_total[i]=np.append(sample_multipole0[i],sample_multipole2[i])

sample_total_avrg=np.zeros(20)
for i in range(nhs):
    sample_total_avrg=sample_total_avrg+sample_total[i]
sample_total_avrg=sample_total_avrg/nhs

##Rijtest
cov_matrix=np.zeros((20,20))
##再读取300个halo样本里的2pcf数据
for i in range(20):
    for j in range(20):
        sum_covij=0
        for k in range(nhs):
            total_data=sample_total[k]
            sum_covij=sum_covij+(total_data[i]-sample_total_avrg[i])*(total_data[j]-sample_total_avrg[j])

        cov_matrix[i][j]=(sum_covij/(nhs-1))

ncov=cov_matrix
invcov=np.linalg.inv(ncov)


    


# In[29]:


def loglike_combine(theta):
    f,b=theta
    flagup = theta>np.array([1, 5])
    flagdown = theta<np.array([0, 1])
    if(sum(flagup) + sum(flagdown) != 0): return -np.inf
    th_vector=np.zeros(20)
    ipara=[f,b]
    for i in range(10):
        th_vector[i]=pcf_disk_multipole(rlist[i],0,ipara)
        th_vector[i+10]=pcf_disk_multipole(rlist[i],2,ipara)
        #th_vector[i+20]=pcf_disk_multipole(rlist[i],4,ipara)



        ##print(np.sum((th_vector-sigma_avrg)**2))
        ##print(np.min((th_vector/sigma_avrg)))
        ##print(np.max((th_vector/sigma_avrg)))
        
        ##matmul
    errorvec=(th_vector-sample_total_avrg)
    print(th_vector/sample_total_avrg)



            
    chi2=np.matmul(errorvec,np.matmul((1-11/119)*invcov,errorvec))


        

    return -0.5*chi2


# In[39]:


import sys

import numpy as np
import pandas as pd
import os
import emcee
import time 
from multiprocessing import Pool

os.environ["OMP_NUM_THREADS"] = "1"
nsteps = 8000

savefile = 'disk.h5'
initial = np.random.random([4, 2]) * 0.01 + np.array([0.4,2.6])
nwalkers, ndim = initial.shape
backend = emcee.backends.HDFBackend(savefile)
backend.reset(nwalkers, ndim)
sampler = emcee.EnsembleSampler(nwalkers, ndim, loglike_combine, pool=Pool(processes=128),backend=backend,moves=emcee.moves.StretchMove()) 

sampler.run_mcmc(initial, nsteps,progress=True)

# In[ ]:




