# BFI for Survival

import numpy as np
import matplotlib.pyplot as plt
import numpy.random as rnd
import scipy.optimize as opt 

n = 100
m = 1
h0 = 0.1
nc = 4
dt =  0.2
S = 1000

def ML(C,T,L,dt):
    #define discrete grid
    nt = int(L/dt)
    tp = np.arange(0,nt,1)*dt
    O = np.empty((nt,2))
    #quicksort alg 
    def qs(x,l,r):
        idx = int(np.floor((l+r)/2))
        if(idx ==0):
            return 0 
        if(idx ==r):
            return int(r)
        if(x>r*dt):
            return int(r)
        if(x>(idx+1)*dt):
            idx = qs(x,(idx+1),r)
        if(x<idx*dt):
            idx = qs(x,l,idx)
        if((x<=(idx+1)*dt) and (x>=idx*dt)):
            return int(idx)
    v_qs = np.vectorize(qs)

    #use quicksort to assign each event time to a time interval 
    d = v_qs(T,0,(nt-1))

    #compute the death and risk indicators
    D = np.zeros((len(T),nt))
    R = np.ones((len(T),nt))
    for i in range(len(T)):
        D[i,d[i]] = 1
        R[i,:] = np.array(tp<T[i],int)

    #compute ML 
    h_ml = np.empty(nt)
    R_p = np.ones(len(T))@R
    #print(np.shape(C))
    #print(np.shape(D))
    D_p = np.ones(len(T))@D
    F_p = C@D
    for k in range(nt):
        if(R_p[k]-D_p[k]>0):
            h_ml[k] = (np.log(R_p[k] + F_p[k]- D_p[k]) - np.log(R_p[k] - D_p[k]))/dt
        else:
            h_ml[k] = 10000#np.inf

    #compute fisher info
    I = np.empty(nt)
    for k in range(nt):
        if(R_p[k]-D_p[k]>0):
            I[k] = 0.25*(dt*dt)*F_p[k]*(1.0-np.tanh(0.5*dt*h_ml[k])**2)
        else:
            I[k] = 0.0
    O[:,0] = h_ml
    O[:,1] = I
    return O

alpha = 3
nt = int(alpha/dt)
BFI = np.zeros((S,nt-1))
MLE = np.zeros((S,nt-1))

for counter in range(S):
    Tc = np.empty((nc,n))
    Cc = np.empty((nc,n))
    T_max = 0 
    for k in range(nc):
        T1 = np.exp(0.5*(np.log(-np.log(rnd.random(size = n)))+np.log(2)))
        T0 = alpha*np.ones(n)
        Tc[k,:] = np.minimum(T1,T0)
        Cc[k,:] = np.array(T1<T0,int)
        #T_max = max(T_max,max(Tc[k,:]))

    h_BFI = np.zeros(nt-1)
    I_BFI = 1.0e-13*np.ones(nt-1)
    for k in range(nc):
        Mk = ML(Cc[k,:],Tc[k,:],alpha,dt)
        #print(Mk[:-1,0])
        #a = Mk[:,0]!=np.inf
        #a = np.array(a,int)
        h_BFI = h_BFI + Mk[:-1,1]*Mk[:-1,0]
        I_BFI = I_BFI + Mk[:-1,1]
        #print(Mk[:-1,0])
        #print(np.shape(I_BFI))
    #print(h_BFI)
    h_BFI = h_BFI/I_BFI
    #print(h_BFI)
    Tk = Tc.reshape(nc*n)
    #print(len(Tk))
    Ck = Cc.reshape(nc*n)
    #print(len(Ck))
    M = ML(Ck,Tk,alpha,dt)
    h = M[:-1,0]
    BFI[counter,:] = h_BFI
    MLE[counter,:] = h

MSE = np.sqrt(np.sum((MLE-BFI)**2,axis = 0)/S)

t = np.linspace(0.0,alpha,nt-1)
MSE_BFI_true = np.sqrt(np.sum((BFI-t)**2,axis = 0)/S)
MSE_MLE_true = np.sqrt(np.sum((MLE-t)**2,axis = 0)/S)


plt.figure()
plt.plot(t,MSE,label = 'MLE-BFI')
plt.plot(t,MSE_BFI_true,label = 'BFI-true')
plt.plot(t,MSE_MLE_true,label = 'MLE-true')
plt.xlabel(r'$t$')
plt.ylabel(r'${\rm MSE}$')
plt.legend()
plt.savefig('MSE.png')
    

plt.figure()
plt.title(r'$h_T$' + '  n ='+str(n)+' k =' + str(nc))
plt.plot(t,BFI.transpose(),'k.',label = 'BFI')
plt.plot(t,MLE.transpose(),'bx',label = 'full data-set')
plt.plot(t,t,'r-',label = 'true')
#plt.legend()
plt.savefig('hazard_rate.png')

H_BFI = np.empty((S,nt-1))
H = np.empty((S,nt-1))
for l in range(S):
    H_BFI[l,:] = np.array([dt*sum(BFI[l,:k]) for k in range(len(h_BFI))])
    H[l,:] = np.array([dt*sum(MLE[l,:k]) for k in range(len(h))]) 

plt.figure()
plt.title(r'$S_T$' + '  n ='+str(n)+' k =' + str(nc))
plt.plot(t,np.exp(-H_BFI.transpose()), 'k.',label = 'BFI')
plt.plot(t,np.exp(-H.transpose()),'bx',label = 'full data-set')
plt.plot(t,np.exp(-0.5*t*t),'r-',label = 'true')
#plt.legend()
plt.savefig('surv.png')
plt.show()
