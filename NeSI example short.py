# -*- coding: utf-8 -*-
"""
Created on Thu May 25 13:51:47 2017

@author: scho849
"""

from pylab import *
import copy
from scipy.stats import *

eta = 0                                            
eta_count = 0
eta_max = 6
eta_steps = len(arange(0,eta_max,0.1)) + 1

t_steps = 100           
run_max = 5                 


phi_store = np.zeros((run_max,eta_steps,t_steps+1))
coh_store = copy.deepcopy(phi_store)



while eta <= eta_max:
    run = 1
    runmax = run_max

    while run <= runmax:                    
        L = 10                              
        N = 40
        radius = 1
        dt = 1
        
        coord = np.random.rand(N,2) * L
        angle = np.random.uniform(-pi,pi,N)
        dimensions = np.array([L,L])
        
        velocity = np.zeros((N,2))
        distance = np.zeros((N,N))
        angle_temp = np.zeros((N,1))
        
        steps = 0.1*eta
        theta_0 = 0
        A, B = 0.1, 0
        
        vector_sum = np.zeros((N,2))
        vec_norm = np.zeros((2,1))
        
        t = 0
        num_steps = t_steps
        
        x = 0
        
        while t <= num_steps:
        ########################## t = t ########################
            dtheta = np.random.uniform(-eta/2,eta/2,N)   
            for i in range(N):
                velocity[i,0] = norm.rvs(scale = steps)
                velocity[i,1] = norm.rvs(scale = steps)
                
                delta = np.abs(coord[i] - coord)
                delta = np.where(delta > 0.5 * dimensions, delta - dimensions, delta)
                dist = np.sqrt((delta**2).sum(axis = -1))
                distance[i] = dist
                
                sin_sum = 0
                cos_sum = 0
                num = 0
                
                for j in range(N):
                    if distance[i,j] < radius:
                        sin_sum += np.sin(angle[j])
                        cos_sum += np.cos(angle[j])
                        num += 1
            
                sin_sum /= num
                cos_sum /= num
\
                vector_sum[i,0] = cos_sum + A*np.cos(theta_0) + B*np.cos(angle[i])
                vector_sum[i,1] = sin_sum + A*np.sin(theta_0) + B*np.sin(angle[i])
                angle_temp[i] = math.atan2(vector_sum[i,1],vector_sum[i,0]) + dtheta[i]

                
            vec_sum = array([np.cos(angle).sum(axis=0),np.sin(angle).sum(axis=0)])
            vec_mod = sqrt((vec_sum**2).sum(axis=0))
            phi = vec_mod/(N)
            
            vec_0 = array([np.cos(theta_0),np.sin(theta_0)])
            coherence = np.dot(vec_0, vec_sum) / N
        
            phi_store[run-1,eta_count,t] = phi
            coh_store[run-1,eta_count,t] = coherence
            
            if x == 100:
                print("eta=",eta,"run=",run,"t=",t,"theta_a=",phi,"theta_a_0=",coherence)
                x = 0        
            
        ################### t = t + 1 ##########################
            if t < num_steps:
                angle = copy.deepcopy(angle_temp)
                coord += velocity
                
                for i in range(N):
                    for j in range(2):
                        if coord[i,j] > L:
                            coord[i,j] += -L
                        if coord[i,j] < L:
                            coord[i,j] += L
        
            
            
            t += 1
            x += 1
        run += 1
    eta += 1
    eta_count += 1
    
np.savetxt('phi N=40 A=0.1.txt',phi_store.sum(axis=0))
np.savetxt('coh N=40 A=0.1.txt',coh_store.sum(axis=0))
