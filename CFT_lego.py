# -*- coding: utf-8 -*-
"""
Created on Fri Jun  7 17:17:32 2024

@author: User

This file contains everything you need to compute crossing equation.
Furthermore, it also can be used to get z_sampling and c_mean and std of z.
"""
import numpy as np
from scipy.special import hyp2f1
from mpmath import matrix as M

class CFT_lego:
    def __init__(self):
        self.d_phi = 1

    def g_spec(self, d, z):
        return z**d * hyp2f1(d,d,2*d,z)
                
    def G(self, d, z):
        return -(z**(2*self.d_phi)) * self.g_spec(d,1-z) + (1-z)**(2*self.d_phi)*self.g_spec(d,z)
    
    def V(self, z):
        return 2*z - 1
    
    def U(self, z, a1, a3):
        return (2*z - 1) - a1*(z-0.5) - a3*(z-0.5)**3
    
    def c_prediction(self, g_mtx, v_mtx):
        c = M(g_mtx)**-1 * M(v_mtx)
        
        return [float(ele) for ele in c]
    
    def z_sampling_idx(self, num_z_partition, num_sampling, num_states, Type='Global'):
        if Type == 'Global':
            choice_mtx = np.random.choice(num_z_partition, num_sampling, replace=False)
            inds = [[ii*10+ele for ii in range(num_states)] for ele in choice_mtx]
            return inds
        elif Type == 'Local':
            ind = [ii*num_states for ii in range(num_z_partition)][np.random.randint(num_z_partition, size=1)[0]]
            inds = [ind+i for i in range(num_states)]
            return inds

    def c_sampling(self, deltas, p_norm=1, sampling_type='Global', **kwargs):
        num_states = len(deltas)
        num_sampling = 10
        num_z_partition = 10
        z_list = np.linspace(0.4, 0.49, num_z_partition*num_states)
        cs=[]
        sys_list = [[self.G(d, z) for d in deltas] for z in z_list]
        z_sampling_sets = self.z_sampling_idx(num_z_partition, num_sampling, num_states, Type=sampling_type)
        for i in range(num_sampling):
            inds = z_sampling_sets[i]
            # print('inds: ', inds)
            G = [sys_list[ind] for ind in inds]
            if kwargs=={}:
                v=[self.V(z_list[ind]) for ind in inds]
            else:
                v=[self.U(z_list[ind], kwargs['a1'], kwargs['a3']) for ind in inds]
            
            cs.append(self.c_prediction(G, v))
        c_mean=np.mean(cs, axis=0)
        c_std=np.mean(np.abs(cs/c_mean-1)**p_norm)
        # print(np.std(cs/c_mean-1))
        return c_mean, c_std
