# -*- coding: utf-8 -*-
"""
Created on Fri Jun  7 16:17:46 2024

@author: User

1a: There is only 1 value in each alphas.
"""

from CFT_lego import CFT_lego
import cvxpy as cp
import numpy as np
from DG_analytic_v1 import DG
from Where_execute import Where_execute

root_dir = Where_execute()
cft = CFT_lego()

class LP:
    
    def __init__(self, low_delta_list, delta_truncation, delta_spacing, delta_max, 
                 rescale_factor, how_many_constraints):
        self.low_delta_list = low_delta_list
        self.delta_truncation = delta_truncation
        self.delta_spacing = delta_spacing
        self.delta_max = delta_max
        self.rescale_factor = rescale_factor
        
        self.num_grid = int((delta_max-max(low_delta_list))/ delta_spacing)
        self.delta_list = low_delta_list+[delta_truncation+delta_spacing*i for i in range(1, self.num_grid+1)]
        
        dn_constraints = [2*i+1 for i in range(how_many_constraints)]
        
        self.derivative_constraints = [DG("LP_data/isingDGn_v5.csv", self.delta_list, 3, dn) for dn in dn_constraints]*np.exp(-rescale_factor*np.array(self.delta_list))
        self.solver = np.array([2]+[0 for _ in range(len(dn_constraints)-1)])
        
        self.a1_eq = np.array(list(DG("LP_data/isingDGn_v5.csv", self.low_delta_list, 3, 1)*np.exp(-rescale_factor*np.array(self.low_delta_list))) + 
                              [0 for _ in range(1, self.num_grid+1)])
    
    def a1_bound(self):
        # Define and solve the CVXPY problem.
        ope = cp.Variable(len(self.delta_list))
        prob_low = cp.Problem(cp.Minimize(self.a1_eq@ope-2),
                               [self.derivative_constraints @ ope == self.solver, 
                                ope[:len(self.low_delta_list)]>=1e-5,
                                ope>=0])
        
        prob_up = cp.Problem(cp.Maximize(self.a1_eq@ope-2),
                               [self.derivative_constraints @ ope == self.solver, 
                                ope[:len(self.low_delta_list)]>=1e-5,
                                ope>=0])
        
        try:
            prob_low.solve()
            prob_up.solve()
            
            if type(prob_low.value) == np.float64 and type(prob_up.value) == np.float64:
                return [prob_low.value, prob_up.value]
            else:
                return None
        except:
            return None
    
    def alpha_criteria(self):
        a1_range = self.a1_bound()
        if a1_range == None:
            return 'No such bound'
        else:
            return 'Exist bound', a1_range

class LP_data:
    
    def __init__(self, low_delta_list, delta_truncation, delta_spacing, delta_max, rescale_factor):
        self.low_delta_list = low_delta_list
        self.delta_truncation = delta_truncation
        self.delta_spacing = delta_spacing
        self.delta_max = delta_max
        self.rescale_factor = rescale_factor
    
    def stored_a_bound(self, how_many_constraints):
        import os
        import json
        file_name = f'LP_data/a_region_constraint_num{how_many_constraints}_with_alphabound.json'
        
        d_property = {'delta':self.low_delta_list, 'd_truncate':self.delta_truncation, 
                      'd_space':self.delta_spacing, 'd_max':self.delta_max, 'd_rf':self.rescale_factor}
        root_dir=Where_execute()
        
        if os.path.isfile(root_dir+file_name):
            with open(os.path.join(root_dir, file_name), 'r') as file:
                combined_data = json.load(file)
            # 将 target_label 转换为 JSON 字符串
            target_label_str = json.dumps(d_property)
            
            # 获取对应的值
            values = combined_data.get(target_label_str)
            
            if values is None:
                print('Building new alpha json element!')
                lp = LP(self.low_delta_list, self.delta_truncation, 
                        self.delta_spacing, self.delta_max, self.rescale_factor, how_many_constraints)
                a_bound = lp.alpha_criteria()
                alpha_container = {'property':a_bound[0], 'alpha bound': a_bound[1]}
                # 创建以 d_property 标签为键的字典，值为 alpha_container
                combined_data[json.dumps(d_property)] = alpha_container
                
                # Write to JSON file
                with open(os.path.join(root_dir, file_name), 'w') as file:
                    json.dump(combined_data, file, indent=4)
                
                return alpha_container
            else:
                return values
        else:
            print('Building new alpha json file!')
            combined_data={}
            lp = LP(self.low_delta_list, self.delta_truncation, 
                    self.delta_spacing, self.delta_max, self.rescale_factor, how_many_constraints)
            a_bound = lp.a_region()
            alpha_container = {'property':a_bound[0], 'alpha bound': a_bound[1]}
            # 创建以 d_property 标签为键的字典，值为 alpha_container
            combined_data[json.dumps(d_property)] = alpha_container
            
            # Write to JSON file
            with open(os.path.join(root_dir, file_name), 'w') as file:
                json.dump(combined_data, file, indent=4)
            return alpha_container
