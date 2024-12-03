# -*- coding: utf-8 -*-
"""
Created on Fri Jun 14 15:48:03 2024

@author: User

Status: Finish

This file contains two main functions: 
    DG: 1. If deltas existed in the file will directly output correspond value.
        2. If deltas don't exist, it will first call mathematica_dG to calc and store the value then output.
        
    DG_new_derivative_constraint: If there're new constraints it will calc all existed deltas and store in a new file.
"""
import numpy as np
import pandas as pd
from Where_execute import Where_execute

root_dir = Where_execute()

def int_to_str(int_num, digit):
    return f'{round(float(int_num), digit)}'

def mathematica_dG(delta, digit, fconstraintnum, lconstraintnum):
    # Call constraintnum by their derivative order
    from wolframclient.evaluation import WolframLanguageSession
    from wolframclient.language import wl, wlexpr
    print('Activate mathematica')
    
    if type(delta) == list:
        delta_m = '{'+str([int_to_str(i, digit) for i in delta]).replace('[', '').replace(']', '').replace("'", "")+'}'
        # 启动一个新的Wolfram Language会话
        with WolframLanguageSession(timeout=300) as session:
            # 载入包含Mathematica函数定义的文件
            session.evaluate(wl.Get(root_dir+"DG(delta)_v1.m"))
            # 在Mathematica中调用函数并传递值
            result = session.evaluate(wlexpr(f"DG[{delta_m}, {fconstraintnum}, {lconstraintnum}]"))
            
            return {int_to_str(k, digit): v for k, v in zip(delta, [[x[i] for x in result] for i in range(len(result[0]))])}
    else:
        delta=int_to_str(delta, digit)
        # 启动一个新的Wolfram Language会话
        with WolframLanguageSession() as session:
            # 载入包含Mathematica函数定义的文件
            session.evaluate(wl.Get(root_dir+"DG(delta)_v1.m"))
            # 在Mathematica中调用函数并传递值
            result = session.evaluate(wlexpr(f"DG[{delta}, {fconstraintnum}, {lconstraintnum}]"))
            
            return delta, list(result)

def find_insert_position(missing_delta, exist_delta_list):
    for i in range(len(exist_delta_list)):
        if exist_delta_list[i] > float(missing_delta):
            return i
    return len(exist_delta_list)

def DG(data, deltas, digit, Dn_max):
    # Read csv
    #print(deltas)
    data = data.iloc[:Dn_max]
    
    fconstraintnum=1
    lconstraintnum=2*Dn_max-1
    derivative_order_list = [(2*i-1) for i in range(fconstraintnum, (lconstraintnum+4)//2)]
    # data.set_index(derivative_order_list, inplace=True)
    data.index = derivative_order_list
    delta_list = [int_to_str(d, digit) for d in data.columns[1:].to_list()]
    data.columns = [data.columns[0]]+delta_list
    
    # Delta preprocess
    deltas_str = [int_to_str(d, digit) for d in deltas]
    false_deltas = [int_to_str(d, digit) for d, is_in in zip(deltas, [d in delta_list for d in deltas_str]) if not is_in]
    
    if len(false_deltas) == 0:
        # print([f'{d}' for d in deltas_str])
        return np.array([data[f'{int_to_str(d, digit)}'].to_list() for d in deltas]).T.tolist()
    else:
        print(false_deltas)
        raise ValueError("data not found")
        
    
  

def DG_new_derivative_constraint(digit):
    # Read csv
    path = root_dir+"isingDGn_v5.csv"
    data = pd.read_csv(path)
    data.index = ['derivative order']+[2*i+1 for i in range(12)]
    
    
    # 定义导数阶数和 delta 的列表
    fconstraintnum=19
    lconstraintnum=25
    derivative_order_list = [fconstraintnum+2*i for i in range(int((lconstraintnum-fconstraintnum)/2)+1)]
    exist_deltas = list(data.columns)[0:5]
    
    # New csv
    new_path = root_dir + "LP_data/isingDGn_new_constraint.csv"
    new_df = pd.DataFrame()
    new_df.index = derivative_order_list
    batch_num = 2
    for i in range(len(exist_deltas)//batch_num+1):
        print('idx: ', i)
        print(exist_deltas[batch_num*i:batch_num*(i+1)])
        new_data = mathematica_dG(exist_deltas[batch_num*i:batch_num*(i+1)], digit, 
                                  fconstraintnum, lconstraintnum)
        new_columns_df = pd.DataFrame(new_data, index=derivative_order_list)
        insert_position = batch_num*i
        new_df = pd.concat([new_df.iloc[:, :insert_position], new_columns_df, new_df.iloc[:, insert_position:]], axis=1)
        
    new_df.to_csv(new_path, index=True)
    
# print(DG_new_derivative_constraint(3))