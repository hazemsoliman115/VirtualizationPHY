# -*- coding: utf-8 -*-
"""
Created on Tue Sep 27 12:25:39 2016

@author: hazem.soliman
"""

# -*- coding: utf-8 -*-
"""
Created on Mon Mar 28 14:50:42 2016

@author: hazem.soliman
"""

# -*- coding: utf-8 -*-
"""
Created on Tue Oct 13 15:18:16 2015

@author: hazem.soliman
"""

import random
import operator
import math
import csv
import numpy as np
import matplotlib.pyplot as plt
from scipy import optimize
from cvxopt import matrix, solvers
from scipy.optimize import minimize

class Resour(object):
    """ A class for the resource node used in the interval scheduling """
    def __init__(self, iden):
        """ constructor with attributes: identity(number), flows that start here and flows that end here """
        self.iden = iden
        self.starting = []
        self.ending = []
        
class RB_element(object):
    """ A class for resource nodes used in SDMA testing """
    def __init__(self, iden):
        """ constructor with attributes: identity """
        self.iden = iden
        self.slice_list = []

class W_Slice(object):
    """ The main node for scheduling entities, each object represents a flow """
    def __init__(self, iden, exp_users_k, QoS_expo, n_antenna, m_antenna):
        """ constructor """
        self.iden = iden
        self.exp_users_k = exp_users_k
        self.QoS_expo = QoS_expo
        self.n_antenna = n_antenna
        self.m_antenna = m_antenna
        self.ReqSize = None
        self.centerf = None
        self.startf = None
        self.endf = None
        self.psi = 0
        
        self.connVNodes = None
        self.NonconnVNodes = None
        self.Sch_Indicator = None
        self.Original_Sch_Indicator = None
        self.Glob_Thru = 0
       
        self.ReqNodes = None
        self.NonReqNodes = None
        self.ReqSize = None
        
        
    def GenResReq(self, NumRes, DepthTree, Allow_DepthTree, MaxDepthTree, TTI):
        """ (WFlow, int, int, int, int) -> NoneTyep Generate a reqsource request as follows: first, randomly select the requested resources size in a uniform fashion considering the depth of the tree to be the lower bound on the size. Second, choose the best set of resources of the chosen size within the allowed tree configuration considering the channel profile"""
        #Generate a request size by genearting a random integer between the maxdepthtree (from the number of global resources) and depthtree(as specified by the running envrionment). Then raise it as a power of 2
        self.req_size = 2**random.randint(DepthTree, MaxDepthTree)
        #Initialize an empty to store the weights of the various resource subsets to be examined
        self.req_weight=[]
        #Get the current channel profile
        col=[row[TTI] for row in self.channel_profile]
        #print(len(col))
        #print(type(col))
        #Get the weight for each subset
        for i in range(0,self.Res_size,max(self.req_size,2**Allow_DepthTree)):
            self.req_weight.append(sum(col[i:i+self.req_size]))
            
        #print(self.req_weight)
        #Find the subset of resources with the largest weight
        self.max_ResReq_index, self.max_ResReq = max(enumerate(self.req_weight), key = operator.itemgetter(1))
        self.Thru = self.max_ResReq
        self.original_Thru = self.max_ResReq
        #print(self.max_ResReq_index)
        #Encode the request, convert the resource subset index into binary, then fill it with zeros from the left according to its position in the tree
        self.ResReq = bin(self.max_ResReq_index)[2:].zfill(int(math.log2(len(self.req_weight))))
        
    def UpdateThru(self):
#        try:
            temp_NonconnVNodes = sorted(self.NonconnVNodes, key = lambda VarNode: len(VarNode.connVNodes), reverse=True)
            norm_Glob_Thru  = []
            for j in temp_NonconnVNodes:
                temp_sub_Thru  = [j.Sch_Indicator]
                for k in j.connVNodes:
                    if k in temp_NonconnVNodes:
                        temp_sub_Thru.append(k.Sch_Indicator)
                        temp_NonconnVNodes.remove(k)
                #norm_Thru = max(temp_sub_Thru)
                norm_Thru = 0
                for k in range(len(temp_sub_Thru)):
                    norm_Thru += np.exp(temp_sub_Thru[k])
                norm_Thru = np.log(norm_Thru)
                norm_Glob_Thru.append(norm_Thru)
            self.Glob_Thru = sum(norm_Glob_Thru)

        
    def Round_Centerf(self, x, NumRes):
        if round(x) >= self.ReqSize/2 and round(x) <= NumRes-self.ReqSize/2:
            self.centerf = round(x)
            self.startf = math.ceil(self.centerf - self.ReqSize/2)
            self.endf = math.ceil(self.centerf + self.ReqSize/2)-1
        elif x > NumRes/2:
            self.centerf = math.floor(x)
            self.startf = math.floor(self.centerf - self.ReqSize/2)
            self.endf = math.floor(self.centerf + self.ReqSize/2)
        else:
            self.centerf = math.ceil(x)
            self.startf = math.ceil(self.centerf - self.ReqSize/2)
            self.endf = math.ceil(self.centerf + self.ReqSize/2)
            
    def Check_QoS(self, total_k, Q_epsilon):
        """ A function to check the QoS is still valid """
        #print((1-total_k/(self.n_antenna-self.m_antenna))/Q_epsilon)
        if ((self.n_antenna-self.m_antenna)-total_k)/Q_epsilon > self.QoS_expo:
            #print("True")
            return True
        else:
            #print("False")
            return False
            
            
class WScheduler(object):
    """The main scheduler class, takes the form of a tree. Then users are sequentially fed to the tree and stored in the appropriate node depending upon their request. A method is used to examine the tree and return the list of scheduled users"""
    
    def __init__(self):
        """ Initialize the tree to an emoty tree """
        
        self.root = WSch_Node()
        
    def CreateConflictGraph(self, List_Flows):
        """ A function that takes the list of flows and finds the connected and no-connected nodes for each flow """
        for fl in List_Flows:
            fl.connVNodes = []
            for c_fl in List_Flows:
                if c_fl is not  fl:
                    if c_fl.ResReq in fl.ResReq[0:len(c_fl.ResReq)] or fl.ResReq in c_fl.ResReq[0:len(fl.ResReq)]:
                        fl.connVNodes.append(c_fl)
                    
        for fl in List_Flows:
            fl.NonconnVNodes = [x for x in List_Flows if x not in fl.connVNodes]
                    
                    
                
    def General_Find_Sch_Flows(self, List_Flows, gain, no_iter):
        """ The general scheduling algorithm """
        for fl in List_Flows:
            fl.UpdateThru()
        for j in range(no_iter):
            for fl in List_Flows:
                #max([x.Sch_Indicator+x.Glob_Thru for x in fl.connVNodes] if len(fl.connVNodes) > 0 else [0])
                norm_thru = 0
                if len(fl.connVNodes) > 0:
                    for x in fl.connVNodes:
                        norm_thru += np.exp(x.Sch_Indicator+x.Glob_Thru)
                norm_thru = np.log(norm_thru)
                fl.Sch_Indicator = fl.original_Thru*(gain + 0.5+0.5*np.tanh((fl.Sch_Indicator+fl.Glob_Thru - norm_thru )))

        
class AccessPoint(object):
    """ A general class for a wireless access point in the scheduling system """
    
    def __init__(self, i, pos_x, pos_y, a = 0):
        self.id = i
        self.pos_x = pos_x
        self.pos_y = pos_y
        self.active_status = a
        self.covered_cells = []
        self.Intf_cells = []
        self.attached_users = []
        self.Sec_attached_users = []
        self.total_Q = None
        self.clustered_cells = []
               
              
            
            
        
class Interference_Graph(object):
    """ A general graph to store the interference graph """
    
    def __init__(self, list_APs, intf_threshold):
        intf_g = [[0 for i in range(len(list_APs))] for i in range(len(list_APs))]
        intf_g_power = [[0 for i in range(len(list_APs))] for i in range(len(list_APs))]
        ap_index = 0
        for ap in list_APs:
            i_ap_index = 0
            for i_ap in list_APs:
                if np.sqrt((ap.pos_x-i_ap.pos_x)**2 + (ap.pos_y-i_ap.pos_y)**2) <= intf_threshold and ap is not i_ap:
                    intf_g[ap_index][i_ap_index] = 1
                    intf_g[i_ap_index][ap_index] = 1
                    intf_g_power[ap_index][i_ap_index] = np.sqrt((ap.pos_x-i_ap.pos_x)**2 + (ap.pos_y-i_ap.pos_y)**2)/intf_threshold
                    intf_g_power[i_ap_index][ap_index] = np.sqrt((ap.pos_x-i_ap.pos_x)**2 + (ap.pos_y-i_ap.pos_y)**2)/intf_threshold
                i_ap_index+=1
            ap_index+=1
            
        self.intf_g = intf_g
        self.intf_g_power = intf_g_power
        self.Max_f = 0
        
    def find_max_freq(self):
        """ A function to find the maximum frequency of occurrence in the graph, for the rounding """
        self.Max_f = 0        
        for l in self.intf_g:
            temp_max = sum(l)+1
            if self.Max_f <= temp_max:
                self.Max_f = temp_max
                
    def find_intf_g_matrix_opt(self):
        """ Prepares the matrix to be used in the optimization problrm """
        intf_g_matrix_opt = [[0 for i in range(len(self.intf_g))] for i in range(len(self.intf_g))]
        for k in range(len(intf_g_matrix_opt)):
            intf_g_matrix_opt[k][k] = -1
            for m in range(len(intf_g_matrix_opt)):
                if self.intf_g[k][m] == 1:
                    intf_g_matrix_opt[k][m] = -1
        return(intf_g_matrix_opt)

            
        
class FDOptimizer(object):
    """ A class to hold the optimization function for frequency spacing """
    def __init__(self):
        self.list_slices = None
        
    def FDOptimize(self, list_slices, NumRes):
        """ A function to decide the center of each slice band """
        
        N = len(list_slices)        

        si = [y.ReqSize for y in list_slices]

        fun2 = lambda x: -x[N]

        cons = []
        for i in range(N):
            #print(i)
            cons.append({'type': 'ineq', 'fun': lambda x, ind = i:  -x[ind] + (NumRes-si[ind]/2)})
            cons.append({'type': 'ineq', 'fun': lambda x, ind = i:  x[ind] - si[ind]/2})
            
        for i in range(N):
            for j in range(i+1,N,1):
                cons.append({'type': 'ineq', 'fun': lambda x, ind1 = i, ind2 = j:  (x[ind1] - x[ind2])**2 -x[N]})
        
        cons = tuple(cons)

        bnds = [(0, NumRes)]*(N) + [(-100,100)]
        
        x_init = [0+i/N*NumRes for i in range(N)]+[4]
        
        res = minimize(fun2, x_init, method='SLSQP', bounds=bnds, constraints=cons, options={'disp': False})
        #print(res)
        return(res)
        
    def Find_Sch_Flows_Interval(self, NumRes, List_Flows):
        """ The general function to schedule in interval graphs """
        List_Res = []
        # Create list of resources and associate to each the requesting flows
        for i in range(NumRes):
            r = Resour(i)
            for v in List_Flows:
                if v.startf == i:
                    r.starting.append(v)
                if v.endf == i:
                    r.ending.append(v)
            List_Res.append(r)
        # Forward Path
        temp_max = 0
        last_interval = None
        for r in List_Res:
            for v in r.starting:
                v.psi = temp_max + v.ReqSize
            for v in r.ending:
                if v.psi > temp_max:
                    temp_max = v.psi
                    last_interval = v
        
        #Backward Path            
        self.Max_Indep_Set = []
        self.Max_Indep_Set.append(last_interval)
        temp_max_reverse = temp_max - last_interval.ReqSize
        for r in reversed(List_Res):
            if r.iden >= last_interval.startf:
                continue
            for v in r.ending:
                if v.psi == temp_max_reverse:
                    last_interval = v
                    temp_max_reverse = temp_max_reverse - v.ReqSize
                    self.Max_Indep_Set.append(v)
        
    def Append_Flows_SDMA(self, list_slices, list_RBs, Q_epsilon):
        """ A function to append more flows using SDMA if the QoS criteria is satsified """
        for x in list_slices:
            if x not in self.Max_Indep_Set:
    #            List_neighbours = []
    #            for i in range(x.startf, x.endf+1,1):
    #                for y in list_RBs[i].slice_list:
    #                    if y not in List_neighbours:
    #                        List_neighbours.append(y)
                init_decision = True
                no_per_RB = []
                for i in range(x.startf, x.endf+1,1):
                    no_per_RB.append(len(list_RBs[i].slice_list))
                    for y in list_RBs[i].slice_list:
                        init_decision = init_decision and y.Check_QoS(sum([z.exp_users_k for z in list_RBs[i].slice_list])+x.exp_users_k - y.exp_users_k, Q_epsilon)
                if init_decision and (max(no_per_RB)+1)*per_slice_antenna <= total_antenna:
                    for i in range(x.startf, x.endf+1,1):
                        list_RBs[i].slice_list.append(x)
                    
        
if __name__ == "__main__":
    no_iter_statistics = 200
    Num_slices_list = [3,5,7,10,12]
    Q_epsilon_list = [i for i in np.arange(0.5,5.5,0.5)]
    total_slice_no = [[0 for i in range(len(Num_slices_list))] for j in range(len(Q_epsilon_list))]
    total_slice_no_SDMA = [[0 for i in range(len(Num_slices_list))] for j in range(len(Q_epsilon_list))]
    for i_num_slice in range(len(Num_slices_list)):
        print(i_num_slice)
        for i_Q_eps in range(len(Q_epsilon_list)):
            for k in range(no_iter_statistics):
                #Num_slices = 5
                Num_slices = Num_slices_list[i_num_slice]
                NumRes = 10
                exp_users = 1
                total_antenna  = 8
                per_slice_antenna = 4
                QoS_expo = 0.95
                sigma = 0.2
                #Q_epsilon = 0.95
                Q_epsilon = Q_epsilon_list[i_Q_eps]
                list_slices = []
                for i in range(Num_slices):
                    list_slices.append(W_Slice(i, np.random.normal(exp_users, sigma, 1), QoS_expo, total_antenna, per_slice_antenna))
                for y in list_slices:
                    #y.GenResReq(NumRes)
                    y.ReqSize = 3
                
                
                # 1. Optimization Step
                FDO = FDOptimizer()
                res = FDO.FDOptimize(list_slices, NumRes)
                for i in range(len(list_slices)):
                    list_slices[i].Round_Centerf(res.x[i], NumRes)
                    #print("Slice number =", list_slices[i].iden)
                    #print("Slice Size =", list_slices[i].ReqSize)
                    #print("Slice center =", list_slices[i].centerf)
                    #print("Slice start = ", list_slices[i].startf)
                    #print("Slice end = ", list_slices[i].endf)
                    
                    
                # 2. Selection through Graph MWIS
                FDO.Find_Sch_Flows_Interval(NumRes, list_slices)
                #print([x.iden for x in FDO.Max_Indep_Set])
                
                list_RBs = []    
                for i in range(NumRes):
                    list_RBs.append(RB_element(i))
                
                for y in list_RBs:
                    for x in FDO.Max_Indep_Set:
                        if x.startf <= y.iden <= x.endf:
                            y.slice_list.append(x)
                #print([[x.iden for x in y.slice_list] for y in list_RBs])
                #print([x.exp_users_k for x in list_slices])
                
                
                # 3. Extension
                FDO.Append_Flows_SDMA(list_slices, list_RBs, Q_epsilon)
                #print([[x.iden for x in y.slice_list] for y in list_RBs])
                
                list_slices_final = []
                temp = 0
                for y in list_RBs:
                    temp += len(y.slice_list)
                    
                    for x in y.slice_list:
                        if x not in list_slices_final:
                            list_slices_final.append(x)
                total_slice_no[i_Q_eps][i_num_slice] += len(list_slices_final)
                
                total_slice_no_SDMA[i_Q_eps][i_num_slice] += temp/len(list_RBs)
                
            total_slice_no[i_Q_eps][i_num_slice] = total_slice_no[i_Q_eps][i_num_slice]/no_iter_statistics
            total_slice_no_SDMA[i_Q_eps][i_num_slice] = total_slice_no_SDMA[i_Q_eps][i_num_slice]/no_iter_statistics
    marker_list='osv<>x'
    for i in range(len(Num_slices_list)):
        plt.plot(Q_epsilon_list,[total_slice_no[j][i] for j in range(len(Q_epsilon_list))], marker = marker_list[i], label = 'S = '+str(Num_slices_list[i]))
        plt.axis([0.5, 6, 2.5, 6.5])
        plt.legend(loc=1)
        plt.xlabel('$\epsilon$')
        plt.ylabel('Number of Selected Slices')
        plt.title('Number of Selected Slices versus QoS')
    plt.savefig('NumberSlices28032016.pdf', bbox_inches='tight')
    plt.savefig('NumberSlices28032016.eps', bbox_inches='tight')
    marker_list='osv<>x'
    for i in range(len(Num_slices_list)):
        plt.plot(Q_epsilon_list,[total_slice_no_SDMA[j][i] for j in range(len(Q_epsilon_list))], marker = marker_list[i], label = 'S = '+str(Num_slices_list[i]))
        plt.axis([0.5, 6, 0.6, 2])
        plt.legend(loc=1)
        plt.xlabel('$\epsilon$')
        plt.ylabel('Number of Slices Per Frequency Resource')
        plt.title('Number of Multiplexed Slices versus QoS')
    plt.savefig('NumberSlices28032016SDMA.pdf', bbox_inches='tight')
    plt.savefig('NumberSlices28032016SDMA.eps', bbox_inches='tight')
    