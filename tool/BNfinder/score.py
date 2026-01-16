# Copyright (c) 2004-2007 Bartek Wilczynski and Norbert Dojer; All Rights Reserved.
#
# This software is distributable under the terms of the GNU
# General Public License (GPL) v2, the text of which can be found at
# http://www.gnu.org/copyleft/gpl-2.0.txt. Installing, importing or otherwise
# using this module constitutes acceptance of the terms of this License.
#
# Disclaimer
# 
# This software is provided "as-is".  There are no expressed or implied
# warranties of any kind, including, but not limited to, the warranties
# of merchantability and fittness for a given application.  In no event
# shall the authors be liable for any direct, indirect, incidental,
# special, exemplary or consequential damages (including, but not limited
# to, loss of use, data or profits, or business interruption) however
# caused and on any theory of liability, whether in contract, strict
# liability or tort (including negligence or otherwise) arising in any way
# out of the use of this software, even if advised of the possibility of
# such damage.
#
import graph
from heapq import heappush, heappop, heapreplace

class score:
    """
    Abstract class implementing a scoring function.

    To obtain a working implementation, one has to define methods graph_score and data_score
    """
    
    def __init__(self,data_factor = 1.0,prior=None):
        self.data_factor=data_factor
        self.sloops = False # do we allow self-loops in the networks ?
        self.epsilon = 10**(-7)
        
    def graph_score(self,n_vert,v,weights,n_data):
        """
        Method for computing graph score factor - ABSTRACT method

        it computes the score given :
        n_vert - number of potential parents of a given gene
        v - given gene
        par_disc - list of weights of parents
        n_data - number of data points
        """
        pass
    def data_min(self,v,data):
        """
        method for computing the lower bound of the data score factor

        parameters:
        v - number of the vertice
        data - list of datapoints :

        each datapoint is a vector containing 2*n_vert values (expression levels for all genes for time points t and t+1
        """
        return 0.0
    
    def data_score(self,v,par,data):
        """
        abstract method for computing the data score factor

        parameters:
        v - number of the vertice
        par - list of numbers of parents
        data - list of datapoints :

        each datapoint is a vector containing 2*n_vert values (expression levels for all genes for time points t and t+1
        """
        
        pass

    def sel_data(self,v,par,data):
        """
        given a vertice, list of his parents and data,
        select all lists of parents values at time t
        and child values at point t+1
        """
        n = len(data[0])/2
        r = []
        for i in data:
            r.append([])
        for i in par:
            for j,d in enumerate(data):
                try:
                    r[j].append(d[i])
                except TypeError: #in case of older (2.4) python not supporting non integer indexes
                    r[j].append(d[i.__index__()])
        for j,d in enumerate(data):
            r[j].append(d[v+n])
        r.sort()
        return r

    def subsets(self,list,k):
        """
        generate all the k-long subsets of a list 
        """

        s = graph.stack()
        s.put((list,k,[]))
        try:
            while True: # exception used to get out of the loop
                list,k,acc = s.get()
                if k==0:
                    yield acc
                elif list==[]:
                    pass
                else:
                    s.put((list[1:],k,acc))
                    s.put((list[1:],k-1,acc+[list[0]]))
        except IndexError:
            pass # we've exhausted all options

    def stats(self,sd,v,par):
        def strec(key,prob,d):
            first,first_disc=d[0]
            rest=d[1:]
            if rest==[]:
                try:
                    stats_par[key]+=prob
                except KeyError:
                    stats_par[key]=prob
                if first_disc:
                    try:
                        stats_all[key+(first,)]+=prob
                    except KeyError:
                        stats_all[key+(first,)]=prob
                else:
                    try:
                        stats_all[key+(0,)]+=prob*(1-first)
                    except KeyError:
                        stats_all[key+(0,)]=prob*(1-first)
                    try:
                        stats_all[key+(1,)]+=prob*first
                    except KeyError:
                        stats_all[key+(1,)]=prob*first
            else:
                if first_disc:
                    strec(key+(first,),prob,rest)
                else:
                    strec(key+(0,),prob*(1-first),rest)
                    strec(key+(1,),prob*first,rest)
        d_disc=map(lambda p: p.n_disc,par+[v])
        if 0 not in d_disc:
            return self.stats_disc(sd)
        stats_all = {}
        stats_par = {}
        for d in sd:
            strec((),1,zip(d,d_disc))
        return stats_all,stats_par
        
    def stats_disc(self,sd):
        stats_all = {}
        stats_par = {}
        for d in sd:
            pars = tuple(d[:-1])
            all = tuple(d)
            try:
                stats_par[pars]+=1
            except KeyError:
                stats_par[pars]=1
            try:
                stats_all[all]+=1
            except KeyError:
                stats_all[all]=1
        return stats_all,stats_par

    def cpd_andor(self,n_par,stats_all,stats_par,prod_in):
        n_all=n_par+1
        eps=self.epsilon
        epsilon=eps*n_all*10
        v_all=0.0
        v_in=0.0
        for all,val in stats_all.items():
            v_all+=val
            if all[-1]==prod_in:
                v_in+=val
        if v_all-v_in<=eps:
            return [0.0]*n_all
        elif v_in<=eps:
            return [1.0]*n_all
        p=[1.0-eps]*n_all
        delta=epsilon+1.0
        while delta>epsilon:
            delta=0.0
            for i in range(n_all):
                p_old=p[i]
                delt=eps+1.0
                while abs(delt)>eps:
                    sumf=0.0
                    sumfprim=0.0
                    for par in stats_par.keys():
                        all=par+(prod_in,)
                        if all[i]==prod_in:
                            if all in stats_all.keys():
                                prod=1.0
                                for j in range(n_all):
                                    if all[j]==prod_in:
                                        prod*=p[j]
                                q=stats_all[all]/(1.0-prod)
                                sumf+=q-stats_par[par]
                                sumfprim+=(q/(1.0-prod))*(prod/p[i])
                            else:
                                sumf-=stats_par[par]
                    delt=sumf/sumfprim
                    p[i]-=delt
                    if p[i]>=1.0:
                        delt=1.0-eps-p[i]-delt
                        p[i]=1.0-eps
                delta+=abs(p[i]-p_old)
        return p
		    
    def learn_1(self,v,vert,data,nd,verbose=None,n_min=1,limit=None,score_max=float("inf"),score_delta=float("inf")):
        
        if verbose:
            print 'Learning parents of', v.name, '...',
        if not self.sloops:
            for (weight,parent) in vert:
                if parent==v:
                    import copy
                    vv = copy.copy(vert)
                    vv.remove((weight,parent))
                    vert = vv
        n = len(vert)
        try:
            lim=int(limit)
        except:
            lim=n
        mindata = self.data_min(v,data)
        
        min_set = minset(n_min,score_max,score_delta,self.data_score(v,[],data)+self.graph_score(n,v,[],nd)) #empty parents set
        if vert: # are there any potential parents?
            if vert[0][0]==vert[-1][0]: # we can use algorithm 2
                weight=vert[0][0]
                parents=map(lambda pair: pair[1],vert)
                size = 1
                mg = self.graph_score(n,v,[weight],nd)
                while min_set.accepts(mg+mindata) and (size<=lim): #we can possibly add (sub-)optimal scores
                    for sub in self.subsets(parents,size):
                        min_set.add(mg+self.data_score(v,sub,data), sub)
                    size+=1
                    mg = self.graph_score(n,v,[weight]*size,nd)
            else: # we have to use algorithm 1
                subsets=[] #successors of considered yet potential parents sets
                for (weight,parent) in vert: #one-element parents sets
                    heappush(subsets,(self.graph_score(n,v,[weight],nd),[weight],[parent]))
                while subsets:
                    mg,weights,sub=heappop(subsets)
                    if not min_set.accepts(mg+mindata): #we cannot improve the score
                        break
                    min_set.add(mg+self.data_score(v,sub,data), sub)
                    #insert sub's successors
                    if len(sub)<lim:
                        last_parent=vert.index((weights[-1],sub[-1]))
                        for (weight,parent) in vert[last_parent+1:]:
                            sub_succ=sub+[parent]
                            weights_succ=weights+[weight]
                            mg_succ=self.graph_score(n,v,weights_succ,nd)
                            heappush(subsets,(mg_succ,weights_succ,sub_succ))
        if verbose:
            print 'done', min_set
        return min_set.optimal, min_set.tolist() 

    def learn_all(self,vertices,data,n_points):
        par = {}
        for v in vertices:
            par[v] = self.learn_1(v,vertices,data,n_points)
        return par
    
    def score_graph(self,g,data):
        s = 0.0
        n_vert = len(g.vertices)
        for i,v in enumerate(g.vertices):
            p = g.parents(v)
            n,d = data.sel_data(i)
            sg = self.graph_score(n_vert,v,map(lambda par:par.n_disc,p),n)
            sd = self.data_score(i,p,d)
            s+=sg+sd
        return s
        
class minset:
    def __init__(self,size,score_max,score_delta,emptyscore):
        self.optimal=([],emptyscore)
        self.escore=emptyscore
        self.free=size
        self.mscore=score_max
        self.dscore=score_delta
        self.mset=[]
        self.add(emptyscore,[])
        
    def add(self,score,parents):
        if score<min(self.optimal[1]+self.dscore,self.mscore):
            if self.free:
                heappush(self.mset,(-score,parents))
                self.free-=1
                if not self.free:
                    self.mscore=-self.mset[0][0]
            elif self.mset:
                heapreplace(self.mset,(-score,parents))
                self.mscore=-self.mset[0][0]
        if score<self.optimal[1]:
            self.optimal=(parents,score)
            while self.mset and -self.mset[0][0]>self.optimal[1]+self.dscore:
                tmp=heappop(self.mset)
                self.free+=1
                
    def __str__(self):
        self.mset.sort()
        minstr=''
        for (score,parents) in reversed(self.mset):
            minstr+='\n Score '+str(-score)+' for parents set:'
            for p in parents:
                minstr+=' '+p.name
        return minstr
        
    def accepts(self,score):
        return score<min(self.optimal[1]+self.dscore,self.mscore)
        
#    def optimal(self):
#        negscore, parset= max(self.mset)
#        return parset, -negscore
        
    def tolist(self):
        self.mset.sort()
        minlist=[]
        for (score,parents) in self.mset:
            minlist.append((score+self.escore,parents))
        minlist.reverse()
        return minlist
        
