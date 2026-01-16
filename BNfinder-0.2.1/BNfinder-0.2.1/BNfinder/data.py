#!/usr/bin/env python2.5
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

import graph,stats,math
from itertools import chain
from BDE import BDE
from MDL import MDL

def param_simple(l):
    mu=stats.mean(l)
    sigma=stats.stdev(l)/2
    mu0=mu-sigma
    mu1=mu+sigma
    p0=p1=0.5
    return mu0, mu1, sigma, p0, p1

def param_2means(l):
    def sumvar(l,m):
        return stats.sum(map(lambda x:(x-m)**2,l))
    med=stats.lmedianscore(l)
    oldmed=med+1
    while med!=oldmed:
        l0=filter(lambda x: x<med,l)
        l1=filter(lambda x: x>med,l)
        lm=filter(lambda x: x==med,l)
        mu0=stats.mean(l0*2+lm)
        mu1=stats.mean(l1*2+lm)
        (oldmed,med)=(med,(mu0+mu1)/2)
    sigma=max(0.1**10,\
             math.sqrt((sumvar(l0+lm,mu0)+sumvar(l1,mu1)+(mu1-med)**2)/(len(l))))
    p0=(len(l0)+len(lm)*0.5)/len(l)
    p1=(len(l1)+len(lm)*0.5)/len(l)
    return mu0, mu1, sigma, p0, p1

def estprob(l,t=1.0,method=param_2means):    
    mu0,mu1,sigma,p0,p1 = method(l)
    trans=lambda x: 1.0/(1.0+(p0/p1)*math.exp((mu1**2-mu0**2+2*(mu0-mu1)*x)/(2*t*sigma**2)))
    return map(trans,l)

epp2m= lambda l : estprob(l,method=param_2means)
epps = lambda l : estprob(l,method=param_simple)

def ratio(sc):
    try:
        prob=str(2**sc)
    except OverflowError:
        sc10=math.log(2,10)*sc
        prob=str(10**(sc10-int(sc10)))+"e+"+str(int(sc10))
    return prob


class experiment:
    def __init__(self,name,data,perturbed,points=[],static=False):
        self._data=data
        self.perturbed=perturbed
        self.data=[]
        self.name=name
        self.points=points
        if static:
            for i,d in enumerate(self._data):
                self.data.append(self._data[i]+d)
        else:
            for i,d in enumerate(self._data[1:]):
                self.data.append(self._data[i]+d)

class gene:
    def __init__(self,name,n_disc,index):
        self.name=name
        self.n_disc=n_disc
        self.parents=[]
        self.values=[]
        self.index=index
        self.cpd=None
        
    def __str__(self):
        return self.name
    
    def disc_vals(self):
        if self.values:
            try:
                return map(int,self.values)
            except:
                return self.values
        else:
            return [0,1]

    def __index__(self):
        return self.index

    def __add__(self,other):
        return other+self.index
        
class dataset:
    def __init__(self,name='no_name'):
        self.name=name
        self.n_series=0
        self.n_gen=0
        self.points=[]
        self.vertices= []
        self.vertice_names=[]
        self.static=False
        self.cont_weight=1.5
        
    def fromFile(self,f,float_transform=epp2m):
        #number of genes in the first line
        self.n_gen=int(f.readline())
        # experiments description
        ln = f.readline().strip().split()
        self.n_series = int(ln[0])
        series = []
        for txt in ln[1:]:
            if txt.find("[")==-1:
                series.append((int(txt),[]))
            else: #perturbed genes
                series.append((int(txt[:txt.find("[")]),
                               [int(txt[1+txt.find("["):txt.find("]")])]))
        lines = []
        self.vertices=[]
        for i in range(self.n_gen):
            ln = f.readline().strip().split()
            n_disc = max(0,int(ln[1]))
            self.vertices.append(gene(ln[0],n_disc,i))
            if n_disc==0:
                if float_transform: #use float_transformn if specified
                    lines.append(float_transform(map(float,ln[2:])))
                else:
                    lines.append(map(float,ln[2:]))
            else:
                lines.append(map(int,ln[2:]))
        self.points=[]
        i=0
        for s,p in series:
            d = []
            for k in range(s): 
                d.append(map(lambda x:x[i],lines))
                i+=1
            name="serie%d"%len(self.points)
            self.points.append(experiment(name,d,p,range(len(d))))

        for v in self.vertices:
            v.parents=map(lambda v: (float(max(self.cont_weight,v.n_disc)),v),self.vertices)
        for e in self.points:
            e.perturbed=map(lambda i: self.vertices[i],e.perturbed)
        return self

    def toNewFile(self,file):
        #perturbations
        for s in self.points:
            for g in s.perturbed:
                file.write("#perturbed\t%s\t%s\n"%(s.name,self.vertices[g]))
        #condition names\
        file.write("genes\conditions")
        for s in self.points:
            for i,d in enumerate(s._data):
                file.write("\t%s:%d"%(s.name,i))
        file.write("\n")
        #expr_data
        for i,v in enumerate(self.vertices):
            file.write(str(v))
            for s in self.points:
                for j,d in enumerate(s._data):
                    if v.n_disc>0: #discrete
                        format="\t%d"
                    else: #continuous
                        format="\t%f"
                    file.write(format%d[i])
            file.write("\n")

    def fromNewFile(self,file,float_transform=epp2m):
        pert=[]
        regul=[]
        parents={}
        disc={}
        prior={}
        cpd={'and':[],'or':[]}
        contin=[]
        self.default=None
        ln=file.readline()
        #################
        #THE PREAMBLE
        ################
        #print 'preamble'
        while ln[0]=="#": #preamble
            rec=ln.strip().split()
            if rec[0]=="#perturbed": #experiment serie name followed by list of its perturbed vertices
                pert.append((rec[1],rec[2:]))
            elif rec[0]=="#regulators": #list of potential regulators of all vertices (except specified in #parents command) not specified in previous or present #regulators command
                reg=filter(lambda x:True not in map(lambda r:x in r,regul),rec[1:])
                for i,r in enumerate(reg):
                    if i!=reg.index(r):
                        reg[i]=None
                reg=filter(None,reg)
                regul.append(reg)
            elif rec[0]=="#parents": #vertice name followed by list of its potential regulators
                reg=rec[2:]
                for i in range(len(reg)):
                    if i!=reg.index(reg[i]):
                        reg[i]=None
                reg=filter(None,reg)
                parents[rec[1]]=reg
            elif rec[0]=="#default": #list of default possible values or 'FLOAT' for default continuous values
                if len(rec)==2 and rec[1]=="FLOAT":
                    self.default=[]
                else:
                    self.default=rec[1:]
            elif rec[0]=="#discrete": #vertice name followed by list of all its possible values
                disc[rec[1]]=rec[2:]
            elif rec[0]=="#continuous": #list of vertices with continuous values
                contin.extend(rec[1:])
            elif rec[0]=="#and": #list of vertices with cpd noisy-and
                cpd['and'].extend(rec[1:])
            elif rec[0]=="#or": #list of vertices with cpd noisy-or
                cpd['or'].extend(rec[1:])
            elif rec[0]=="#priorvert": #prior weight of all outgoing edges (except specified in #prioredge command) followed by list of related vertices
                val=float(rec[1])
                for name in rec[2:]:
                    prior[name]=val
            elif rec[0]=="#prioredge": #prior weight of edges preceded by edge target (common for all specified edges) and followed by list of edge sources
                target_name=rec[1]
                val=float(rec[2])
                for name in rec[3:]:
                    prior[(target_name,name)]=val
            else:
                pass #ignoring other stuff
            ln=file.readline()
        self.regulators=regul
        
        ################
        #condition names
        ################
        #print 'conditions'
        rec=ln.strip().split()
        #dynamic or static data?
        if len(rec[1].split(":"))==1:
            self.static=True
        #else:
        #    self.static=False
        if self.static:
            conds=map(lambda x : [x],rec[1:])
            names=conds
        else:
            serie=rec[1].split(":")[0]
            cond=[rec[1].split(":")[1]]
            conds=[cond]
            names=[serie]
            for c in rec[2:]:
                r=c.split(":")
                if r[0]==serie:
                    cond.append(r[1])
                else:
                    serie=r[0]
                    names.append(serie)
                    cond=[r[1]]
                    conds.append(cond)
        ln=file.readline()
        ###################
        # EXPRESSION DATA #
        ###################
        #print 'expression'
        self.vertices=[]
        self.vertice_names=[]
        #prepare the lists for expression data
        n_cond=sum(map(len,conds))
        exp_data=map(lambda l: [ [] for el in l],conds)
        geneindex=0
        while ln: #expression data
            rec=ln.strip().split()
            name=rec[0]
            vals=[]
            #are the values discrete or not?
            if name in contin:#did we say it's continuous?
                discrete=False
            elif name in disc.keys(): #did we say is's discrete?
                discrete=True
                vals=disc[name]
            elif self.default!=None:# did we specify a default
                if self.default==[]:
                    discrete=False
                else:
                    discrete=True
                    vals=self.default
            else: #We did not specify the type of values for this gene, try to guess it
                try:
                    i = int(rec[1])
                except ValueError: #it's not an integer
                    try:
                        i=float(rec[1])
                    except ValueError: #it's not a float -> it's an alphanumerical value
                        discrete=True
                    else:
                        discrete=False # it's a float
                else:
                    discrete=True #it's an int, let's map the strings to ints
            if discrete:
                line=rec[1:]
                if vals==[]:#we don't know what are the possible values
                    for v in line:
                        if v not in vals:
                            vals.append(v)
                    #let's sort the values
                    vals.sort()
            else: #not discrete
                line=map(float,rec[1:])
                if float_transform: #use float_transformn if specified
                    line=float_transform(line)
            
            self.vertices.append(gene(name,len(vals),geneindex))
            self.vertices[-1].values=vals
            self.vertice_names.append(name)
            #append this line to the expression data
            for el,l in zip(line,chain(*exp_data)):
                if discrete:
                    l.append(vals.index(el))
                else:
                    l.append(el)
            ln=file.readline()
            geneindex+=1
        ################
        #POSTPROCESSING#
        ################
        #print 'postprocessing'
        #process regulatory constraints
        for v in self.vertices:
            #select potential parents
            if parents.has_key(v.name):
                vparents=map(lambda n: self.vertices[self.vertice_names.index(n)],parents[v.name])
            else:
                #maybe there are regulators specified
                if regul==[]:
                    vparents=self.vertices
                else:#there is a list of regulators
                    vparents=[]
                    for r in regul:
                        if v.name in r:
                            break
                        else:
                            vparents.extend(map(lambda n: self.vertices[self.vertice_names.index(n)],r))
            #assign weights to parents
            for par in vparents:
                weight=float(max(self.cont_weight,par.n_disc))
                if (v.name,par.name) in prior.keys():
                    weight=weight**prior[(v.name,par.name)]
                elif par.name in prior.keys():
                    weight=weight**prior[par.name]
                v.parents.append((weight,par))
            v.parents.sort()
            #specify cpd
            if v.name in cpd['and']:
                v.cpd='and'
            elif v.name in cpd['or']:
                v.cpd='or'
        #process the lists of perturbed genes and make the list of experiments
        for i,(c,n) in enumerate(zip(conds,names)):
            pg = []
            for x,p in pert:
                if x==n:
                    pg.extend(map(lambda n: self.vertices[self.vertice_names.index(n)],p))# map gene names to indexes
            
            self.points.append(experiment(n,exp_data[i],pg,c,self.static))

        #set the number of genes
        self.n_gen=geneindex
        return self

    def sel_data(self,vert_num):
        data = []
        n_points = 0
        for e in self.points:
            if vert_num not in e.perturbed:
                data +=e.data
                n_points+=len(e.data)
        return n_points,data

    def orientation(self,par,child):
        nd,d = self.sel_data(child)
        pl = []
        cl = []
        n = len(d[0])/2
        for i in range(nd):
            try:
                pl.append(d[i][par])
            except TypeError: #python 2.4 compatibility wit __index__
                pl.append(d[i][par.__index__()])
            cl.append(d[i][child+n])
            
        sign = stats.pearsonr(pl,cl)[0]
        if sign >=0:
            return "+"
        else:
            return "-"

    
    def learn(self,score,data_factor,prior=None,distrib=False,\
            verbose=None,n_min=1,limit=None,min_empty=None,min_optim=None):
        scr=eval(score)(data_factor=data_factor,prior=prior)
        
        if distrib:
            from dispatch import CacheDispatch
            dispatch = CacheDispatch(scr)
            verts ={}
            from Queue import Empty
        score = 0.0
        pars ={}
        subpars ={}
        
        for x in self.vertices:
            n_points,data = self.sel_data(x)
            if distrib: #round robin dispatching
                dispatch.put_job(x,"learn_1",x,x.parents,data,n_points)
            else:
                if min_empty:
                    score_empty=scr.data_score(x,[],data)+scr.graph_score(len(x.parents),x,[],n_points)
                    score_max=score_empty-math.log(min_empty,2)
                else:
                    score_max=float("inf")
                if min_optim:
                    score_delta=-math.log(min_optim,2)
                else:
                    score_delta=float("inf")
                (par,sc),minlist = scr.learn_1(x,x.parents,data,n_points,\
                    verbose,n_min,limit,score_max,score_delta)
                score+=sc
                pars[x]=par
                subpars[x]=minlist
                
        if distrib: #clean up
            dispatch.clean()
            try:
                while True:
                    k,v = dispatch.get_result()
                    (par,sc),minlist = v
                    pars[k]=par
                    subpars[k]=minlist
                    score+=sc
            except Empty:
                pass
        
        #move parents from a dictionary to a list
        par_list = []
        for v in self.vertices:
            par_list.append(pars[v])
        pars=par_list

        g = graph.graph()
        g.fromParents(self.vertices ,pars)
        for v,par in zip(self.vertices,pars):
            g.vertice_labelling[v]=v.name
            for p in par:
                g.edge_labelling[p,v]=self.orientation(p,v)
        return score,g,subpars
   
    def get_stats(self,g):
        from score import score
        scr=score()
        stats={}
        for v in self.vertices:
            n_points,data = self.sel_data(v)
            sd = scr.sel_data(v,g.parents(v),data)
            stats_all,stats_par = scr.stats(sd,v,g.parents(v))
            stats[v]=stats_all, stats_par
        return stats
        
    def write_txt(self,subpars,file_name):
        f=open(file_name,"w")
        for v,minlist in subpars.items():
            f.write('\n%s'% v.name)
            for prob, pars in minlist:
                prob_s=ratio(prob)
                f.write('\n %s '% prob_s)
                for p in pars:
                    f.write(' %s'% p.name)
        f.close()
        
    def write_cpd(self,g,file_name,dirichlet=None):
        if dirichlet==None:
            dirichlet=1.0
        else:
            dirichlet=float(dirichlet)
        from score import score
        scr=score()
        f=open(file_name,"w")
        f.write('{\n')
        for v in self.vertices:
            n_disc=max(2,v.n_disc)
            n_points,data = self.sel_data(v)
            sd = scr.sel_data(v,g.parents(v),data)
            stats_all,stats_par = scr.stats(sd,v,g.parents(v))
            f.write('\''+str(v)+'\' : {\n')
            f.write('  \'vals\' : '+str(v.values)+' ,\n')
            f.write('  \'pars\' : '+str(map(str,g.parents(v)))+' ,\n')
            if v.cpd=='and':
                f.write('  \'type\' : \'and\' ,\n')
            elif v.cpd=='or':
                f.write('  \'type\' : \'or\' ,\n')
            f.write('  \'cpds\' : {\n')
            if v.cpd:
                p=scr.cpd_andor(len(g.parents(v)),stats_all,stats_par,v.cpd=='or')
                for (i,pv) in enumerate(p[:-1]):
                    f.write('     '+str(g.parents(v)[i])+' : '+str(pv)+'\n')
                f.write('     '+str(v)+' : '+str(p[-1])+' } } ,\n')
            else:
                stats_values={}
                for par_val in stats_par.keys():
                    stats_values[par_val]=[]
                for all_val in stats_all.keys():
                    stats_values[all_val[:-1]].append(all_val[-1])
                for par_val in stats_par.keys():
                    f.write('     '+str(par_val)+' : { ')
                    for val in stats_values[par_val]:
                        f.write(str(val)+' : '+str((stats_all[par_val+(val,)]+dirichlet)/(stats_par[par_val]+dirichlet*n_disc))+' , ')
                    f.write('None : '+str(dirichlet/(stats_par[par_val]+dirichlet*n_disc))+' } ,\n')
                f.write('     None : '+str(1.0/n_disc)+' } } ,\n')
        f.write('}\n')
        f.close()

    def write_bif(self,g,file_name,dirichlet=None,comments=[]):
        if dirichlet==None:
            dirichlet=1.0
        else:
            dirichlet=float(dirichlet)
        from score import score
        scr=score()
        f=open(file_name,"w")
        f.write('\\\\ File generated by BNfinder\n')
        import time
        f.write('\\\\ '+time.strftime('%x %X')+'\n')
        for com in comments:
            f.write('\\\\ '+com+'\n')
        f.write('\\\\ Conditional probability distributions generated with total pseudocounts number %f\n' % dirichlet)
        f.write('\n')
        f.write('network \"%s\" {}\n' % self.name)
        f.write('\n')
        for v in self.vertices:
            n_disc=max(2,v.n_disc)
            n_points,data = self.sel_data(v)
            sd = scr.sel_data(v,g.parents(v),data)
            stats_all,stats_par = scr.stats(sd,v,g.parents(v))
            f.write('variable \"%s\" {\n' % str(v))
            f.write('   type discrete[%d] { \"%s\" }\n'% (n_disc,'\" \"'.join(map(str,v.disc_vals()))))
            f.write('    }\n')
            if g.parents(v):
                f.write('probability ( \"%s\" | \"%s\" ) {\n'% (str(v),'\" \"'.join(map(str,g.parents(v)))))
                f.write('     default %s ;\n' % ' '.join([str(1.0/n_disc)]*n_disc))
            else:
                f.write('probability ( \"%s\" ) {\n'% str(v))
            stats_values={}
            for par_val in stats_par.keys():
                stats_values[par_val]=[]
            for all_val in stats_all.keys():
                stats_values[all_val[:-1]].append(all_val[-1])
            for par_val in stats_par.keys():
                if par_val:
                    f.write('     ( \"%s\" ) ' % '\" \"'.join(map(str,par_val)))
                else:
                    f.write('     table ')
                for val in v.disc_vals():
                    if val in stats_values[par_val]:
                        f.write(str((stats_all[par_val+(val,)]+dirichlet)/(stats_par[par_val]+dirichlet*n_disc))+' ')
                    else:
                        f.write(str(dirichlet/(stats_par[par_val]+dirichlet*n_disc))+' ')
                f.write(';\n')
            f.write('    }\n\n')
        f.close()

def rand_exp(n_genes,n_exp,n_disc):
    import random
    res = []
    for i in range(n_exp):
        exp = []
        for j in range(n_genes*2):
            exp.append(random.randint(0,n_disc-1))
        res.append(exp)
    return res
    
if __name__=="__main__":
    import sys
    new=open(sys.argv[2],"w")
    old=open(sys.argv[1])
    dataset().fromFile(old,None).toNewFile(new)
    new.close()
    old.close()
    
