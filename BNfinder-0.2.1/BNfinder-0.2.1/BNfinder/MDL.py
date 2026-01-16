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

import math
from score import score

class MDL(score):
    
    def graph_score(self,n_vert,v,weights,n_data):
        graph_struct=(1+len(weights))*math.log(n_vert+1,2)
        graph_cpd=0.5*(max(1.1,v.n_disc)-1)*math.log(n_data+1,2)
        if v.cpd in ['or','and','preproc']:
            graph_cpd*=1+sum(weights)/2.0
        else:
            for weight in weights:
                graph_cpd*=weight
        return graph_struct+graph_cpd
    
    
    def data_min(self,v,data):
        if v.n_disc:
            return 0.0
        sd = self.sel_data(v,[],data)
        s=0.0
        for d in sd:
            try:
                s -= d[0]*math.log(d[0],2) + (1-d[0])*math.log(1-d[0],2)
            except:
                pass
        return s*self.data_factor

    def data_score(self,v,par,data):
        sd = self.sel_data(v,par,data)
        stats_all,stats_par = self.stats(sd,v,par)
        s = 0.0
        if v.cpd in ['or','and']:
            prod_in=(v.cpd=='or')
            prod_out=1-prod_in
            p=self.cpd_andor(len(par),stats_all,stats_par,prod_in)
            for a in stats_par.keys():
                prob_a=p[-1]
                for i,av in enumerate(a):
                    if av==prod_in:
                        prob_a*=p[i]
                try:
                    s-=stats_all[a+(prod_out,)]*math.log(prob_a,2)
                except:
                    pass
                try:
                    s-=stats_all[a+(prod_in,)]*math.log(1-prob_a,2)
                except:
                    pass
        else:
            for a,cv in stats_all.items():
                par = a[:-1]
                cp = stats_par[par]
                try:
                    sc = math.log(1.0*cp/cv,2)
                    s+=cv*sc
                except ZeroDivisionError:
                    pass
        return s*self.data_factor


