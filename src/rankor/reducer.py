"""
Copyright 2020 RICHARD TJÖRNHAMMAR

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

    http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.
"""
import pandas as pd
import numpy as np
from scipy.stats import rankdata
from impetuous.quantification import qvalues, permuter
from rankor.quantification import pi0
from rankor.contrasts import contrast

def svd_reduced_mean ( x,axis=0,keep=[0] ) :
    if True :
        sk = set ( keep )
        if len ( np.shape(x) ) > 1 :
            u , s , vt = np .linalg .svd( x , full_matrices=False )
            xred = np.mean( np.dot(u*[s[i_] if i_ in sk else 0 for i_ in range(len(s))],vt) , axis)
            if 'pandas' in str(type(x)) :
                if not 'series' in str(type(x)) :
                    xname = x.index.values[0]
                    return ( pd.DataFrame( [xred] , index=[xname] , columns=x.columns ) )
                else :
                    xname = x.name
                    return ( pd.Series( xred , name=xname , index=x.columns ) )
            else :
                return ( xred )
    return ( x )

def hyper_params ( df_ , label = 'generic' , sep = ',', power=1. ):
    #
    idx_ = df_.index.values
    N_s = len ( df_.columns )
    u,s,vt = np.linalg.svd ( df_.values**power , full_matrices=False )
    #
    rdf_ = pd.Series ( np.sum(u**2,1) , index=idx_ , name = label+sep+"u" )
    rdf_ = pd.concat ( [ pd.DataFrame(rdf_) ,
                         pd.DataFrame( pd.Series( np.mean( df_.values,1 ) ,
                                                  index = idx_ , name=label+sep+"m") ) ] , axis=1 )
    w_ = rdf_ .loc[ :,label+sep+"u" ].values
    r_ = rankdata ( [ w for w in w_ ] ,'average' )
    N = len ( r_ )
    #
    df0_ = pd.DataFrame( [ [a for a in range(N)],w_,r_ ],index=['idx','w_','r_'], columns=idx_ ).T
    #
    from scipy.special import erf as erf_
    loc_pval = lambda X , mean , stdev : [ 1. - 0.5*( 1. + erf_( ( x - mean )/stdev/np.sqrt(2.) ) ) for x in X ]
    lvs      = np.log( df0_.loc[ :,'w_'].values )
    #
    return ( np.mean(lvs) , np.std(lvs) )


def hyper_rdf ( df_ , label = 'generic' , sep = ',' , power=1. ,
                diagnostic_output = False , MEAN=None , STD=None ) :
    #
    idx_= df_.index.values
    N_s = len ( df_.columns )
    u,s,vt = np.linalg.svd ( df_.values**power , full_matrices=False )
    #
    rdf_ = pd.Series ( np.sum(u**2,1) , index=idx_ , name = label+sep+"u" )
    rdf_ = pd.concat ( [ pd.DataFrame(rdf_) ,
                         pd.DataFrame( pd.Series( np.mean( df_.values,1 ) ,
                                                  index = idx_ , name=label+sep+"m") ) ] , axis=1 )
    w_ = rdf_ .loc[ :,label+sep+"u" ].values
    r_ = rankdata ( [ w for w in w_] ,'average' )
    N  = len ( r_ )
    #
    # HERE WE CONSTRUCT A DIANGOSTICS AND INTERMITTENT CALCULATION
    # DATAFRAME FOR EASIER DEBUGGING AND INSPECTION
    df0_ = pd.DataFrame( [ [a for a in range(N)],w_,r_ ],index=['idx','w_','r_'], columns=idx_ )
    df0_ = df0_.T.sort_values ( by='r_' )
    df0_ .loc[ : , 'd_'  ] = [ v for v in ( df0_.loc[:, 'w_' ] * df0_.loc[:,'r_'] ) ]
    df0_ .loc[ : , 'da_' ] = np.cumsum ( df0_ .loc[ : , 'd_' ].values  )
    #
    # HOW MANY EQUALLY OR MORE EXTREME VALUES ARE THERE? ( RANK BASED )
    df0_ .loc[ : , 'dt_' ] = np.cumsum ( df0_ .loc[ : , 'd_' ].values[::-1] )[::-1]
    df0_ .loc[ : , 'rank_ps' ] = df0_ .loc[ :,'dt_' ] / np.sum( df0_ .loc[ :,'d_' ] )
    #
    # DRAW HISTOGRAM TO SEE DISTRIBUTION OF THE DISTANCES
    from scipy.special import erf as erf_
    loc_pval = lambda X , mean , stdev : [ 1. - 0.5*( 1. + erf_( ( x - mean )/stdev/np.sqrt(2.) ) ) for x in X ]
    lvs      = np.log( df0_.loc[ :,'w_'].values )
    if MEAN is None or STD is None:
        if len(lvs)>3 :
            ps = loc_pval( lvs , np.mean( lvs ) , np.std( lvs ) )
        else :
            ps = df0_ .loc[ : , 'rank_ps' ]
    else :
        ps = loc_pval( lvs , MEAN , STD )
    #
    if diagnostic_output :
        import scipy.stats       as scs
        NB    = 100
        lv    = np.log(  rdf_.loc[:,[label+sep+"u"]].values  )
        y , x = np.histogram( lv , bins=NB , density=True )
        skw = scs.skew(lv)[0]
        kur = scs.kurtosis(lv)[0]
        shape_stats = "kurtosis: " + "{:.2f} ".format( kur ) + "skewness: "+ "{:.2f}".format( skw )
        locd  = lambda x,M,s : (1./s/np.sqrt(2.*np.pi))*np.exp(-0.5*((x-M)/s)**2 )
        lin_x = 0.5 * ( x[1:] + x[:-1] )
        his_y = y
        rsy = sorted( [ (y_,x_) for (y_,x_) in zip(y,lin_x) ] )[::-1]
        hm , hs = np.mean([rsy[i][1] for i in range(5)]) , np.mean([rsy[i][0] for i in range(5)])
        h_mod_y = locd( lin_x , hm , 1.0/(hs*np.sqrt(2.*np.pi)) )
        d_mod_y = locd( lv , np.mean(lv) , np.std(lv) )
        rem_y = [ (y_-m_) for (y_,m_) in zip(y, locd(0.5*(x[1:]+x[:-1]),np.mean(lv),np.std(lv))) ]
        prc_y = [ 100.*np.abs( contrast(y_,m_) ) for (y_,m_) in zip(y, locd(0.5*(x[1:]+x[:-1]),np.mean(lv),np.std(lv))) ]
        RMSD = np.sqrt(np.sum([ ry**2 for ry in rem_y ]))
        PMAE = np.mean(prc_y)
    #
    df0_ .loc[ : , 'pvalues' ] = ps
    #
    # ASSIGN THE P VALUES
    rdf_.loc[df0_.index.values,label+sep+"p"] = df0_.loc[:,'pvalues']
    rdf_.loc[df0_.index.values,label+sep+"q"] = [ qvs[0] for qvs in qvalues( df0_.loc[:,'pvalues'].values , pi0 = pi0(df0_.loc[:,'pvalues'].values) ) ]
    rdf_.loc[df0_.index.values,label+sep+"r"] = df0_ .loc[ : , 'rank_ps' ]
    #
    # AND RETURN THE RESULTS
    if diagnostic_output :
        return ( rdf_ , RMSD , PMAE , kur , skw , rem_y )
    else :
        return ( rdf_ )


if __name__ == '__main__' :
    #
    print ( 'REDUCER :: TESTS ' )
    a   = 2*np.random.rand(10)
    b   = 4*np.random.rand(10)
    X   = [ [*(a[:5]+1),*a[5:]],[*(b[:5]+3),*(b[5:])] ]
    Xdf = pd.DataFrame( X , columns=['a','b','s','c','d','e','f','g','h','i'] , index=['i0','i1'])
    #
    print ( Xdf )
    print ( svd_reduced_mean ( X ) )
    print ( svd_reduced_mean ( Xdf ) )
