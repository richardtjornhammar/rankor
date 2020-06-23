import pandas as pd
import numpy as np

def aggregate_duplicate_indices( mdat , aggregate ,
                                merge_and_split=False ,
                                keep_old_cols=True ) :
    rdat = None
    if merge_and_split :
        mdat.index = [ "#".join(idx) for idx in mdat.index.values ]
    all_unique_indices = sorted( list(set(mdat.index.values)) )
    Ntot = len ( all_unique_indices )
    for idx_ in range( Ntot ) :
        print ( np.floor(idx_/Ntot*100),'%',end='\r' )
        idx    = all_unique_indices[ idx_ ]
        loc_df = mdat.loc[ idx,: ]
        loc_df = loc_df.groupby(loc_df.index).apply(aggregate)
        if len ( loc_df ) > 0 :
            ldf = loc_df
            if 'serie' in str(type(loc_df)).lower() :
                ldf = pd.DataFrame( loc_df.values , columns=[loc_df.name] ,
                                   index=[idx[0] if 'tuple' in str(type(idx)) else idx for idx in loc_df.index.values] ).T
            if rdat is None :
                rdat = ldf  
            else :
                rdat = pd.concat( [rdat,ldf] )
    if merge_and_split :
        rdat .index = [ tuple(a.split('#')) if '#' in a else a for a in rdat.index ]
    return ( rdat )


def low_missing_value_imputation ( fdf , fraction = 0.9 , absolute = 'True' ) :
    import numpy as np
    #
    # fdf is a dataframe with NaN values
    # fraction is the fraction of information that should be kept
    # absolute is used if the data is positive
    #
    V = fdf.apply(pd.to_numeric).fillna(0).values
    u,s,vt = np.linalg.svd(V,full_matrices=False)
    s =  np.array( [ s[i_] if i_<np.floor(len(s)*fraction) else 0 for i_ in range(len(s)) ] )
    nan_values = np.dot(np.dot(u,np.diag(s)),vt)
    if absolute :
        nan_values = np.abs(nan_values)
    #
    # THIS CAN BE DONE BETTER
    for j in range(len(fdf.columns.values)):
        for i in range(len(fdf.index.values)):
            if 'nan' in str(fdf.iloc[i,j]).lower():
                fdf.iloc[i,j] = nan_values[i,j]
    return ( fdf )

def convert ( filename,to='csv' ) :
    sep_d = { 'tsv':'\t' , 'csv':',' }
    aggregate_duplicate_indices ( pd.read_csv( filename,sep_d[filename.split('.')[-1]],index_col=0 ) , np.mean )\
        .to_csv( filename.replace( '.'+filename.split('.')[-1],'.'+to )\
            .replace(filename.split('.')[0],filename.split('.')[0]+'_mod') , sep_d[to] )

def differential ( expression,level=2 ,use_contrast=False ):
    edf = expression.copy()
    if 'tuple' in str( type(edf.columns.values[0]) ) :
        pt   = [ v[level] for v in edf.columns.values ]
    else :
        pt   = [ v.split('_')[level] for v in edf.columns.values ]
    keep_complete = edf.T.groupby(pt).apply( lambda x:len(x.values)==2 )
    ptk = [ p for (p,k) in zip(keep_complete.index.values,keep_complete) if k ]
    edf = edf.iloc[:, [c[level] in set(ptk) for c in edf.columns.values ] ]
    
    differential = lambda x : (x.values[1]-x.values[0]) # FORWARD DIFF (OGTT - Basal) 
    if use_contrast :
            differential = lambda x : (x.values[1]-x.values[0])/(x.values[0]+x.values[1])
    ddf = edf.T.groupby(level=level).apply( \
            lambda x: pd.Series( differential(x) ,
                          index = x.columns ) ).T
    return ( ddf )

def output_separate_tsv(dedf_,which,suffix,status_d,sex_d):

    dedf_  = dedf_.loc[:,[c for c in dedf_.columns if which in sex_d[c]] ]

    pdf = pd.DataFrame( [ status_d[c] for c in dedf_.columns ] , columns=['Status'] )
    pdf .index = [ c for c in dedf_.columns.values ]
    pdf .to_csv( which.lower()+'_desc.tsv','\t' )

    dedf_ .index.name = 'Symbol'
    dedf_ .columns = pdf.index.values 
    dedf_ .to_csv( which.lower()+suffix,'\t' )

    return ( dedf_ )

def make_cpm( df , which = [0,'Female'] ,
             ledger={1:'status',-1:'time',2:'replicate'} ,
             sep='_' ,suffix ='wp1_data.tsv' ) :
    #
    #   FOR CREATING R EQUIVALENT STRUCTURES
    #
    df = df.loc[ :,[c for c in df.columns if which[1] in c.split(sep)[which[0]] ] ]
    sdf = None
    for key , value in ledger.items() :
        tdf = pd.DataFrame( [ c.split(sep)[key] for c in df.columns ] , columns=[value] )
        if sdf is None :
            sdf = tdf.T
        else :
            sdf = pd.concat([tdf.T,sdf])
    sdf.columns = df.columns
    
    sdf = sdf.T
    CPM = df.copy()
    CNT = df.apply( lambda x:[ int(v) for v in np.round(2**x) ] )
    sdf.to_csv( which[1].lower()+'_sam_'+suffix,'\t' )
    CPM.to_csv( which[1].lower()+'_cpm_'+suffix,'\t' )
    CNT.to_csv( which[1].lower()+'_cnt_'+suffix,'\t' )
    
    
if __name__ == '__main__' :

    print ( 'HERE WE MAKE THINGS UP!' )
