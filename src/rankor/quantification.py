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
from impetuous.quantification import group_significance , qvalues
from impetuous.convert import *
pi0 = lambda pvs : 1.
#
# THIS IS SOMEWHAT REDUNDANT AND A DUPLICATION OF
# SOME OF MY CODE THAT IS PART OF THE IMPETUOUS REPO
#
def pathway_frame_from_file ( filename , 
        delimiter = '\t' , item_sep = ',' ) :
    pdf = None
    with open( filename,'r' ) as input :
        for line in input :
            lspl = line.replace('\n','').split(delimiter)
            analytes_ = lspl[2:]
            desc = lspl[1]
            iden = lspl[0]
            ps = pd.Series( [desc , item_sep.join(analytes_) ] , 
                        name = iden , index = ['description','analytes'] )
            pdf = pd.concat([pdf,pd.DataFrame(ps).T])
    return ( pdf )

def create_dag_representation_df ( pathway_file = '../data/GROUPDEFINITIONS.gmt' , 
                                   pcfile = '../data/PCLIST.txt' 
                                 ) :
    pc_list_file = pcfile
    tree , ance , desc = parent_child_to_dag ( pc_list_file )
    pdf = make_pathway_ancestor_data_frame ( ance )
    pdf_ = pathway_frame_from_file( pathway_file )
    pdf.index = [v.replace(' ','') for v in  pdf.index.values]
    pdf_.index= [v.replace(' ','') for v in pdf_.index.values]
    dag_df = pd.concat([pdf.T,pdf_.T]).T
    return ( dag_df , tree )


def HierarchicalEnrichment (
            analyte_df , dag_df = None , dag_level_label = 'DAG,l' ,
            ancestors_id_label = 'aid' , id_name = None , threshold = 0.05 ,
            enrich_label = None , analyte_name_label = 'analytes' ,
            item_delimiter = ',' , alexa_elim = False , alternative = 'two-sided' ,
            remove_analytes = set(['']) , group_defs = None , pc_file = None ,
            fallback_sig_type = ',r', power = 2., use_global_params = False,
            centered = -1 , local_bh_correction = False 
        ) :
    from rankor.reducer import hyper_rdf, hyper_params
    #rankor.
    #
    # NEEDS AN ANALYTE FRAME :
    #     INCLUDING P VALUES OF ANALYTES (IF ENRICH LABEL IS SET
    #     OR EXPRESSION VALUES IF ENRICH LABEL IS NOT SET
    #
    # DAG GRAPH DESCRIPTION FRAME :
    #     INCLUDING NODE ID, NODE ANALYTES FIELD (SEPERATED BY ITEM DELIMITER)
    #     INCLUDING ANCESTORS FIELD (SEPERATED BY ITEM DELIMITER)
    #     DAG LEVEL OF EACH NODE OR A PARENT CHILD LOOKUP AND A GROUP DEFINITION
    #     FILE
    #
    global_mean , global_std = None , None
    sig_type = fallback_sig_type
    if enrich_label is None :
        #
        # IF TRUE DO NOT USE LOCAL PARAMETERS FOR THE P VALUES
        if use_global_params and not fallback_sig_type == ',r' :
            global_mean , global_std = hyper_params( analyte_df , power=power , label = 'global' )
        
    tolerance = threshold
    AllAnalytes = set( analyte_df.index.values ) ; nidx = len( AllAnalytes )
    SigAnalytes = None
    if not enrich_label is None :
        if enrich_label+sig_type in set( analyte_df.columns.values ) :
            SigAnalytes = set( analyte_df.iloc[ (analyte_df.loc[:,enrich_label+sig_type].values < tolerance), : ].index.values )
        else :
            print ( 'ERROR : COULD NOT FIND QUANTIFICATION LABEL' )
            exit ( 1 )
            
        if len( AllAnalytes ) == len( SigAnalytes ) :
            print ( 'WARNING : THIS STATISTICAL TEST WILL BE NONSENSE' )
            print ( 'WARNING : TRY A DIFFERENT THRESHOLD' )
    if len( AllAnalytes ) <= 1 :
        print ( ' WARNING : NOTHING TO TEST ' )
    #
    if dag_df is None :
        if group_defs is None or pc_file is None:
            print ( 'ERROR : MISSING DAG ' )
            exit(1)
        dag_df , tree = create_dag_representation_df (
            pathway_file    = group_defs ,
            pcfile          = pc_file 
        )
        dag_df = dag_df .rename( columns = { 'DAG,level' : 'DAG,l' ,
                                         'DAG,ancestors' : 'aid' } )
    df = dag_df ; dag_depth = np.max( df[dag_level_label].values )
    marked_analytes = {} ; used_analytes = {} ; node_sig = {}
    for d in range( dag_depth, 0, -1 ) : 
        # ROOT IS NOT INCLUDED
        filter_ = df [ dag_level_label ] == d
        nodes = df.iloc[ [i for i in np.where(filter_)[ 0 ]] ,:].index.values
        for node in nodes :
            if 'nan' in str(df.loc[node,analyte_name_label]).lower() :
                continue
            analytes_ = df.loc[node,analyte_name_label].replace('\n','').replace(' ','').split(item_delimiter)
            try :
                group = analyte_df.loc[[a for a in analytes_ if a in AllAnalytes] ].dropna( axis=0, how='any', thresh=analyte_df.shape[1]/2 ).drop_duplicates()
            except KeyError as e :
                continue
            if node in marked_analytes :
                unused_group = group.loc[ list( set(group.index.values) - marked_analytes[node] ) ]
                group = unused_group
            L_ = len( group ) ; str_analytes = ','.join(group.index.values)
            if L_ > 0 :
                used_analytes[node] = ',' .join ( group.index.values )
                SigAnalytes_ = SigAnalytes
                if SigAnalytes is None :
                     hdf_ = hyper_rdf ( analyte_df.loc[group.index.values,:] ,
                                        power = power , label = 'generic',
                                        MEAN = global_mean , STD = global_std ,
                                        centered = centered )
                loc_tol = tolerance
                if local_bh_correction :
                    loc_tol = loc_tol / float(L_)
                    if not enrich_label is None :
                        bSel = [ v < loc_tol for v in analyte_df.loc[group.index.values,enrich_label+sig_type].values ]
                        SigAnalytes_ = set ( group.index.values [ [i for i in range(L_) if bSel[i]] ] )
                    else :
                        SigAnalytes_ = set ( hdf_.index.values[ [ v < loc_tol for v in hdf_.loc[:,'generic'+sig_type].values ] ] )
                     
                pv , odds = group_significance ( group , AllAnalytes = AllAnalytes,
                                             SigAnalytes = SigAnalytes_ , tolerance = tolerance ,
                                             alternative = alternative )
                
                node_sig[node] = pv ; marked_ = set( group.index.values )
                ancestors = df.loc[node,ancestors_id_label].replace('\n','').replace(' ','').split(item_delimiter)
                #
                # IF TRUE THEN USE ALEXAS ELIM ALGORITHM : IT IS NOT DEFAULT
                # IN WHICH CASE ONLY THE SIGNIFICANT ANALYTES BELONGING TO
                # SIGNIFICANT GROUPS ARE REMOVED. ALEXAS ALGO FAVOURS LARGER GROUPS
                # WITH MORE ANALYTES
                #
                if alexa_elim and pv > threshold :
                    continue
                for u in ancestors :
                    if u in marked_analytes :
                        us = marked_analytes[u]
                        marked_analytes[u] = us | marked_
                    else :
                        marked_analytes[u] = marked_
    df['Hierarchical,p'] = [ node_sig[idx] if idx in node_sig else 1. for idx in df.index.values ]
    df['Hierarchical,q'] = [ qvs[0] for qvs in qvalues( df['Hierarchical,p'].values , pi0 = pi0(df.loc[:,'Hierarchical,p'].values) ) ]
    df['Included analytes,ids'] = [ used_analytes[idx] if idx in used_analytes else '' for idx in df.index.values ]
    df = df.dropna()
    return ( df )



if __name__ == '__main__' :
    print ( 'QUANT' )

