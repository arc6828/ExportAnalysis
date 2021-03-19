import warnings
import dynamo as dyn
import pandas as pd
import numpy as np
import anndata


def config() :
    warnings.filterwarnings('ignore')
    dyn.get_all_dependencies_version()

def removeOutlierObservation(a,b, obs):
    intersect_index = list( set(list(a.index)) & set(list(b.index)) & set(list(obs.index)) )
    intersect_index.sort()
    intersect_index

    a = a.loc[intersect_index]
    b = b.loc[intersect_index]
    obs = obs.loc[intersect_index]
    return [a,b,obs]

def removeOutlierObservationExtra(a, obs):
    intersect_index = list( set(list(a.index)) & set(list(obs.index)) )
    intersect_index.sort()
    intersect_index

    a = a.loc[intersect_index]
    obs = obs.loc[intersect_index]
    return [a,obs]

def createAnnotationData(a,b, obs, var, filter=True) :
    # intersect_index = list( set(list(a.index)) & set(list(b.index)) & set(list(obs.index)) )
    # intersect_index.sort()
    # intersect_index

    # a = a.loc[intersect_index]
    # b = b.loc[intersect_index]
    # obs = obs.loc[intersect_index]
    if filter :
        [a,b,obs] = removeOutlierObservation(a,b, obs)
    
    layers = {
        "spliced" : a.to_numpy(),
        "unspliced" : b.to_numpy(),
    }
    adata = anndata.AnnData(X=a.to_numpy(), obs=obs,var=var, layers=layers)
    return adata

    # print("555")
    # print("888")

def dynamoProcess(adata) :
    # print("555")
    # print("hello")
    dyn.pp.recipe_monocle(adata)
    dyn.tl.dynamics(adata, model='stochastic', cores=8)
    dyn.tl.reduceDimension(adata)
    dyn.tl.cell_velocities(adata)
    return adata

def dynamoProcess2(adata) :
    # print("555")
    # print("hello")
    dyn.pp.recipe_monocle(adata)
    dyn.tl.dynamics(adata, model='stochastic', cores=8)
    dyn.tl.reduceDimension(adata)
    dyn.tl.cell_velocities(adata)
    return adata

def dynamoPlot(adata, group, filter_key="", filter_values=[]) :
    temp_adata = adata.copy()
    if len(filter_values) > 0 :
        print()
        # List Region
        # temp_adata.obs['parentName'].drop_duplicates().to_numpy()
        temp_adata.obs['custom_'+group] = temp_adata.obs[group] 
        temp_adata.obs['custom_'+group] = temp_adata.obs.apply(lambda x: x[group]*(x[filter_key] in filter_values), axis=1)
        group = 'custom_'+group



    dyn.pl.streamline_plot(temp_adata, color=[group], basis='umap', show_legend='on data', show_arrowed_spines=True,cut_off_velocity=False)
    

def dynamoLargePlot(adata, group, filter_key="", filter_values=[]) :
    temp_adata = adata.copy()
    if len(filter_values) > 0 :
        print()
        # List Region
        # temp_adata.obs['parentName'].drop_duplicates().to_numpy()
        temp_adata.obs['custom_'+group] = temp_adata.obs[group] 
        temp_adata.obs['custom_'+group] = temp_adata.obs.apply(lambda x: x[group]*(x[filter_key] in filter_values), axis=1)
        group = 'custom_'+group
    dyn.pl.streamline_plot(temp_adata, color=[group], basis='umap', show_legend='on data', show_arrowed_spines=True,cut_off_velocity=False,figsize=[12,8])

def dynamoVeryLargePlot(adata, group, filter_key="", filter_values=[]) :
    temp_adata = adata.copy()
    if len(filter_values) > 0 :
        print()
        # List Region
        # temp_adata.obs['parentName'].drop_duplicates().to_numpy()
        temp_adata.obs['custom_'+group] = temp_adata.obs[group] 
        temp_adata.obs['custom_'+group] = temp_adata.obs.apply(lambda x: x[group]*(x[filter_key] in filter_values), axis=1)
        group = 'custom_'+group
    dyn.pl.streamline_plot(temp_adata, color=[group], basis='umap', show_legend='on data', show_arrowed_spines=True,cut_off_velocity=False,figsize=[24,16])
 

def concatDataFrame(df1, df2):
    pieces = {'latest': df1, 'oldest': df2}
    df_piece = pd.concat(pieces)
    # df_piece = df_piece.reset_index()
    # df_piece["key"] = df_piece["level_0"] + "_" + df_piece["exporter"]
    # df_piece = df_piece.drop(['level_0'], axis=1)
    # df_piece = df_piece.drop(['exporter'], axis=1)
    # df_piece = df_piece.set_index(['key'])
    return df_piece

def concatDataFrameExtra(dfs):
    # pieces = {'latest': df1, 'oldest': df2}
    df_piece = pd.concat(dfs)
    # df_piece = df_piece.reset_index()
    # df_piece["key"] = df_piece["level_0"] + "_" + df_piece["exporter"]
    # df_piece = df_piece.drop(['level_0'], axis=1)
    # df_piece = df_piece.drop(['exporter'], axis=1)
    # df_piece = df_piece.set_index(['key'])
    return df_piece

def createVectorField(adata) :    
    print("Prepare to Create Vector Field")

    a = pd.DataFrame.sparse.from_spmatrix(adata.layers['spliced'])
    b = pd.DataFrame.sparse.from_spmatrix(adata.layers['unspliced'])

    df_main = concatDataFrame(a, b)
    
    df_obs = concatDataFrame(adata.obs, adata.obs)

    layers_obs = {
        "spliced" :  df_main.to_numpy(),
        "unspliced" : df_main.to_numpy(),
    }
    adata2 = anndata.AnnData(X=df_main.to_numpy(), obs=df_obs, var=adata.var, layers=layers_obs)

    adata2 = dynamoProcess(adata2)
    
    print("Create Vector Field ...")

    v = np.split(adata2.obsm['X_umap'], 2)
    velocity_umap = v[0] - v[1]
    
    adata3 = adata.copy()
    adata3.obsm['X_umap'] = v[0]
    adata3.obsm['velocity_umap'] = velocity_umap
    adata3.uns['neighbors']['indices'] = np.split(adata2.uns['neighbors']['indices'],2)[0]
    nrows = len(adata3.uns['neighbors']['indices'])
    for i in range(nrows) : 
        # adata3.uns['neighbors']['indices'][i] = adata3.uns['neighbors']['indices'][i][adata3.uns['neighbors']['indices'][i] >= nrows]
        for j in range(len(adata3.uns['neighbors']['indices'][i])-2,-1,-1) :
            if adata3.uns['neighbors']['indices'][i][j] >= nrows : 
                adata3.uns['neighbors']['indices'][i][j] = adata3.uns['neighbors']['indices'][i][j+1]

    return adata3

def knnToVectorField(adata) :   
    print("KNN")
    adata4 = adata.copy()
    velocity_umap = adata.obsm['velocity_umap']
    neighbors = adata.uns['neighbors']['indices']
    for indices in adata.uns['neighbors']['indices']:
        i = indices[0]
        adata4.obsm['velocity_umap'][i] = averageVectorField(velocity_umap,indices)
    return adata4

def averageVectorField(data,indices ):
    indices = indices[indices<len(data)]
    filter_v = np.take(data, indices, axis=0)    
    return  np.average(filter_v,axis=0)

