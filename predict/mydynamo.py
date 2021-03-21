import warnings
import dynamo as dyn
import pandas as pd
import numpy as np
import anndata
from sklearn.neighbors import KDTree
from sklearn.metrics import mean_squared_error



def config() :
    warnings.filterwarnings('ignore')
    dyn.get_all_dependencies_version()

def removeOutlier(X_list):
    #CHECK FIRST
    if len(X_list) == 0:
        return []

    #REMOVE DUPLICATE
    intersect_index = set(list(X_list[0].index))
    for df in X_list :
        intersect_index = intersect_index & set(list(df.index))      
    intersect_index = list(intersect_index)
    intersect_index.sort()
    intersect_index

    #PACK ARRAY
    X_list_new = []
    for df in X_list :
        X_list_new.append(df.loc[intersect_index])
    return X_list_new

def removeOutlierObservation(a,b, obs):
    # DEPLICATED
    
    intersect_index = list( set(list(a.index)) & set(list(b.index)) & set(list(obs.index)) )
    intersect_index.sort()
    intersect_index

    a = a.loc[intersect_index]
    b = b.loc[intersect_index]
    obs = obs.loc[intersect_index]
    return [a,b,obs]

def removeOutlierObservationExtra(a, obs):
    # DEPLICATED
    
    intersect_index = list( set(list(a.index)) & set(list(obs.index)) )
    intersect_index.sort()
    intersect_index

    a = a.loc[intersect_index]
    obs = obs.loc[intersect_index]
    return [a,obs]

def createAnnotationDataForVectorField(X_list, obs_list, var) :
    X_all = pd.concat(X_list)
    obs_all = pd.concat(obs_list)

    layers = {
        "spliced" : X_all.to_numpy(),
        "unspliced" : X_all.to_numpy(),
    }
    adata = anndata.AnnData(X=X_all.to_numpy(), obs=obs_all,var=var, layers=layers)
    return adata

    # print("555")
    # print("888")

def createAnnotationData(a,b, obs, var, filter=True) :
    # DEPLICATED
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
    # DEPLICATED
    
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

def createVectorFieldForEvaluation(adata) :    
    print("Create Vector Field ...")

    v = np.split(adata.obsm['X_umap'], 3)
    # velocity_umap = v[1] - v[0]
    velocity_umap = adata.obsm['X_umap'] - np.concatenate((v[0], v[0], v[2]), axis=0)

    # PREPARE
    adata3 = adata.copy()
    # adata3.obsm['X_umap'] = v[0]
    adata3.obsm['velocity_umap'] = velocity_umap
    # adata3.uns['neighbors']['indices'] = np.split(adata2.uns['neighbors']['indices'],2)[0]
    # nrows = len(adata3.uns['neighbors']['indices'])
    # for i in range(nrows) : 
    #     # adata3.uns['neighbors']['indices'][i] = adata3.uns['neighbors']['indices'][i][adata3.uns['neighbors']['indices'][i] >= nrows]
    #     for j in range(len(adata3.uns['neighbors']['indices'][i])-2,-1,-1) :
    #         if adata3.uns['neighbors']['indices'][i][j] >= nrows : 
    #             adata3.uns['neighbors']['indices'][i][j] = adata3.uns['neighbors']['indices'][i][j+1]

    return adata3

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

def knnToVectorField(adata, k=30) :   
    print("KNN")
    adata4 = adata.copy()
    # KNN
    # X = np.array([[-1, -1], [-2, -1], [-3, -2], [1, 1], [2, 1], [3, 2]])
    X = adata4.obsm['X_umap']
    kdt = KDTree(X, leaf_size=30, metric='euclidean')
    adata4.uns['neighbors']['indices'] = kdt.query(X, k=len(X), return_distance=False)

    # adata4 = adata.copy()
    velocity_umap = adata4.obsm['velocity_umap']
    neighbors = adata4.uns['neighbors']['indices']
    
    # FILTER CURRENT ONLY + LIMIT with K 
    new_neighbors = []
    nrows = len(neighbors)
    # k = 60
    # GENERATE NEIGHBOR
    for i in range(nrows):
        indices = neighbors[i]
        temp = neighbors[i][ np.logical_and(indices >= nrows/3 , indices < nrows/3*2) ]
        new_neighbors.append( temp[:k] )
    #     adata4.obsm['velocity_umap'][i] = averageVectorField(velocity_umap,indices)
        # UPDATE velocity_umap
        if i >= nrows/3 and i < nrows/3*2 :     
            adata4.obsm['velocity_umap'][i] = averageVectorField(velocity_umap, np.array(new_neighbors) )
        else :
            adata4.obsm['velocity_umap'][i] = np.array([0,0])
    # UPDATE NEIGHBOR
    adata4.uns['neighbors']['indices'] = np.array(new_neighbors)
    
    return adata4

def averageVectorField(data,indices ):
    indices = indices[indices<len(data)]
    filter_v = np.take(data, indices, axis=0)    
    return  np.average(filter_v,axis=0)

def evaluate(adata) :
    X = adata.obsm['X_umap']
    V = adata.obsm['velocity_umap']
    
    Xs = np.split(X,3)
    Vs = np.split(V,3)

    #[0] : HISTORY
    #[1] : PRESENT
    #[2] : FUTURE
    mean_squared_errors = []
    X_predict = []    
    X_true = []
    for i in range(len(Xs[1])) :
        X_predict.append(Xs[1][i] + Vs[1][i])
        X_true.append(Xs[2][i])
    return mean_squared_error(X_true, X_predict)
        
