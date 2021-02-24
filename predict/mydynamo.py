import warnings
import dynamo as dyn
import pandas as pd
import numpy as np
import anndata


def config() :
    warnings.filterwarnings('ignore')
    dyn.get_all_dependencies_version()

def createAnnotationData(a,b, obs, var) :
    intersect_index = list( set(list(a.index)) & set(list(b.index)) & set(list(obs.index)) )
    intersect_index.sort()
    intersect_index

    a = a.loc[intersect_index]
    b = b.loc[intersect_index]
    obs = obs.loc[intersect_index]
    
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

def dynamoPlot(adata) :
    dyn.pl.streamline_plot(adata, color=['parentName'], basis='umap', show_legend='on data', show_arrowed_spines=True)
   

def hello2() :
    # print("555")
    print("hello2")

