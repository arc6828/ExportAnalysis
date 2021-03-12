import pandas as pd
import os
import numpy as np

def create() :
    big_list = []
    for filename in os.listdir(folder_name) :
        #print(filename)
        path = folder_name+"/"+filename       
        print(path, end='\r')
        df = pd.read_json(path, orient='values')
    #     print("HELLO : ",df[0][3])
        df = pd.read_html(df[0][3])[0].replace(np.nan,0)    
        
        
        for i in range(len(df)) : 
            if len(df) == 1 :  continue
            for j in range(2,len(df.columns)) : 
                if df.iloc[i, 1] == "World" :
                    continue
                item = [getYearFromString(df.columns[j]), df.iloc[i, 1], "world", filename.replace(".json",""), df.iloc[i, j] ] 
                #print( item )
                big_list.append(item)
        
        #break
    # df
    # 