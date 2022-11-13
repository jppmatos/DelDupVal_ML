#!/usr/bin/env python
# coding: utf-8
from sys import argv

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
from scipy.stats import norm
from sklearn.preprocessing import StandardScaler
from scipy import stats
#import warnings
#warnings.filterwarnings('ignore')

# # Load data:


#df_Centromeres = pd.read_csv('dataset3/Centromero&Telomero/centromeres_bands.txt', sep= "\t" ,header=0)

#df_Centromeres

# ## Functions:

# ### Suplementar functions:


def filter_chrName(v):
    """Filter string like Yp11.2Yq11.223 to Y """
    v = str(v)
    #print(v)
    chromossome = ''
    if 'q' in v:
        chromossome = v.split("q")[0]
        if 'p' in chromossome:
            chromossome = chromossome.split("p")[0]
            
    elif 'p' in v:
        chromossome = v.split("p")[0]
        if 'q' in chromossome:
            chromossome = chromossome.split("q")[0]
            
    elif 'chr' in v:
        chromossome = v.split("chr")[1]
        if '_' in chromossome:
            chromossome = chromossome.split("_")[0]
    else:
        chromossome = v
 
    
    return str(chromossome)


# In[6]:


def get_Overlap(a,b):
    ovl=max(0,min(a[1],b[1])-max(a[0],b[0]))#calcula o overlap
    #sz_a=a[1]-a[0] #calcula o tamanho da região definida pelo agrupamento
    #sz_b=b[1]-b[0] #calcula o tamanho da deleção
    
    #if sz_a <= 0 or sz_b <= 0:
    #    return False
    
    #ta=ovl/sz_a # calcula a proporção da região definida pelo cluster que se sobrepõe à deleção
    #tb=ovl/sz_b # calcula a proporção da deleção que se sobrepõe à região definida pelo cluster
    #if ta>=0.9 and tb>=0.9:#verifica se ambas as proporções são superiores a 90%
    #        return True
    #else:
    #        return False
    if ovl != 0:
        return True
    return False


# ### Main functions:
# ### With Telomeres:
def lookfor_telomere_dist(CNV_start,CNV_end,sub_df_Telo,position_item):
    """ Look for the distance of the closest telomere.
    """
    if position_item > 0:
        selected_telomere = sub_df_Telo.iloc[1]
        distance = selected_telomere['chromStart'] - CNV_end

    else:
        selected_telomere = sub_df_Telo.iloc[0]
        distance = CNV_start - selected_telomere['chromEnd']

    #
    #
    dist_telo = float(distance) / 1000000 # to Mb
    dist_telo = np.log10(dist_telo +1)
  
    return dist_telo
    

# #### With Centromers:

# In[7]:


def closest_to_CNV(df_CNVs,df_Centro,df_Telo):
    """Associates each CNV_ID (key) to a  Centromere or Telomeres/Gap
        Requieres: two dataframes, one is the CNV and the other is Centromeres or Telomeres/Gaps (df_Centro)
        Ensures: retruns a dict that associates each CNV_ID (key) to a  Centromere or Telomeres/Gap
    """
    #filter chrom column:
    df_Centro['#chrom'] = df_Centro['#chrom'].transform(filter_chrName)
    df_Telo['chrom'] = df_Telo['chrom'].transform(filter_chrName)
    #Start by each CNV of df_CNVs:
    
    dict_items = {} # key is a CNV id
    
    for CNV_id in df_CNVs.ID:
        sub_df_CNV = df_CNVs.loc[df_CNVs['ID'] == CNV_id]
        #
        selected_chrom = sub_df_CNV['chr'].iloc[0]
        sub_df_Centro = df_Centro.loc[df_Centro['#chrom'] == selected_chrom]
        sub_df_Telo = df_Telo.loc[df_Telo['chrom'] == selected_chrom]
        #
        #print("DEBUG2",sub_df_Centro)#DEBUG
        #
        #find_dist(selected_chrom,sub_df_Centro)
        dist=999**99
        position_item = None
        a = [sub_df_Centro['chromStart'].iloc[0].item() ,sub_df_Centro['chromEnd'].iloc[1].item() ] #Centomere
        b = [sub_df_CNV['start'].item() , sub_df_CNV['end'].item() ] #CNV
        if get_Overlap(a,b) == True:
            position_item = -1 # if 0 is 5', if 1 is 3'
            dist = 0
            #
        elif get_Overlap(a,b) == False and dist != 0:
            posA = (sub_df_CNV['start'].item() - a[1])
            posB = (sub_df_CNV['end'].item() - a[0])
            posC = (sub_df_CNV['start'].item() - a[0])
            posD = (sub_df_CNV['end'].item() - a[1])

            orientation = posA # if > 0 , (sub_df_Telo['chromStart'].iloc[i] - sub_df_CNV['end']).item()
            # if < 0, (sub_df_CNV['start'] - sub_df_Telo['chromEnd'].iloc[i]).item()
                
            dist = min([abs(posA), abs(posB), abs(posC), abs(posD)])
         

        dist_centro = float(dist) / 1000000 # to Mb
        dist_centro = np.log10(dist_centro +1)
        #print(dist_centro, dist, b, CNV_id) # DEBUG
        #add_list.append( dist_centro) #log10 +1 aplied!
        
        telomere_dist = lookfor_telomere_dist(b[0],b[1],sub_df_Telo,orientation)


        dict_items[CNV_id]=[CNV_id,dist_centro,telomere_dist]
        
    return dict_items
    
 
    


# ## Load CNVs and start: 

'''PATH_CNVs = argv[1] 
sheetName = argv[2]
headerline = int(argv[3])

df_CNVs = pd.read_excel(PATH_CNVs, header=headerline, sheet_name=sheetName)
df_CNVs['chr'] = df_CNVs['chr'].transform(filter_chrName)

'''
region_path = argv[1]
df_regions = pd.read_csv(region_path, sep= ";" ,header=0,index_col=0)
#df_regions = pd.read_csv('./dataset3/rep_regions/match_CNVs_regions_TP.csv', sep= ";" ,header=0,index_col=0)

df_regions['chr_A'] = df_regions['chr_A'].astype(str)
df_regions['chr_B'] = df_regions['chr_B'].astype(str)

df_CNVs = df_regions[['ID','chr_A','start_position_max','end_position_max']].copy()
df_CNVs.rename(columns={'chr_A': 'chr', 'start_position_max': 'start', 'end_position_max':'end'}, inplace=True) 

#-----------------------------------------------------------------------------------------------------------------
df_Centromeres = pd.read_csv('./dataset3/Centromero&Telomero/centromeres_bands.txt', sep= "\t" ,header=0)
df_TelomeresGaps = pd.read_csv('./dataset3/Centromero&Telomero/Telomeres_just.txt', sep= "\t" ,header=0)
#Labels: #bin	chrom	chromStart	chromEnd	name
#Labels:#bin	chrom	chromStart	chromEnd	ix	n	size	type	bridge

centromers = closest_to_CNV(df_CNVs,df_Centromeres,df_TelomeresGaps)




#print(list(centromers.items())[:3], '\n \n',list(telomeres.items())[:3])




df_selec_centromers = pd.DataFrame.from_dict(centromers, orient='index',columns=['CNV_id','dist_Centromere','dist_Telomere'])




#Save data in .csv
outPATH = argv[2]
#outPATH = "'./dataset3/Centromero&Telomero/"
df_selec_centromers.to_csv(outPATH + 'centromere_telomere_matched.csv', sep = ';' ,index_label=False)
