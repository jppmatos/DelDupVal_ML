#!/usr/bin/env python
# coding: utf-8

# In[1]:
from sys import argv

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
from scipy.stats import norm
from sklearn.preprocessing import StandardScaler
from scipy import stats
import warnings
warnings.filterwarnings('ignore')


# # Axiliar Functions:

# In[2]:


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
        # add if for 'chr22_KI270739v1_random'
        chromossome = v.split("chr")[1]
        if '_' in chromossome:
            chromossome = chromossome.split("_")[0]
 
    
    return str(chromossome)


# In[3]:


def get_Overlap(a,b):
    """
        Requires: a= region start and end (list), b=  chromStart and chromEnd from rR (list)
    """
    #region_start, region_end = a
    #rR_start, rR_end = b
    #if region_start < rR_start and region_end < rR_start:
    #    return False
    #elif region_start > rR_end and region_end > rR_end:
    #    return False
    #return True
    ovl=max(0,min(a[1],b[1])-max(a[0],b[0]))
    if ovl != 0:
        return True
    return False


# # Function:

# In[4]:
def get_RecombRate(start,end,chrom,df_rR_info):
    """
        Ensures: Returns the row (info) form df_rR_info thats overlaps with the regions.
    """
    target = []#,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan] 
    #
    df_rR_info = df_rR_info.loc[df_rR_info["chrom"] == ('chr' + str(chrom))]
    #display(df_rR_info)#DEBUG
    is_in=False #new
    for index, row in df_rR_info.iterrows():
        rR_Start = row['chromStart']
        rR_End = row['chromEnd']
        b = [rR_Start,rR_End]
        #
        a_r = [start,end]
        #print(a_r,' - ',b)#DEBUG
        #print(is_in)#DEBUG
        if get_Overlap(a_r,b) == True:
            is_in=True #new
            #target = [row['name'],row['decodeAvg'],row['decodeFemale'],row['decodeMale'],row['marshfieldAvg'],row['marshfieldFemale'],row['marshfieldMale'],row['genethonAvg'],row['genethonFemale'],row['genethonMale']]
            target.append(row['decodeAvg'])
            #print('targeted',target,is_in)#DEBUG
            #break
        elif is_in == True: #new
            #print(target)#DEBUG
            #decodeAvg_mean = np.mean(target)
            break

    if target == []:
        decodeAvg_mean = np.nan#[]
    else:
        decodeAvg_mean = np.mean(target)
    #if len(target) > 1:#DEBUG
    #    print(target,'->',decodeAvg_mean)#DEBUG
    return decodeAvg_mean
#old
def get_RecombRate_old(start,end,chrom,df_rR_info):
    """
        Ensures: Returns the row (info) form df_rR_info thats overlaps with the regions.
    """
    target = [np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan] 
    df_rR_info = df_rR_info.loc[df_rR_info["chrom"] == ('chr' + str(chrom))]
    #display(df_rR_info)#DEBUG
    for index, row in df_rR_info.iterrows():
        rR_Start = row['chromStart']
        rR_End = row['chromEnd']
        b = [rR_Start,rR_End]
        #
        a_r = [start,end]
        #print(a_r,' - ',b)#DEBUG
        if get_Overlap(a_r,b) == True:
            target = [row['name'],row['decodeAvg'],row['decodeFemale'],row['decodeMale'],row['marshfieldAvg'],row['marshfieldFemale'],row['marshfieldMale'],row['genethonAvg'],row['genethonFemale'],row['genethonMale']]
            break
            
    return target
#/old


# #test
# #Recombination Rate data:
# df_rR_info = pd.read_csv('./dataset3/RecombinationRate/'+'hg38/'+'RecombinationRate_hg_38.csv',
#                          names=['chrom','chromStart','chromEnd','name','decodeAvg','decodeFemale','decodeMale','marshfieldAvg','marshfieldFemale','marshfieldMale','genethonAvg','genethonFemale','genethonMale'])
# 
# get_RecombRate(58877520,58878767,4,df_rR_info)

# In[6]:


def main_RecombRate(df_regions,df_rR_info):
    """
    """
    list_recombRates = []
    for index, row in df_regions.iterrows():
        #A
        rR_A = get_RecombRate(row['regionA_stat'],row['regionA_end'],row['chr_A'],df_rR_info)
        #B
        rR_B = get_RecombRate(row['regionB_stat'],row['regionB_end'],row['chr_B'],df_rR_info)
        #
        list_add = [row['ID'],row['Cluster_id'],rR_A, rR_B]
        #list_add = [row['ID'],row['Cluster_id']] + rR_A + rR_B
        #print(list_add)#DEBUG
        list_recombRates.append(list_add)
    return list_recombRates


# # Load Data

# In[7]:

region_path = argv[1]
df_regions = pd.read_csv(region_path, sep= ";" ,header=0,index_col=0)
#df_regions = pd.read_csv('./dataset3/rep_regions/match_CNVs_regions_TP.csv', sep= ";" ,header=0,index_col=0)
df_regions['chr_A'] = df_regions['chr_A'].astype(str)
df_regions['chr_B'] = df_regions['chr_B'].astype(str)


# In[8]:


#df_regions


# In[9]:


#Recombination Rate data:
df_rR_info = pd.read_csv('./dataset3/RecombinationRate/'+'hg38/'+'RecombinationRate_hg_38.csv',
                         names=['chrom','chromStart','chromEnd','name','decodeAvg','decodeFemale','decodeMale','marshfieldAvg','marshfieldFemale','marshfieldMale','genethonAvg','genethonFemale','genethonMale'])


# In[10]:


#df_rR_info


# ## Run

# In[11]:


path_RR_content = './dataset3/RecombinationRate/'

list_rR = main_RecombRate(df_regions,df_rR_info)


# In[12]:


#list_rR[0]


# In[13]:


#To a dataframe:

#old
#df_rR = pd.DataFrame(list_rR, columns=['CNV_ID','Cluster_ID',
#                                       #'A_name_rR','A_decodeAvg','A_decodeFemale','A_decodeMale','A_marshfieldAvg','A_marshfieldFemale','A_marshfieldMale','A_genethonAvg','A_genethonFemale','A_genethonMale',
                                     # 'B_name_rR','B_decodeAvg','B_decodeFemale','B_decodeMale','B_marshfieldAvg','B_marshfieldFemale','B_marshfieldMale','B_genethonAvg','B_genethonFemale','B_genethonMale'])
#/old
    
df_rR = pd.DataFrame(list_rR, columns=['CNV_ID','Cluster_ID',
                                       'A_decodeAvg',
                                     'B_decodeAvg'])

# In[14]:


#df_rR


# ### Quick analisis:
#Note: Keep it like that, some data was lost in liftover, see liftover error log.


# ## Save to a csv:

# In[20]:


#Save to a csv:
output_path = argv[2]
df_rR.to_csv(output_path + 'RecombinationRate_regionsAB.csv', sep=';',index_label=False)
#df_rR.to_csv('./dataset3/RecombinationRate/'+'RecombinationRate_regionsAB.csv', sep=';',index_label=False)

