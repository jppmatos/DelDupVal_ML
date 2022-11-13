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


# # Function:

def get_Overlap_plus(a,b):
    """
        Requires: a= region start and end (list), b=  chromStart and chromEnd from rR (list)
        it count as overlap when the end of interval_A matches with the start of inverval_B, nad vice-versa.
    """
    if a[0] == b[1] or b[0] == a[1]: # 
        return True
    
    ovl=max(0,min(a[1],b[1])-max(a[0],b[0]))#calcula o overlap
    if ovl != 0:
        return True
    return False #ovl
# In[3]:


def ReplicationTime(star_pos,end_pos,chr_id,path_RT_content):#,celllines_list):
    """
        Requires: path_RT_content is replication time scores 
    """
    #paths folders
    #celllines_list = ['BG02ES/','GM06990/','H7-hESC/','H9-hESC/','HeLa-S3/','IMR90/','iPS_hFib2_iPS4/','iPS_hFib2_iPS5/']
    #path_RT_content = './dataset3/ReplicationTime/'
    #
    CDs_RTimes_list = []
    #
    #
    df_RepT = pd.read_csv(path_RT_content, sep='\t',names=['chrom','chromStart','chromEnd','score']) #header: '#chrom','chromStar','chromEnd','score'
    #df_RepT["chrom"] = df_RepT["chrom"].astype(str)# Needed because of chr x and Y.
    df_RepT["chrom"] = df_RepT["chrom"].transform(filter_chrName)
    df_RepT = df_RepT.loc[df_RepT['chrom'] == chr_id] #reduce
    #new: ------
    #reduce size hack:
    df_RepT = df_RepT.loc[df_RepT['chromStart'] <= end_pos]#hack
    for index, row_RepT in df_RepT.iterrows():
        b = [row_RepT['chromStart'],row_RepT['chromEnd']]
        overlap  = get_Overlap_plus([star_pos,end_pos],b) #region overlap
        if overlap == True:# overlap != 0:
            #seleted_row_RepT = row_RepT
            rep_Time = row_RepT['score']
            CDs_RTimes_list.append(rep_Time)
            
    if CDs_RTimes_list != []:
        RepTime_mean = np.mean(CDs_RTimes_list) 
        CDs_RTimes_list = RepTime_mean 
    else:
        CDs_RTimes_list = np.nan
        
    return CDs_RTimes_list        
        

# In[5]:


def main_ReplicationTime(df_regions,Path):#,celllines_list):
    """
        Requires: ex: path = './dataset3/ReplicationTime/ReplicationTimes_mean.csv' , and regions file
    """
    list_Reg_Dust = []
    for index, row in df_regions.iterrows():
        #ReplicationTimes:
        #A
        RT_A = ReplicationTime(row['regionA_stat'],row['regionA_end'],row['chr_A'],Path)#,celllines_list)
        #B
        RT_B = ReplicationTime(row['regionB_stat'],row['regionB_end'],row['chr_B'],Path)#,celllines_list)
        #
        #print(RT_A,'\n',RT_B,'\n')#DEBUG
        list_Reg_Dust.append([row['ID'],row['Cluster_id'],RT_A,RT_B])
    return list_Reg_Dust
    
        


# # Load Data

# In[6]:

region_path = argv[1]
df_regions = pd.read_csv(region_path, sep= ";" ,header=0,index_col=0)
#df_regions = pd.read_csv('./dataset3/rep_regions/match_CNVs_regions_TP.csv', sep= ";" ,header=0,index_col=0)
df_regions['chr_A'] = df_regions['chr_A'].astype(str)
df_regions['chr_B'] = df_regions['chr_B'].astype(str)


# In[7]:


#df_regions


# ## Run

# In[8]:


#path_RT_content = './dataset3/ReplicationTime/ReplicationTimes_mean.csv'
path_RT_content = './dataset3/ReplicationTime/hg38_mean_repTime_SCORES'
#celllines_list = ['BG02ES/','GM06990/','H7-hESC/','H9-hESC/','HeLa-S3/','IMR90/','iPS_hFib2_iPS4/','iPS_hFib2_iPS5/']


# In[9]:


list_ReplicationTimes = main_ReplicationTime(df_regions,path_RT_content)#,celllines_list)


# ### Make dataframe:

# In[10]:


df_RepTimes = pd.DataFrame(list_ReplicationTimes, columns=['CNV_ID','Cluster_ID','ReplicationTime_A','ReplicationTime_B'])


# In[11]:


#df_RepTimes


# #### Save df to file:
# In[12]:
output_path = argv[2]
df_RepTimes.to_csv(output_path + 'ReplicationTimes_regionsAB.csv', sep=';',index_label=False)
#df_RepTimes.to_csv(path_RT_content+'ReplicationTimes_regionsAB.csv', sep=';',index_label=False)

