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
def get_Overlap(a,b):
    """
        Requires: a= region start and end (list), b=  chromStart and chromEnd from rR (list)
    """
    ovl=max(0,min(a[1],b[1])-max(a[0],b[0]))#calcula o overlap
    if ovl != 0:
        return True
    return False #ovl
# In[3]:


def calc_GC_content_val(list_data):
    """
        faz o calculo do GC content para a região 
        GC data percentage is in intervals of 5 bps.
    """
    GC_content = sum(list_data)/(5*len(list_data))  
    return GC_content #/100



# In[4]:
def GC_match(path_GCcontent,start,end):
    """Match and Calculation of the GC percentage data
        Requires: start and end positition int values of the region \
            path str to the GC Content file
        Ensures: value int of the GC % of the region
    """
    list_data = []
    GC_content_val = np.NaN # PODE ACABAR SEM ENCONTRAR % GC ????
    with open(path_GCcontent,'r') as GC_content_file:
        #
        is_in=False
        for line in GC_content_file:
            line = line.split(",")
            el = int(line[0]) #position_start
            GCper = float(line[1]) #GC percentage in 5bp interval
            #DEBUG
            #print(line,type(line))
            #DEBUG
            ##if el >= start and el < end:
            if get_Overlap([start,end],[el,el+5])==True:
                is_in=True
                #adiciona para fazer a percentagem do gc content
                GCbp = (GCper*5)/100 #!
                list_data.append(GCbp)
            elif is_in == True:
                #faz o calculo do GC content para a região - %gc_20/(2*100)
                #guarda o resultado com o ID da variante
                GC_content_val = calc_GC_content_val(list_data)
                break
    return GC_content_val
#
def GC_match_old(path_GCcontent,start,end):
    """Match and Calculation of the GC percentage data
        Requires: start and end positition int values of the region \
            path str to the GC Content file
        Ensures: value int of the GC % of the region
    """
    list_data = []
    GC_content_val = np.NaN # PODE ACABAR SEM ENCONTRAR % GC ????
    with open(path_GCcontent,'r') as GC_content_file:
        #
        is_in=False
        for line in GC_content_file:
            line = line.split(",")
            el = int(line[0]) #position_start
            GCper = float(line[1]) #GC percentage in 5bp interval
            #DEBUG
            #print(line,type(line))
            #DEBUG
            if el >= start and el < end:
                is_in=True
                #adiciona para fazer a percentagem do gc content
                list_data.append(GCper)
            elif is_in == True:
                #faz o calculo do GC content para a região - %gc_20/(2*100)
                #guarda o resultado com o ID da variante
                GC_content_val = calc_GC_content_val(list_data)
                break
    return GC_content_val


# In[5]:


def main_regionsAB(df_regions,path_GC_content):
    """
    """
    GC_regions = []
    n = 0#DEBUG
    for index, row in df_regions.iterrows():
        n += 1#DEBUG
 
        #A
        chr_ID = row["chr_A"]
        path_GCcontent = path_GC_content + 'GCcontent_chr{}.csv'.format(chr_ID)
        #print(row["regionA_stat"],type(row["regionA_stat"]))#DEBUG
        GC_A = GC_match(path_GCcontent,row["regionA_stat"],row["regionA_end"])
        
        #B
        chr_ID = row["chr_B"]
        path_GCcontent = path_GC_content + 'GCcontent_chr{}.csv'.format(chr_ID)
        GC_B = GC_match(path_GCcontent,int(row["regionB_stat"]),int(row["regionB_end"]))
        
        #
        regionA = [row['ID'],row['Cluster_id'], 'A' ,GC_A] #TO ADD??? row['Cluster_id'],row['region_stat'],row['region_end']
        regionB = [row['ID'],row['Cluster_id'], 'B' ,GC_B] #TO ADD??? row['Cluster_id'],row['region_stat'],row['region_end']
        #
        GC_regions.append(regionA)
        GC_regions.append(regionB) 
            
        #if n > 2:#DEBUG
        #    return GC_regions#DEBUG
        #    break#DEBUG
    return GC_regions


# # Load Data

# In[6]:

region_path = argv[1]
#'./dataset3/rep_regions/match_CNVs_regions_TP.csv'

df_regions = pd.read_csv(region_path, sep= ";" ,header=0,index_col=0)


# In[7]:


#df_regions.head()


# In[8]:


#
path_GC_content = './dataset3/GCcontent/'

list_regionsAB = main_regionsAB(df_regions,path_GC_content) #Tudo menos telemoero e centromerlos.


# print(list_regionsAB)
# """
# [['DGRC_DEL_31', 'A', 0.2970464135021097], ['DGRC_DEL_31', 'B', 0.33774834437086093], ['DGRC_DEL_666', 'A', 0.425], ['DGRC_DEL_666', 'B', 0.4942857142857143], ['DGRC_DEL_1163', 'A', 0.605586592178771], ['DGRC_DEL_1163', 'B', 0.5761421319796954]]"""

# ### Make dataframe:

# In[9]:


df_regionsAB = pd.DataFrame(list_regionsAB, columns=['CNV_ID','Cluster_ID','Region','GC_%'])


# In[10]:


#df_regionsAB


# #### Save df to file:

# In[11]:
output_path = argv[2]
#output_path = path_GC_content+'GCcontent_regionsAB.csv'
df_regionsAB.to_csv(output_path + 'GCcontent_regionsAB.csv' , sep=';',index_label=False)

