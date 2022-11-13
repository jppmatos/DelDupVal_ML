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
    ovl=max(0,min(a[1],b[1])-max(a[0],b[0]))#calcula o overlap
    return ovl


# # Function:

# In[4]:


def get_Lamina_scores(region_stat,region_end,chrom,df_Lamina):
    """
    """
    result = []
    #a=[region_stat,region_end]#region
    df_Lamina.loc[df_Lamina["chrom"] == str(chrom)]
    
    #Line by line form df_Lamina:
    for index, row in df_Lamina.iterrows():
  
        b = [row['start'],row['end']]#gene
        #------------------------------------------------------------------------------------
        #associated genes
        if get_Overlap([region_stat,region_end],b) != 0:
            #df_exons = df_genes.loc[(df_genes["gene_name"] == row['gene_name']) & (df_genes["gene_name"] == 'exon')]
            #cov_Exons, cov_Introns = coverage(region_stat,region_end,df_exons) 
            #if cov_Exons != [] and cov_Introns != []:#DEBUG
            #    print("|region start & end:",region_stat,region_end,'\n|chrom start & end:',b)#DEBUG
            #    print("\nCov_Exons:",cov_Exons, "\ncov_Introns:",cov_Introns) #DEBUG
            #    print("\n _________________________________________") #DEBUG   
            #result += [[row['gene_name'],cov_Exons, cov_Introns]]
            result += [row['score']]
    result = np.mean(result)
    return result


# In[5]:


def main_Lamina_Dom_means(df_regions,df_Lamina):
    """
    """
    list_Lamina_means = []
    #Line by line from the A nad B regions:
    for index, row in df_regions.iterrows():
        #A
        rR_A = get_Lamina_scores(row['regionA_stat'],row['regionA_end'],row['chr_A'],df_Lamina)
        #B
        rR_B = get_Lamina_scores(row['regionB_stat'],row['regionB_end'],row['chr_B'],df_Lamina)
        #Prep info by line:
        
        list_add = [row['ID'],row['Cluster_id'],rR_A, rR_B] 
        list_Lamina_means.append(list_add)
    return list_Lamina_means


# # Load Data

# In[6]:

region_path = argv[1]
df_regions = pd.read_csv(region_path, sep= ";" ,header=0,index_col=0)
#df_regions = pd.read_csv('./dataset3/rep_regions/match_CNVs_regions_TP.csv', sep= ";" ,header=0,index_col=0)

df_regions['chr_A'] = df_regions['chr_A'].astype(str)
df_regions['chr_B'] = df_regions['chr_B'].astype(str)


# In[7]:


#df_regions


# ### Lamina File:

# In[8]:


df_lamina_domains= pd.read_csv('./dataset3/Lamina_associated_domains/'+'41586_2008_BFnature06947_MOESM252_ESM.txt', sep= "\t",  
                            names = ['chrom','name_paper','name','start','end','score','|','||','|||'])


# In[9]:


df_lamina_domains["chrom"] = df_lamina_domains["chrom"].transform(filter_chrName)
#df_lamina_domains


# ### Run

# In[10]:


path_Lamina_data = './dataset3/Lamina_associated_domains/'

Lamina_Data = main_Lamina_Dom_means(df_regions,df_lamina_domains)


# In[11]:


#Lamina_Data[0] #cov_Exons, cov_Introns
#:> ['DGRC_DEL_31', 1139, 0.9191106952443581, 0.9191106952443581]


# In[16]:


#To a dataframe:
df_Lamina_scores= pd.DataFrame(Lamina_Data, columns=['CNV_ID','Cluster_ID','A_score','B_score'])


# In[17]:


#df_Lamina_scores


# ## Save to a csv:

# In[18]:


#Save to a csv:
output_path = argv[2]
df_Lamina_scores.to_csv(output_path + 'Lamina_scores_regionsAB.csv', sep=';',index_label=False)
#df_Lamina_scores.to_csv(path_Lamina_data+'Lamina_scores_regionsAB.csv', sep=';',index_label=False)

