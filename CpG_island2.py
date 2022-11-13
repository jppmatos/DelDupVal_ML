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
        #if '_' in chromossome:
        #    chromossome = chromossome.split("_")[0]
 
    
    return str(chromossome)


# In[3]:


def get_Overlap(a,b):
    """
        Requires: a= region start and end (list), b=  chromStart and chromEnd from rR (list)
    """
    ovl=max(0,min(a[1],b[1])-max(a[0],b[0]))#calcula o overlap
    if ovl != 0:
        return True
    return False #ovl


# # Function:

# In[4]:


def closest_to_region(region_stat,region_end,chrom,df_CpG_Islands):
    """Associates each region to a CpG_Island
        Requieres: stat and end int value of the region, chromossome ,and the CpG_Islands dataframe
        
        Ensures: a list, row, with the data fromthe closest CpG island 
    """
    selected_island = [np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan] #DEBUG #This'll cause missing data!
    #Start by each region of df_regions:
    
       
    df_CpG_Islands = df_CpG_Islands.loc[df_CpG_Islands['chrom'] == chrom]
    df_CpG_Islands = df_CpG_Islands.sort_values(by=['chromStart'])
    #display(df_CpG_Islands) #DEBUG

    dist=999**99 #!!!
    closest_island_dist = np.nan #!!!
        
    for index, row in df_CpG_Islands.iterrows():
        b = [row['chromStart'],row['chromEnd']]
        overlap  = get_Overlap([region_stat,region_end],b)
        #print([region_stat,region_end],'\n',b,'\n',overlap)#DEBUG
        
        if overlap == True:# overlap != 0:
            #print(overlap)#DEBUG
            #new_dist = 0
            #dist = overlap
            selected_island = row
            dist = 0
            #Done -----
        elif overlap == False and dist !=0:
            posA = (region_stat - row['chromEnd'])  # item() ?
            #posB = (row['chromStart'] - region_end)  # item() ?
            
            if posA >= 0:
                new_dist = posA

            #elif posB >= 0:
            #    new_dist = posB

            #elif posA >= 0 & posB >= 0:
            #    new_dist = min([posA, posB])

                if new_dist < dist:
                    dist = new_dist
                    selected_island = row
            else:
                break
                #break here!!!
        
    if type(selected_island) != list:    
        selected_island = selected_island.tolist() #selecteted_item = sub_df_Centro.iloc[position_item]
        
    
    #distance_log
    dist = float(dist) / 1000000 # to Mb
    distance_log = np.log10(dist +1) #was log e
    
    selected_island.append(distance_log)
    return selected_island #add_list       
#____Region_B:______________________________________________________
def closest_to_region_B(region_stat,region_end,chrom,df_CpG_Islands):
    """Associates each region to a CpG_Island
        Requieres: stat and end int value of the region, chromossome ,and the CpG_Islands dataframe
        
        Ensures: a list, row, with the data fromthe closest CpG island 
    """
    selected_island = [np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan] #DEBUG #This'll cause missing data!
    #Start by each region of df_regions:
    
       
    df_CpG_Islands = df_CpG_Islands.loc[df_CpG_Islands['chrom'] == chrom]
    df_CpG_Islands = df_CpG_Islands.sort_values(by=['chromStart'])
    #display(df_CpG_Islands) #DEBUG
    
    #Since for Region_B its considered the closest from its right, so invert order from 5'-3' to 3'-5', "degrowing_order" :
    df_CpG_Islands = df_CpG_Islands.iloc[::-1] #B!
    #
    dist=999**99 #!!!
    #closest_island_dist = np.nan #!!!
        
    for index, row in df_CpG_Islands.iterrows():
        b = [row['chromStart'],row['chromEnd']]
        overlap  = get_Overlap([region_stat,region_end],b)
        #print([region_stat,region_end],'\n',b,'\n',overlap)#DEBUG
        
        if overlap == True:# overlap != 0:
            #print(overlap)#DEBUG
            #new_dist = 0
            #dist = overlap
            selected_island = row
            dist = 0
            #Done -----
        elif overlap == False and dist !=0:
            #posA = (region_stat - row['chromEnd']) # start_CNV_region - repRegion_end  #B?
            posB = (row['chromStart'] - region_end) # repRegion_start - end_CNV_region #B?
            
            if posB >= 0:
                new_dist = posB

            #elif posA >= 0:
            #    new_dist = posA

            #elif posA >= 0 & posB >= 0:
            #    new_dist = min([posA, posB])

                if new_dist < dist:
                    dist = new_dist
                    selected_island = row
            else:
                break
                #break here!!!
        
    if type(selected_island) != list:    
        selected_island = selected_island.tolist() #selecteted_item = sub_df_Centro.iloc[position_item]       
    #selected_island = selected_island.tolist() #selecteted_item = sub_df_Centro.iloc[position_item]
    
    #distance_log
    dist = float(dist) / 1000000 # to Mb
    distance_log = np.log10(dist +1) #was log e
    
    selected_island.append(distance_log)
    return selected_island #add_list
#___________________________________________________________________
# In[5]:


def main_get_CpG_Island(df_regions,df_CpG_Islands):
    """
    """
    list_Genes = []
    #Line by line from the A nad B regions:
    for index, row in df_regions.iterrows():
        #A
        rR_A = closest_to_region(row['regionA_stat'],row['regionA_end'],row['chr_A'],df_CpG_Islands)
        #B
        rR_B = closest_to_region_B(row['regionB_stat'],row['regionB_end'],row['chr_B'],df_CpG_Islands) #B!
        
        #Prep info by line:
        list_add = [row['ID'],row['Cluster_id']] + rR_A + rR_B
        list_Genes.append(list_add)
    return list_Genes


# # Load Data

# In[6]:

region_path = argv[1]
df_regions = pd.read_csv(region_path, sep= ";" ,header=0,index_col=0)
#df_regions = pd.read_csv('./dataset3/rep_regions/match_CNVs_regions_TP.csv', sep= ";" ,header=0,index_col=0)

df_regions['chr_A'] = df_regions['chr_A'].astype(str)
df_regions['chr_B'] = df_regions['chr_B'].astype(str)


# In[7]:


#df_regions


# ### CpG islands:

# In[8]:


df_CpG_Islands= pd.read_csv('./dataset3/CpG_islands/'+'CpG_Islands_info.txt', sep= "\t")


# In[9]:


df_CpG_Islands["chrom"] = df_CpG_Islands["chrom"].transform(filter_chrName)
#df_CpG_Islands


# In[10]:


#remove contigs and others alike (ex: 22_KI270738v1_random) :
df_CpG_Islands = df_CpG_Islands[~df_CpG_Islands.chrom.str.contains('_')]
#checking
#print(df_CpG_Islands["chrom"].unique())


# In[11]:


#df_CpG_Islands


# ### Run

# In[12]:


path_CpG_Island = './dataset3/CpG_islands/'

CpG_Island_Data = main_get_CpG_Island(df_regions,df_CpG_Islands)


# In[13]:


#print(CpG_Island_Data[:10]) #cov_Exons, cov_Introns


# In[14]:


#To a dataframe:
df_CpG_Island_AB_data= pd.DataFrame(CpG_Island_Data, columns=['CNV_ID','Cluster_ID',
                                                          'A_#bin','A_chrom','A_chromStart','A_chromEnd','A_name', 'A_length', 'A_cpgNum', 'A_gcNum', 'A_perCpg', 'A_perGc', 'A_obsExp', 'A_dist',
                                                          'B_#bin','B_chrom','B_chromStart','B_chromEnd','B_name', 'B_length', 'B_cpgNum', 'B_gcNum', 'B_perCpg', 'B_perGc', 'B_obsExp', 'B_dist'
                                                         ])


# In[15]:


#df_CpG_Island_AB_data


# ## Save to a csv:

# In[19]:


#path_CpG_Island = './dataset3/CpG_islands/'


# In[20]:


#Save to a csv:
output_path = argv[2]
df_CpG_Island_AB_data.to_csv(output_path + 'CpG_Island_AB.csv', sep=';',index_label=False)
#df_CpG_Island_AB_data.to_csv(path_CpG_Island +'CpG_Island_AB.csv', sep=';',index_label=False)

