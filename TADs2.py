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


# # Function:
def get_Overlap(a,b):
    """
        Requires: a= region start and end (list), b=  chromStart and chromEnd from rR (list)
    """
    if a[0] == b[1] or b[0] == a[1]:
        return True
    
    ovl=max(0,min(a[1],b[1])-max(a[0],b[0]))#calcula o overlap
    if ovl != 0:
        return True
    return False #ovl
# In[3]:

def closest_to_region(region_start,region_end,chrom,df_TADs):
    """Associates each region to a CpG_Island
        Requieres: stat and end int value of the region, chromossome ,and the CpG_Islands dataframe
        
        Ensures: a list, row, with the data fromthe closest CpG island 
    """
    selected_island = [np.nan,np.nan,np.nan] #DEBUG #This'll cause missing data!
    #Start by each region of df_regions:
    
    list_results = [] # key is a CNV id
       
    df_TADs = df_TADs.loc[df_TADs['chrom'] == chrom]
    df_TADs = df_TADs.sort_values(by=['start'])
    #display(df_TADs) #DEBUG

    dist=999**99 #!!!
    #closest_island_dist = np.nan #!!!
        
    for index, row in df_TADs.iterrows():
        b = [row['start'],row['end']]
        overlap  = get_Overlap([region_start,region_end],b) #region overlap
        #print([region_start,region_end],'\n',b,'\n',overlap)#DEBUG
        
        if overlap == True:# overlap != 0:
            #print(overlap)#DEBUG
            #print([region_start,region_end],b)#DEBUG
            #dist = 0
            if (region_end >= row['start'] and region_start <= row['start']):
                dist = 0
                selected_island = row
                break         
            elif (region_start <= row['end'] and region_end >= row['end']):
                dist = 0
                selected_island = row
                break
                
            else:
                dist =  region_start - row['start']
                selected_island = row
                #print(dist,region_start,row['start'])#DEBUG
                #print(region_start,region_end)#DEBUG
                #print(row['end',row['end']])#DEBUG
                break #?
            #Done -----
        elif overlap == False and dist !=0:
            posA = (region_start - row['end'])  
            #posB = (row['start'] - region_end)  # 
            
            if posA >= 0:
                new_dist = posA

                if new_dist < dist:
                    dist = new_dist
                    selected_island = row
            else:
                break
                #break here!!!
     
    if type(selected_island) != list:    
        selected_island = selected_island.tolist()
        
    #selected_island = selected_island.tolist() #selecteted_item = sub_df_Centro.iloc[position_item]
    
    #distance_log
    selected_island.append(dist)
    
    dist_Mb = float(dist) / 1000000 # to Mb
    distance_log = np.log10(dist_Mb +1)
    

    selected_island.append(distance_log)
    return selected_island #add_list       
#____Region_B:______________________________________________________
def closest_to_region_B(region_start,region_end,chrom,df_TADs):
    """Associates each region to a CpG_Island
        Requieres: stat and end int value of the region, chromossome ,and the CpG_Islands dataframe
        
        Ensures: a list, row, with the data fromthe closest CpG island 
    """
    selected_island = [np.nan,np.nan,np.nan] #DEBUG #This'll cause missing data!
    #Start by each region of df_regions:
    
    list_results = [] # key is a CNV id
       
    df_TADs = df_TADs.loc[df_TADs['chrom'] == chrom]
    df_TADs = df_TADs.sort_values(by=['start'])
    #display(df_TADs) #DEBUG
    
    #Since for Region_B its considered the closest from its right, so invert order from 5'-3' to 3'-5', "degrowing_order" :
    df_TADs = df_TADs.iloc[::-1] #B!
    #
    dist=999**99 #!!!
    #closest_island_dist = np.nan #!!!
        
    for index, row in df_TADs.iterrows():
        b = [row['start'],row['end']]
        overlap  = get_Overlap([region_start,region_end],b)
        #print([region_start,region_end],'\n',b,'\n',overlap)#DEBUG
        
        if overlap == True:# overlap != 0:
            #print(overlap)#DEBUG
            #dist = 0
            if (region_end >= row['start'] and region_start <= row['start']):
                dist = 0
                selected_island = row
                break
            elif (region_start <= row['end'] and region_end >= row['end']):
                dist = 0
                selected_island = row
                break
                
            else:
                dist =  row['end'] - region_end
                selected_island = row
                break #?
            #Done -----
        elif overlap == False and dist !=0:
            #posA = (region_start - row['end']) # start_CNV_region - repRegion_end  #B?
            posB = (row['start'] - region_end) # repRegion_start - end_CNV_region #B?
            
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
        selected_island = selected_island.tolist()     
    #selected_island = selected_island.tolist() #selecteted_item = sub_df_Centro.iloc[position_item]
    
    #distance_log
    selected_island.append(dist)
    
    dist_Mb = float(dist) / 1000000 # to Mb
    distance_log = np.log10(dist_Mb +1) #was log e
    
    
    selected_island.append(distance_log)
    return selected_island #add_list
#___________________________________________________________________
#___________________________________________________________________
# In[4]:


def main_get_TADs(df_regions,df_TADs):
    """
    """
    list_TADs = []
    #Line by line from the A nad B regions:
    for index, row in df_regions.iterrows():
        #A
        rR_A = closest_to_region(row['regionA_stat'],row['regionA_end'],row['chr_A'],df_TADs)
        #B
        rR_B = closest_to_region_B(row['regionB_stat'],row['regionB_end'],row['chr_B'],df_TADs) #B!
        
        #Prep info by line:
        list_add = [row['ID'],row['Cluster_id']] + rR_A + rR_B
        list_TADs.append(list_add)
    return list_TADs    


# # Load Data

# In[5]:

region_path = argv[1]
df_regions = pd.read_csv(region_path, sep= ";" ,header=0,index_col=0)
#df_regions = pd.read_csv('./dataset3/rep_regions/match_CNVs_regions_TP.csv', sep= ";" ,header=0,index_col=0)

df_regions['chr_A'] = df_regions['chr_A'].astype(str)
df_regions['chr_B'] = df_regions['chr_B'].astype(str)


# In[6]:


#df_regions


# ## Tags

# In[7]:


df_TADs= pd.read_csv('./dataset3/TAD_boundaries/'+'IMR90_Rao_2014-raw_TADs.txt', sep= "\t",
                    names=['chrom','start','end'])


# In[8]:


df_TADs['chrom'] = df_TADs['chrom'].astype(str)
df_TADs['chrom'] = df_TADs['chrom'].replace(['23','24'],['X','Y'])


# In[9]:


#df_TADs


# ### Run

# In[10]:


path_CpG_Island = './dataset3/TAD_boundaries/'

TADs_Data = main_get_TADs(df_regions,df_TADs)


# In[11]:


#print(TADs_Data[:10])


# ### To a dataframe:

# In[12]:


#To a dataframe:
df_TADs_AB_data= pd.DataFrame(TADs_Data, columns=['CNV_ID','Cluster_ID',
                                                          'A_chr','A_start','A_end','A_dist_Bp','A_dist_Log',
                                                          'B_chr','B_start','B_end','B_dist_Bp','B_dist_Log'
                                                         ])


# In[13]:


#df_TADs_AB_data


# ## Save to a csv:

# In[14]:


#path_CpG_Island = './dataset3/CpG_islands/'


# In[15]:


#Save to a csv:
output_path = argv[2]
df_TADs_AB_data.to_csv(output_path + 'TAD_boundary_AB_test.csv', sep=';',index_label=False)
#df_TADs_AB_data.to_csv(path_CpG_Island +'TAD_boundary_AB.csv', sep=';',index_label=False)

