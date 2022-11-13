#coding: utf-8
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

# ## Functions:
# ### Suplementar functions:


def filter_chrName(v):
    """Filter string like Yp11.2Yq11.223 to Y """
    v = str(v)
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
 
    
    return str(chromossome)


# In[6]:

def get_Overlap(a,b):
    #print("AB:",a,b)
    ovl=max(0,min(int(a[1]),int(b[1]))-max(int(a[0]),int(b[0])))#calcula o overlap
    #print("olv:",ovl)
    if ovl != 0:
        return True
    return False


# ### Main functions:

# #### With Centromers:

# In[7]:

def closest_to_CNV(df_raw_regions,dict_Centro):
    """Associates each CNV_ID (key) to a region
        Requieres: a dataframe that is the CNV and the dict is a region (dict_Centro)
        Ensures: retruns a dict that associates each CNV_ID (key) to a region
    """
    #filter chrom column:
    #Start by each CNV of df_CNVs:
    
    list_items = []
    
    #for CNV_id in df_CNVs.ID:
    for index, row in df_raw_regions.iterrows():
            #sub_df_CNV = df_CNVs.loc[df_CNVs['ID'] == CNV_id]
            #sub_df_CNV = row # regions A and B  per line
            #print("ID>516")
        sub_df_region_A = row[['regionA_stat','regionA_end','chr_A']]
        sub_df_region_B = row[['regionB_stat','regionB_end','chr_B']]
            #
        regionA_data = assigned_region_data(sub_df_region_A,dict_Centro)#left
        regionB_data = assigned_region_data_B(sub_df_region_B,dict_Centro)#right #B!
        
        ##row[['Case_id','ID',''Cluster_id','star_position_max','end_position_max']]
        #print("_______")
        #print(row.iloc[0:5].to_list())#DEBUG
        #row[['regionA_stat','regionA_end','chr_A']]
        #print(sub_df_region_A.to_list())# DEBUG
        #
        #print(regionA_data[['region_name','region_family','region_class','region_strand','r_chrom']].to_list())# DEBUG
        #print(regionA_data[2:9])# DEBUG
        #
        #print(sub_df_region_B.to_list())#DEBUG
        #
        #print(regionB_data[['region_name','region_family','region_class','region_strand','r_chrom']].to_list())#DEBUG
        #print(regionB_data[2:9])# DEBUG
        #print('____ \n \n')
        #
        data_regions = row.iloc[0:5].to_list() + sub_df_region_A.to_list() + regionA_data[2:9] + sub_df_region_B.to_list() + regionB_data[2:9]
        #
        list_items.append(data_regions)
        #print("DR:",data_regions)

    return list_items 

def assigned_region_data(sub_df_CNV,dict_Centro):
    """Finds region name, family and class
        sub_df_CNV is a region!!
    """
    selected_chrom = sub_df_CNV.iloc[-1]
    #selected_chrom = sub_df_CNV['chr'].iloc[0]
    sub_df_Centro = pd.DataFrame.from_dict(dict_Centro['chr' + str(selected_chrom)])
    sub_df_Centro.columns=['chromStart','chromEnd','region_name','region_family','region_class','region_strand','r_chrom'] #possible mistake in the columns order of 'region_name','region_family','region_class'. 
    #Remove segmental duplications:
    sub_df_Centro = sub_df_Centro.loc[sub_df_Centro['region_name'] != 'Dup']
    dist=999**99
    position_item = None
        
        
    for i in range(len(sub_df_Centro)):
        a = [sub_df_Centro['chromStart'].iloc[i].item(),sub_df_Centro['chromEnd'].iloc[i].item()]
        b = [sub_df_CNV.iloc[0] , sub_df_CNV.iloc[1]]#start,end #!!! .item()
        #print(a,'\n',b,'\n')#DEBUG
        if get_Overlap(a,b) == True:
            dist = 0
            position_item = i
            break
        elif get_Overlap(a,b) == False and dist != 0:
            posA = (sub_df_CNV.iloc[0] - sub_df_Centro['chromEnd'].iloc[i]).item() # start_CNV_region - repRegion_end 
            #posB = (sub_df_Centro['chromStart'].iloc[i] - sub_df_CNV.iloc[1]).item() # repRegion_start - end_CNV_region
            
            if posA >= 0:
                new_dist = posA
            
            #elif posB >= 0: 
            #    new_dist = posB
                
            #elif posA >= 0 & posB >= 0:
            #    new_dist = min([posA, posB])
            
                if new_dist < dist:
                    dist = new_dist
                    position_item = i
            else:
                break
                #break here!!!
    if position_item != None:
        selecteted_item = sub_df_Centro.iloc[position_item]
        add_list=selecteted_item.tolist()
    else:
        add_list = [np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan]
            
    dist = float(dist) / 1000000 # to Mb
    add_list.append( np.log10(dist +1) ) #log10 +1 aplied!
    return add_list

#____________________for_Region_B_______________________________

def assigned_region_data_B(sub_df_CNV,dict_Centro): #check it if it is OKKAY!!!!
    """Finds region name, family and class
        sub_df_CNV is a region!!
    """
    selected_chrom = sub_df_CNV.iloc[-1]
    #selected_chrom = sub_df_CNV['chr'].iloc[0]
    sub_df_Centro = pd.DataFrame.from_dict(dict_Centro['chr' + str(selected_chrom)])
    sub_df_Centro.columns=['chromStart','chromEnd','region_name','region_family','region_class','region_strand','r_chrom'] #possible mistake in the columns order of 'region_name','region_family','region_class'. 
    #Remove segmental duplications:
    sub_df_Centro = sub_df_Centro.loc[sub_df_Centro['region_name'] != 'Dup']
    #reverce sub_df_Centro order since it is interesested the rep region that's on the right of the Region_B:
    sub_df_Centro = sub_df_Centro.iloc[::-1] #B!
    #
    dist=999**99
    position_item = None
    #
    #sub_df_Centro
    #
    for i in range(len(sub_df_Centro)):
        a = [sub_df_Centro['chromStart'].iloc[i].item(),sub_df_Centro['chromEnd'].iloc[i].item()]
        b = [sub_df_CNV.iloc[0] , sub_df_CNV.iloc[1]]#start,end #!!! .item() #B?
        #print(a,'\n',b,'\n')#DEBUG
        if get_Overlap(b,a) == True: #B!
            dist = 0
            position_item = i
            break
        elif get_Overlap(b,a) == False and dist != 0: #B!
            #posA = (sub_df_CNV.iloc[0] - sub_df_Centro['chromEnd'].iloc[i]).item() # start_CNV_region - repRegion_end  #B?
            posB = (sub_df_Centro['chromStart'].iloc[i] - sub_df_CNV.iloc[1]).item() # repRegion_start - end_CNV_region #B?
            
            if posB >= 0: #B!
                new_dist = posB #B!
            
            #elif posA >= 0: #B!
            #    new_dist = posA #B!
                
            #elif posA >= 0 & posB >= 0:
            #    new_dist = min([posA, posB])
            #    print("DEBUG")#DEBUG
            
                if new_dist < dist:
                    dist = new_dist
                    position_item = i
            else:
                break
                #break here!!!
    if position_item != None:
        selecteted_item = sub_df_Centro.iloc[position_item]
        add_list=selecteted_item.tolist()
    else:
        add_list = [np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan]
                        
    dist = float(dist) / 1000000 # to Mb
    add_list.append( np.log10(dist +1) ) #log10 +1 aplied!
    return add_list
#__________________________________________________________________________________
# ## Load CNVs and start: 

# In[8]:
# In[9]:
# In[10]:

PATH_CNVs = argv[1] 

#sheetName = argv[2]

#df_CNVs = pd.read_excel(PATH_CNVs, header=3, sheet_name=sheetName)
#df_CNVs['chr'] = df_CNVs['chr'].transform(filter_chrName)
df_CNVs = pd.read_csv(PATH_CNVs, sep= ";" ,header=0)
#df_CNVs = df_vars_TP 

import json
with open("./dataset3/rep_regions/dit_RepDups.json", "r") as outfile:
    dict_regions = json.load(outfile)
#Labels: #bin	chrom	chromStart	chromEnd	name

regions = closest_to_CNV(df_CNVs,dict_regions)

# In[11]:
# In[12]:

#df_selec_centromers
#df_selec_regions = pd.DataFrame.from_dict(regions, orient='index').reset_index()
#df_selec_regions.columns=['CNV_id','chromStart','chromEnd','region_name','region_family','region_class','region_strand','r_chrom','region_dist']
df_selec_regions = pd.DataFrame(regions)
#print(df_selec_regions)
df_selec_regions.columns=['Case_id', 'ID', 'Cluster_id', 'star_position_max', 'end_position_max','regionA_stat','regionA_end','chr_A','A_region_name','A_region_family','A_region_class','A_region_strand', 'A_r_chrom','A_r_dist','regionB_stat','regionB_end','chr_B','B_region_name','B_region_family','B_region_class','B_region_strand', 'B_r_chrom','B_r_dist']


#Save data in .csv
outPATH = argv[2]
#outPATH = "'./dataset3/"

#df_selec_centromers.to_csv(outPATH + '.csv', sep = ';' ,index_label=False) #main
df_selec_regions.to_csv(outPATH + 'info_regions_matched_test.csv', sep = ';' ,index_label=False)
