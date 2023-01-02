#!/usr/bin/env python
# coding: utf-8

# In[1]:


import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
from scipy.stats import norm
from sklearn.preprocessing import StandardScaler
from scipy import stats
import warnings
warnings.filterwarnings('ignore')
#get_ipython().run_line_magic('matplotlib', 'inline')


# ## Load data:

# In[2]:


#load_folder
load_folder = './results_CSV/'


# In[3]:


#df_results_open = pd.read_csv('result2.csv', sep= ";" ,header=0,index_col=0)
#df_results_open


# ## CNVs:

# In[4]:

df_vars_TP = pd.read_excel('./dataset2/SVs_V2_5_2020.xlsx', header=3, sheet_name='Sheet_DEL_DUP')
#remove "*_50overlap" when not for 50 %overlap!
#___
#df_50_TP = pd.read_csv('list_tp_50over',names=['ID'])
#df_vars_TP = df_vars_TP.loc[df_vars_TP['ID'].isin(df_50_TP['ID'].tolist())]
#__________________________________________________________________________

df_vars_FP = pd.read_excel('./dataset2/fp_01_2020.xlsx', header=3, sheet_name='Folha3')
#remove "*_50overlap" when not for 50 %overlap!
#___
#df_50_FP = pd.read_csv('list_tp_50over',names=['ID'])
#df_vars_FP = df_vars_FP.loc[df_vars_FP['ID'].isin(df_50_FP['ID'].tolist())]
#__________________________________________________________________________

# In[5]:


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
    #print(chromossome)
    return str(chromossome)


# In[6]:


#df_vars_TP['chr'].transform(filter_chrName)


# In[7]:


df_vars_TP['chr'] = df_vars_TP['chr'].transform(filter_chrName) #keep just the chromossome identifier number
#display(df_vars_TP)
df_vars_FP['chr'] = df_vars_FP['chr'].transform(filter_chrName)
#display(df_vars_FP)


# ## Start:

# In[8]:


def select_clustered_CNVs(dataframe,CNV_type):
    dataframe['Type'] = dataframe['Type'].str.lower()
    df_SplitRead = dataframe.loc[dataframe['Type'] == CNV_type]
    
    #df_SplitRead['Cluster'] = df_SplitRead['Cluster'].str.lower()
    df_SplitRead = df_SplitRead.loc[(df_SplitRead['Cluster'] != 'No') & (df_SplitRead['Cluster'] != 'no')]
    return df_SplitRead


# In[9]:


##Other solution 

def get_Overlap(a,b):
    ovl=max(0,min(a[1],b[1])-max(a[0],b[0]))#calcula o overlap
    sz_a=a[1]-a[0] #calcula o tamanho da região definida pelo agrupamento
    sz_b=b[1]-b[0] #calcula o tamanho da deleção
    
    if sz_a <= 0 or sz_a <= 0:
        return False
    
    ta=ovl/sz_a # calcula a proporção da região definida pelo cluster que se sobrepõe à deleção
    tb=ovl/sz_b # calcula a proporção da deleção que se sobrepõe à região definida pelo cluster
    if ta>=0.7 and tb>=0.7:#verifica se ambas as proporções são superiores a 90%
            return True
    else:
            return False
        
def match_CNVs_regions(clusters_dict,df_SplitRead,list_CNVsIDS):
    """
    """
    #display(df_SplitRead)#BEBUG
    clusters = []
    for c in clusters_dict.keys():
        #DEBUG
        #print("Cluster_ID: ", c,'\n',clusters_dict[c])
        #\DEBUG
        
        df_r= df_SplitRead.loc[((df_SplitRead['chr']) == str(clusters_dict[c][6]))
                            & ((df_SplitRead['chr']) == str(clusters_dict[c][7]))]
        #print(df_r)#DEBUG ______________________________________________________________________
        #clusters_dict
        Max_Start_pos = clusters_dict[c][2]#[0] #!!!!
        Min_end_pos = clusters_dict[c][5]#[1] #!!!!
        a = [Max_Start_pos, Min_end_pos]#limites do agrupamento que identificam a deleção
    
    
        #df_r
        for index, row in df_r.iterrows(): #CNV
            if row['ID'] not in list_CNVsIDS:
                #df_r
                #ID      Type chr      start        end                       Region  ... Coverage Array Sanger Reference Observations Talkowski cases
                #if row['ID'] == 'DGRC_DEL_255': print(row['ID'])#DEBUG !!!!
                #DEBUG

                Start_del = row['start']
                end_del = row['end']
                b=[Start_del, end_del]#limites da deleçao    
                #print(row['ID'], a,b,'\n')#DEBUG

                #DEBUG
                #if Max_Start_pos <Min_end_pos is False:
                #    print('\n error!!!! \n', c ,'\n', index, ' - ',row)

                #if a[1]-a[0] <= 0: #THE CLUSTER HAS REGIONS OVERFIT!!! 
                #    print("Cluster: ",c,"|",Max_Start_pos,Min_end_pos)
                #\DEBUG

                if get_Overlap(a,b) is True:
                    #print('TRUE',row.iloc[0]) #DEBUG--------------------------------------------------------------------------------------------------------------
                    #clusters.append([c,row['ID']]) #c is the cluster id and ID is the CNV id.
                    to_add = [row['ID'],c]
                    to_add.extend(clusters_dict[c]) 
                    #clusters.append([row['ID'],c,clusters_dict[c]]) 
                    clusters.append(to_add)
                    #
                    list_CNVsIDS.append(row['ID'])
                    #print(list_CNVsIDS)

                    #if row['ID'] == 'DGRC_DEL_255': #DEBUG
                    #print(['ID', '____>>__>>>_>_>'])#DEBUG

            #drop row:
            #df = df.drop()    
                    break #!!!
        #
    return clusters,list_CNVsIDS 


# In[10]:
# In[11]:
def study_cases_main(case_list,df_CNVs,CNV_type):
    """
    """
    #
    df_CNVs['Type'] = df_CNVs['Type'].str.lower() #!!!
    #display(df_CNVs)
    #
    matrix_list_clusters = []
    list_CNVsIDS = []
    #Path to the case_clusters folders:
    repo3 = 'dataset3/'
    rep_regions_path = "rep_regions/"
    
    df_SplitRead = select_clustered_CNVs(df_CNVs,CNV_type) #reduced CNV dataframe with only the CNV witn the <"Cluster" = Yes>
    #print(df_SplitRead)#DEBUG
    #df_SplitRead = df_CNVs
    for selected_case in case_list:
        ##DEBUG
        #print('\n | %s \n' % selected_case)
        #if selected_case == 'DGRC0027':
        #    print("\n selected_case: \n",selected_case)
        #\DEBUG
        
        #case_study = df_SplitRead.loc[df_SplitRead['Cases with the alteration'].str.contains(selected_case)] ## WHY I USE IT OR MADE IT??!??! DON'T REMEMBER...
        #----------------------------------------------1
        CNV_type.lower()
        if CNV_type == 'deletion':
            CNV_type = 'del'
        elif CNV_type == 'duplication':
            CNV_type = 'dup'
            
        case_path = repo3 + rep_regions_path + selected_case + '/' + CNV_type + '.clusters.txt'
        #df_case = case_path
        df_case= pd.read_csv(case_path, sep= "\t" ,
                                names=['id','2','3','chr_A','chr_B','star_position','end_position','1_readpair_quality','2_readpair_quality','1_readpair_size','2_readpair_size','SAMflag_1_readpair','SAMflag_2_readpair','14','15','16','17','18'],
                                header=None)
        #----------------------------------------------2
        list_regions = df_case['id'].unique()
        
        
        clusters_dict = {}
        
        for i in list_regions:
            #print(i, 'item_regions')#DEBUG
            
            df_select_region = df_case.loc[df_case['id'] == i] #df for each cluster
            #print(df_select_region, 'df_select_region')#DEBUG
            #asdsdas
            star_position_max = df_select_region['star_position'].max() #star dos Dels  #ex: 106130914 
            end_position_min = df_select_region['end_position'].min() #end dos Dels
    
            #region A 
            regionA_stat= df_select_region['star_position'].min() #df_case['star_position'].min()   
            regionA_end = star_position_max #df_case['star_position'].max()

            #region B
            regionB_stat =end_position_min #df_case['end_position'].min()
            regionB_end = df_select_region['end_position'].max() #df_case['end_position'].max()

            cluster_n = i
    
            clusters_dict[i] = [star_position_max,end_position_min,
                                 regionA_stat,regionA_end,
                                 regionB_stat,regionB_end,
                                df_select_region['chr_A'].iloc[0], df_select_region['chr_B'].iloc[0]]
        #----------------------------------------------
        #list_CNVsIDS of the already matched CNVs:
        list_clusters, list_CNVsIDS = match_CNVs_regions(clusters_dict,df_SplitRead,list_CNVsIDS)
        #print(list_CNVsIDS)#DEBUG
        #----------------------------------------------
        #print(clusters_dict)
        #print(list_clusters) #[['DGRC_DEL_31', 1139, 49686901, 49736702, 49685715, 49686901, 49736702, 49738215, '11', '11'], [...], ... ]
        matrix_list_clusters.append([selected_case,list_clusters])
        
    return matrix_list_clusters
        
        
        


# ### Start of aplication:

# # With TPs:

# In[12]:

#Case selection TPs:
#df_vars_TP = pd.read_excel('./dataset2/SVs_V2_5_2020.xlsx', header=3, sheet_name='Sheet_DEL_DUP')
#case_list = ['DGRC0017','DGRC0015','DGRC0014','DGRC0023','DGRC0013','DGRC0024'] 
#case_list = ['DGRC0013']
case_list = ['DGRC0017','DGRC0015','DGRC0014','DGRC0023','DGRC0013','DGRC0024','DGRC0022','DGRC0027','DGRC0037','DGRC0050']
CNV_type = 'duplication'#'deletion'#'duplication'
df_CNVs = df_vars_TP 
#df_CNVs = df_vars_FP


match_CNVs_regionsC = study_cases_main(case_list,df_CNVs,CNV_type)


# In[13]:


match_CNVs_regionsC = [(x[0],y) for x in match_CNVs_regionsC for y in x[1] ] #reorganize, it's stupid but works


# In[14]:


#it's stupid but works
match_CNVs_regionsC_2 = []

for i in match_CNVs_regionsC: #('DGRC0017', ['DGRC_DEL_31', 1139, 49686901, 49736702, 49685715, 49686901, 49736702, 49738215, '11', '11']),
    sub_list = []
    #print(i[0],i[1])
    sub_list = [i[0]]
    sub_list.extend(i[1])
    match_CNVs_regionsC_2.append(sub_list)
    
#print(match_CNVs_regionsC_2)

df_match_CNVs_regionsC_TP = pd.DataFrame.from_records(match_CNVs_regionsC_2,columns=['Case_id','ID','Cluster_id','star_position_max','end_position_max','regionA_stat','regionA_end','regionB_stat','regionB_end','chr_A','chr_B'])
#df_match_CNVs_regionsC
#CNV_id = ID


# In[15]:


#df_match_CNVs_regionsC_TP


# In[16]:


#SAVE DATAFRAME df_match_CNVs_regionsC_TP ...
#revmove "*_50overlap" when not for 50 %overlap!
df_match_CNVs_regionsC_TP.to_csv('./dataset6/rep_regions/match_CNVs_regions_TP_dup.csv', sep=';',index_label=False)


# In[17]:


#


# df_match_CNVs_regionsC.loc[df_match_CNVs_regionsC_TP['ID'] == 'DGRC_DEL_1140']

# df_match_CNVs_regionsC.loc[df_match_CNVs_regionsC['Cluster_id'] == 1694]

# In[18]:


df_CNVs.merge(df_match_CNVs_regionsC_TP, on="ID") #REPLACE MERGE WITH JOINTOR SOMETHING SIMILAR !!!!!


# # With FPs:

# In[19]:

#--------------------------------------------------------------------------------------------------------------------------------
#Case selection FPs:
'''
#df_vars_FP = pd.read_excel('./dataset2/fp_01_2020.xlsx', header=3, sheet_name='Folha3')
case_list = ['DGRC0017','DGRC0015','DGRC0014','DGRC0023','DGRC0013','DGRC0024','DGRC0022','DGRC0027','DGRC0037','DGRC0050']
CNV_type = 'deletion'
df_CNVs = df_vars_FP 
#df_CNVs = df_vars_FP


match_CNVs_regionsC = study_cases_main(case_list,df_CNVs,CNV_type)


# In[20]:


match_CNVs_regionsC = [(x[0],y) for x in match_CNVs_regionsC for y in x[1] ] #reorganize, it's stupid but works


# In[21]:


#it's stupid but works
match_CNVs_regionsC_2 = []

for i in match_CNVs_regionsC: #('DGRC0017', ['DGRC_DEL_31', 1139, 49686901, 49736702, 49685715, 49686901, 49736702, 49738215, '11', '11']),
    sub_list = []
    #print(i[0],i[1])
    sub_list = [i[0]]
    sub_list.extend(i[1])
    match_CNVs_regionsC_2.append(sub_list)
    
#print(match_CNVs_regionsC_2)

df_match_CNVs_regionsC_FP = pd.DataFrame.from_records(match_CNVs_regionsC_2,columns=['Case_id','ID','Cluster_id','star_position_max','end_position_max','regionA_stat','regionA_end','regionB_stat','regionB_end','chr_A','chr_B'])
#df_match_CNVs_regionsC
#CNV_id = ID


# In[22]:


df_match_CNVs_regionsC_FP


# In[23]:


#SAVE DATAFRAME df_match_CNVs_regionsC_FP ...
df_match_CNVs_regionsC_FP.to_csv('./dataset4/rep_regions/df_match_CNVs_regions_FP.csv', sep=';',index_label=False)

# In[24]:


#


# In[25]:


df_CNVs.merge(df_match_CNVs_regionsC_FP, on="ID") #REPLACE MERGE WITH JOINTOR SOMETHING SIMILAR !!!!!
'''

#----------------------------------------------------------------------------------------------------------------------

# In[ ]:





# # With FP and invalid cluster arent removed:

# In[26]:


# With FP and invalid cluster arent removed:
'''
def get_Overlap(a,b): #unfiltered
    ovl=max(0,min(a[1],b[1])-max(a[0],b[0]))#calcula o overlap
    sz_a=a[1]-a[0] #calcula o tamanho da região definida pelo agrupamento
    sz_b=b[1]-b[0] #calcula o tamanho da deleção
    
    if sz_a <= 0 or sz_a <= 0:
        return False
    
    ta=ovl/sz_a # calcula a proporção da região definida pelo cluster que se sobrepõe à deleção
    tb=ovl/sz_b # calcula a proporção da deleção que se sobrepõe à região definida pelo cluster
    if ta>=0.9 and tb>=0.9:#verifica se ambas as proporções são superiores a 90%
            return True
    else:
            return False
'''

# # With FPs:

# In[27]:


#Case selection FPs:
#case_list = ['DGRC0017','DGRC0015','DGRC0014','DGRC0023','DGRC0013','DGRC0024','DGRC0022','DGRC0027','DGRC0037','DGRC0050']
case_list = ['DGRC0017','DGRC0015','DGRC0014','DGRC0023','DGRC0013','DGRC0024'] #case_list = ['DGRC0017','DGRC0015']

CNV_type = 'duplication'
df_CNVs = df_vars_FP 

#print("\n ___________________________________----______________________________----__________________________________________ \n")

match_CNVs_regionsC = study_cases_main(case_list,df_CNVs,CNV_type)


# In[28]:


match_CNVs_regionsC = [(x[0],y) for x in match_CNVs_regionsC for y in x[1] ] #reorganize, it's stupid but works


# In[29]:


#it's stupid but works
match_CNVs_regionsC_2 = []

for i in match_CNVs_regionsC: #('DGRC0017', ['DGRC_DEL_31', 1139, 49686901, 49736702, 49685715, 49686901, 49736702, 49738215, '11', '11']),
    sub_list = []
    #print(i[0],i[1])
    sub_list = [i[0]]
    sub_list.extend(i[1])
    match_CNVs_regionsC_2.append(sub_list)
    
#print(match_CNVs_regionsC_2)

df_match_CNVs_regionsC_FP_unfiltered = pd.DataFrame.from_records(match_CNVs_regionsC_2,columns=['Case_id','ID','Cluster_id','star_position_max','end_position_max','regionA_stat','regionA_end','regionB_stat','regionB_end','chr_A','chr_B'])
#df_match_CNVs_regionsC
#CNV_id = ID


# In[30]:


#df_match_CNVs_regionsC_FP_unfiltered


# In[31]:


#SAVE DATAFRAME df_match_CNVs_regionsC_FP_unfiltered ...
#revmove "*_50overlap" when not for 50 %overlap!
df_match_CNVs_regionsC_FP_unfiltered.to_csv('./dataset6/rep_regions/df_match_CNVs_regions_FP_unfiltered_dup.csv', sep=';',index_label=False)


# In[32]:


#


# In[33]:


df_CNVs.merge(df_match_CNVs_regionsC_FP_unfiltered, on="ID") #REPLACE MERGE WITH JOINTOR SOMETHING SIMILAR !!!!!

