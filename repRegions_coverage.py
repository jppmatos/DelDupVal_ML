from sys import argv
import json
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
from scipy.stats import norm
from sklearn.preprocessing import StandardScaler
from scipy import stats
import warnings
warnings.filterwarnings('ignore')

## Axiliar Functions:
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
# In[]:
def get_Overlap(a,b):
    """
        Requires: a= region start and end (list), b=  chromStart and chromEnd from rR (list)
    """
    ovl=max(0,min(a[1],b[1])-max(a[0],b[0]))#calcula o overlap
    return ovl

'''
def coverage(region_stat,region_end,df_exons): #def coverage(region_stat,region_end,gene_data):
    """ 
    """
    cov_percent_list_EXONS = 0

    #
    block_start_exon = region_stat#!!!
    for index, row in df_exons.iterrows():
        # exon, intron ,exon , ...

        #exon:
        block_start_exon = row["start"]
        block_end_exon = row["end"]
        #
        overlap  = get_Overlap([region_stat,region_end],[block_start_exon,block_end_exon])
        #print("EXON",overlap)#DEBUG
        #save overlap exon!!!
        cov_percent_list_EXONS += overlap

    return cov_percent_list_EXONS, 
'''       
#...

def get_segmentalDUPS(region_stat,region_end,chrom,dict_Centro):#,df_genes):
    """
        Finds overlaped reptitive regions with the selected region(A or B) then gets the coverage of the overlaps. \
        INFO: As the same script used for segmental duplications is usable to reptitive regions,\
        the variable with names as **segDUPS** will refere to *reptitive regions data* !!!
    """
    result = []
    a=[region_stat,region_end]#region
    #df_genes = df_genes.loc[df_genes["chrom"] == str(chrom)] #replaced
    sub_df_segDUPS = pd.DataFrame.from_dict(dict_Centro['chr' + str(chrom)])
    #sub_df_segDUPS.columns=['start','end','region_name','region_family','region_class','region_strand','r_chrom']
    sub_df_segDUPS.columns=['start','end','region_class','region_name','region_family','region_strand','r_chrom']
    #Just segmental duplications:
    sub_df_segDUPS = sub_df_segDUPS.loc[sub_df_segDUPS['region_class'] != 'Dup'] #FOR ONLY REPETITIVE REGIONS !!!
    #---    
    #Line by line form sub_df_segDUPS (Gene by gene):
    for index, row in sub_df_segDUPS.iterrows():
        b = [row['start'],row['end']]#segmental duplication
 #------------------------------------------------------------------------------------
        overlap = get_Overlap([region_stat,region_end],b)
        cov_segDUPS = 0
        if overlap != 0:
            cov_segDUPS += overlap
            #print("OVERLAP")#DEBUG
            #print(row['region_name'])
            #add info:
            info = (row['region_name'],row['region_family'],row['region_class'])
            #print("\n _________________________________________") #DEBUG
            #print([tuple(b),cov_segDUPS,info])#DEBUG
            result += [[tuple(b),cov_segDUPS,info]]
            
    return result
#
def deal_multi_segDUPS(rR_segDUPS, region_size,region_start ,region_end ,count=0):
    """check a pair os segmental duplications, refered as intervals, does overlap and unite then to calqulate the coverage.
    
    """
    #print("III",rR_segDUPS)#DEBUG
    rR_intervals = [ i[0] for i in rR_segDUPS ] #adpator for merge_intervals
    rR_segDUPS = merge_intervals(rR_intervals)#old -,region_start ,region_end)
    coverage = 0
    #print("xxx", rR_segDUPS)#DEBUG
    for segDup in rR_segDUPS:
        b = [ segDup[0], segDup[1] ] #  b = [ segDup[0][0], segDup[0][1] ]
        overlap = get_Overlap([region_start,region_end],b)
        coverage += overlap
    #print(coverage)
    return coverage


def merge_intervals(rR_intervals):#,region_start ,region_end):
	#rR_intervals = [ i[0] for i in pre_rR_intervals ]
	new_inters=[]
	count=0
	new_row_a=[]
	start_b = 0
	while count <= len(rR_intervals)-1:
		if new_row_a == []:
			new_row_a = list(rR_intervals[count])
		else:
			start_b, end_b= rR_intervals[count][0:3]
		if new_row_a[0] != 0 and start_b != 0:
			if get_Overlap((new_row_a[0], new_row_a[1]),(start_b,end_b)) != 0:
				new_row_a = (min(new_row_a[0],start_b),max(new_row_a[1],end_b))#,new_row_a[2])
			elif new_row_a[1] == start_b:
				new_row_a=(new_row_a[0],end_b)#,new_row_a[2])
			else:
				new_inters.append(new_row_a)
				new_row_a = [start_b,end_b]
		count += 1
	if new_row_a != []:
		new_inters.append((new_row_a[0],new_row_a[1]))#,new_row_a[2]))
	#print(new_inters)#DEBUG
	return new_inters
#        
def define_interval(rR_segDUPS, region_end,region_start):
    """trash???
    """
    #if count +1 = len(rR_segDUPS):
    #    return new
    #    
    #start_seg, end_seg = rR_segDUPS[count][0] # tuple
    #start_new, end_new = (0,0)
    #
    #if region_start <= start_seg:
          

#
def merge(rR_segDUPS, region_end, region_stat):
    """
        ex:(rR_A,(row['regionA_end'] - row['regionA_stat']))
        ex:([[(start,end),cov_segDUPS],[(start,end),cov_segDUPS]],(row['regionA_end'] - row['regionA_stat']))
    """
    region_size = region_end - region_stat
    cov_segDUPS = 0
    cov = 0
    #print(rR_segDUPS)#DEBUG

    if rR_segDUPS != []:
        if len(rR_segDUPS) > 1:
            #print("||",rR_segDUPS)#DEBUG
            cov = deal_multi_segDUPS(rR_segDUPS, region_size, region_stat, region_end)# when two or more segDuplication overlap it self!
            #print(cov)#DEBUG
            cov_segDUPS += cov
        else:
            #print("|",rR_segDUPS)#DEBUG
            cov= rR_segDUPS[0][1]
            #print(cov)#DEBUG
            cov_segDUPS += cov
    #____________________________
    cov_segDUPS /= region_size

    return cov_segDUPS


def find_info_maxcov(items):
    """Gets the name, family and class of the reptitiveRegion/segmentalDuplication of the one with the biggest overlap.
        Requires:list of lists of all overlaped items, like: [[(125384121, 125385230), 807], ('L1MEf', 'L1', 'Rep_LINE'), ...] .
        Ensures: return the list with the reptitiveRegion/segmentalDuplication with the biggest overlap. 
    """
    #print(items)#DEBUG
    if items != []:
        cov_list = sorted(items,reverse=True,key=lambda x: x[1])
        return cov_list[0]
        #return max(items,key=lambda x:x[1]) #ERROR: TypeError: '>' not supported between instances of 'str' and 'int'
    return items

def get_info_maxcov(items):
    """ return the name, family and class of the of reptitiveRegion/segmentalDuplication with the biggest overlap"""
    name_, family_, class_ = "no_overlap","no_overlap","no_overlap"
    found_item = find_info_maxcov(items)
    if found_item != []:
        name_, family_, class_ = found_item[2]
    return name_, family_, class_ 

# In[]:
def main_segDuplications(df_regions,dict_segDUPS):
    """
    """
    list_segDUPS= []
    #Line by line from the A nad B regions:
    for index, row in df_regions.iterrows():
        #rR_A = ['','','']
        #rR_B = []
        #A
        #print("|----A----|\n")#DEBUG
        rR_A = get_segmentalDUPS(row['regionA_stat'],row['regionA_end'],row['chr_A'],dict_segDUPS)
        #Find best overlaped segDup/repRegion: 
        A_covMax_A_name, A_covMax_family, A_covMax_class = get_info_maxcov(rR_A)
        #merge data from the genes:
        #print(rR_A,'\n')#DEBUG
        rR_A = merge(rR_A,row['regionA_end'], row['regionA_stat'])
        #
        #print(rR_A,'\n')#DEBUG
        #B
        #print("|----B----|")#DEBUG
        rR_B = get_segmentalDUPS(row['regionB_stat'],row['regionB_end'],row['chr_B'],dict_segDUPS)
        #Find best overlaped segDup/repRegion:
        B_covMax_name, B_covMax_family, B_covMax_class = get_info_maxcov(rR_B)
        #print(rR_B,'\n')#DEBUG
        #merge data from the genes:
        rR_B = merge(rR_B,row['regionB_end'], row['regionB_stat'])
        #
        #print(rR_B,'\n')#DEBUG
        #Prep info by line:
        list_add = [row['ID'],row['Cluster_id']] + [rR_A] + [rR_B] + [A_covMax_A_name] + [A_covMax_family] + [A_covMax_class] + [B_covMax_name] + [B_covMax_family] + [B_covMax_class]
        list_segDUPS.append(list_add)
    return list_segDUPS

# # Load Data

region_path = argv[1]
df_regions = pd.read_csv(region_path, sep= ";" ,header=0,index_col=0)
#df_regions = pd.read_csv('./dataset3/rep_regions/match_CNVs_regions_TP.csv', sep= ";" ,header=0,index_col=0)

df_regions['chr_A'] = df_regions['chr_A'].astype(str)
df_regions['chr_B'] = df_regions['chr_B'].astype(str)
#
with open("./dataset3/rep_regions/dit_RepDups.json", "r") as outfile:
    dict_segmentalduplications = json.load(outfile)


path_Genes = './dataset3/DUPs_regions/'

Genes_Data = main_segDuplications(df_regions,dict_segmentalduplications)
# In[16]:
#To a dataframe:
df_Genes = pd.DataFrame(Genes_Data, columns=['CNV_ID','Cluster_ID',
                                       'A_cov_repRegion','B_cov_repRegion',
                                             'Name_A','Family_A','Class_A',
                                                 'Name_B','Family_B','Class_B'])

#Save to a csv:
output_path = argv[2]
df_Genes.to_csv(output_path + 'repRegions_coverage_regionsAB.csv', sep=';',index_label=False)
#df_Genes.to_csv('./dataset3/Genes/'+'Genes_regionsAB.csv', sep=';',index_label=False)

