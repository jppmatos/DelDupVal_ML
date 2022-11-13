#!/usr/bin/env python
# coding: utf-8
from sys import argv
import pandas as pd
#import matplotlib.pyplot as plt
#import seaborn as sns
import numpy as np
from scipy.stats import norm
#from sklearn.preprocessing import StandardScaler
from scipy import stats
#import warnings
#warnings.filterwarnings('ignore')
import pyfaidx

## Based of the DUST algorithm.
## check: https://www.ncbi.nlm.nih.gov/IEB/ToolBox/CPP_DOC/lxr/source/src/app/dustmasker/


# Axiliar Functions:
def n_masked(string):
    """Execpt 'n' char (masked)(no sequnece info, likey a repeat?), and soft-masked elements (lower case letters)? """
    return sum(1 for c in string if c.islower() and c != 'n')


#Function:
def dust(trinuc_motif_counts_in_pixel,nonN_pixel_width):
    """Classifies the seleted region with a score that higher represents a region/sequence complexity. \n (Sum of squares of trinuc. Motif counts in pixel, divided by square of nn-N pixel width)
        Requires:trinuc_motif_counts_in_pixel, total of bases/pixels excpt the masked regions ('n'); \n nonN_pixel_width, total leght of the sequence - total of 'masks'.
        Ensures: return a score than represents the sequence complexity, higher the score lower the seq complexity.
    """
    low_complexity_val = (trinuc_motif_counts_in_pixel/3)**2 / (nonN_pixel_width)**2
    return low_complexity_val


# path = './dataset3/Sequence_Complexity/'
# genes = pyfaidx.Fasta(path + 'Homo_sapiens.GRCh38.dna_sm.chromosome.1.fa')
# 
# to_analise = genes['1'][9980:13000].seq
# 
# print(to_analise.count('n'))

def collect_dust(star_pos,end_pos,chr_id,Path):
    """In each region
    """
    chr_id = str(chr_id)
    star_pos = round(star_pos)
    end_pos = round(end_pos)
    #Path = './dataset3/Sequence_Complexity/'
    genes = pyfaidx.Fasta(Path + 'Homo_sapiens.GRCh38.dna_sm.chromosome.' + str(chr_id) +'.fa')
    to_analise = genes[chr_id][star_pos:end_pos +1].seq
    #
    trinuc_motif_counts_in_pixel = n_masked(to_analise) # = total of bases/pixels excpt the masked regions ('n')
    nonN_pixel_width = len(to_analise) - to_analise.count('n') # = total lenght - total of masked out low complexity parts of a genome
    #
    return dust(trinuc_motif_counts_in_pixel,nonN_pixel_width)
    


# In[6]:


def Seq_Complexity(df_regions,Path):
    """
    """
    list_Reg_Dust = []
    for index, row in df_regions.iterrows():
        #A
        dust_A = collect_dust(row['regionA_stat'],row['regionA_end'],row['chr_A'],Path)
        #B
        dust_B = collect_dust(row['regionB_stat'],row['regionB_end'],row['chr_B'],Path)
        #CREATE NEW COLS FOR Dust_A and Dust_B
        # New df: ID Cluster_id Dust_A  Dust_B
        
        list_Reg_Dust.append([row['ID'],row['Cluster_id'],dust_A,dust_B])
    
    return list_Reg_Dust
     


# # Load Data

# In[7]:

region_path = argv[1]
df_regions = pd.read_csv(region_path, sep= ";" ,header=0,index_col=0)
#df_regions = pd.read_csv('./dataset3/rep_regions/match_CNVs_regions_TP.csv', sep= ";" ,header=0,index_col=0)
#
#path_GC_content = './dataset3/GCcontent/'
Path_Seq_Comp = './dataset3/Sequence_Complexity/'

data_Seq_Complexity = Seq_Complexity(df_regions,Path_Seq_Comp)
df_Seq_Complexity = pd.DataFrame(data_Seq_Complexity, columns=['CNV_ID','Cluster_ID','Dust_A','Dust_B'])

#Save to a csv:
output_path = argv[2]
df_Seq_Complexity.to_csv(output_path + 'Regs_Seq_Complexity.csv' , sep=';',index_label=False)
# In[ ]:
#df_Seq_Complexity.loc[df_Seq_Complexity['CNV_ID']=='DGRC_DEL_1281']

