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

#---
def get_Overlap(a,b):
    """
        Requires: a= region start and end (list), b=  chromStart and chromEnd from rR (list)
    """
    ovl=max(0,min(a[1],b[1])-max(a[0],b[0]))#calcula o overlap
    return ovl

# # Functions:
#---new--- 
def match_overlaped_exons(df_exons):
    """ Match df_exons with df_intron, and find possible overlap between diff genes.
        If there is overlap: remove the portion between position of the overlap.
        Ensures: a list with (chrom, start, end) of the exons/introns without the genes overlap.
    """
    new_list_exons = []
    chr_list = df_exons['chrom'].unique().tolist()
    #chr_list = [1,2] # test
    for id_chr in chr_list:
        df_chr = df_exons.loc[df_exons['chrom']==id_chr]
        #sort by start position:
        df_chr.sort_values(by=['start'],inplace=True)
        #into a list
        pos_l_exons = [[row['start'],row['end'],row['chrom']] for index, row in df_chr.iterrows()]
        #merge them:
        new_exons_chr = merge_intervals(pos_l_exons)
        #print(new_exons_chr)#DEBUG
        new_list_exons.extend(new_exons_chr)
        
    return new_list_exons
#--/new---
#--new2---
def merge_genes_intervals(df_genes):
    """Merge diferent genes interval, exons are the limits, into one interval when overlaped.
    """
    list_genes = []
    #
    chr_list = df_genes['chrom'].unique().tolist()
    df_exons = df_genes.loc[(df_genes["type"] == 'gene')]  #df_exons is genes
    df_exons.drop_duplicates(inplace = True)
    #
    for chrom_id in chr_list:
        df_chr = df_exons.loc[(df_exons["chrom"] == chrom_id)]
        df_chr.sort_values(by=['start'],inplace=True)   
        genes_exons_chr_intervals = df_chr[['start','end','chrom','gene_name']].values.tolist()
        #:
        interval_chr = merge_intervals(genes_exons_chr_intervals) #it's a list of list, of each chr
        list_genes += interval_chr
    return list_genes # ex: [['start','end','chrom','gene_name'],...]

def merge_intervals(rR_intervals):#,region_start ,region_end):
	new_inters=[]
	count=0
	new_row_a=[]
	start_b = 0
	while count <= len(rR_intervals)-1:
		if new_row_a == []:
			new_row_a = list(rR_intervals[count][0:3])
		else:
			start_b, end_b, chr_b= rR_intervals[count][0:3]
		if new_row_a[0] != 0 and start_b != 0:
			if get_Overlap((new_row_a[0], new_row_a[1]),(start_b,end_b)) != 0:
				new_row_a = (min(new_row_a[0],start_b),max(new_row_a[1],end_b),new_row_a[2])
			elif new_row_a[1] == start_b:
				new_row_a=(new_row_a[0],end_b,new_row_a[2])
			else:
				new_inters.append(new_row_a)
				new_row_a = [start_b,end_b,chr_b]
		count += 1
	if new_row_a != []:
		new_inters.append((new_row_a[0],new_row_a[1],new_row_a[2]))
	#print(new_inters) #DEBUG
	return new_inters

#--/new2--
#--new3---
def find_introns(df_exons,df_genes_intervals):
    """ Find all introns by identifying the exons that belongs to a gene interval,
        and by the 'gaps' between exons identify the introns.
    """
    introns_list = []
    for index, row in df_genes_intervals.iterrows():
        gene_start = row['start']
        gene_end = row['end']
        gene_chr = row['chrom']
        df_selec_exons = df_exons.loc[
                            (df_exons['chrom']==gene_chr)
                                &(df_exons['start']>=gene_start)
                                &(df_exons['end']<=gene_end)]
        list_introns_chr = get_introns(df_selec_exons)
        if list_introns_chr != []:
            introns_list += list_introns_chr
        #DEBUG
        #print("start: %s , end: %s , chrom: %s " % (gene_start,gene_end,gene_chr))
        #print("exons:", df_selec_exons)
        #print(list_introns_chr ,'\n')
        #/DEBUG
    return introns_list #-> (['start'],['end'],['chrom']
    
#--/new3--
def get_introns(df_exons):
    """ Get introns from the interval between exons.
    """
    #select by gene:
    pos_l_introns = []
    df_gene = df_exons
    #use info of star and end of the first and last exon.
    pos_l_exons = [(row['start'],row['end'],row['chrom']) for index, row in df_gene.iterrows()] #[(1, 1), (1, 2), (2, 2)]
    for pos in range(len(pos_l_exons)-1):
        # save the gap between the end of the first exon to the star of the secound exon, in loop.
        start_intron = pos_l_exons[pos][1]
        end_intron = pos_l_exons[pos+1][0]
        pos_l_introns.append((start_intron,end_intron,pos_l_exons[pos][2]))
    #save data in df_introns:
    #df_introns = pd.DataFrame(pos_l_introns,columns =['start','end','chrom'])
    ##pos_l_introns -> (['start'],['end'],['chrom'])
    return pos_l_introns 

def prep_genes_data_main(df_genes):
    """ Prepareand get the necessary data from 'genes' data, which are exons and introns. 
    """
    #gene overlap matching:
    df_exons = df_genes.loc[(df_genes["type"] == 'exon')] # removes type genes
    df_exons.drop_duplicates(inplace = True) # remove duplicate rows
    #Not Fine VVVV
    list_exons = match_overlaped_exons(df_exons) # merge all overlaped exons
    df_exons = pd.DataFrame(list_exons,columns =['start','end','chrom'])
    #Fine VVVV-------------------------------------------------------------
    #gene overlap matching: ex: start_a,end_b;
    genes_intervals = merge_genes_intervals(df_genes) #reuse: merge_intervals()
    #print(genes_intervals)#DEBUG
    df_genes_intervals = pd.DataFrame(genes_intervals,columns =['start','end','chrom'])#,'gene_name'])
    #Fine VVVV-------------------------------------------------------------
    list_introns = find_introns(df_exons,df_genes_intervals)
    df_introns = pd.DataFrame(list_introns,columns =['start','end','chrom'])
    ###df_introns.drop_duplicates(inplace = True) # may be necessary due to possivle duplications !!!
    #DEBUG
    #print(df_exons.shape, df_introns.shape,df_genes_intervals.shape)#DEBUG
    #display(df_exons)
    #display(df_genes_intervals.sort_values(['start']))
    #/DEBUG
    
    return df_exons, df_introns
#---
def check_overlap_cov(region_stat,region_end,df_item):
    """Check possible items, exons or intros, than overlaps with the selected region"""
    results_exons = 0
    cov_item = 0#desnecessary???
    #print(region_stat,region_end,df_item.shape)#DEBUG
    #Line by line form df_item (item by item):
    for index, row in df_item.iterrows():
        b = [row['start'],row['end']]#exon
        #associated exons
        cov_item = get_Overlap([region_stat,region_end],b)
        #if cov_item != 0: #DEBUG
        #    print(b,'-> ', cov_item)#DEBUG
        results_exons += cov_item
    region_size = region_end - region_stat
    #
    try:
        cov_perc = results_exons/region_size
    except ZeroDivisionError:
        print('ZeroDivisionError: division by zero')    
        print(region_stat,region_end,b)
        print(results_exons,region_size)
    return cov_perc

def get_Genes(region_stat,region_end,chrom,df_exons, df_introns):
    """ exons or introns
    """
    #print(df_exons.shape,df_introns.shape)
    df_exons = df_exons.loc[df_exons["chrom"] == chrom]
    df_introns = df_introns.loc[df_introns["chrom"] == chrom]
    #Line by line form df_exons (Exon by Exon):
    cov_Exons = check_overlap_cov(region_stat,region_end,df_exons)
    #Line by line form df_introns (Intron by Intron):
    cov_Introns = check_overlap_cov(region_stat,region_end,df_introns)
    return cov_Exons, cov_Introns

#---main---:
def main_genes(df_regions,df_genes):
    """
    """
    list_Genes = []
    #transform row['chrom'] items to strings:
    df_genes['chrom'] = df_genes['chrom'].astype(str)
    df_regions['chr_A'] = df_regions['chr_A'].astype(str)
    df_regions['chr_B'] = df_regions['chr_B'].astype(str)
    #Get exons and introns:
    df_exons, df_introns = prep_genes_data_main(df_genes)
    #Line by line from the A nad B regions:
    for index, row in df_regions.iterrows():
        #A
        #print("|----A----|\n")#DEBUG
        rR_A = get_Genes(row['regionA_stat'],row['regionA_end'],row['chr_A'],df_exons, df_introns)
        #B
        #print("|----B----|")#DEBUG
        rR_B = get_Genes(row['regionB_stat'],row['regionB_end'],row['chr_B'],df_exons, df_introns)
        #Prep info by line:
        ##print([row['ID'],row['Cluster_id'],rR_A[0],rR_A[1],rR_B[0],rR_B[1]])#DEBUG
        list_add = [row['ID'],row['Cluster_id'],rR_A[0],rR_A[1],rR_B[0],rR_B[1]]
        #list_add = [row['ID'],row['Cluster_id']] + rR_A + rR_B
        list_Genes.append(list_add)
    return list_Genes

    
##---Load Data---:
region_path = argv[1]
df_regions = pd.read_csv(region_path, sep= ";" ,header=0,index_col=0)
#df_regions = pd.read_csv('./dataset3/rep_regions/match_CNVs_regions_TP.csv', sep= ";" ,header=0,index_col=0)

df_regions['chr_A'] = df_regions['chr_A'].astype(str)
df_regions['chr_B'] = df_regions['chr_B'].astype(str)


#df_genes_INFO = df_genes = pd.read_csv('./dataset3/Genes/'+'hg38_Map50.txt', sep= "\t") #xxxxxxx


df_genes_Hg38 = pd.read_csv('./dataset3/Genes/'+'genes_hg38_processed.csv',  
                            names = ['chrom','type','start','end','strand','gene_name'])
#df_genes_Hg38

##--Run--:
nes = './dataset3/Genes/'

Genes_Data = main_genes(df_regions,df_genes_Hg38)


##Get dataframe:
df_Genes = pd.DataFrame(Genes_Data, columns=['CNV_ID','Cluster_ID','A_cov_Exons','A_cov_Introns','B_cov_Exons','B_cov_Introns'])





###--Save to a csv--:
output_path = argv[2]
df_Genes.to_csv(output_path + 'Genes_regionsAB.csv', sep=';',index_label=False)

#df_Genes.to_csv('./dataset3/Genes/'+'Genes_regionsAB.csv', sep=';',index_label=False)

