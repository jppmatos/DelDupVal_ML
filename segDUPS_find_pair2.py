import pandas as pd
from sys import argv

def check_overlap(a,b):
    """ checks for overlap between two itervals.
        Requires: a= region start and end (list), b=  chromStart and chromEnd from rR (list)
    """
    ovl=max(0,min(a[1],b[1])-max(a[0],b[0]))#calcula o overlap
    #print(ovl)#DEBUG
    return ovl
    


def get_sD(CNV_row,df_sD):
    """
        Search for overlaping segmental duplication with region A of the CNV and then check if its pairs does overlap with region B.
        Requires: CNV_row as the CNV region A and B information, df_sD is the dataframe with the segmental duplication pairs.
        Ensures: returns a boolen as True when the pair has overlap match, else it is retuned as False.
    """
    start_A = CNV_row['regionA_stat']
    end_A = CNV_row['regionA_end']
    chrom = CNV_row['chr_A']
    start_B = CNV_row['regionB_stat']
    end_B = CNV_row['regionB_end']
    #
    l_IDs_sD = []
    #Reduce sD dataset to a selected chromossome:
    df_chr_sD = df_sD.loc[df_sD['chrom']==('chr'+str(chrom))]
    #print(chrom,df_sD.shape)#DEBUG
    for index, row in df_chr_sD.iterrows():
        overlap = check_overlap( [start_A,end_A] ,[ row['chromStart'],row['chromEnd'] ] )
        if overlap != 0:
            overlap = check_overlap( [ row['otherStart'],row['otherEnd'] ], [start_B,end_B] )
            if overlap != 0:
                return True
    return False

def match_segDups(df_regions,df_sD):
    """
        In each CNV there'll be found a matching, overlap, segmental duplication to the region A or B of the CNV, \
        then same'll be none to the pair of the segmatal duplication that overlaps with the other region of the CNV.
    """
    # Just intrachromossomal sD pairs:
    df_sD = df_sD.loc[(df_sD['chrom']==df_sD['otherChrom'])]
    #
    list_result=[]
    for index, row in df_regions.iterrows():
        rR_pair = get_sD(row,df_sD) # is there any overlap match sD pair? (True,False)
        #
        result_item = [row['Case_id'],row['ID'], row['Cluster_id'],
                        row['regionA_stat'],row['regionA_end'],row['chr_A'],
                        row['regionB_stat'],row['regionB_end'],row['chr_B'],
                        rR_pair]
                        #pair_A,pair_B]
        list_result.append(result_item)   
    return list_result
#____________________________________________________________________________
input_segDUPS_pairs = "./dataset3/rep_regions/segDuplications_pairs.csv"
df_segDups = pd.read_csv(input_segDUPS_pairs, sep= ';')

region_path = argv[1] #"./dataset3/rep_regions/match_CNVs_regions_TP.csv" 
output_path = argv[2] #"./output/tp/"
#
df_CNV_regions = pd.read_csv(region_path, sep= ";" ,header=0,index_col=0)

df_CNV_regions['chr_A'] = df_CNV_regions['chr_A'].astype(str)
df_CNV_regions['chr_B'] = df_CNV_regions['chr_B'].astype(str)
#
l_sD_pairs = match_segDups(df_CNV_regions,df_segDups)

df_sD_pairs = pd.DataFrame(l_sD_pairs, columns=['Case_id','CNV_ID','Cluster_id',
                                                'regionA_stat','regionA_end','chr_A',
                                                'regionB_stat','regionB_end','chr_B',
                                                'sD_pair_found',])
df_sD_pairs.to_csv(output_path + 'segDups_pairs_match.csv' , sep=';',index_label=False)

