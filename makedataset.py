import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from sys import argv
import warnings
warnings.filterwarnings('ignore')

#functions:
def check_for_cols_nans(df_main):
	print("checking:")
	for cols in df_main.columns:
		num_nans = df_main[cols].isna().sum()
		if num_nans != 0:
			print(cols,' -> ',num_nans)
	print("|")
#/

#input and output:
inputf_path = argv[1] # "./output/duplications/"
output_path = argv[2] # "./output/duplication/duplication_"

#Open data:
##Reptitive Regions data:
df_info_regions = pd.read_csv(inputf_path + 'info_regions_matched_test.csv', sep= "\;" ,header=0)
df_info_regions = df_info_regions.rename(columns={'A_region_name':'A_region_class','A_region_family':'A_region_name','A_region_class':'A_region_family','B_region_name':'B_region_class','B_region_family':'B_region_name','B_region_class':'B_region_family'}) #it's due to some lost typo
df_info_regions = df_info_regions.rename(columns={'A_r_dist':'A_repitiveRegion_dist','B_r_dist':'B_repitiveRegion_dist'})
#### Region A and B size:
df_info_regions['A_size_bp'] = df_info_regions['regionA_end'] - df_info_regions['regionA_stat']
df_info_regions['B_size_bp'] = df_info_regions['regionB_end'] - df_info_regions['regionB_stat']
#### Region A and B size in log(Mb +1):
df_info_regions['A_size'] = np.log10((df_info_regions['A_size_bp'] / 1000000) +1) #log10(Mbp +1)
df_info_regions['B_size'] = np.log10((df_info_regions['B_size_bp'] / 1000000) +1)

##Reptitive Regions coverage:
df_coverage_regions = pd.read_csv(inputf_path + 'repRegions_coverage_regionsAB.csv', sep="\;" ,header=0)
df_coverage_regions = df_coverage_regions.rename(columns={'A_cov_segDup':'A_cov_repRegion','B_cov_segDup':'B_cov_repRegion'})

##Segmental duplicaiotns distance:
df_info_segDUPS = pd.read_csv(inputf_path + 'info_segmentalDups_matched_test.csv', sep="\;" ,header=0)
df_info_segDUPS = df_info_segDUPS.rename(columns={'A_r_dist':'A_SegDup_dist','B_r_dist':'B_SegDup_dist'})
##Segmental duplications coverage:
df_coverage_segDUPS = pd.read_csv(inputf_path + 'segDuplications_coverage_regionsAB.csv', sep="\;" ,header=0)
##Segmental duplication Pairs:
df_pairs_segDUPS = pd.read_csv(inputf_path + 'segDups_pairs_match.csv', sep="\;" ,header=0)
df_pairs_segDUPS = df_pairs_segDUPS.rename(columns={'sD_pair_found':'sD_pair'})

##Centromeres and Telomeres distance:
df_centro_telo = pd.read_csv(inputf_path + 'centromere_telomere_matched.csv', sep="\;" ,header=0)

##GC content:
df_GC = pd.read_csv(inputf_path + 'GCcontent_regionsAB.csv', sep="\;" ,header=0)
###formating:
df_GC_A = df_GC.loc[df_GC['Region']=='A']
df_GC_A = df_GC_A.rename(columns={'GC_%':'GC_perc_A'})
df_GC_B = df_GC.loc[df_GC['Region']=='B']
df_GC_B = df_GC_B.rename(columns={'GC_%':'GC_perc_B'})
###Converting to 0->1, because of other used percentages:
df_GC_A['GC_perc_A'] = df_GC_A['GC_perc_A'].apply(lambda x: x/100)
df_GC_B['GC_perc_B'] = df_GC_B['GC_perc_B'].apply(lambda x: x/100)

##Sequence complexity:
df_seqC = pd.read_csv(inputf_path + 'Regs_Seq_Complexity.csv', sep="\;" ,header=0)

##Genes, intron and exon coverage:
df_genes = pd.read_csv(inputf_path + 'Genes_regionsAB.csv', sep="\;" ,header=0)

##Lamina associated domais:
df_lamina = pd.read_csv(inputf_path + 'Lamina_scores_regionsAB.csv', sep="\;" ,header=0)
df_lamina = df_lamina.rename(columns={'A_score':'A_Lamina_score','B_score':'B_Lamina_score'})

##CpG islands:
df_CpGi = pd.read_csv(inputf_path + 'CpG_Island_AB.csv', sep="\;" ,header=0)
df_CpGi = df_CpGi.rename(columns={'A_dist':'A_CpGi_dist','B_dist':'B_CpGi_dist'})

##TAD boundaries:
df_TADs = pd.read_csv(inputf_path + 'TAD_boundary_AB_test.csv', sep="\;" ,header=0)
df_TADs = df_TADs.rename(columns={'A_dist_Log':'A_TAD_dist_Log','B_dist_Log':'B_TAD_dist_Log'})

#Merge data in one dataset:
df_main = pd.concat([
	df_info_regions.set_index('ID').Case_id,
	df_info_regions.set_index('ID').regionA_stat,df_info_regions.set_index('ID').regionA_end,df_info_regions.set_index('ID').chr_A,
	df_info_regions.set_index('ID').regionB_stat,df_info_regions.set_index('ID').regionB_end,df_info_regions.set_index('ID').chr_B,
        df_info_regions.set_index('ID').A_size_bp,df_info_regions.set_index('ID').B_size_bp,
	df_info_regions.set_index('ID').A_size,df_info_regions.set_index('ID').B_size,
	df_info_regions.set_index('ID').A_region_name,df_info_regions.set_index('ID').A_region_family,df_info_regions.set_index('ID').A_region_class,
	df_info_regions.set_index('ID').B_region_name,df_info_regions.set_index('ID').B_region_family,df_info_regions.set_index('ID').B_region_class,
	df_coverage_regions.set_index('CNV_ID').Name_A,df_coverage_regions.set_index('CNV_ID').Family_A,df_coverage_regions.set_index('CNV_ID').Class_A, 
	df_coverage_regions.set_index('CNV_ID').Name_B,df_coverage_regions.set_index('CNV_ID').Family_B,df_coverage_regions.set_index('CNV_ID').Class_B,
	df_info_regions.set_index('ID').A_repitiveRegion_dist,df_info_regions.set_index('ID').B_repitiveRegion_dist,
	df_coverage_regions.set_index('CNV_ID').A_cov_repRegion,df_coverage_regions.set_index('CNV_ID').B_cov_repRegion,
	df_info_segDUPS.set_index('ID').A_SegDup_dist,df_info_segDUPS.set_index('ID').B_SegDup_dist,
	df_coverage_segDUPS.set_index('CNV_ID').A_cov_segDup,df_coverage_segDUPS.set_index('CNV_ID').B_cov_segDup,
	df_pairs_segDUPS.set_index('CNV_ID').sD_pair,
	df_centro_telo.set_index('CNV_id').dist_Centromere,df_centro_telo.set_index('CNV_id').dist_Telomere,
	df_GC_A.set_index('CNV_ID')['GC_perc_A'],df_GC_B.set_index('CNV_ID')['GC_perc_B'],
	df_seqC.set_index('CNV_ID').Dust_A,df_seqC.set_index('CNV_ID').Dust_B,
	df_genes.set_index('CNV_ID').A_cov_Exons,df_genes.set_index('CNV_ID').A_cov_Introns,
	df_genes.set_index('CNV_ID').B_cov_Exons,df_genes.set_index('CNV_ID').B_cov_Introns,
	df_lamina.set_index('CNV_ID').A_Lamina_score,df_lamina.set_index('CNV_ID').B_Lamina_score,
	df_CpGi.set_index('CNV_ID').A_CpGi_dist,df_CpGi.set_index('CNV_ID').B_CpGi_dist,
	df_TADs.set_index('CNV_ID').A_TAD_dist_Log,df_TADs.set_index('CNV_ID').B_TAD_dist_Log],
			ignore_index=False,axis=1)
###
df_main.reset_index(inplace=True)
df_main.rename(columns={"index":"CNV_ID"},inplace=True)

##Solve  possible missing data in Lamina Associated Domains data:
###Replace Lamina scores nan with 0 !!!
df_main["A_Lamina_score"].fillna(0, inplace=True)
df_main["B_Lamina_score"].fillna(0, inplace=True)

#Check if remains some missing data:
if df_main.isnull().sum().sum() != 0:
	check_for_cols_nans(df_main)

### Reptitive regions name, family and class by overlaping and/or closest distance:
### **Replace "no_overlap" with the data of the closest repRegions:**
for index, row in (df_main.loc[(df_main['Name_A'] == 'no_overlap') | (df_main['Name_B'] == 'no_overlap')]).iterrows():
	r_id = row.CNV_ID
    
	if row['Name_A'] == 'no_overlap':
		df_main["Name_A"].loc[df_main['CNV_ID'] == r_id] = row['A_region_name']
		df_main["Family_A"].loc[df_main['CNV_ID'] == r_id] = row['A_region_family']
		df_main["Class_A"].loc[df_main['CNV_ID'] == r_id] = row['A_region_class']
        
	if row['Name_B'] == 'no_overlap':
		df_main["Name_B"].loc[df_main['CNV_ID'] == r_id] = row['B_region_name']
		df_main["Family_B"].loc[df_main['CNV_ID'] == r_id] = row['B_region_family']
		df_main["Class_B"].loc[df_main['CNV_ID'] == r_id] = row['B_region_class']
### **Drop columns of the closest repRegions name,family,class:**
df_main.drop(['A_region_name', 'A_region_family', 'A_region_class','B_region_name', 'B_region_family', 'B_region_class'], axis=1, inplace=True)
### **Rename repRegions name,family,class columns:**
df_main.rename({'Name_A': 'A_region_name', 'Family_A': 'A_region_family', 'Class_A':'A_region_class',
		'Name_B': 'B_region_name', 'Family_B': 'B_region_family', 'Class_B':'B_region_class'}, axis=1, inplace = True)

# Save dataset:
df_main.to_csv(output_path + 'dataset.csv', sep = ';' ,index_label=False)
