import pandas as pd
#import matplotlib.pyplot as plt
#import seaborn as sns
import numpy as np
#from scipy.stats import norm
#from sklearn.preprocessing import StandardScaler
#from scipy import stats
import os
import csv
def main_1(PATH,file_name,new_file_name):
	"""
	"""
	with open(PATH + new_file_name, 'w', newline='') as file:
		writer = csv.writer(file) 
		with open(PATH + file_name, 'r', newline='') as r_file:
			for line in r_file:
				if '#!' not in line:
					line = line.split('\t')
					if (line[2] == 'gene' or line[2] == 'exon') and (line[0] != 'MT'):
						to_add = [line[0],line[2],line[3],line[4],line[6],line[8].strip('gene_id "')]
						#print(line)#DEBUG
						#print(to_add)#DEBUG
						#break#DEBUG
						writer.writerow(to_add)




PATH = './dataset3/Genes/'
file_name = 'genes_hg38.txt'
new_file_name = 'genes_hg38_processed.csv'

main_1(PATH,file_name,new_file_name)
