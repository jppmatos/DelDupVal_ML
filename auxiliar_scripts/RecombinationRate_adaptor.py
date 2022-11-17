import pandas as pd
#import matplotlib.pyplot as plt
#import seaborn as sns
import numpy as np
#from scipy.stats import norm
#from sklearn.preprocessing import StandardScaler
#from scipy import stats
import os
import csv


def adapt(PATH, file,new_name):
	"""
	"""
	counter = 0
	with open(PATH + file,'r', newline='') as file_r:
		with open(PATH + new_name, 'w', newline='') as file_w:
			writer = csv.writer(file_w)
			for line in file_r:
				if "#chrom" not in line:
					counter += 1
					#print(line)
					line = line.split()
					
					id_rr = str(counter)
					line[3] = ('rR_'+ id_rr)
					writer.writerow(line)

def adapt_bed(PATH, file,new_name):
	"""
	"""
	counter = 0
	with open(PATH + file,'r', newline='') as file_r:
		fw = open(PATH + new_name, 'w')d
		#writer = csv.writer(file_w)
		for line in file_r:
			if "#chrom" not in line:
				counter += 1
				id_rr = str(counter)
				line = line.replace("recombRate","rR_" + id_rr)
				fw.write(line)
		fw.close()



def make_match(list_csv, list_bed):
	"""
	"""
	list_new_data = []
	start = 0
	for item in list_bed:
		for y in range(start,len(list_csv)):
			if item[3] == list_csv[y][3]:
				start = y
				data_bed = item#[0:3]
				data_csv = list_csv[y][4:]
				data_new = data_bed + data_csv
				list_new_data.append(data_new)
	return list_new_data



def make_csv_with_bed(bed_file,csv_file,new_csv_file):
	"""
	"""
	with open(csv_file,'r', newline='') as file_csv_info:
		reader = csv.reader(file_csv_info)
		list_csv_info = list(reader)

	list_bed = []
	with open(bed_file,'r', newline='') as file_bed:
		for line_bed in file_bed:
			line_bed = line_bed.split()
			#Beliving that it is in order!!!
			list_bed.append(line_bed)
	#print(line_bed)#DEBUG OK!
	list_match = make_match(list_csv_info,list_bed)

	#print(list_match[:10])#DEBUG

	with open(new_csv_file, "w", newline="") as new_f:
		writer = csv.writer(new_f)
		writer.writerows(list_match)





path_RR_content_hg19 = './dataset3/RecombinationRate/hg19/'
name_txt = 'recombRate' 
name_txt_new = 'recombRate_info_hg19.csv'
#adapt(path_RR_content_hg19,name_txt,name_txt_new)


##name_bed = 'recombRate.bed' 
##name_bed_new = 'recombRate_bed_hg19.bed'
#adapt_bed(path_RR_content_hg19,name_bed,name_bed_new)

#Convert 'recombRate_bed_hg19.bed' to Hg38 with LiftOver, http://genome.ucsc.edu/cgi-bin/hgLiftOver

path_RR_content_hg38 = './dataset3/RecombinationRate/hg38/'
bed_file_hg_38 = 'hglft_genome_target.bed'
file_info = 'RecombinationRate_hg_38.csv'
make_csv_with_bed(path_RR_content_hg38+bed_file_hg_38,path_RR_content_hg19+name_txt_new,path_RR_content_hg38+file_info)