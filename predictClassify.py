import pandas as pd
import numpy as np
import autosklearn
from sys import argv
import pickle

def make_prediction(loaded_classifier,df_testset):
	""" Use a selected classifier model to predict result from a selected dataframe.
	"""
	#Defining X and Y adata:
	X = df_testset.drop(["CNV_ID","Case_id","regionA_stat","regionA_end","chr_A",
			"regionB_stat","regionB_end","chr_B","A_size_bp","B_size_bp"],axis=1)
	#
	desired_categorical_columns = ['A_region_name','A_region_family','A_region_class',
			'B_region_name','B_region_family','B_region_class','sD_pair']
	for column in X.columns:
		if column in desired_categorical_columns:
			X[column] = X[column].astype('category')
	#target
	#y = df_testset["target"].astype('bool')
	#Start auto-sklearn:
	#X_test, y_test = X, y
	X_test = X

	# predict
	#y_true = y_test
	y_pred = loaded_classifier.predict(X_test)
	#	
	#
	return y_pred

def select_classify(df_main,CNV_type):
	"""	Selects model for deletion or duplication and makes predidtion if it is True or False.
		Requires: dataset with features of the CNVs for classification/predition, model in .pkl files
		Ensures: returns a .csv table with the CNVs ids and classification result, True or False. 
	"""
	CNV_type = CNV_type.lower()
	if CNV_type == "del" or CNV_type == "deletion":
		model_path = "./models/ensembled_model_deletions_balanced.pkl"
	elif CNV_type == "dup" or CNV_type == "duplication":
		model_path = "./models/ensembled_model_duplications_balanced.pkl"
	else:
		print("Please insert CNV type as deletion/del or duplication/dup.")
		return "None"
	#
	with open(model_path, 'rb') as f:
		loaded_model = pickle.load(f)
	#Results dataframe:
	df_results = df_main[['CNV_ID']].copy()
	#Classify:
	predicted = make_prediction(loaded_model,df_main)
	df_results['predicted_CNV']= predicted
	df_results['predicted_CNV'] = df_results['predicted_CNV'].astype(bool)
	#
	return df_results

def main_process(dataset_path,CNV_type,output_path):
	"""Start process of prediction the avaliable deletions or duplications in a dataset if they are True or False.
	"""
	df_main = pd.read_csv(dataset_path ,sep=';')
	df_results = select_classify(df_main,CNV_type)
	df_results.to_csv(output_path + CNV_type + "_prediction_results.csv", sep=';',index_label=False)
#___
input_dataset = argv[1] + "dataset.csv"
CNV_type = argv[2]
output_path = argv[3]

main_process(input_dataset,CNV_type,output_path)
