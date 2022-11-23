import pandas as pd
import numpy as np
import autosklearn
import os
import pickle
from sklearn.model_selection import train_test_split
from sklearn.metrics import accuracy_score
from autosklearn.metrics import balanced_accuracy, precision, recall, f1, roc_auc
from autosklearn.classification import AutoSklearnClassifier
from sklearn.metrics import confusion_matrix
from imblearn.under_sampling import RandomUnderSampler

from sklearn.model_selection import StratifiedKFold

#Total of threads in use:
TOTAL_JOBS=5
#File settings:
out_path = '../auto_sklearn_data'
path = '../dataset3/data_main/'


###file_name = 'refined_data_DUPS_Talkowski-autoSK.csv'

file_name = 'DUPS_12_04/DGRC_liGS-refined_data_DUPS-autoSK_train.csv'

df_main = pd.read_csv(path + file_name,sep=';')
#Make directory instructions:
directory_name = "DUPS_Talkowski_4hour-22-04_balanced_noSize" #name of output folder!
storage_folder='/'+ directory_name

#Do datasampling? (pseudo_sampling): (20% extra)
pbalancing_switch = False #!!!!!!!!!!!!!!!
undersamp_balancing_switch = True
#:
#path = os.path.join(out_path, directory_name)
#os.mkdir(path)
#
logfilename = out_path + storage_folder + '/autosklearn_Logfile.txt'
resultsfilename = out_path + storage_folder + '/autosklearn_cv_results_mean_test_score.txt'
paramsfilename = out_path + storage_folder + '/autosklearn_cv_results_params.txt'
temp_folder = out_path + storage_folder + '/tmp' 
infoTable = out_path + storage_folder + '/info_models.csv'
LeaderBoard = out_path + storage_folder + '/learderboard.csv'
models_folder = out_path + storage_folder
#auto-sklearn settings:
time_run = 14400 #10800 ##7200 #3600 #86400 #  #57600 ##86400 #7200 #36000 #10*60 #7200#21600#(30*60)
time_lim_per = 5400 #2400 #43200 # #30*60 #50*60 #3*60 #5*60 #2300 #3600#30*60

#folds/splits for cross valitation method:
n_folds = 10 

#why 10 fold?: https://machinelearningmastery.com/cross-validation-for-imbalanced-classification/
#— Page 205, Imbalanced Learning: Foundations, Algorithms, and Applications, 2013.

#____________________________________________________
def error(solution, prediction):
	#custom function defining error
	return np.mean(solution != prediction)

def FalsePositive_rate(solution, prediction):
	#probability of false alarm
	tp, fp, fn, tn = confusion_matrix(solution,prediction).ravel()
	return (fp) / (tp + tn)

def FalseNegative_rate(solution, prediction):
        #miss rate
        tp, fp, fn, tn = confusion_matrix(solution,prediction).ravel()
        return (fn) / (tp + tn)
#___________________________________________________
def get_metric_result(cv_results):
	results = pd.DataFrame.from_dict(cv_results)
	results = results[results['status'] == "Success"]
	cols = ['rank_test_scores', 'param_classifier:__choice__', 'mean_test_score']
	cols.extend([key for key in cv_results.keys() if key.startswith('metric_')])
	#Make info_table:
	info = cols.copy()
	info.extend(('metric_balanced_accuracy','params'))
	df_infomatrix = results[info]
	#df_infomatrix = results[['rank_test_scores','param_classifier:__choice__', 'mean_test_score','metric_balanced_accuracy','params']]
	df_infomatrix.to_csv(infoTable, sep = ';' ,index_label=False)
	#
	return results[cols], df_infomatrix
#
#
def find_best_log_models(df_resultsLOG,model_autoSK,n_lim,save_path):
	"""Look for the n best models in between *.cv_results* and *.automL_.models_.items()*.
		Requires: df_resultsLOG, model.cv_results dataframe; model_autoSK, ensembled models; n_lim, the n first from df_resultsLOG; save_path, path to save the models.
		Ensures: found models in model_autoSK.automL_.models_.items() are saved in *save_path*.
	"""
	df_bestLOG = df_resultsLOG.reset_index()
	for model_place in range(n_lim):
		selected_model_LOG_params = df_bestLOG['params'].loc[model_place]
		#selected_model:
		searchfor_model(selected_model_LOG_params,model_autoSK,save_path,model_place)

def searchfor_model(wanted_model_params,model,save_path,n_rank):
	"""Search of a model identified in models.cs_results() and save it."""
	#wanted_model_params = df_wanted_model['params'].loc[model_place]
	wanted_model_params = str(wanted_model_params)
	for (seed, model_id, budget), model in model.automl_.models_.items():
		models_params = model.__str__()
		search_params = models_params.split("SimpleClassificationPipeline")[1].strip("()")
		search_params = search_params.split(",\ndataset_properties")[0]
		#print(type(search_params),'\n',search_params)
		#print(search_params,'\n')
		if search_params == wanted_model_params:
			model_file = '/classifier_model_' + str(n_rank) +'_best.pkl'
			print(model_file,'\n \n',model)
			with open(save_path + model_file, 'wb') as f:
				pickle.dump(model, f)

def save_ensembled_models(models,save_path):
	""" Save models variable that has all enseble models by auto-SKlean-
	"""
	with open(save_path + "/all_ensembled_models.pkl", 'wb') as f:
		pickle.dump(models,f)
        
#----------------------------------------------------------------
def pseudo_balancing(df_main_T,perc_keep_extra = .2):
	""" Removes a part of the data of the majority class
		Requires: df_main_T, perc_keep_extra: percentage of data to keep extra over the minority class (in decimal)
		Ensures: returns the samples dataset
	"""
	#check majo class :
	if df_main_T["target"].value_counts()[0] > df_main_T["target"].value_counts()[1]:
		maj = 0
		mino = 1
	else:
		maj = 1
		mino = 0
	#_________________________________________________________________________________________
	maj_val = df_main_T["target"].value_counts()[maj]
	mino_val = df_main_T["target"].value_counts()[mino]
    #
	to_add = mino_val + (mino_val * perc_keep_extra)
	#_________________________________________________________________________________________  
	perc_remove = 1 - ((1 * (to_add))/ maj_val)
	df_sampled = df_main_T.drop(df_main_T.loc[df_main_T['target'] == maj].sample(frac= perc_remove, random_state=1 ).index)
    #---
	count_original = df_main_T["target"].value_counts()
	count_new = df_sampled["target"].value_counts()
    #---
	return df_sampled, count_original, count_new
#----------------------------------------------------------------
#_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|
#Pseudo_balancing undersampling:
count_original = df_main["target"].value_counts()
count_new = count_original
#
if pbalancing_switch == True:
	df_main,count_original,count_new  = pseudo_balancing(df_main, 0.2) 
#_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|

#Defining X and Y adata:
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
##!X = df_main.drop(["CNV_ID","target"],axis=1)
X = df_main.drop(["CNV_ID",'A_size','B_size',"target"],axis=1) #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#Categorical_columns:
desired_categorical_columns = ['A_region_name','A_region_family','A_region_class','B_region_name','B_region_family','B_region_class','sD_pair']
for column in X.columns:
    if column in desired_categorical_columns:
        X[column] = X[column].astype('category')
#    else:
#        X[column] = pd.to_numeric(X[column])

#Target
y = df_main["target"].astype('bool')

#Start auto-sklearn:
###X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.33, random_state=1)
X_train, y_train = X, y

#_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|
###Data sampling:
if undersamp_balancing_switch == True:
    rus = RandomUnderSampler(random_state=0)
    X_train, y_train = rus.fit_resample(X_train, y_train) #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|

#models_run:__________________________________________________________________________________________________________________
if __name__ == '__main__':
	#Making direcory:
	path = os.path.join(out_path, directory_name)
	os.mkdir(path)
	print("running ...")
	print('-|-|-|-'*10)
    #Save a sample from train data:
	X_train.head(5).to_csv(out_path + storage_folder + '/trainset_sample.csv', sep = ';' ,index_label=False)
	print('X_train size: ', X_train.shape)
	print('Trainset info: \n')
	#print("Dataset size: \n",count_original)
	print("Trainset size: \n",count_new)
	print('_|_|_|_'*10)
	#Custom scores setup:
	error_rate = autosklearn.metrics.make_scorer(
		name='custom_error',
		score_func=error,
		optimum=0,
		greater_is_better=False,
		needs_proba=False,
		needs_threshold=False)

	FP_rate = autosklearn.metrics.make_scorer(
		name='FP_Rate',
		score_func=FalsePositive_rate,
		optimum=0,
		greater_is_better=False,
		needs_proba=False,
		needs_threshold=False)

	FN_rate = autosklearn.metrics.make_scorer(
		name='FN_Rate',
		score_func=FalseNegative_rate,
		optimum=0,
		greater_is_better=False,
		needs_proba=False,
		needs_threshold=False)

	#Stratified Cross validadtion: (25/03)
	skf = StratifiedKFold(n_splits=n_folds)

	#models_run:
	model_autoSK = AutoSklearnClassifier(time_left_for_this_task=time_run, per_run_time_limit=time_lim_per,
						n_jobs=TOTAL_JOBS,
						resampling_strategy=skf,
						ensemble_size=1000,
						ensemble_nbest=200,
						max_models_on_disc=500,
						tmp_folder= temp_folder,
						delete_tmp_folder_after_terminate= True,
						scoring_functions=[balanced_accuracy,precision, recall,f1,roc_auc])
						#scoring_functions=[balanced_accuracy,f1,roc_auc,FP_rate])
						#scoring_functions=[balanced_accuracy,precision, recall,f1,roc_auc])

	#preform the seach:
	model_autoSK.fit(X_train,y_train)
	#
	print("Used dataset: %s \n" % (file_name),file=open(logfilename,'a'))
	#
	#print("X_train size: %s \n" % (X_train.shape),file=open(logfilename,'a'))
	print("Trainset values:",file=open(logfilename,'a'))
	print("%s \n" % (count_new),file=open(logfilename,'a'))
	#
	#
	print("Used_jobs: %s " % (TOTAL_JOBS),file=open(logfilename,'a'))
	print("Time_run: %s " % (time_run),file=open(logfilename,'a'))
	print("Time_lim_per: %s " % (time_lim_per),file=open(logfilename,'a'))
	print("Time_run: %s " % (time_run))
	print("Time_lim_per: %s " % (time_lim_per))
	#results, logging:
	print("time_left_for_this_task: %s \n per_run_time_limit: %s \n" % (time_run,time_lim_per))
	print("time_left_for_this_task: %s \n per_run_time_limit: %s \n" % (time_run,time_lim_per),file=open(logfilename,'a'))
	print(model_autoSK.sprint_statistics(),file=open(logfilename,'a'))

	#Show runing models:
	print(model_autoSK.show_models(),file=open(logfilename,'a'))
	
	#
	#

	#### evaluate best model
	###y_hat = model_autoSK.predict(X_test)
	###print('Before re-fit',file=open(logfilename,'a'))
	###acc = accuracy_score(y_test, y_hat)
	
	#Save params data:
	print(np.argmax(model_autoSK.cv_results_['mean_test_score']),file=open(resultsfilename,'a'))
	print(model_autoSK.cv_results_['params'],file=open(paramsfilename,'a'))

	###print("Accuracy: %.3f" % acc,file=open(logfilename,'a'))
	
        #Confusion Matrix:
	###print(confusion_matrix(y_test,y_hat))
	###tp, fp, fn, tn = confusion_matrix(y_test,y_hat).ravel()
	###print("|TP: %s,FP: %s| \n|FN: %s ,TN: %s| \n \n " % (tp, fp, fn, tn),file=open(logfilename,'a'))

	#Refit
	print('After re-fit',file=open(logfilename,'a'))
	
	#Show runing models:
	print(model_autoSK.show_models(),file=open(logfilename,'a'))

	#Save params data:
	print(np.argmax(model_autoSK.cv_results_['mean_test_score']),file=open(resultsfilename,'a'))
	print(model_autoSK.cv_results_['params'],file=open(paramsfilename,'a'))

	
	#Refit, use all avaliable data: -----------------------------------------------------------------------------
	model_autoSK.refit(X_train.copy(), y_train.copy())
	###predictions = model_autoSK.predict(X_test)
	###print("Accuracy score CV:", accuracy_score(y_test, predictions),file=open(logfilename,'a'))
		
	#Confusion Matrix:
	##print(confusion_matrix(y_test,predictions))
	##r_tp, r_fp, r_fn, r_tn = confusion_matrix(y_test,predictions).ravel()
	#print(confusion_matrix(y_test,y_test).ravel())#DEBUG
	##print("|TP: %s,FP: %s| \n|FN: %s ,TN: %s| \n \n " % (r_tp, r_fp, r_fn, r_tn),file=open(logfilename,'a'))
	
	#*.cv_results_*
	print('#'* 80,file=open(logfilename,'a'))
	print("Metrics results",file=open(logfilename,'a'))
	data_metrics, df_metrics_results = get_metric_result(model_autoSK.cv_results_)
	print(data_metrics.to_string(index=False),file=open(logfilename,'a'))
	#start extrating the models:
	df_metrics_results = df_metrics_results.loc[:,~df_metrics_results.columns.duplicated()].sort_values(['metric_balanced_accuracy','metric_recall'],axis=0, ascending=False)
	find_best_log_models(df_metrics_results,model_autoSK,12,models_folder) #grab models comun to metrics_results
	save_ensembled_models(model_autoSK,models_folder) # save all ensembled models
	#Leaderboard:
	df_LBoard = model_autoSK.leaderboard()
	df_LBoard.to_csv(LeaderBoard, sep = ';' ,index_label=False)
