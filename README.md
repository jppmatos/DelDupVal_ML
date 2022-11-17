# DDval_ML
DDval_ML is a ML tool for validation of pos-mapping deletion and duplications copy number variations (CNVs). It takes use of clusters of identified improper read-pairs to define a flanking region of the CNV it self (Figure 1), then genomic information is overlaped against the flacking regions of it is provided enough information to chracterize the CNV. 

![image](https://user-images.githubusercontent.com/44948470/201972334-130a94aa-25a1-41cb-b4fc-18ff0cb761c8.png)
![image](https://user-images.githubusercontent.com/44948470/201972350-422dce3d-42d7-4fc3-96f1-1c73b04da057.png)

Figure 1: Representation of a (A) Deletion and a (B) Tandem Duplication with their respectively flanking regions. Improper pairs are depicted as green arrows connected by a dashed line. X’ and Y’ represent the CNV’s breakpoints and X and Y mark the limits of the improper pair cluster.


The CNVs used for the model's training were reads mapped with large-insert genomic sequencing (liGS) libraries. 


## How to run it
1. Download DDval_ML repository
2. Download the [subdatasets](https://www.dropbox.com/s/lvbga9cnay5dwq5/dataset3.zip?dl=0) data
3. Uncompress the dataset3.zip file, it has to be in the same diretory within the repository
4. Select the bash script for the designed task, RUN_validate_DELS-DGRC0005.sh for deletions or RUN_validate_DUPS-DGRC0005.sh for duplications;
5. Modify the selected bash script by inserting the input and output path;
6. Run the bash script;
7. At the output directory there'll be a cvs report file with the classification results, **DEL_prediction_results** (deletions) or **DUP_prediction_results.csv** (duplications).


### Input:
. The input has to be a csv file the CNVs, each CNVs has to have their mapping infomation such as start and end position in basepairs of the CNV it self and the flanking regions made by the improper read pairs:
```
Case_id;ID;Cluster_id;start_position_max;end_position_max;regionA_stat;regionA_end;regionB_stat;regionB_end;chr_A;chr_B;libraries;A_size;B_size
0;DGRC0005;fp_del_56;56;4064277;4069587;4062862;4064277;4069587;4071328;1;1;liGS;1415;1741
1;DGRC0005;fp_del_126;126;7630994;7636797;7629784;7630994;7636797;7637618;1;1;liGS;1210;821

```
### Output:
. The output will be a csv file with the CNV id and the classification result as True or False:
```
CNV_ID;predicted_CNV
0;DGRC0005_dup_1;True
1;DGRC0005_dup_2;False

```

## Example:
. Edit and run RUN_validate_DELS-DGRC0005_SAMPLE.sh for deletions or RUN_validate_DUPS-DGRC0005_SAMPLE.sh for duplications.
```
sh RUN_validate_DELS-DGRC0005.sh 
```
### How the script works:
1. The variable **input_TP_4** it's the input were the CNV are inserted;
```bash
input_TP_4=./DGRC0005_sampled/DGRC0005_dels_SAMPLE.csv
```
3. **output_path_TP_4** variable is the defined output, were the results and other content will be avaliable;
```bash
mkdir ./DGRC0005_sampled/output/
mkdir ./DGRC0005_sampled/output/deletions/
output_path_TP_4=./DGRC0005_sampled/output/deletions/
```
5. Several python scripts will process a series of sub-datasets by matching genomic data (present in *./dataset3* folder) to the given CNVs. These are: 
	 - **look_repRegions2.py**
	 - **look_segDuplications2.py**
	 - **segDuplication_coverage.py**
	 - **repRegions_coverage.py**
	 - **segDUPS_find_pair2.py**
	 - **Centro_Telo_match2.py**
	 - **GCcontent2.py**
	 - **SeqComplexity.py**
	 - **Recombination_Rate2.py** (exclusive for deletions)
	 - **Genes_V3.py**
	 - **Lamina_Associ_Dom.py**
	 - **CpG_island2.py**
	 - **TADs2.py**
6. All the sub-datasets will be merged and processed by **makedataset-deletions.py** (if deletions, or **makedatase.py** if duplications);
7. Then after having the complete dataset (**dataset.csv**), **preditClassify.py** will classify the CNVs as possible True CNV or False CNV. The classification script need to know if the CNVs are deletions (DEL) or duplication (DUPS);
8. A *.csv* report is made at the defined output path in **output_path_TP_4** named as **DEL_prediction_results** (deletions) or **DUP_prediction_results.csv** (duplications).
  
 ## Requirements:
 - Linux operative system
 - python3  (>=3.9.7)
   - auto-sklearn (==0.13.0) 
   - numpy (>=1.22.3)
   - pyfaidx (>=0.5.9.5)
   - scikit-learn (>=1.0.2)
   - pandas (>=1.1.3)

---

## Adicional Notes:
Most of the genomic data used as features for the dataset came from [Genome Browser](https://genome.ucsc.edu/). Nearly all data add to be previouslly processed by script in *auxiliar_scripts* folder before being process by the script for that match the genomic data the CNVs. Just the script of **SeqComplexity.py**, **Lamina_Associ_Dom.py** and **CpG_island2.py** use directly the extrated raw data. While the other ones had to use auxliliar scripts:
- **look_repRegions2.py**, **look_segDuplications2.py** ,**segDuplication_coverage.py** and **segDUPS_find_pair2.py** uses the **adaptData_repRegions.py**
- Centromeres and Telomeres raw data had processed by only keeping the row that where *P11.1* and *q11* (centromeres), and *telomere* in the *type* column  (**Centro_Telo_match2.py**) 
- Recombination Rate's raw data had been process from the hg19 format to hg38 by the [LiftOver](https://genome.ucsc.edu/cgi-bin/hgLiftOver) software (**Recombination_Rate2.py**)
- Gene raw data came from [Ensembl](https://www.ensembl.org/index.html), process by **extrat_genes_data.sh** and **Genes_data_processing.py** (**Genes_V3.py**)
- TAD's raw data where from IMR90_Rao_2014's data (**TADs2.py**)
