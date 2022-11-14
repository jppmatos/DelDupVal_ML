# DDval_ML
A ML tool for validation of pos-mapping CNVs.

### How to run it
1. Start with RUN_validate_DUPS-DGRC0005.sh for duplications.
2. ...
3. ...
```
sh RUN_validate_DUPS-DGRC0005.sh
```
4. the output will be avaliable as:
```
./output/...
```

## Example:
. Start with RUN_validate_DELS-DGRC0005_SAMPLE.sh for deletions or RUN_validate_DUPS-DGRC0005_SAMPLE.sh for duplications.
```
sh RUN_validate_DELS-DGRC0005.sh 
```
### How the script works:
1. The variable **input_TP_4** it's the input were the CNV are inserted;
2. **output_path_TP_4** variable is the defined output, were the results and other content will be avaliable;
3. Several python scripts will process a series of sub-datasets by matching genomic data ver the given CNVs. These are: 
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
4. All the sub-datasets will be merged and processed by **makedataset-deletions.py** (if deletions, or **makedatase.py** if duplications);
5. Then after having the complete dataset (**dataset.csv**), **preditClassify.py** will classify the CNVs as possible True CNV or False CNV. The classification script need to know if the CNVs are deletions (DEL) or duplication (DUPS);
6. A *.csv* report is made at the defined output path in **output_path_TP_4** named as **DEL_prediction_results** (deletions) or **DUP_prediction_results.csv** (duplications).
  
