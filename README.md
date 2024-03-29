# DelDupVal_ML
DelDupVal_ML is a ML tool for validation of pos-mapping deletions and duplications, copy number variations (CNVs), sequenced with large-insert genomic sequencing libraries that could been mislabeled as an actual CNV due to possible artifacts. 

## How it works?
DelDupVal_ML makes use of clusters of identified improper read-pairs to define a flanking region of the CNV it self (Figure 1), then genomic information is overlaped against the flacking regions of it is provided enough information to chracterize the CNV. 

![image](https://user-images.githubusercontent.com/44948470/201972334-130a94aa-25a1-41cb-b4fc-18ff0cb761c8.png)
![image](https://user-images.githubusercontent.com/44948470/201972350-422dce3d-42d7-4fc3-96f1-1c73b04da057.png)

Figure 1: Representation of a (A) Deletion and a (B) Tandem Duplication with their respectively flanking regions. Improper pairs are depicted as green arrows connected by a dashed line. X’ and Y’ represent the CNV’s breakpoints and X and Y mark the limits of the improper pair cluster.

To characterize the genomic regions X and Y (refered as regions A and B in the entire code), genomic annotation data based on a set of genomic parameters according to *Li et al., 2020* was retrieved for those regions and organized into a dataset as features:
1. **Regions Size**: Regions X and Y size, calculated as log10(Mb of distance +1).
2. **Repetitive Regions** (look_repRegions2.py; repRegions_coverage.py): Distance between regions X/Y and Repetitive Region (RRegions) was calculated as log10(Mb of distance +1), where the value 0 was given when overlap was detected. In case of overlap, the percentage of the X/Y regions covered by RRegions was included. Also, name, family and class of RRegions was included on the dataset. Data was obtained from [genome browser’s RepeatMasker UCSC genome browser track](http://genome-euro.ucsc.edu/cgi-bin/hgTrackUi?hgsid=258916625_zlDCrzMw0jagpUpAa8G5uTCgOjX4&c=chrY&g=joinedRmsk).
3. **Segmental Duplications** (look_segDuplications2.py; segDuplications_coverage.py; segDUPS_find_pair2.py): Distance and overlap data were retrieved as indicated for 2 – RRegions. Additionally, pair segmental duplications (SegDups) overlapping both X/Y regions, were included on the dataset. Data was retrieved from [UCSC genome browser’s Segmental duplications track](http://genome-euro.ucsc.edu/cgi-bin/hgTrackUi?hgsid=258916625_zlDCrzMw0jagpUpAa8G5uTCgOjX4&c=chrY&g=joinedRmsk).
4. **Centromere** (Centro_Telo_match2.py): Data about the distance and overlap of regions X/Y with centromeres was retrieved as indicated in 2 - RRegions. Centromere position was retrieved from [UCSC genome browser’s Centromere Locations track](https://genome.ucsc.edu/cgi-bin/hgTrackUi?hgsid=1061569099_etPvsC1wDoYejDPDAUYzvBRt8AAL&g=centromeres).
5. **Telomere** (Centro_Telo_match2.py): Data about the distance and overlap of regions X/Y with telomers was retrieved as indicated in 2 - RRegions. Data was obtained from [UCSC genome browser’s Gap Locations track](https://genome.ucsc.edu/cgi-bin/hgTrackUi?hgsid=1061569099_etPvsC1wDoYejDPDAUYzvBRt8AAL&g=gap).
6. **GC content** (GCcontent2.py): The GC content ratio was calculated for regions X and Y using [UCSC genome browser’s gc content track](http://genome-euro.ucsc.edu/cgi-bin/hgTrackUi?hgsid=258917648_O9QaAJaEVNiYWFygclXVJcge5p8O&g=gc5BaseBw).
7. **Sequence Complexity** (SeqComplexity.py): Sequence Complexity (SeqComp) calculation was based in the  Dust algorithm (Morgulis et al., 2006). This algorithm calculates a score that represents the complexity of the genomic region based on the total of repeated nucleotides in a sequence fragment: higher the score, less complex is the region. 
8. **Replication time** (Replication_Time2.py): The average Replication time (RTime) of all available tissues was extracted for regions X and Y using [UCSC browser’s Replication Time track by Repli-chip from ENCODE/FSU](http://genome.ucsc.edu/cgi-bin/hgTrackUi?db=hg19&g=wgEncodeFsuRepliChip). As RTime data was only available for hg19 reference genome assembly, the coordinates were converted to GRCh38/hg38 using the LiftOver tool (Haeussler et al., 2019). 
9. **Recombination Rate** (Recombination_Rate2.py): The average of Recombination rate (RRate) values that overlaps with regions X and Y was calculated with [UCSC genome browser’s recombination rate track, "deCODE Sex-Averaged Rate score”](http://genome.ucsc.edu/cgi-bin/hgTrackUi?hgsid=1061587053_WeaKCkkGIjm9bXWbZAVn3ei9VGLx&g=recombRate). The RRate data was converted to GRCh38/hg38 genome version using the LiftOver tool (Haeussler et al., 2019).
10. **Genes** (Genes_V3.py): Gene data was extracted as the percentage of overlap of regions X and Y with exons and introns using [ensembl‘s gff anotation file for GRCh38/hg38](http://ftp.ensembl.org/pub/release-104/gtf/homo_sapiens/Homo_sapiens.GRCh38.104.chr.gtf.gz).
11. **Lamina associated domains** (Lamina_Associ_Dom.py): The Lamina associated domains (LADs) data was extracted from Guelen et al., 2018 publication. The average LAD was calculated for regions X and Y.
12. **CpG Islands** (CpG_island2.py): The distance and overlap with CpG islands (CpGIs) were calculated as described in 2 - RRegions with the [UCSC genome browser’s CpG islands track](http://genome.ucsc.edu/cgi-bin/hgTrackUi?hgsid=1061596481_5AATrBam7ZLtUehsvZmGAPZrucwZ&g=cpgIslandExt).
13. **TAD boundaries** (TADs2.py): The distance to [TAD boundaries for cell line hESC](http://3dgenome.fsm.northwestern.edu/publications.html) (Dixon et al., 2012) was measured as log10(Mb of distance +1) between regions X/Y and the closest boundary.
 
The models were trained with CNVs that were previously identified through improper pair clustering and coverage analysis (David et al., 2020), from large-insert libraries (Talkowski et al., 2011) sequencing data, using as reference the human genome version **GRCh38/hg38**. 

## How to run it
1. Download DelDupVal_ML repository
2. Download the [subdatasets](https://www.dropbox.com/s/lvbga9cnay5dwq5/dataset3.zip?dl=0) data
3. Uncompress the dataset3.zip file, it has to be in the same diretory within the repository
4. Select the bash script for the designed task, RUN_validate_DELS.sh for deletions or RUN_validate_DUPS.sh for duplications;
```
sh RUN_validate_DELS.sh <input_file> <output_path>
```
5. Run the selected bash script by inserting the input file and output path;
6. At the output directory there'll be a cvs report file with the classification results, **DEL_prediction_results** (deletions) or **DUP_prediction_results.csv** (duplications).


### Input:
. The input has to be a csv file the CNVs, each CNVs has to have their mapping infomation such as start and end position in basepairs of the CNV it self and the flanking regions made by the improper read pairs. The CNV's position *start_position_max* and *end_position_max* correspond to *regionA_end* and *regionB_stat* in case of deletions, *regionA_stat* and *regionB_end* in case of duplications. All columns of region position and size are in base paris.
```
Case_id;ID;Cluster_id;start_position_max;end_position_max;regionA_stat;regionA_end;regionB_stat;regionB_end;chr_A;chr_B;libraries;A_size;B_size
0;case_01;DEL_TP_01;null;1092351;1095174;1089917;1092351;1095174;1097820;2;2;liGS;2434;2646
1;case_01;DEL_FP_01;null;6649534;6651327;6645792;6649534;6651327;6655100;4;4;liGS;3742;3773
```

### Output:
. The output will be a csv file with the CNV id and the classification result as True or False:
```
CNV_ID;predicted_CNV
0;DEL_TP_01;True
1;DEL_FP_01;False

```

## Example:
. Run RUN_validate_DELS.sh for deletions or RUN_validate_DUPS.sh for duplications.
```
#Deletions:
sh RUN_validate_DELS.sh ./test_sample/deletions_SAMPLE.csv ./test_sample/

#Duplications:
sh RUN_validate_DUPS.sh ./test_sample/duplications_SAMPLE.csv ./test_sample/
```
### How the script works:
1. The variables **input_TP_4** and **output_path_TP_4** are the input and output variables, at the output there'll be created a folder with all output content.
2. Several python scripts will process a series of sub-datasets by matching genomic data (present in *./dataset3* folder) to the given CNVs. These are: 
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
3. All the sub-datasets will be merged and processed by **makedataset-deletions.py** (if deletions, or **makedatase.py** if duplications);4. Then after having the complete dataset (**dataset.csv**), **preditClassify.py** will classify the CNVs as possible True CNV or False CNV. The classification script need to know if the CNVs are deletions (DEL) or duplication (DUPS);
5. A *.csv* report is made at the defined output path in **output_path_TP_4** named as **DEL_prediction_results** (deletions) or **DUP_prediction_results.csv** (duplications).
  
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

In folder *train_scripts* there're the used script to train the classification models avaliable in the folder *models*.

### More documentation:
This project was part of my Master Thesis (dissertation) at the Faculty of Sciences of the University of Lisbon (FCUL).
The dissertation is available at: http://hdl.handle.net/10451/58954

## References:
- David, D., Freixo, J. P., Fino, J., Carvalho, I., Marques, M., Cardoso, M., Piña-Aguilar, R. E., & Morton, C. C. (2020). Comprehensive clinically oriented workflow for nucleotide level resolution and interpretation in prenatal diagnosis of de novo apparently balanced chromosomal translocations in their genomic landscape. Human Genetics, 139(4), 531–543. https://doi.org/10.1007/s00439-020-02121-x
- Dixon, J. R., Selvaraj, S., Yue, F., Kim, A., Li, Y., Shen, Y., Hu, M., Liu, J. S., & Ren, B. (2012). Topological Domains in Mammalian Genomes Identified by Analysis of Chromatin Interactions. Nature, 485(7398), 376. https://doi.org/10.1038/NATURE11082
- Du, Q., Bert, S. A., Armstrong, N. J., Caldon, C. E., Song, J. Z., Nair, S. S., Gould, C. M., Luu, P. L., Peters, T., Khoury, A., Qu, W., Zotenko, E., Stirzaker, C., & Clark, S. J. (2019). Replication timing and epigenome remodelling are associated with the nature of chromosomal rearrangements in cancer. Nature Communications 2019 10:1, 10(1), 1–15. https://doi.org/10.1038/s41467-019-08302-1
- Guelen, L., Pagie, L., Brasset, E., Meuleman, W., Faza, M. B., Talhout, W., Eussen, B. H., De Klein, A., Wessels, L., De Laat, W., & Van Steensel, B. (2008). Domain organization of human chromosomes revealed by mapping of nuclear lamina interactions. Nature, 453(7197), 948–951. https://doi.org/10.1038/nature06947
- Haeussler, M., Zweig, A. S., Tyner, C., Speir, M. L., Rosenbloom, K. R., Raney, B. J., Lee, C. M., Lee, B. T., Hinrichs, A. S., Gonzalez, J. N., Gibson, D., Diekhans, M., Clawson, H., Casper, J., Barber, G. P., Haussler, D., Kuhn, R. M., & Kent, W. J. (2019). The UCSC Genome Browser database: 2019 update. Nucleic Acids Research, 47(D1), D853–D858. https://doi.org/10.1093/NAR/GKY1095
- Li, Y., Roberts, N. D., Wala, J. A., Shapira, O., Schumacher, S. E., Kumar, K., Khurana, E., Waszak, S., Korbel, J. O., Haber, J. E., Imielinski, M., Akdemir, K. C., Alvarez, E. G., Baez-Ortega, A., Beroukhim, R., Boutros, P. C., Bowtell, D. D. L., Brors, B., Burns, K. H., … Campbell, P. J. (2020). Patterns of somatic structural variation in human cancer genomes. Nature, 578(7793), 112–121. https://doi.org/10.1038/s41586-019-1913-9
- Morgulis, A., Gertz, E. M., & Schäffer, A. A. (2006). A Fast and Symmetric DUST Implementation to Mask Low-Complexity DNA Sequences. JOURNAL OF COMPUTATIONAL BIOLOGY, 13(5).
- Talkowski, M. E., Ernst, C., Heilbut, A., Chiang, C., Hanscom, C., Lindgren, A., Kirby, A., Liu, S., Muddukrishna, B., Ohsumi, T. K., Shen, Y., Borowsky, M., Daly, M. J., Morton, C. C., & Gusella, J. F. (2011). Next-Generation Sequencing Strategies Enable Routine Detection of Balanced Chromosome Rearrangements for Clinical Diagnostics and Genetic Research. American Journal of Human Genetics, 88(4), 469. https://doi.org/10.1016/J.AJHG.2011.03.013

