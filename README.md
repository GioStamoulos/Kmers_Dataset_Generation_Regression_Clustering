## About
The generation of a kmers dataset that is associated with multiple genome sequences
and the further manipulation of this generated dataset are the main contents of the
current project. Specifically, psbH photosystem II protein H [1] that constitutes the
organism Coffea arabica was utilized and was set as a query in the BLAST tool-
algorithm [2]. Then, the output of the search, namely all the genes and proteins that
were defined as most relevant by BLAST, was downloaded into the fasta complete
file.

Every gene/protein was split into its unique fasta file using python scripts that
are described in Section 2, and fasta files are given as input to the count kmers
perl script to calculate the frequency of occurence of every possible combination
of k nucleotide letters. Every output of the perl count_kmers script [3] is written
in a csv that corresponds to a unique gene/protein, and finally, all csv files are
concatenated into the complete kmers dataset csv file. This complete csv file has as
headers the gene name and all the possible combinations of k nucleotide letters.

The summary kmers csv file is further utilized for regression and clustering im
plementations. Firstly, normalization of the dataset is applied in a range of 0 to
1 value. Then the k-means unsupervised learning algorithm is utilized for dataset
clustering, and then the Isomap dimensionality reduction algorithm is performed for
the 2D visualizations of the clustered dataset.

Furthermore, 2D plots of dataset labeled samples are performed using genes’
accession length and the BLAST search metrics, namely E-value, and total score,
which are included in the description file that can be downloaded from the BLAST
result output page. Previous plots are being implemented to investigate if the gene
kmers feature-based clustering could approximate the relationship between genes
depending on BLAST metrics. Examining whether the labeled samples are gath
ered in some 2D BLAST metrics dimensional space is one way to evaluate this
approximation.

The Support Vector Regression algorithm is utilized to predict the e-value and
total score values in order to have a fast algorithm for extracting useful information
from large datasets. The workflow described above is shown in Figure 1.
  
  
![pipeline_bio_project](https://user-images.githubusercontent.com/60938391/148859532-d7368803-d3e5-4e12-9105-e7485d03e6af.png)

<p align="center">
<img src=  (https://user-images.githubusercontent.com/60938391/148859841-1e0dbe60-c8d8-45a5-a984-4fcc55dd01ca.png)>
 </p>

```
Figure 1: Workflow.
```                                                   
                                                        

## Data harvesting

Two classes was created, namely create_kmers_output and
kmers_dataset for this purpose. All the implementations were developed utilizing
python language. Additionally, all the the function of these classes were utilized
through main python script. The pattern (way) to run main script through com
mand line is :

cd scripts_path

python main.py < requirement_inputs.txt

```
Algorithm 1: Bash commands to run main script for generating kmers dataset.
```
### Create_kmers_output Class

The create_kmers_output class requires the following to be defined: k_mers (inte
ger k number), script_name (perl script that counts kmers), genome_file (BLAST
fasta output), genome_folder, input_folder, and output_folder. It is composed of
two functions, namely, the fasta_generator and the folder_loop. In the fasta_generator
function,

- the BLAST search result file (its directory path) is given as input.


- The input file is added to the string type variable.
- Every record (different genome sequence) that is contained in a string variable
    is spitted out by using char ’>’ and inserted into a new list.
- List elements are looped through and written to different files whose names
    are associated with the names of current genome elements.
- The input result file contains all genome sequences that are defined after a
    BLAST search as related to the genome of interest.

In the folder_loop function,

- An output folder is created.
- All fasta files that are contained in the input folder are inserted into the
    count_kmers perl script.
- Every perl script output is written in a csv file in the output folder’s directory.

### Kmers_dataset Class

The kmers_dataset class requires the following to be defined: k_mers (integer k
number), output_folder an genome_folder. It is composed of two functions, namely,
create_csv_headers and create_dataset. In the create_csv_headers function, a
list with all possible combinations of the 4 nucleotide letters (A, C, G, and T) that
could create a string with a length equal to the k number of characters (k_mers) is
generated. In the create_dataset function,

- A list of headers is required to be defined.
- An initially empty dictionary is created that has as its keys the elements of
    the header list.
- Then, the output folder is looped through and the k_mers information is
    inserted into a temporary dictionary.
- All the temporary dictionary’s values are appended to the initially created
    dictionary.
- After the folder loop, the finally formed dictionary is inserted into the DataFrame.
- Finally, the DataFrame is saved to a csv file.

### Pipeline

The format of the file that is comprised of the nucleotide sequence of every gene
that is most related to the genome query of BLAST. The previous file is split into
fasta files that correspond to unique genes. The format of a fasta file is shown in
the figure 2.

![fasta](https://user-images.githubusercontent.com/60938391/148865934-123f677e-3140-4e41-abd4-81adddedacda.png)
```
Figure 2: A gene’s fasta file example.
```
Then every fasta file is given as input to the count_kmers script, and a csv file
is given as output with the format that is captured in the figure. Then, all output
files are concatenated and finally compiled into the kmers dataset csv file (figure 3)
with the first dimension equal to number of genes and the second equal to k^4 .

![Screenshot_4](https://user-images.githubusercontent.com/60938391/148866180-539e86bd-bdc0-4341-b6c4-444f44353cdb.png)
```
Figure 3: A kmers dataset example.
```

## Data manipulation

### Overview

After the data harvesting, the kmers dataset was utilized for regression and cluster
ing tasks. Specifically, the dataset was used for the prediction of BLAST metrics
for every gene and protein that constitutes the dataset. The clustering task focuses
on the evaluation of a k-mer analysis as a gene discrimination method that could
approximate the BLAST metrics.

### Regression

After kmers dataset generation, the values of the BLAST metrics, namely E-value
and total score, were utilized as targets in different regression implementations.The
Pandas library and the SVG algorithm from the Sklearn library are used to process
data.

Then the initial dataset is split into two parts: 60% for the train set and 40% for
the final test set.The GridSearch technique is used with five repetitions per combi
nation of parameters to estimate the best regression model for the current dataset.
Finally, the SVR model with the best grid search hyperparameters is applied to the
kmers dataset and the Root Mean Square Error (RMSE) metric is calculated for
the evaluation of the regression.


### Clustering

A function, namely dataset_clustering, was created for the clustering task. The
dataset_clustering function takes as input the path of a kmers dataset, the path to
save the clustering implementation plots, and the path of a description file containing
the BLAST metrics for each gene. In the dataset_clustering function,

- the Kmers dataset was imported.
- the dataset’s values were normalized with a range of 0 to 1.
- the sk-learn library’s Kmean algorithm is used to cluster datasets into three
    groups.
- the isomap dimensional reduction algorithm from the SK-Learn library is ap-
    plied to eliminate three dataset dimensions.
- a plot of a 2D dimensionally reduced dataset with labels of clustering is im-
    plemented utilizing the matplotlib library.
- a plot of the clustered dataset using BLAST metrics (Total score, E-value, and
    Accuracy length) as the axes is implemented.
- all plots are saved to the output path as png files.


##  Results

### Overview

The complete BLAST fasta file is constituted of 1000 unique genes. For the k value
equal to 3, 4, 5, and 6, different kmers datasets were generated, respectively. From
the BLAST description file, only the headers E-value, total score, and accession
length were used for further implementations.

### Regression

Regression tasks were applied in four generated kmers datasets, namely 3mers,
4mers, 5mers, and 6mers. For every task the regression task is focusing on pre
diction of BLAST metrics, namely either E-value or total score.

The main metric that was utilized for evaluation of every regression progress is
Root Square Mean Error(RSME). In table 1 are captured all the RSME of a every
regression experiment.


```
Table 1: RMSE regression error for different BLAST metrics targets and different kmers datasets.
```

|K-mers        |Total Score     |   E-value     |
| :---         |     :---:      |          ---: |
|3 | 0.8374 | 3.793e-96|
|4 | 0.7224 | 3.7938e-96|
|5 | 0.9948 | 3.7938e-96|
|6 | 1.1007 | 3.79383e-96|



### Clustering

In this section, all the output-plots that are associated with clustering implementa
tion are presented. There are four types of plots that are implemented. The first
type shows the distibution of the kmers dataset after kmeans clustering and 2D
isomap dimensionally reduction, using as label the output labels of clustering. The
other three plots capture the kmers clustered dataset’s distribution in dimensions
that correspond to BLAST metrics (E-value or/and total score) or/and the accession
length of every gene. The figures 4 to 7 are associated with the 3mers dataset’s clustering.

![clustering_isomap](https://user-images.githubusercontent.com/60938391/148864884-ff916aab-b8be-44e7-a5e7-2ae39ddf85f1.png)
```
Figure 4: Isomap 2D plot for 3-mers clustered dataset.
```
![Acc_length_Acc_length](https://user-images.githubusercontent.com/60938391/148864907-38d9c8cd-9edc-47e3-a4c1-4483943bd318.png)
```
Figure 5: Accession_lenght - E_value plot for 3-mers clustered dataset.
```
![E_value_Acc_length](https://user-images.githubusercontent.com/60938391/148864961-85b5617f-3adc-43ab-8f82-b7278f48bf46.png)
```
Figure 6: E_value - Total_score plot for 3-mers clustered dataset.
```
![Total_score_Acc_length](https://user-images.githubusercontent.com/60938391/148864982-26057031-db3c-4720-a507-e97a799bd35f.png)
```
Figure 7: Total_score - Accession_length plot for 3-mers clustered dataset.
```


## Conclusion

### Regression

A very low RMSE (3.7938343428779805e-96) for the E-value regression is observed,
which is encouraging. However, it is observed that the SVR result for each value
is zero (0). This leads to the conclusion that due to the nature of E-value values,
which not only get very low prices but at the same time their values are very different
from each other, it is not easy to compare them. Therefore, SVR model utilization
is impossible without proper dataset alteration.

On the other hand, for the use of the algorithm to find the total score value,
a more normal result is observed with a rather low RMSE value despite the large
initial values and their range observed in the total score. Finally, using an initially
small dataset, one could utilize the SVR algorithm, as described above, to extract
faster results for a much larger dataset.

### Clustering

For all the datasets that were utilized, their clustering coclusions are almost al
ways the same. The plots of the Isomap dimensionally reduced kmers datasets
after kmeans clustering demonstrate that clusters can be easily distinguished in this
dimensional space.One of the three clusters is remote from the other two and is
composed of fewer samples than the other two.

In the plots of clustered datasets in the BLAST metrics (E-value and Total score)
2D dimensional space, the samples are shuffled and can not be discriminate. This
fact is reasonable because the datasets were not clustered with these features. Pre
vious plots were performed to investigate the case that kmers clustering way could
approximate the hierarchical discrimination that is provided by BLAST metrics.

For plots where one axis corresponds to accession length, clusters are distin
guished and don’t overlap. This fact reveals that genes with the same accuracy
length scale are in the same cluster. Thus, k-mers dataset clustering is indirectly
influenced by the lengths of genes’ nucleotide sequences.


## References

[1] https://www.ncbi.nlm.nih.gov/gene/?term=4421757.

[2] https://blast.ncbi.nlm.nih.gov/Blast.cgi.

[3] https://github.com/lindberg-m/count_kmers.


