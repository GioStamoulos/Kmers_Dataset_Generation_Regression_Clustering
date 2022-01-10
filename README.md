# Kmers_Dataset_Generation-Clustering
 
**Overview**

The generation of a k-mer dataset that is associated with multiple genome sequences
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
  The summary kmers csv file is further utilized for regression and clustering im-
plementations. Firstly, normalization of the dataset is applied in a range of 0 to
1 value. Then the k-means unsupervised learning algorithm is utilized for dataset
clustering, and then the Isomap dimensionality reduction algorithm is performed for
the 2D visualizations of the clustered dataset.
  Furthermore, 2D plots of dataset labeled samples are performed using genes’
accession length and the BLAST search metrics, namely E-value, and total score,
which are included in the description file that can be downloaded from the BLAST
result output page. Previous plots are being implemented to investigate if the gene
kmers feature-based clustering could approximate the relationship between genes depending on BLAST metrics. Examining whether the labeled samples are gath-
ered in some 2D BLAST metrics dimensional space is one way to evaluate this
approximation.
  The Support Vector Regression algorithm is utilized to predict the e-value and
total score values in order to have a fast algorithm for extracting useful information
from large datasets. The workflow described above is shown in Figure 1.
