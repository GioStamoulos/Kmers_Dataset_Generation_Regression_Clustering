import matplotlib.pyplot as plt
from sklearn.cluster import KMeans
import pandas as pd
from sklearn.manifold import Isomap
import os
from sklearn import preprocessing

def dataset_clustering(dataset_path, output_file, description_file):
  #import dataset
  df1 = pd.read_csv(dataset_path)
  df_n_genome = list(df1['Genome_Name'])
  df1.pop('Genome_Name')

  #normalize data
  
  x = df1.values #returns a numpy array
  min_max_scaler = preprocessing.MinMaxScaler()
  x_scaled = min_max_scaler.fit_transform(x)
  df1 = pd.DataFrame(x_scaled)
  print(df1)

  #apply k-means clustering - 3 number of clusters
  kmeans = KMeans(n_clusters=3).fit(df1)
  kmeans_labels = kmeans.labels_

  # Dimensionality reduction (2D) using Isomap
  embedding = Isomap(n_components=2)
  X_transformed = embedding.fit_transform(df1)
  X_transformed = pd.DataFrame(X_transformed)

  #normalize
  X_transformed[0] = X_transformed[0]/max(X_transformed[0])
  X_transformed[1] = X_transformed[1]/max(X_transformed[1])
  
  #2D clustering plot
  plt.figure(figsize=(10,10))
  plt.scatter(X_transformed[0][kmeans_labels== 0], X_transformed[1][kmeans_labels==0], s=10, c='red')
  plt.scatter(X_transformed[0][kmeans_labels== 1], X_transformed[1][kmeans_labels==1], s=10, c='green')
  plt.scatter(X_transformed[0][kmeans_labels== 2], X_transformed[1][kmeans_labels==2], s=10, c='blue')
  plt.savefig('{}/clustering_isomap.png'.format(output_file))
  plt.show()

  #Match gene products with blast components (total score and E-value) by utilizing Accession column
  df2 = pd.read_csv(description_file, header=0)
  total_score = []
  E_value = []
  acc_length = []
  for i in df_n_genome:
    accession_index = df2.index[df2['Accession  '] == i].tolist()
    total_score.append(df2['Total Score'][accession_index[0]])
    E_value.append(df2['E value'][accession_index[0]])
    acc_length.append(df2['Acc. Len'][accession_index[0]])
  
  #insert to dataframe
  data = {'Total_score':total_score, 'E_value':E_value, 'Acc_length': acc_length}
  data2d = pd.DataFrame(data)

  #2D plot using blast metrics 
  blast_metrics= ['Total_score', 'E_value', 'Acc_length']

  for i in range(3):
    plt.figure(figsize=(11,11))
    plt.scatter(y=data2d[blast_metrics[i]][kmeans_labels== 0], x=data2d[blast_metrics[i-1]][kmeans_labels==0], s=10, c='red',  label = 'Cluster 1')
    plt.scatter(y=data2d[blast_metrics[i]][kmeans_labels==1], x=data2d[blast_metrics[i-1]][kmeans_labels==1],  c='green', s=10, label = 'Cluster 2')
    plt.scatter(y=data2d[blast_metrics[i]][kmeans_labels==2], x=data2d[blast_metrics[i-1]][kmeans_labels==2],  c='blue', s=10, label = 'Cluster 3')
    plt.ylabel(blast_metrics[i], fontsize=25)
    plt.xlabel(blast_metrics[i-1], fontsize=25) 
    plt.savefig('{0}/{1}_{2}'.format(output_file, blast_metrics[i], blast_metrics[-1]))


def main():

    #6_mers dataset
    dataset_clustering('/content/drive/MyDrive/Bio_project/NC_008535.1/6/kmers_dataset.csv', '/content/drive/MyDrive/Bio_project/NC_008535.1/6', '/content/drive/MyDrive/Bio_project/NC_008535.1/X33Z0KDD01R-Alignment-Descriptions.csv')
    
    #5_mers dataset
    dataset_clustering('/content/drive/MyDrive/Bio_project/NC_008535.1/5/kmers_dataset.csv', '/content/drive/MyDrive/Bio_project/NC_008535.1/5', '/content/drive/MyDrive/Bio_project/NC_008535.1/X33Z0KDD01R-Alignment-Descriptions.csv' )
    
    #4_mers dataset
    dataset_clustering('/content/drive/MyDrive/Bio_project/NC_008535.1/4/kmers_dataset.csv', '/content/drive/MyDrive/Bio_project/NC_008535.1/4', '/content/drive/MyDrive/Bio_project/NC_008535.1/X33Z0KDD01R-Alignment-Descriptions.csv')
    
    #3_mers dataset
    dataset_clustering('/content/drive/MyDrive/Bio_project/NC_008535.1/3/kmers_dataset.csv', '/content/drive/MyDrive/Bio_project/NC_008535.1/3', '/content/drive/MyDrive/Bio_project/NC_008535.1/X33Z0KDD01R-Alignment-Descriptions.csv')

if __name__ == "__main__":
    main()