import pandas as pd
from sklearn.metrics import mean_squared_error
import math
from sklearn.model_selection import GridSearchCV
import sklearn
from sklearn.svm import SVR, LinearSVR
from sklearn.linear_model import SGDRegressor
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import MinMaxScaler

def blast_metric_regression(X, blast_metric):
    # Split data using train_test_split
    x_train, x_test, y_train, y_test = train_test_split(X, blast_metric, test_size = 0.4, random_state = 0)  

    #normalize data
    scaler = MinMaxScaler()
    x_train, x_test = scaler.fit_transform(x_train), scaler.transform(x_test)

    #grid search
    param_grid = {'C': [0.1, 1, 10, 100], 'gamma': [0.1, 1, 10, 100],'kernel': ['rbf', 'linear']}
    grid = GridSearchCV(sklearn.svm.SVR(),param_grid,refit=True,verbose=2)
    grid.fit(x_train,y_train) 
    Best = grid.best_estimator_
    print(grid.best_estimator_)

    #SVR regression
    clf = grid.best_estimator_
    rbf_pred2 = clf.fit(x_test,y_test)
    y_predrbf2 = clf.predict(x_test)

    #calculate Root Mean Square Error
    MSE = mean_squared_error(y_test, y_predrbf2)
    RMSE = math.sqrt(MSE)
    print("Root Mean Square Error:\n")
    print(RMSE)  
  

def kmers_regression(dataset_path, description_path):
    # import dataset
    df1 = pd.read_csv(dataset_path)
    df_n_genome = list(df1['Genome_Name'])
    df1.pop('Genome_Name')

    # import description file
    df2 = pd.read_csv(description_path, header=0)
    
    total_score = []
    E_value = []
    for i in df_n_genome:
      accession_index = df2.index[df2['Accession  '] == i].tolist()
      total_score.append(df2['Total Score'][accession_index[0]])
      E_value.append(df2['E value'][accession_index[0]])

    X = df1
    blast_metric_regression(X, total_score)
    blast_metric_regression(X, E_value)


def main():

    kmers_regression('./NC_008535.1/3/kmers_dataset.csv','./NC_008535.1/X33Z0KDD01R-Alignment-Descriptions.csv')
    kmers_regression('./NC_008535.1/4/kmers_dataset.csv','./NC_008535.1/X33Z0KDD01R-Alignment-Descriptions.csv')
    kmers_regression('./NC_008535.1/5/kmers_dataset.csv','./NC_008535.1/X33Z0KDD01R-Alignment-Descriptions.csv')
    kmers_regression('./NC_008535.1/6/kmers_dataset.csv','./NC_008535.1/X33Z0KDD01R-Alignment-Descriptions.csv')

if __name__ == "__main__":
    main()