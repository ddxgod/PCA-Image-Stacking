# import functions needed for GBM model
from sklearn.model_selection import KFold
from sklearn.linear_model import LogisticRegression
from sklearn.ensemble import RandomForestClassifier
import pandas as pd
import numpy as np
from sklearn.ensemble import GradientBoostingClassifier  
from sklearn import metrics
from sklearn.model_selection import cross_validate
from sklearn.model_selection import cross_val_score
from sklearn.metrics import accuracy_score
from sklearn.grid_search import GridSearchCV   #Performing grid search
import matplotlib.pyplot as plt

companies = ['AAPL', 'MSFT', 'AMZN', 'INTC', 'TSLA']
# all_data = dict.fromkeys(companies, pd.DataFrame())

#first two dataframes are the X_train and y_train
main_data = dict.fromkeys(companies, pd.DataFrame)
X_split_data = dict.fromkeys(companies, pd.DataFrame)
X_train_data = dict.fromkeys(companies, pd.DataFrame)
y_train_data = dict.fromkeys(companies, [])
X_test_data = dict.fromkeys(companies, pd.DataFrame)
y_test_data = dict.fromkeys(companies, [])


split_percentagesplit_p  = 0.8
def init_data(companies):
#     pd.DataFrame(),pd.DataFrame(), pd.DataFrame(),[],pd.DataFrame(),
    for comp in companies:
        main_data[comp] = pd.read_csv(f'./rsc/{comp}_data.csv', sep = ',')
        main_data[comp].dropna()
        y = np.where(main_data[comp]['close'].shift(-1) > main_data[comp]['close'],1,-1)
        main_data[comp]['Open-Close'] = main_data[comp].open - main_data[comp].close
        main_data[comp]['High-Low'] = main_data[comp].high - main_data[comp].low
        X = main_data[comp][main_data[comp].columns[1:]]
        split = int(split_percentage*len(main_data[comp]['date']))
        X_split_data[comp] = main_data[comp][split:]
        X_train_data[comp] = X[:split]
        y_train_data[comp] = y[:split]
        X_test_data[comp] = X[split:]
        y_test_data[comp] = y[split:]
        
init_data(companies)


# optimum parameters found from ML_Model_Training
parameter_dict = {
    'max_depth': 7,
    'min_samples_split': 420,
    'min_samples_leaf': 30,
    'subsample': 0.8,
    'random_state': 10,
    'max_features': 11
}
predictors = X_train_data['AAPL'].columns

print(predictors)
