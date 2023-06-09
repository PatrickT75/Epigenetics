# -*- coding: utf-8 -*-
"""Pacemaker 4.ipynb

Automatically generated by Colaboratory.

Original file is located at
    https://colab.research.google.com/drive/1Y6EjqDRdquybH-4Qqnmr-WtaQUrIcQ1M
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import math
import statistics

from sklearn.model_selection import train_test_split
from sklearn.model_selection import LeaveOneOut
from sklearn.linear_model import LassoCV
from sklearn.linear_model import LassoLars
from sklearn.linear_model import LassoLarsCV

from matplotlib import rc
from scipy import optimize
import scipy.stats as stats

#### Data Cleanup

CG_table = pd.read_csv("methylation_used.csv", header=None)
CG_table = CG_table.drop(columns = 0)

isna = CG_table.isna().any()
ones = np.ones((CG_table.shape[1]))
keep = ones-isna
keep = keep.to_numpy()
keep = keep.astype(bool)

CG_table = CG_table.dropna(axis = 'columns')


#### Split features into dogs and wolves
UCLA_CG_table = CG_table[0:147]
Wolves_CG_table = CG_table[147:]


#### Process traits
traits_short = pd.read_csv("traits_used.csv", header=None)

UCLA_traits = traits_short[0:147]
Wolves_traits = traits_short[147:]

name_temp = UCLA_traits.loc[:,0]
age_temp = UCLA_traits.loc[:,3]
sex_temp = UCLA_traits.loc[:,4]
spayed_temp = UCLA_traits.loc[:,6]
weight_temp = UCLA_traits.loc[:,7]
hetero_temp = UCLA_traits.loc[:,9]
FOH_temp = UCLA_traits.loc[:,10]
pred_ages_temp = UCLA_traits.loc[:,12]
wolves_age_temp = Wolves_traits.loc[:,3]

names = []
ages = []
sexes = []
spayeds = []
weights = []
heteros = []
FOHs = []
wolves_ages = []

for i in range(len(name_temp)):
    names.append(name_temp[i])
    ages.append(age_temp[i])
    sexes.append(sex_temp[i])
    spayeds.append(spayed_temp[i])
    weights.append(weight_temp[i])
    heteros.append(hetero_temp[i])
    FOHs.append(FOH_temp[i])
    
for i in range(len(wolves_age_temp)):
    wolves_ages.append(wolves_age_temp[i+147])

names = np.array(names)
ages = np.array(ages)
sexes = np.array(sexes)
weights = np.array(weights)
spayeds = np.array(spayeds)
heteros = np.array(heteros)
FOHs = np.array(FOHs)
wolves_ages = np.array(wolves_ages)

#### Keep information on samples that are part of triplicates
triplicates = [[25,26,54],[27,28,55],[29,30,56],[31,32,57],\
               [33,34,58],[35,36,59],[37,38,60],[39,40,61],\
               [41,42,62],[43,44,63],[45,46,64],[47,48,65],\
               [49,50,66],[52,53,146]]


#### If coverage data provided:
# coverages = pd.read_csv("dog_counts_217.csv", header=None)
# coverages_short = coverages[91:]
# coverages_short = coverages_short.drop(columns = 0)
# coverages_short = coverages_short.to_numpy()

# coverages_short_keep = coverages_short[:,keep]

# coverages_mean = coverages_short_keep.mean(axis = 0)

#### Parser function to select only proteins in the top percentile of counts
def parser(coverage_means, quantile):
    coverage_means_q = statistics.quantiles(coverage_means, n=quantile)[quantile-2]
    best = coverage_means > coverage_means_q
    return best

# best = parser(coverages_mean, 3)

methylation = UCLA_CG_table.to_numpy()
IDs = methylation[:,0]
methylation = np.delete(methylation, 0, 1)
# methylation = methylation[:,best]


########################################
#### LARS Prediction

def get_predictions(X_train, Y_train):
    loo = LeaveOneOut()
    predicted_ages = []

    for train_index, test_index in loo.split(X_train):
        X_retrain, X_retest = X_train[train_index], X_train[test_index]
        
        IDcheck = np.ones(len(IDs))
        dupIDs = IDs
        for i in range(len(IDcheck)):
            if IDs[test_index[0]] == dupIDs[i]:
                IDcheck[i] = 0
        IDcheck = np.delete(IDcheck, test_index[0])
        
        X_retrain2 = X_retrain[IDcheck.astype(bool)]
        print(X_retrain2.shape)
        
        Y_retrain, Y_retest = Y_train[train_index], Y_train[test_index]
        Y_retrain2 = Y_retrain[IDcheck.astype(bool)]
        
        lars = LassoLars(alpha = 0.001)
        lars.fit(X_retrain2, Y_retrain2)
        predicted_ages.append(lars.predict(X_retest)[0])
    
    return predicted_ages

#### Try predicting either ages or sqrt-transformed ages
# lars_predict = get_predictions(methylation, ages)
lars_predict_sqrt = get_predictions(methylation, np.sqrt(ages))

#### Functions to calculate prediction quality
def get_median_abs_error(prediction, actual):
    errors = []
    for i in range(len(prediction)):
        difference = np.abs(prediction[i] - actual[i])
        errors.append(difference)
    return statistics.median(errors)

def get_r_value(predicted_ages_list, actual_ages_list):
    r_value = np.corrcoef(actual_ages_list, predicted_ages_list)[0,1]
    return r_value

def get_r_squared(predicted_ages_list, actual_ages_list):
    return (get_r_value(predicted_ages_list, actual_ages_list)**2)

#### Plot predictions
def plot_predictions(predicted_ages_list, actual_ages_list):
    plt.plot(actual_ages_list, predicted_ages_list, 'ko')
    plt.ylabel("Predicted Age")
    plt.xlabel("Actual Age")
    plt.title("Lasso LARS Model of Predicted Age (Square Root Transform)")
    plt.text(8, 2, "Median Absolute Error: 0.6522")
    plt.text(8, 1, "R squared: 0.8939")
    plt.plot(actual_ages_list, actual_ages_list, 'r--')
    plt.savefig('Predicted_vs_actual_lars_sqrt_unduplicated.png', dpi = 300)
    plt.show()
    print("R squared: ", get_r_squared(predicted_ages_list, actual_ages_list))
    print("R value: ", get_r_value(predicted_ages_list, actual_ages_list))
    print("Median Absolute Error: ", get_median_abs_error(predicted_ages_list, actual_ages_list))

plot_predictions(np.square(lars_predict_sqrt), ages)


####################################
#### Predict wolves (i.e. test set) using LARS model

wolves_methylation = Wolves_CG_table

wolves_predict = LassoLars(alpha = 0.02)
wolves_predict.fit(methylation, ages)
predictions = wolves_predict.predict(wolves_methylation)

plt.plot(wolves_ages, predictions, 'ko')
plt.plot(wolves_ages, wolves_ages, 'r--')
plt.title('Predicting Wolves Using Lasso LARS Model')
plt.xlabel('Actual Age')
plt.ylabel('Predicted Age')
plt.savefig('Wolves_predictions.png', dpi=300)

print("R squared: ", get_r_squared(predictions, wolves_ages))
print("R value: ", get_r_value(predictions, wolves_ages))
print("Median Absolute Error: ", get_median_abs_error(predictions, wolves_ages))
#############################
#### Show boxplot of predicted ages by triplicates

triplicate_ranges = {1:[], 2:[], 3:[], 4:[], 5:[], 6:[], 7:[],\
                     8:[], 9:[], 10:[], 11:[], 12:[], 13:[]}
for x in triplicate_ranges:
    temp = []
    for y in triplicates[x]:
        temp.append(np.square(lars_predict_sqrt[y]))
    triplicate_ranges[x] = temp
    
fig, ax = plt.subplots()
ax.boxplot(triplicate_ranges.values())
ax.set_xticklabels(triplicate_ranges.keys())
ax.title.set_text('Prediction Ranges per Individual')
ax.set_ylabel('Predicted Age')
ax.set_xlabel('Individual')
# plt.savefig('Prediction_ranges.png', dpi = 300)

triplicate_range_for_median = []
for x in triplicate_ranges:
    triplicate_range_for_median.append(np.std(triplicate_ranges[x]))
    
triplicate_range_for_median = np.array(triplicate_range_for_median)
print("Median of Standard Deviation for Triplicates:", statistics.median(triplicate_range_for_median))


###############################
#### Different approach: Epigenetic Pacemaker Prediction

from EpigeneticPacemaker.EpigeneticPacemaker import EpigeneticPacemaker

def pearson_correlation(meth_matrix: np.array, phenotype: np.array) -> np.array:
    """calculate pearson correlation coefficient between rows of input matrix and phenotype"""
    # calculate mean for each row and phenotype mean
    matrix_means = np.mean(meth_matrix, axis=0)
    phenotype_mean = np.mean(phenotype)

    # subtract means from observed values
    transformed_matrix = np.transpose(meth_matrix) - matrix_means.reshape([-1,1])
    transformed_phenotype = phenotype - phenotype_mean

    # calculate covariance
    covariance = np.sum(transformed_matrix * transformed_phenotype, axis=1)
    sums = np.sum(transformed_matrix ** 2, axis=1)
    variance_meth = []
    for i in sums:
        variance_meth.append(np.sqrt(i))
    variance_meth = np.array(variance_meth)
    variance_phenotype = np.sqrt(np.sum(transformed_phenotype ** 2))

    return covariance / (variance_meth * variance_phenotype)

# get the absolute value of the correlation coefficient
abs_pcc_coefficients = abs(pearson_correlation(methylation, ages))

# return list of site indices with a high absolute correlation coefficient
training_sites = np.where(abs_pcc_coefficients > .6)[0]

print(training_sites)

from matplotlib import rc
from scipy import optimize
import scipy.stats as stats

#### Plot graph for pacemaker prediction using 3 different functions (square root, linear, log)
def plot_known_predicted_ages(known_ages, predicted_ages, label=None):
    # define optimization function
    def func(x, a, b, c):
        return a * np.asarray(x)**0.5 + c
    def func_2(x, a, b, c):
        return a * np.asarray(x) + c
    def func_3(x, a, b, c):
        return a * np.log(np.asarray(x)) + c
    
    # fit trend line
    popt, pcov = optimize.curve_fit(func, [1 + x for x in known_ages], predicted_ages)
    popt_2, pcov_2 = optimize.curve_fit(func_2, [1 + x for x in known_ages], predicted_ages)
    popt_3, pcov_3 = optimize.curve_fit(func_3, [1 + x for x in known_ages], predicted_ages)
    
    # get r squared
    rsquared = r2(predicted_ages, func([1 + x for x in known_ages], *popt))
    rsquared_2 = r2(predicted_ages, func_2([1 + x for x in known_ages], *popt_2))
    rsquared_3 = r2(predicted_ages, func_3([1 + x for x in known_ages], *popt_3))
    
    # format plot label
    plot_label = f'$f(x)={popt[0]:.2f}x^{{1/2}} {popt[2]:.2f}, R^{{2}}={rsquared:.2f}$'
    plot_label_2 = f'$f(x)={popt_2[0]:.2f}x {popt_2[2]:.2f}, R^{{2}}={rsquared_2:.2f}$'
    plot_label_3 = f'$f(x)={popt_3[0]:.2f}ln(x) {popt_3[2]:.2f}, R^{{2}}={rsquared_3:.2f}$'
    
    # initialize plt plot
    fig, ax = plt.subplots(figsize=(12,12))
    # plot trend line
    ax.plot(sorted(known_ages), func(sorted([1 + x for x in known_ages]), *popt), 'r--', label=plot_label)
    ax.plot(sorted(known_ages), func_2(sorted([1 + x for x in known_ages]), *popt_2), 'b--', label=plot_label_2)
    ax.plot(sorted(known_ages), func_3(sorted([1 + x for x in known_ages]), *popt_3), 'g--', label=plot_label_3)
    
    # scatter plot
    ax.scatter(known_ages, predicted_ages, marker='o', alpha=0.8, color='k')
    ax.set_title(label, fontsize=18)
    ax.set_xlabel('Chronological Age', fontsize=16)
    ax.set_ylabel('Epigenetic State', fontsize=16)
    ax.tick_params(axis='both', which='major', labelsize=16)
    ax.legend(fontsize=16)
    plt.show()
    
# use latex formatting for plots
rc('text', usetex=False)

def r2(x,y):
    # return r squared
    return stats.pearsonr(x,y)[0] **2

# initialize the EPM model 
epm = EpigeneticPacemaker(iter_limit=100, error_tolerance=0.00001)

table = methylation
short_table = np.transpose(table[:,training_sites])
print("Short Table: ", short_table.shape)

# fit the model using the training data
epm.fit(short_table, ages)

# generate predicted ages using the test data
pacemaker_state_predict = epm.predict(short_table)

plot_known_predicted_ages(ages, pacemaker_state_predict, 'Epigenetic Pacemaker')
# Predict square root of age
# Put in table - both predictions + actual


####################################
#### Testing for Moderation and Interaction 

import statsmodels.api as sm

sexes_int = []
weights_int = []

male_num = 0
female_num = 0
other_num = 0

for i in range(len(sexes)):
    if sexes[i][0] == 'M':
        sexes_int.append(0)
        male_num += 1
    elif sexes[i][0] == 'F':
        sexes_int.append(1)
        female_num += 1
    else:
        print(i)
        other_num += 1
        
for weight in weights:
    weights_int.append(float(weight))

sexes_int = np.array(sexes_int)
print(sexes_int.shape)
print("Males: ", male_num)
print("Females: ", female_num)
print("Others: ", other_num)

spayeds_int = []

for j in range(len(spayeds)):
    if spayeds[j] == 'nan':
        spayeds_int.append(0)
    elif spayeds[j][0] == 'Y' or spayeds[j][0] == 'y':
        spayeds_int.append(1)
    elif spayeds[j][1] == 'E' or spayeds[j][1] == 'e':
        spayeds_int.append(1)
    elif spayeds[j][0] == 'S' or spayeds[j][0] == 's':
        spayeds_int.append(1)
    else:
        spayeds_int.append(0)
        
spayeds_int = np.array(spayeds_int)

# age_sexes = ages * sexes_int
age_weights = ages * weights_int
# age_spayed = ages * spayeds_int
# age_hetero = ages * heteros
# age_FOH = ages * FOHs

stack = np.stack((ages, age_weights), axis = 1)

# Instead of age-squared, use square root of age
stack_df = pd.DataFrame(stack, columns = ['Age', 'Age*Weight'])

#### Generate plots for moderation analysis

lars_sqrt_model = sm.OLS(lars_predict, stack_df)
lars_sqrt_fit = lars_sqrt_model.fit()
lars_sqrt_fit.summary()

plt.rc('figure', figsize=(12, 7))
#plt.text(0.01, 0.05, str(model.summary()), {'fontsize': 12}) old approach
plt.text(0.01, 0.05, str(lars_sqrt_fit.summary()), {'fontsize': 10}, fontproperties = 'monospace') # approach improved by OP -> monospace!
plt.axis('off')
plt.tight_layout()
# plt.savefig('Lars_Transform_Statsmodels_only.png', dpi=300)

# Put this in docs (export figure), and fix the labeling

lars_model = sm.OLS(lars_predict, stack_df)
lars_fit = lars_model.fit()
lars_fit.summary()

plt.rc('figure', figsize=(12, 7))
#plt.text(0.01, 0.05, str(model.summary()), {'fontsize': 12}) old approach
plt.text(0.01, 0.05, str(lars_fit.summary()), {'fontsize': 10}, fontproperties = 'monospace') # approach improved by OP -> monospace!
plt.axis('off')
plt.tight_layout()
plt.savefig('Lars_Statsmodels_no_sqrt.png')

pacemaker_model = sm.OLS(pacemaker_state_predict, stack_df)
pacemaker_fit = pacemaker_model.fit()
pacemaker_fit.summary()

plt.rc('figure', figsize=(12, 7))
#plt.text(0.01, 0.05, str(model.summary()), {'fontsize': 12}) old approach
plt.text(0.01, 0.05, str(pacemaker_fit.summary()), {'fontsize': 10}, fontproperties = 'monospace') # approach improved by OP -> monospace!
plt.axis('off')
plt.tight_layout()
plt.savefig('Pacemaker_Statsmodels.png')
