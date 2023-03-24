# -*- coding: utf-8 -*-
"""Smoking Moderation.ipynb

Automatically generated by Colaboratory.

Original file is located at
    https://colab.research.google.com/drive/13eY5OsLkxC1ONWwYXKnC9aE9gcjKnY_i
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import math
import statistics
import statsmodels.api as sm

file = pd.read_csv("MetadataHuman2021_554.csv")

file

civic = file[224:286]

civic_ages = civic['age'].to_numpy()
civic_smoker = civic['smoker'].to_numpy()
civic_ages_smoker = civic_ages * civic_smoker

civic_stack = np.stack((civic_ages, civic_smoker, civic_ages_smoker), axis = 1)
civic_stack_df = pd.DataFrame(stack, columns = ['Age', 'Smoker', 'Age*Smoker'])

print(civic_stack_df)
print(civic_ages)

civic_smoker_model = sm.OLS(civic_ages, civic_stack_df)
civic_smoker_fit = civic_smoker_model.fit()
civic_smoker_fit.summary()

saliva = file[365:461]

saliva_ages = saliva['age'].to_numpy()
saliva_smoker = saliva['smoker'].to_numpy()
saliva_ages_smoker = saliva_ages * saliva_smoker

print(len(saliva_smoker))
print(len(saliva_ages))
print(len(saliva_ages_smoker))

saliva_stack = np.stack((saliva_ages, saliva_smoker, saliva_ages_smoker), axis = 1)
saliva_stack_df = pd.DataFrame(saliva_stack, columns = ['Age', 'Smoker', 'Age*Smoker'])

print(saliva_stack_df)
print(len(saliva_ages))

saliva_smoker_model = sm.OLS(saliva_ages, saliva_stack_df)
saliva_smoker_fit = saliva_smoker_model.fit()
saliva_smoker_fit.summary()

n_groups = 2
smoker = (-4.619e-14, -7.105e-15)
smokerp = (0.462, 0.575)
smokerstars = (' ',' ')
smokermod = (-5.551e-17, 4.441e-16)
smokerpmod = (0.968, 0.016)
smokermodstars = (' ', '*')

figsmok, axsmok = plt.subplots()
index = np.arange(n_groups)
bar_width = 0.35
opacity = 0.8

smok1 = plt.bar(index, smoker, bar_width,
alpha=opacity,
color='b',
label='Smoker')

smok2 = plt.bar(index + bar_width, smokermod, bar_width,
alpha=opacity,
color='g',
label='Smoker * Age')

axsmok.set_ylabel('Coefficient')
axsmok.set_title('Smoker')
axsmok.set_xticks(index + 0.35/2)
axsmok.set_xticklabels(('Civic', 'Saliva'))
axsmok.legend()

axsmok.bar_label(smok1, labels = smokerstars)
axsmok.bar_label(smok2, labels = smokermodstars)

# figbuc.savefig('Buccal Moderators.png', dpi=300)



