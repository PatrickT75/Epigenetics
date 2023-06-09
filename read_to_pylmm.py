# -*- coding: utf-8 -*-
"""Read to PyLMM.ipynb

Automatically generated by Colaboratory.

Original file is located at
    https://colab.research.google.com/drive/1t5ZQ0CmyIvt6X1zkmSae9QvNMdk-oxvg
"""

import numpy as np
import matplotlib.pyplot as plt
from lmm import *

def get_table(file):
    m_list = []
    f = open(file, 'r')
    f.readline()
    for line in f:
        temp = []
        line_split = line.split(' ')
        for term in line_split:
            temp.append(float(term))
        m_list.append(temp)

    np_list = np.array(m_list)
    return np_list

def get_phenotypes():
    # Fill as necessary...
    phenotypes = [25.76, 28.09, 23.01, 25.78, 21.18, 24.76, 25.93, 22.73, 26.73, 22.19, 26.23, 30.57, 22.63, 23.39, 25.91, \
                  27.6, 25.51, 25.09, 30.64, 22.6, 20.31, 27.8, 25.49, 28.88, 30, 24.87, 30.91, 22.18, 27.02, 24.59, \
                 24.46, 25.14, 27.37, 27.76, 30.06, 36.22, 32.77, 29.41, 22.13, 27.62, 26.64, 25.56, 20.95, 30.12, 28.67, \
                 21.74, 28.57, 21.67]
    phenotypes = np.array(phenotypes)
    return phenotypes

def transpose(np_list):
    transposed = np.transpose(np_list)
    return transposed

def transform_p(p, threshold):
    probable_dict = {}
    for i in range(len(p)):
        p[i] = np.abs(p[i])
        if p[i] < threshold:
            probable_dict[i] = p[i]

    print(probable_dict)
    neg_log_p = np.log10(np.array(p))*-1
    
    return neg_log_p

def plot_results(transformed, threshold):
    thresh = threshold
    boundary = np.full((len(transformed), 1), thresh)
    plt.plot(range(len(transformed)), transformed, 'bo')
    plt.plot(range(len(transformed)), boundary, 'r-')
    plt.title("Chromosome 21 with 5 kb Regions")
    plt.savefig("Chr_21_5kb.png")

"""**Running the Code**"""

transpose_table = transpose(get_table("Chr_21_full_5kb.txt"))
phenotypes = get_phenotypes()
kinship = calculateKinship(transpose_table)
CG_table = get_table("Chr_21_full_5kb.txt")

print(transpose(CG_table))

###

t_vals, p_vals = GWAS(phenotypes, transpose(CG_table), kinship)

transformed = transform_p(p_vals, 10**-(3.0))

plot_results(transformed, 3.0)

