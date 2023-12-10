# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .py
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.16.0
#   kernelspec:
#     display_name: oktoberfest
#     language: python
#     name: python3
# ---

import pandas as pd
import matplotlib.pyplot as plt
from matplotlib_venn import venn2
from matplotlib_venn import venn3, venn3_circles
import re
import numpy as np

#import files
A1101_psm_OK_MSF = pd.read_csv('/home/sara/HLA/OK_MSfragger_filtered/filtered_target_psm_A1101.csv')
A1101_psm_MSB = pd.read_csv('/home/sara/HLA/MSB Concat_and_Filtered/MSB_A1101_filtered.csv')

A1101_psm_OK_MSF.head(15)

A1101_psm_MSB.tail(15)

# +
venn2(
      [set(A1101_psm_OK_MSF['modi_peptide'].to_list()), 
       set(A1101_psm_MSB['modi_peptide_2'].to_list())],
       set_labels=('Oktoberfest', 'MSB'),
       set_colors = ('yellow', 'skyblue'),
       alpha = 0.5
     )

plt.title("Oktoberfest vs MSBooster: MSFragger Search Engine (1% FDR on the psm level)")
# -

# # Stacked bar plot

#SE = pd.read_csv ('/home/sara/Metaproteomics/forKarim/sprot_all/sprot_all_Oktoberfest_MSFragger/msms.txt')
percolator_output_andro = pd.read_csv('/home/sara/Metaproteomics/forKarim/sprot_all/sprot_all_Oktoberfest_MSFragger/results/percolator/original_target.psms', sep= '\t')
percolator_output_prosit = pd.read_csv("/home/sara/Metaproteomics/forKarim/sprot_all/sprot_all_Oktoberfest_MSFragger/results/percolator/rescore_target.psms", sep= '\t')


percolator_output_andro

# +
#add a new column where you will clean the peptide column in it / leave the original as it is 
percolator_output_andro["modi_peptide"] = percolator_output_andro['peptide']

#rename the dataframe 
filt_andro = percolator_output_andro

#clean the modi_peptide column
filt_andro['modi_peptide'] = filt_andro['modi_peptide'].str.replace(r'[_,.:,\d]+', '', regex=True)
filt_andro['modi_peptide'] = filt_andro['modi_peptide'].str.replace(r'\[UNIMOD\]', '', regex= True)

filt_andro
# -

percolator_output_prosit

# +
#add new column to clean peptides in it
percolator_output_prosit["modi_peptide"] = percolator_output_prosit['peptide']

#rename the dataframe 
filt_prosit = percolator_output_prosit
filt_prosit
# -

#clean the modi_peptide column
filt_prosit['modi_peptide'] = filt_prosit['modi_peptide'].str.replace(r'[_,.:,\d]+', '', regex=True)
filt_prosit['modi_peptide'] = filt_prosit['modi_peptide'].str.replace(r'\[UNIMOD\]', '', regex= True)
filt_prosit['pep_length'] = filt_prosit['modi_peptide'].str.len()
filt_prosit

# +
#filtered to 1% FDR
filt_prosit = filt_prosit[filt_prosit['q-value'] <= 0.01]

filt_andro = filt_andro[filt_andro['q-value']<=0.01]
# -

# unique peptides to the andromeda_target (percolator target output for rescoring without peptide property prediction) 
lost = ~filt_andro['modi_peptide'].isin(filt_prosit['modi_peptide'])
lost

# +
lost_count = lost.value_counts()

# Access the counts
true_count = lost_count.get(True, 0)
false_count = lost_count.get(False, 0)

print("True count:", true_count)
print("False count:", false_count)
# -

#unique peptides to the prosit_target (percolator target output for rescoring with peptide property prediction)
gained = ~filt_prosit['modi_peptide'].isin(filt_andro['modi_peptide'])
gained

# +
gained_count = gained.value_counts()

# Access the counts
true_count = gained_count.get(True, 0)
false_count = gained_count.get(False, 0)

print("True count:", true_count)
print("False count:", false_count)
# -

# common peptides
common = filt_prosit.merge(filt_andro, how="inner", on=['modi_peptide','PSMId'])
common



shared = len(common.index)
gained = len(gained.index) - shared
lost = len(lost.index) - shared

shared

gained

lost

# +
# geeks for geeks
# # create data
x = ['A']
y1 = lost
y2 = shared
y3 = gained
 
# plot bars in stack manner
plt.bar(x, y1, color='orange')
plt.bar(x, y2, bottom=y1, color='b')
plt.bar(x, y3, bottom=y1+y2, color='g')


plt.xlabel("")
plt.ylabel("Percentage")
plt.legend(["lost", "common", "gain"])
plt.title("PSMs 1% FDR")
plt.show()
# -




