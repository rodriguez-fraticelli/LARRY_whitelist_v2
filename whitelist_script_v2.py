
# coding: utf-8

# In[16]:


from umi_tools import network as nk
import umi_tools as umi
import numpy as np
import sys
import pandas as pd
import numpy as np, matplotlib.pyplot as plt, networkx as nx, pickle, json, gzip
import sys


# In[2]:


f1 = sys.argv[1] 
f2 = sys.argv[2]
HD = sys.argv[3]
N_READS = sys.argv[4] 


# In[57]:


library1 = np.loadtxt(f1, dtype=str)  
library2 = np.loadtxt(f2, dtype=str)  


# In[58]:


library1.shape
library2.shape


# In[59]:


def is_valid(bc):
	return bc[4:6]=='TG' and bc[10:12]=='CA' and bc[16:18]=='AC' and bc[22:24]=='GA' 


# In[60]:


out = open('filtered_'+f1, 'wb')

library1_clean = []
for xx in library1:
    if is_valid(xx):
        library1_clean.append(xx)
        out.write((xx+'\n').encode('utf-8'))
        out.write('\n'.encode('utf-8'))


# In[61]:


out = open('filtered_'+f2, 'wb')

library2_clean = []
for xx in library2:
    if is_valid(xx):
        library2_clean.append(xx)
        out.write((xx+'\n').encode('utf-8'))
        out.write('\n'.encode('utf-8'))


# In[62]:


barcodes_unique1 = np.unique(library1_clean)
barcodes_unique2 = np.unique(library2_clean)


# In[63]:


print(barcodes_unique1.shape)
print(barcodes_unique2.shape)


# In[64]:


###Here we generate a list of barcodes filtered by minimum number of reads, and print it in csv form

barcodes1 , barcodes_counts1 = np.unique(library1_clean, return_counts=True)
cbc1 = dict(zip(barcodes1, barcodes_counts1))
cbc1_filtered = {k:v for k,v in cbc1.items() if v >= N_READS}
print('Replicate 1 retaining '+repr(len(cbc1_filtered))+ ' out of '+repr(len(cbc1))+' barcodes')
library1_df = pd.DataFrame.from_dict(cbc1_filtered, orient='index')
library1_df.to_csv('library1_min_'+repr(N_READS)+'_reads.csv', '\t')

barcodes2 , barcodes_counts2 = np.unique(library2_clean, return_counts=True)
cbc2 = dict(zip(barcodes2, barcodes_counts2))
cbc2_filtered = {k:v for k,v in cbc2.items() if v >= N_READS}
print('Replicate 2 retaining '+repr(len(cbc2_filtered))+ ' out of '+repr(len(cbc2))+' barcodes')
library2_df = pd.DataFrame.from_dict(cbc2_filtered, orient='index')
library2_df.to_csv('library2_min_'+repr(N_READS)+'_reads.csv', '\t')


# In[65]:


###Here we apply UMIClusterer to collapse barcodes and generate a whitelist

uc = nk.UMIClusterer()
CBclusters = uc(cbc1_filtered,threshold=5)
cbFinal = dict()
for l in CBclusters: 
        cbFinal[l[0]] = 0
        for x in l: 
                cbFinal[l[0]] += cbc1_filtered[x]
                
whitelist1 = pd.DataFrame.from_dict(cbFinal, orient='index')
whitelist1.to_csv('whitelist1_test.csv','\t')

uc = nk.UMIClusterer()
CBclusters = uc(cbc2_filtered,threshold=5)
cbFinal = dict()
for l in CBclusters: 
        cbFinal[l[0]] = 0
        for x in l: 
                cbFinal[l[0]] += cbc2_filtered[x]
                
whitelist2 = pd.DataFrame.from_dict(cbFinal, orient='index')
whitelist2.to_csv('whitelist2_test.csv','\t')


# In[80]:


def hamming(bc1,bc2): return np.sum([x1 != x2 for x1,x2 in zip(bc1,bc2)])
all_bcs_rep1 = sorted(cbc1_filtered)
all_bcs_rep2 = sorted(cbc2_filtered)


# In[82]:


all_bcs_rep2


# In[91]:


reproducible_bcs = []
notgood_bcs = []
bc_map = {}
for i,bc1 in enumerate(all_bcs_rep1):
    if i > 0 and i % 500 == 0: print('Mapped '+repr(i)+' out of '+repr(len(all_bcs_rep1))+' barcodes')
    mapped = False
    for bc2 in all_bcs_rep2:
        if hamming(bc1,bc2) <= HD:
            mapped = True
            reproducible_bcs.append(bc1)
            bc_map[bc1] = bc2
            break
    if not mapped:
        notgood_bcs.append(bc1)

print('\nCollapsed '+repr(len(bc_map))+' barcodes')
for bc1 in reproducible_bcs: bc_map[bc1] = bc1


# In[92]:


collapsed_gfp_bcs = []
bc_map = {}
for i,bc1 in enumerate(reproducible_bcs):
    if i > 0 and i % 500 == 0: print('Mapped '+repr(i)+' out of '+repr(len(reproducible_bcs))+' barcodes')
    mapped = False
    for bc2 in reproducible_bcs:
        if hamming(bc1,bc2) <= HD:
            mapped = True
            bc_map[bc1] = bc2
            break
    if not mapped:
        collapsed_gfp_bcs.append(bc1)

print('\nCollapsed '+repr(len(bc_map))+' barcodes')
for bc1 in collapsed_gfp_bcs: bc_map[bc1] = bc1


# In[87]:


lista = whitelist1.loc[collapsed_gfp_bcs, :]


# In[89]:


lista.to_csv('whitelist_merged.csv','\t')

