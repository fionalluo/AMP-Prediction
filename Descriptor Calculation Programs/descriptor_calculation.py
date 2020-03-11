# this code caculates 20 physical properties related desciptors for QSAR based
# ANN model (equations from Fjell et al. J. Med. Chem. 2009, 52,2006-2015)

#!/usr/bin/python

import os, time, sys

# Covalent atomic radius input
def radis(element):
    if element == 'C':
        radi=0.77
        chi = 2.55
    if element == 'H':
        radi = 0.37
        chi = 2.20
    if element == 'O':
        radi = 0.73
        chi = 3.44
    if element == 'N':
        radi = 0.75
        chi = 3.04
    if element == 'S':
        radi = 1.02
        chi = 2.58
    return radi

# Electronegativity for atom input
def part(element):
    if element == 'C':
        chi = 2.55
    if element == 'H':
        chi = 2.20
    if element == 'O':
        chi = 3.44
    if element == 'N':
        chi = 3.04
    if element == 'S':
        chi = 2.58
    return chi

# Assign the position, radius, and eo information to atoms on peptides
coord = [x.strip().decode('utf-8') for x in open('peptide.pdb')]
length0 = len(coord)

ctestx =[]
ctesty =[]
ctestz =[]
cid = []
rid =[]
parc =[]

length = length0-2
for num in range(0,length):

    ctestx.append(num)
    ctesty.append(num)
    ctestz.append(num)
    cid.append(num)
    rid.append(num)
    parc.append(num)

    ctest = coord[num]
    cid[num] = ctest[13]

    ctestx0 = float(ctest[27:38])
    ctesty0 = float(ctest[39:46])
    ctestz0 = float(ctest[47:54])    
    ctestx[num] = ctestx0
    ctesty[num] = ctesty0
    ctestz[num] = ctestz0
    cele = str(cid[num])
    rid[num] = radis(cele)
    parc[num] = part(cele)

# Descriptor calculations based on hardness/softness/electrogativity of a peptide

smol = 0
smol_pos = 0.0
smol_neg = 0.0

hmol = 0
hmol_pos = 0.0
hmol_neg = 0.0

eo = 0
eo_pos = 0.0
eo_neg = 0.0

len_pos = 0
len_neg = 0

hmol_pos_smallest = 100000
hmol_pos_largest= -100000
hmol_neg_smallest = 100000
hmol_neg_largest= -100000

sigma = 0

for n1 in range(0,length):
    smol1 = 0
    hmol1 = 0
    eo1 = 0

    for n2 in range(0,length):
        if n1 != n2 :
            dr2 = (ctestx[n1]-ctestx[n2])**2.0+(ctesty[n1]-ctesty[n2])**2.0+(ctestz[n1]-ctestz[n2])**2.0
            if dr2 != 0:
                smol = smol + (rid[n1]**2.0+rid[n2]**2.0)/dr2
                smol1 = smol1 + (rid[n1]**2.0+rid[n2]**2.0)/dr2
                eo = eo + (parc[n1]-parc[n2])*(rid[n1]**2.0+rid[n2]**2.0)/dr2
                eo1 = eo1 + (parc[n1]-parc[n2])*(rid[n1]**2.0+rid[n2]**2.0)/dr2
    if(eo1 >=0):
        len_pos = len_pos + 1
        smol_pos = smol_pos+smol1
        hmol_pos = hmol_pos+1/smol1
        eo_pos = eo_pos + eo1
        if(1/smol1 >= hmol_pos_largest):
            hmol_pos_largest = 1/smol1
        if(1/smol1 < hmol_pos_smallest):
            hmol_pos_smallest = 1/smol1
    if(eo1 <0):
        len_neg = len_neg + 1
        smol_neg = smol_neg +smol1
        hmol_neg = hmol_neg +1/smol1
        eo_neg = eo_neg + eo1
        if(1/smol1 >= hmol_neg_largest):
            hmol_neg_largest = 1/smol1
        if(1/smol1 < hmol_neg_smallest):
            hmol_neg_smallest = 1/smol1

smol = smol/2
eo = eo/2

for n1 in range(0,length-2):
    s1 = 0
    s2 = 0
    for n2 in range(0,length-2):
        if n2 != n1:
            dr2 = (ctestx[n1]-ctestx[n2])**2.0+(ctesty[n1]-ctesty[n2])**2.0+(ctestz[n1]-ctestz[n2])**2.0
            if dr2 != 0:
                s1 = s1 +  (rid[n1]**2.0+rid[n2]**2.0)/dr2
                s2 = s2 + parc[n2]*(rid[n1]**2.0+rid[n2]**2.0)/dr2
    sigma = sigma + s2/s1

ster_mol_smallest = 100000
ster_mol_largest= -100000
ster_atm_smallest = 100000
ster_atm_largest= -100000

for n1 in range(0,length-2):
    s1 = 0
    s2 = 0
    for n2 in range(0,length-2):
        if n2 != n1:
            dr2 = (ctestx[n1]-ctestx[n2])**2.0+(ctesty[n1]-ctesty[n2])**2.0+(ctestz[n1]-ctestz[n2])**2.0
            if dr2 != 0:
                s1 = s1 + rid[n2]**2.0/dr2
                s2 = s2 + rid[n1]**2.0/dr2
    
    if(s1 >= ster_mol_largest):
        ster_mol_largest = s1
    if(s1 < ster_mol_smallest):
        ster_mol_smallest = s1
    if(s2 >= ster_atm_largest):
        ster_atm_largest = s2
    if(s2 < ster_atm_smallest):
        ster_atm_smallest = s2

print str(1/smol)+ ",",str(hmol_pos)+ ",",str(hmol_neg)+ ",",str(smol_pos)+ ",",str(smol_neg)+ ",",str(1/smol/length)+ ",",str(hmol_pos/len_pos)+ ",",str(hmol_neg/len_neg)+ ",",str(smol_pos/len_pos)+ ",",str(smol_neg/len_neg)+ ",",str(hmol_pos_largest)+ ",",str(hmol_neg_largest)+ ",",str(hmol_pos_smallest)+ ",",str(hmol_neg_smallest)+ ",",str(eo_pos)+ ",",str(sigma)+ ",",str(ster_mol_largest)+ ",",str(ster_mol_smallest)+ ",",str(ster_atm_largest)+ ",",str(ster_atm_smallest)
