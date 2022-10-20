#!/usr/bin/env python
# coding: utf-8

# # FINDING SHORT TANDEM REPEATS

# In[1]:


#fasta file was read with UGENE
#content from fasta was cleaned using sequence_cleaner (http://www.cellbiol.com/scripts/cleaner/dna_protein_sequence_cleaner.php)
#read the text file containing the sequence and print it as a string
with open('mtDapma3clones_MH683636_1.txt', 'r') as sequence:
    seq = sequence.readlines()
    seq = ''.join(seq)
    #print(seq)


# In[ ]:


#sample for debugging
#seq = 'ATGATGATGATGATATAT'


# In[ ]:





# In[2]:


#create a dictionary for combinations of nucleotides that could appear as a repeat; di and tri
repeatsDict = {
    "dinucleotides": ["AT","TA", "CG", "GC", "AA", "TT",
    "CC", "GG", "AC", "AG", "TC", "TG", "CA", "CT", "GA", "GT"],
    
    "trinucleotides": ["AAA","AAT","AAC","AAG","ATG","ATC","ATA",
    "ATT","ACA","ACT","ACC","ACG","AGA","AGT","AGC","AGG","TTT",
    "TTC","TTA","TTG","TCT","TCC", "TCA","TCG","TAT","TAC","TGT", 
    "TGC","TGG","CTT","CTC","CTA","CTG","CCT", "CCC", "CCA","CCG",
    "CAT","CAC","CAA", "CAG", "CGT","CGC","CGA","CGG","GTT","GTC",
    "GTA","GTG", "GCT","GCC","GCA","GCG","GAT","GAC","GAA","GAG",
    "GGT","GGC","GGA","GGG"],
    
    "stop": ["TAG","TAA","TGA"]
    }

#print(repeatsDict["stop"])


# In[ ]:





# In[3]:


#import library with built-in function of finding mutliple substrings in a string
import re


# In[4]:


#import library with built-in function of calculating binomial distribution
from scipy import stats


# In[ ]:





# In[5]:


#BINOMIAL DISTRIBUTION
# binom(k, n, p)
#Define the number of successes (k), number of trials (n), and expected probability success (p).
# k = total STRs; total
# n = len of seq / 2 or 3; len(seq)/len(diN/tri/stC)... enter the diN/triN/stC as n
# p = freqX of each possible STR; freqX(diN/triN/stC)  


# In[6]:


#frequencey of nucleotides (Pb)
def Pb(nucl):
    return seq.count(nucl)/len(seq)

fA = Pb('A')
fT = Pb('T')
fG = Pb('G')
fC = Pb('C')
#print(fA, fT, fG, fC)
#print(fA+fT+fG+fC)

def freqX(STRs): 
    x = 1
    for i in STRs:
        #print(i)
        if i == 'A':
            x *= fA
        elif i == 'T':
            x *= fT
        elif i == 'G':
            x *= fG
        elif i == 'C':
            x *= fC
    #print(STRs, "expected freq:", x)
    return x

#freqX('TT')


# In[7]:


def nTrials(STRs):
    n = len(seq)/len(STRs)
    return int(n)

#nTrials('TGA')


# In[8]:


def binom(STRs, k, n, p):
    statBD = stats.binom.pmf(k, n, p)
    print( STRs, "binomial distribution, Pb:", statBD)
    
#binom('TGA', total, nTrials('TGA'), freqX('TGA'))


# In[ ]:





# # RESULTS

# In[9]:


#create a dictionary for all matching dinucleotides found and print position of all occurences of each dinucleotides
diNucl = {}
for i in repeatsDict["dinucleotides"]:
    match = re.finditer(pattern=i, string=seq)
    diNucl[i] = [index.start() for index in match]
#print(diNucl)


# In[10]:


#for each item in the dinucleotides, find the positions showing at least 3 repeats 
for dn, diN in zip(diNucl.values(), diNucl.keys()):
    #print(diN)  
    
    total = 0
    count = 0
    x = 0
    # extend length by 2 to allow entire element to run, this is just a random number
    dn.extend((-7, -7))
    # iterate through the length of the selected dinucleotide
    while x < len(dn)-2:
        x += 1
        
        #for each item in the list of the selected dinucleotide
        for i in dn:
            # select for next 2, so that it can be identified as at least 3 repeats
            if i == dn[x]-2 and i == dn[x+1]-4:
                count += 1
                #print(i, count)

            # else check next 1 and preceeding 1   
            elif i == dn[x]-2 and i == dn[x-2]+2:
                count += 1
                #print(i, count)

            # else check preceeding 2   
            elif i == dn[x-2]+2 and i == dn[x-3]+4:
                count += 1
                #print(i, count)
                print(count, 'STRs ending at pos.', i+2)
                
                # add count to total, then reset counter 
                total += count
                count = 0
     
    # print only nucleotides with repeats
    if total > 0:
        print(diN, 'total STRs in the sequence :', total)
        
        # print binomial probability
        binom(diN, total, nTrials(diN), freqX(diN))
        print()


# In[ ]:





# In[11]:


#create a dictionary for all matching trinucleotides found and print position of all occurences of each trinucleotides
triNucl = {}
for i in repeatsDict["trinucleotides"]:
    match = re.finditer(pattern=i, string=seq)
    triNucl[i] = [index.start() for index in match]
#print(triNucl)


# In[12]:


#for each item in the trinucleotides,  find the positions showing at least 3 repeats
for tn, triN in zip(triNucl.values(), triNucl.keys()):
    #print(triN)  
    
    total = 0
    count = 0
    x = 0
    # extend length by 2 values to allow entire element to run, this is just a random number
    tn.extend((-7, -7))
    # iterate through the length of the selected trinucleotides
    while x < len(tn)-2:
        x += 1
                
        #for each item in the list of the selected trinucleotides
        for i in tn:
            # select for next 2, so that it can be identified as at least 3 repeats
            if i == tn[x]-3 and i == tn[x+1]-6:
                count += 1
                #print(i, count)
            
            # else check next 1, and preceeding 1   
            elif i == tn[x]-3 and i == tn[x-2]+3:
                count += 1
                #print(i, count)
                
            # else check preceeding 2
            elif i == tn[x-2]+3 and i == tn[x-3]+6:
                count += 1
                #print(i, count)
                print(count, 'STRs ending at pos.', i+3)
                
                # add count to total, then reset counter 
                total += count
                count = 0
    
    # print only nucleotides with repeats
    if total > 0:
        print(triN, 'total STRs in the sequence :', total, )
        
        # print binomial probability
        binom(triN, total, nTrials(triN), freqX(triN))
        print()


# In[ ]:





# In[13]:


#create a dictionary for all matching stop codons found and print position of all occurences of each stop codons
stopCodons = {}
for i in repeatsDict["stop"]:
    match = re.finditer(pattern=i, string=seq)
    stopCodons[i] = [index.start() for index in match]
#print(stopCodons)


# In[14]:


#for each item in the stop codons,  find the positions showing at least 3 repeats
for sc, stC in zip(stopCodons.values(), stopCodons.keys()):
    #print(stC)  
    
    total = 0
    count = 0
    x = 0
    
    # extend length by 2 values to allow entire element to run, this is just a random number
    sc.extend((-7,-7))
    # iterate through the length of the selected stop codon
    while x < len(sc)-2:
        x += 1
                
        #for each item in the list of the selected stop codon
        for i in sc:
            # select for next 2, so that it can be identified as at least 3 repeats
            if i == sc[x]-3 and i == sc[x+1]-6:
                count += 1
                #print(i, count)
            
            # else check next 1, and preceeding 1   
            elif i == sc[x]-3 and i == sc[x-2]+3:
                count += 1
                #print(i, count)
                
            # else check preceeding 2
            elif i == sc[x-2]+3 and i == sc[x-3]+6:
                count += 1
                #print(i, count)
                print(count, 'STRs ending at pos.', i+3)
                
                # add count to total, then reset counter 
                total += count
                count = 0
                
    # print only nucleotides with repeats
    if total > 0:
        print(stC, 'total STRs in the sequence :', total)
        
        # print binomial probability
        binom(stC, total, nTrials(stC), freqX(stC))
        print()


# In[ ]:




