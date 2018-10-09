# Question 3
# this function implementation relies on a biopython library which contains the blosum62 library.
from Bio.SubsMat import MatrixInfo
blosum = MatrixInfo.blosum62
def affine_align_score(s, t):
    m, n = len(s), len(t)
    L = [[0] * (n+1) for k in range(m+1)] # an m+1 by n+1 matrix
    best, best_loc = -(2**100), None
    for j in range(1,n+1):
        for i in range(1,m+1):
            if j == 0:
                L[i][j] = i*-11
            elif i == 0:
                L[i][j] = j*-11
            else: #match = blosum[(s[i-1], t[j-1])] #non_match = -11 - 1*(gap_length-1)
                try:
                    match = blosum[(s[i-1], t[j-1])]
                except:
                    match = blosum[(t[j-1], s[i-1])]
                
                L[i][j] = max(L[i-1][j-1] + match, 
                              max([L[k][j] + (-11 - 1*((i-k)-1)) for k in range(i)]) , 
                              max([L[i][k] + (-11 - 1*((j-k)-1)) for k in range(j)]))
    return L

def alignment(s, t, L):
    # start to trace back start from the highest score
    suffix_s, prefix_t = '', ''
    i, j = len(s), len(t)
    while i > 0 and j > 0:
        try:
            match = blosum[(s[i-1], t[j-1])]
        except:
            match = blosum[(t[j-1], s[i-1])]
            
        if L[i][j] == L[i-1][j-1] + match:
            suffix_s = s[i-1] + suffix_s 
            prefix_t = t[j-1] + prefix_t
            i -= 1
            j -= 1
        elif L[i][j] == max([L[i][k] + (-11 - 1*((j-k)-1)) for k in range(j)]):
            suffix_s = '-' + suffix_s
            prefix_t = t[j-1] + prefix_t
            j -= 1
        elif L[i][j] == max([L[k][j] + (-11 - 1*((i-k)-1)) for k in range(i)]):
            suffix_s = s[i-1] + suffix_s
            prefix_t = '-' + prefix_t
            i -= 1
        else:
            print('error')
    return suffix_s, prefix_t

def global_alignment(s, t):
    L = affine_align_score(s, t)
    return L[-1][-1], alignment(s, t, L) 

# sample test
# s, t = 'DVKVDDRQHGRINCPCNSRPKPPLVLLPKWQAKGLFRPFPDPNHRPKDWSFGCFEFIRFRRWNRHTDYAIGSNLMHSYYIHMAWI', \
# 'DVKVDDRQHGRINCAEYHTFCNSRPKPPLVLLPKWQAFLSLFRPFPWSFGCFEFIRFRRWNGSYYIHMAMI'
# global_alignment(s, t)
