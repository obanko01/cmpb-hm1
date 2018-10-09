def matrix_length(X, Y):
    m, n = len(X), len(Y)
    L = [[0] * (n+1) for k in range(m+1)] # an m+1 by n+1 matrix
    for j in range(m+1):
        for k in range(n+1):
            if j == 0:
                L[j][k] = k
            elif k == 0:
                L[j][k] = j
            elif X[j-1] == Y[k-1]:
                L[j][k] = 1 + L[j-1][k-1]
            else:
                L[j][k] = 1 + min(L[j-1][k], L[j][k-1])
    return L

# ouputting the shortest common subsequence
def sequence(X, Y, L):
    seq, j, k = [], len(X), len(Y)
    while L[j][k]>0:
        if j == 0:
            seq.append(Y[0:k])
            k=0
        elif k == 0:
            seq.append(X[0:j])
            j=0
        elif X[j-1] == Y[k-1]:
            seq.append(X[j-1])
            j=j-1
            k=k-1
        elif L[j-1][k] < L[j][k-1]:
            seq.append(X[j-1])
            j=j-1
        else:
            seq.append(Y[k-1])
            k=k-1
    return ''.join(seq[::-1])

def scs(X,Y):
    L = matrix_length(X, Y)
    value = sequence(X, Y, L)
    return value

# sample run
# X, Y = 'GTCAGTATCGCCGCCTCCACTAACCCGTAGGGACCTTGTACAGGAACAGCGCGTTGCTTACTAGTAGTAATGGGCACATTTAGTAACGCTT', \
# 'TGCAAACACTACATACCTTTCAAACCTACCCGCTTCTGTGTGTGTTATCGATATTAGGTACGAGTGGGATGAGCTGGATTAAATCTTATCGAAACTCTCG'
# scs(X, Y)
    
