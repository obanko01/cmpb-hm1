# Question 2
match, non_match = 1, -2

def align_score(s, t):
    m, n = len(s), len(t)
    L = [[0] * (n+1) for k in range(m+1)] # an m+1 by n+1 matrix
    best, best_loc = -(2**100), None
    for j in range(n+1):
        for i in range(m+1):
            if j == 0:
                L[i][j] = 0
            elif i == 0:
                L[i][j] = j*non_match
            else:
                match = 1 if s[i-1] == t[j-1] else -2
                L[i][j] = max(L[i-1][j-1] + match, 
                              L[i-1][j] + non_match, 
                              L[i][j-1] + non_match)
            if L[i][j] >= best and i == m:
                best = L[i][j]
                best_loc = i,j
    return best, best_loc, L

def alignment(s, t, best_loc, L):
# def alignment(s, t, L):
    # start to trace back start from the highest score
    i, j = best_loc
    suffix_s, prefix_t = '', ''
    while j > 0:
        match = 1 if s[i-1] == t[j-1] else -2
        if L[i][j] == L[i][j-1] + non_match:
            suffix_s = '-' + suffix_s
            prefix_t = t[j-1] + prefix_t
            j -= 1
        elif L[i][j] == L[i-1][j-1] + match:
            suffix_s = s[i-1] + suffix_s 
            prefix_t = t[j-1] + prefix_t
            i -= 1
            j -= 1
        elif L[i][j] == L[i-1][j] + non_match:
            suffix_s = s[i-1] + suffix_s
            prefix_t = '-' + prefix_t
            i -= 1
    return suffix_s, prefix_t
        
def local_alignment(s, t):
    best, best_loc, L = align_score(s, t)
    return best, alignment(s, t, best_loc, L)
    
# sample run
# s, t = 'CTAAGGGATTCCGGTAATTAGACAG','ATAGACCATATGTCAGTGACTGTGTAA'
# local_alignment(s, t)
