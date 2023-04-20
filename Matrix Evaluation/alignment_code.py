import sys

def transpose_matrix(input_file):
    matrix_f = open("ScoringMatrix.txt", "r")
    matrix_lines = matrix_f.readlines()
    
def tester():
    files_f = open("95_file_list.txt")
    files = [f.strip() for f in files_f.readlines()]
    output = open("BLOSUM100_95_results.txt","a")
    total = 0
    counter = 0
    for f in files:
        temp = open("Inputs_Astral_95/"+f)
        if len(temp.readlines()) < 2:
            continue
        print (f)
        scores = LAMultipleCaller ("Inputs_Astral_95/"+f, "BLOSUM100ScoringMatrix.txt")
        print (scores)
        for s in scores:
            output.write(str(s)+"\n")
            total += s
        counter += len(scores)
    print ("Average: " + str(total/counter))

def LAMultipleCaller (input_file, matrix_file):
    f = open (input_file)
    lines = f.readlines()

    scores = []

    for i in range(len(lines)):
        for j in range(i + 1, len(lines)):
            v = lines[i].strip().upper()
            w = lines[j].strip().upper()
            scores.append(LACaller(v,w,matrix_file))
    
    return scores

def LACaller (v, w, matrix_file):
    matrix_f = open(matrix_file)
    matrix_lines = matrix_f.readlines()

    scoring_matrix = [[0 for i in range(23)] for j in range (23)]

    for i in range (23):
        scoring_matrix[i] = [int(x) for x in (matrix_lines[1 + i].strip().split())[1:]]

    for i in range (23):
        for j in range(i):
            scoring_matrix[j].append (scoring_matrix[i][j])


    scoring_dict = dict()

    letters = matrix_lines[0].strip().split()

    for i in range (23):
        scoring_dict[letters[i]] = i

    backtrack, scores = LABackTrack (v, w, scoring_matrix, scoring_dict)
    max_i, max_j = findMaxIndices(scores)

    s1, s2, score = OutputLA (backtrack, scores, v, w, max_i, max_j, scoring_matrix, scoring_dict)

    return score

    print (score)
    print (s1)
    print (s2, end = "")

def LABackTrack (v,w, scoring_matrix, scoring_dict):
    indel_pen = 5
    s = [[0 for i in range(len(w)+1)] for j in range(len(v)+1)]
    backtrack = [["" for i in range(len(w)+1)] for j in range(len(v)+1)]
    
    for i in range (0, len(v)+1):
        backtrack[i][0] = "stop"
    for j in range (0, len(w)+1):
        backtrack[0][j] = "stop"

    for i in range (1, len(v)+1):
        for j in range (1, len(w)+1):
            s[i][j] = max(s[i-1][j] - indel_pen, max(s[i][j-1] - indel_pen, max (s[i-1][j-1] + scoring_matrix[scoring_dict[v[i-1]]][scoring_dict[w[j-1]]],0)))
            if s[i][j] == s[i-1][j] - indel_pen:
                backtrack[i][j] = "down"
            elif s[i][j] == s[i][j-1] - indel_pen:
                backtrack[i][j] = "right"
            elif s[i][j] == s[i-1][j-1] + scoring_matrix[scoring_dict[v[i-1]]][scoring_dict[w[j-1]]]:
                backtrack[i][j] = "diag"
            elif s[i][j] == 0:
                backtrack[i][j] = "stop"
    return backtrack,s

def findMaxIndices (scores):
    max = 0
    max_indices = 0,0
    for i in range (len(scores)):
        for j in range(len(scores[0])):
            if scores[i][j] > max:
                max_indices = i,j
                max = scores[i][j]
    
    return max_indices

def OutputLA(backtrack, scores, v, w, i, j, scoring_matrix, scoring_dict):
    indel_pen = 5
    if i == 0:
        return "", "", 0
    elif j == 0:
        return "", "", 0
    if backtrack[i][j] == "down":
        temp1, temp2, score = OutputLA(backtrack, scores, v, w, i - 1, j, scoring_matrix, scoring_dict)
        temp1 += v[i-1]
        temp2 += "-"
        score -= indel_pen
        return temp1, temp2, score
    elif backtrack[i][j] == "right":
        temp1, temp2, score = OutputLA(backtrack, scores, v, w, i, j-1, scoring_matrix, scoring_dict)
        temp1 += "-"
        temp2 += w[j-1]
        score -= indel_pen
        return temp1, temp2, score
    elif backtrack[i][j] == "diag":
        temp1, temp2, score = OutputLA(backtrack,scores, v, w, i-1, j-1, scoring_matrix, scoring_dict) 
        temp1 += v[i-1]
        temp2 += w[j-1]
        score += scoring_matrix[scoring_dict[v[i-1]]][scoring_dict[w[j-1]]]
        return temp1, temp2, score
    else:
        return "", "", 0


def GABackTrack (v,w, match_rew, mis_pen, open_pen, ext_pen):
    lower = [[0 for i in range(len(w)+1)] for j in range(len(v)+1)]
    middle = [[0 for i in range(len(w)+1)] for j in range(len(v)+1)]
    upper = [[0 for i in range(len(w)+1)] for j in range(len(v)+1)]
    backtrack = [["" for i in range(len(w)+1)] for j in range(len(v)+1)]

    for i in range (1, len(v) + 1):
        middle[i][0] = - 1000000 #- open_pen - ext_pen * (i - 1)
        lower[i][0] = - open_pen - ext_pen * (i - 1)
        upper[i][0] = - 1000000
        backtrack[i][0] = "down"
    for i in range (1, len(w) + 1):
        middle[0][i] = - 1000000 #- open_pen - ext_pen * (i - 1)
        lower[0][i] = - 1000000
        upper[0][i] = - open_pen - ext_pen * (i - 1)
        backtrack[0][i] = "right"

    for i in range (1, len(v)+1):
        for j in range (1, len(w)+1):

            match = -mis_pen
            if v[i-1] == w[j-1]:
                match = match_rew
            
            lower[i][j] = max (lower[i-1][j] - ext_pen, middle[i-1][j] - open_pen)

            upper[i][j] = max (upper[i][j-1] - ext_pen, middle[i][j-1] - open_pen)

            middle[i][j] = max (max(lower[i][j], middle[i-1][j-1] + match),upper[i][j])


          #  print (middle[3][2])

            if (middle[i][j] == upper[i][j]):
                backtrack[i][j] = "right"
            elif (middle[i][j] == middle[i-1][j-1] + match):
                backtrack[i][j] = "diag"
            elif (middle[i][j] == lower[i][j]):
                backtrack[i][j] = "down"
            

    return backtrack, lower, middle, upper

def OutputGA(backtrack, v, w, i, j, match_rew, mis_pen, open_pen, ext_pen):
    if (i == 0 and j == 0):
        return "", ""
    if (i == 0): #implement base case
        return "-" * j, w[0:j]
    if (j == 0): #implement other base case
        return v[0:i], "-" * i
   
    if backtrack[i][j] == "diag":
        temp1, temp2 = OutputGA(backtrack, v, w, i-1, j-1, match_rew, mis_pen, open_pen, ext_pen) 
        temp1 += v[i-1]
        temp2 += w[j-1]
        return temp1, temp2
    elif backtrack[i][j] == "down":
        s1, s2 = OutputGA(backtrack, v, w, i-1, j, match_rew, mis_pen, open_pen, ext_pen) 
        s1 += v[i-1]
        s2 += "-"
        return s1, s2
    elif backtrack[i][j] == "right":
        temp1, temp2 = OutputGA(backtrack, v, w, i, j-1, match_rew, mis_pen, open_pen, ext_pen) 
        temp1 += "-"
        temp2 += w[j-1]
        return temp1, temp2
    else:
        print ("oops")
    

def GACaller (input_file):
    f = open (input_file)
    lines = f.readlines()
    numbers = lines[0].strip().split()
    match_rew = int(numbers[0])
    mis_pen = int(numbers[1])
    open_pen = int(numbers[2])
    ext_pen = int(numbers[3])

    #print ("match:",match_rew," miss: ", mis_pen, " open: ", open_pen, " ext: ",ext_pen)

    s = lines[1].strip()
    t = lines[2].strip()

    backtrack, lower, middle, upper = GABackTrack(s,t, match_rew, mis_pen, open_pen, ext_pen)

    score = 0

    result1, result2 = OutputGA(backtrack, s, t, len(s), len(t), match_rew, mis_pen, open_pen, ext_pen)

    open_gap1 = False
    open_gap2 = False
    for x in range (0, len(result1)):
        if result1[x] == result2[x]:
            open_gap1 = False
            open_gap2 = False
            score += match_rew
        elif result1[x] == "-":
            if open_gap1:
                score -= ext_pen
            else:
                open_gap1 = True
                open_gap2 = False
                score -= open_pen
        elif result2[x] == "-":
            if open_gap2:
                score -= ext_pen
            else:
                open_gap2 = True
                open_gap1 = False
                score -= open_pen
        else:
            score -= mis_pen
            open_gap1 = False
            open_gap2 = False

    score = middle[len(s)][len(t)]

    print (score)
    print (result1)
    print (result2, end="")


def MABackTrack (v,w,x):
    s = [[[0 for i in range(len(x)+1)] for j in range(len(w)+1)] for l in range(len(v)+1)]
    backtrack = [[["" for i in range(len(x)+1)] for j in range(len(w)+1)] for l in range(len(v) + 1)]

    for i in range (1, len(v) + 1):
        backtrack[i][0][0] = "a"
    for i in range (1, len(w) + 1):
        backtrack[0][i][0] = "b"
    for i in range (1, len(x) + 1):
        backtrack[0][0][i] = "c"

    for i in range (0, len(v)+1):
        for j in range (0, len(w)+1):
            for k in range (0, len(x) + 1):

                if (i == 0 and j == 0 and k == 0):
                    s[i][j][k] = 0

                match = 0
                if i > 0 and j > 0 and k > 0 and v[i-1] == w[j-1] and v[i-1] == x[k-1]:
                    match = 1
                
                if i > 0 and j > 0 and k > 0:
                    s[i][j][k] = max(s[i-1][j][k], s[i][j-1][k], s[i][j][k-1], s[i-1][j-1][k], s[i-1][j][k-1], s[i][j-1][k-1], s[i-1][j-1][k-1] + match)
                    if (s[i][j][k] == s[i-1][j][k]):
                        backtrack[i][j][k] = "a"
                    elif (s[i][j][k] == s[i][j-1][k]):
                        backtrack[i][j][k] = "b"
                    elif (s[i][j][k] == s[i][j][k-1]):
                        backtrack[i][j][k] = "c"
                    elif (s[i][j][k] == s[i-1][j-1][k]):
                        backtrack[i][j][k] = "d"
                    elif (s[i][j][k] == s[i-1][j][k-1]):
                        backtrack[i][j][k] = "e"
                    elif (s[i][j][k] == s[i][j-1][k-1]):
                        backtrack[i][j][k] = "f"
                    elif (s[i][j][k] == s[i-1][j-1][k-1] + match):
                        backtrack[i][j][k] = "g"
                elif i > 0 and j > 0:
                    s[i][j][k] = max(s[i-1][j][k], s[i][j-1][k], s[i-1][j-1][k])
                    if (s[i][j][k] == s[i-1][j][k]):
                        backtrack[i][j][k] = "a"
                    elif (s[i][j][k] == s[i][j-1][k]):
                        backtrack[i][j][k] = "b"
                    elif (s[i][j][k] == s[i-1][j-1][k]):
                        backtrack[i][j][k] = "d"
                elif i > 0 and k > 0:
                    s[i][j][k] = max(s[i-1][j][k], s[i][j][k-1], s[i-1][j][k-1])
                    if (s[i][j][k] == s[i-1][j][k]):
                        backtrack[i][j][k] = "a"
                    elif (s[i][j][k] == s[i][j][k-1]):
                        backtrack[i][j][k] = "c"
                    elif (s[i][j][k] == s[i-1][j][k-1]):
                        backtrack[i][j][k] = "e"
                elif j > 0 and k > 0:
                    s[i][j][k] = max(s[i][j-1][k], s[i][j][k-1], s[i][j-1][k-1])
                    if (s[i][j][k] == s[i][j-1][k]):
                        backtrack[i][j][k] = "b"
                    elif (s[i][j][k] == s[i][j][k-1]):
                        backtrack[i][j][k] = "c"
                    elif (s[i][j][k] == s[i][j-1][k-1]):
                        backtrack[i][j][k] = "f"

    return s, backtrack

def OutputMA(backtrack, v, w, x, i, j, k):
    if i == 0 and j == 0 and k == 0:
        return "", "", ""
    if backtrack[i][j][k] == "a":
        temp1, temp2, temp3 = OutputMA(backtrack, v, w, x, i - 1, j, k)
        temp1 += v[i-1]
        temp2 += "-"
        temp3 += "-"
        return temp1, temp2, temp3
    elif backtrack[i][j][k] == "b":
        temp1, temp2, temp3 = OutputMA(backtrack, v, w, x, i, j - 1, k)
        temp1 += "-"
        temp2 += w[j-1]
        temp3 += "-"
        return temp1, temp2, temp3
    elif backtrack[i][j][k] == "c":
        temp1, temp2, temp3 = OutputMA(backtrack, v, w, x, i, j, k - 1)
        temp1 += "-"
        temp2 += "-"
        temp3 += x[k-1]
        return temp1, temp2, temp3
    elif backtrack[i][j][k] == "d":
        temp1, temp2, temp3 = OutputMA(backtrack, v, w, x, i - 1, j - 1, k)
        temp1 += v[i-1]
        temp2 += w[j -1]
        temp3 += "-"
        return temp1, temp2, temp3
    elif backtrack[i][j][k] == "e":
        temp1, temp2, temp3 = OutputMA(backtrack, v, w, x, i -1 , j, k -1)
        temp1 += v[i-1]
        temp2 += "-"
        temp3 += x[k-1]
        return temp1, temp2, temp3
    elif backtrack[i][j][k] == "f":
        temp1, temp2, temp3 = OutputMA(backtrack, v, w, x, i, j -1, k - 1)
        temp1 += "-"
        temp2 += w[j-1]
        temp3 += x[k-1]
        return temp1, temp2, temp3
    else:
        temp1, temp2, temp3 = OutputMA(backtrack, v, w, x, i - 1, j -1, k - 1)
        temp1 += v[i-1]
        temp2 += w[j-1]
        temp3 += x[k-1]
        return temp1, temp2, temp3

def MACaller (input_file):
    f = open (input_file)
    lines = f.readlines()

    s = lines[0].strip()
    t = lines[1].strip()
    u = lines[2].strip()

    scores,backtrack = MABackTrack(s,t,u)

    score = scores[-1][-1][-1]

    result1, result2, result3 = OutputMA(backtrack, s, t, u, len(s), len(t), len(u))

    print (score)
    print (result1)
    print (result2)
    print (result3, end ="")


sys.setrecursionlimit(10000)
tester()



