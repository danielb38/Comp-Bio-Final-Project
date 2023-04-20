import csv
import matplotlib.pyplot

def plot_matrices(matrix2, matrix1):
    print (matrix1)
    matrix_good = [[0 for k in matrix1.keys()] for x in matrix1.keys()]
    keys = list(matrix1.keys())
    print (keys)
    for i in range(len(matrix1.keys())):
        for j in range(i+1):
            matrix_good[i][j] = int(matrix1[keys[i]][keys[j]])
    
    matplotlib.pyplot.matshow(matrix_good,cmap=matplotlib.cm.Blues)
    matplotlib.pyplot.savefig("plot.jpg")
    

def matrix_from_csv(input_file):
    with open(input_file, newline="") as csvfile:
        lines = csvfile.readlines()
        matrix = dict()
        keys = [l.strip() for l in (lines[1].strip().split(","))[1:]]
        for k in keys:
            matrix[k] = dict()
        for i in range(0,len(keys)):
            vals = lines[2 + i].strip().split(",")
            for x in range(i+1):
                matrix[keys[i]][keys[x]] = vals[1 + i]
    
    return matrix

meme_noff = matrix_from_csv("Matrices - MEME - No Offset Matrix.csv")
meme_pos = matrix_from_csv("Matrices - MEME - Positive Matrix.csv")
meme_log = matrix_from_csv("Matrices - MEME - 10 Times Log Base 10.csv") 
gibbs_noff = matrix_from_csv("Matrices - Gibbs - No Offset Matrix.csv") #use
gibbs_pos = matrix_from_csv("Matrices - Gibbs - Positive Matrix.csv")
gibbs_log = matrix_from_csv("Matrices - Gibbs - 10 Times Log Base 10.csv")
blosum62 = matrix_from_csv("Matrices - Sheet2.csv") #use

plot_matrices(blosum62,gibbs_noff)