from collections import defaultdict
import random
from random import choice

def returnProbableKmer(string, k, dictionary):
  # print("Enter Probable Kmer")
  n = len(string)
  mostProbableKmer = string[0:k]
  mostProbableValue = 0
  for i in range(n - k - 1):
    currString = string[i:i+k]
    m = len(currString)
    currProbableValue = 1
    for j in range(m):
      currProbableValue *= dictionary[currString[j]][j] 
    if (currProbableValue > mostProbableValue):
      mostProbableValue = currProbableValue
      mostProbableKmer = currString
  # print("Exit Probable Kmer")
  return mostProbableKmer

def Score(BestMotifs, Profile):
  #generate best string using profile
  AAstring = 'ARNDCQEGHILKMFPSTWYVXBZ'
  AA = []
  for char in AAstring:
    AA.append(char)

  bestString = ''
  bestStringLength = len(BestMotifs[0])

  for i in range(bestStringLength): 
    tempDict = {}
    
    for c in AA:
      tempDict[c] = Profile[c][i]
    
    bestChar = max(tempDict, key=tempDict.get) 
    bestString += bestChar
  score = 0
  for j in range(len(BestMotifs)):
    currMotif = BestMotifs[j]
    for z in range(len(currMotif)):
      if currMotif[z] != bestString[z]:
        score += 1
  # print("Exit Score")
  return score

def ProfileCreator(Motifs):
  m = len(Motifs)
  n = len(Motifs[0])
  denominator = 2*m  ##to account for pseudocounts
  AAstring = 'ARNDCQEGHILKMFPSTWYVXBZ'
  AA = []
  for char in AAstring:
    AA.append(char)
  profile = {}
  for c in AA:
    profile[c] = [0] * n

##go through each motif, and update the counts of the number of nucleotides in each position.
  for motif in Motifs:
    for i in range(n):
      profile[motif[i]][i] +=1

  keyList = profile.keys()
  for j in keyList:
    currList = profile[j]
    profile[j] = list(map(lambda x : (x+1) / denominator, currList)) ##pseudocount

  return profile

def MotifsCreator(Profile, stringList, k):
  newMotifs = []
  for i in range(len(stringList)):
     currString = stringList[i]
     updatedString = returnProbableKmer(currString, k, Profile)
     newMotifs.append(updatedString)
  return newMotifs

def parseInput2_4(fileName):
  fileObject = open(fileName, 'r')
  fileData = fileObject.read()
  listOfLines = fileData.splitlines()
  lenListLines = [len(l) for l in listOfLines]
  k = min(36, min(lenListLines))
  t = len(listOfLines)-1
  N = 2000
  stringList = listOfLines[:-1]
  finalList = []
  for s in stringList:
    x = s.split()
    finalList += x
  return (k, t, N, finalList)

def returnDistributionKmer(string, k, dictionary):
  n = len(string)
  listOfKmers = []
  listOfProbabilities = []
  sumProb = 0
  for i in range(n - k):
    currString = string[i:i+k]
    listOfKmers.append(currString)
    m = len(currString)
    currProbableValue = 1
    for j in range(m):
      currProbableValue *= dictionary[currString[j]][j] 
    listOfProbabilities.append(currProbableValue)
    sumProb += currProbableValue

  if (n == k):
    listOfProbabilities = [1]
    sumProb = 1
    listOfKmers = [string]
  
  listOfProbabilities = list(map(lambda x: x / sumProb, listOfProbabilities))
  returnKmerList = random.choices(listOfKmers, weights = listOfProbabilities, k = 1)

  return returnKmerList[0]

def GibbsSampler(k, t, N, stringList):
  randMotifArray = [''] * t

  for i in range(t):
    n = len(stringList[i])
    if (n <= (k + 2)):
      starting = 0
    else:
      starting = random.randint(0, n-k-2)
    ending = starting + k
    randMotifArray[i] = stringList[i][starting:ending]

  Motifs = randMotifArray
  BestMotifs = randMotifArray


  for j in range(N):
    randNum = random.randint(0,t-1)
    Profile = ProfileCreator(Motifs) ##need to update Motif_i 
    newString = returnDistributionKmer(stringList[randNum], k, Profile)
    Motifs[randNum] = newString
    Motifs = MotifsCreator(Profile, stringList, k)
    if (Score(Motifs, Profile) < Score(BestMotifs, Profile)): 
      BestMotifs = Motifs
    else:
      return (BestMotifs, Score(BestMotifs, Profile))

def main2_4(fileName):
  (k, t, N, stringList) = parseInput2_4(fileName)
  # print(k)

  (bestMotifs, score) = GibbsSampler(k, t, N, stringList)

  for i in range(50):
    (bestMotifs2, score2) = GibbsSampler(k, t, N, stringList)
    if (score2 < score):
      score = score2
      bestMotifs = bestMotifs2
  
  retString = ''
  for s in bestMotifs:
    retString += (s+'\n')
  return retString

def mainFolder(fileWithNames):
  fileObject = open(fileWithNames, 'r')
  fileData = fileObject.read()
  listOfLines = fileData.splitlines()
  for l in listOfLines:
    print(l+'STR\n')
    s = main2_4(l)
    f = open('OUTPUT2'+l, 'a')
    print(l+'WRT\n')
    f.write(s)
    f.close()
    print(l+'END\n')

mainFolder('outputNames.txt')