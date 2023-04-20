from collections import defaultdict

d1 = 12
d2 = 14

dRange = 12

numSignificantMatches = 0
numInternalPatternRepeats = 0
maxD1D2 = max(d1,d2)

(s,r,d) = (numSignificantMatches, numInternalPatternRepeats, maxD1D2)

def smithAlgorithm(sequenceList):
  motifSet = set()
  aa1Pos = defaultdict(list)
  motifPos = {} 
  #denotes the location of each motif (first index is which sequence, second is where it starts in sequence)
  #this is for the second attempts to filter out motifs within 10 of each other

  sequenceCounter = 0
  for sequence in sequenceList:
    currSequenceDict = {}
    n = len(sequence)
    for i in range(n - dRange):
      for d1 in range(i+1,i+1+dRange):
        for d2 in range(d1+1,d1+1+dRange):

          aa1 = sequence[i]
          aa2 = sequence[d1]
          aa3 = sequence[d2]
          motif = aa1+aa2+aa3
          indexList = aa1Pos[aa1]

          if i not in indexList: #filtering out duplicates
            motifSet.add(motif)
            indexList.append(i)
            aa1Pos[aa1] = indexList
            currSequenceDict[i] = motif

    motifPos[sequenceCounter] = currSequenceDict
    sequenceCounter += 1

  return motifSet

def score(motif):
  return 0

def allClosePairs(distList, dist):
  tupleList = []
  for i in range(len(distList)):
    for j in range(i+1, len(distList)):
      if (distList[j]-distList[i] <= dist):
        tupleList.append((i,j))
      else:
        j = len(distList)
  
  return tupleList



def filterScoreFunction(motifPos):
  for seqNum in motifPos:
    motifToPos = motifPos[seqNum]
    keyList = list(motifToPos.keys())
    removeTuples = allClosePairs(keyList, 10)
    removedIdx = set()
    for (i,j) in removeTuples: #properly removing tuples?
      aminoAcidI = motifPos[i]
      aminoAcidJ = motifPos[j]
      if score(aminoAcidI) < score(aminoAcidJ):
        removedIdx.add(i)
    #filter out all tuples in the above list
    #by removing the lower of the two scoring ones


  return True

def scoringFunction(L):
  return True

