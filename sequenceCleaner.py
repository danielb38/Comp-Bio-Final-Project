
hardCodedFile = 'uniprotTestDirty.txt'

def fileCleaner(file):
  fileObject = open(file, 'r')
  fileData = fileObject.read()
  listOfLines = fileData.splitlines()

  startAA = False
  sSTART:str = "SQ   SEQUENCE "
  sEND:str = "//"
  currString = ''
  listSequences = []
  for l in listOfLines:
    # print("reached3\n")
    if (startAA == True):
      if (l[0:2] == sEND):
        startAA = False
        listSequences.append(currString)
        currString = ''
      else:
        l = l.replace("     ","") #remove the original blankspace in front
        l = l.replace(" ","") #remove all spaces together
        currString += l
    else:
      lineLen = len(l)
      startingLine = ''
      if (lineLen > 14):
        startingLine = l[0:14]
      if (startingLine == sSTART):
        startAA = True
  
  return listSequences

def main(file):
  finSequence = fileCleaner(file)
  # finSequence = fileCleaner(file)
  for s in finSequence:
    print(s)

# main(hardCodedFile)

def multFiles(fileWithNames):
  fileObject = open(fileWithNames, 'r')
  fileData = fileObject.read()
  listOfLines = fileData.splitlines()
  for f in listOfLines:
    print(f)
    main(f)
    print('*****')
  
multFiles('fileNames.txt')