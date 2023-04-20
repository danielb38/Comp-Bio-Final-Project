def main(file):
  fileObject = open(file, 'r')
  fileData = fileObject.read()
  listOfLines = fileData.splitlines()
  prevEnd = True
  fileName = ''
  for l in listOfLines:
    if (prevEnd):
      fileName = l
      prevEnd = False
    elif (l == '*****'):
      prevEnd = True
      f = open('OUTPUT'+fileName, 'a')
      f.close()
      continue
    else:
      f = open('OUTPUT'+fileName, 'a')
      f.write(l+'\n')

main('cleanedSeq.txt')