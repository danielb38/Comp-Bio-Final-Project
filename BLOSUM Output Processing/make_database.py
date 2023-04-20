
def make_file_list():
    f = open ("file_list.txt")
    lines = f.readlines()
    lines = [l.strip() for l in lines]
    return lines

def format_block(input_file):
    output = open("MEME_blocks_output.txt", "a")
    input = open("MEME_Output_Original/"+input_file)
    lines = input.readlines()
    print (input_file)
    print (input_file.index("2023"))
    print (input_file.index(".fasta.result.txt"))
    identifier = input_file[input_file.index("2023"):input_file.index(".fasta.result.txt")]

    for i in range(len(lines)):
        if "BLOCKS format" in lines[i]:
            #block found!
            output.write("ID   "+identifier+"; BLOCK\nAC   AAA; distance from previous block=(1)\nDE   AAA\n")
            index_M = lines[i+2].index("MOTIF")
            output.write(lines[i+2])
            line = lines[i+3]
            x = 3
            while (line != "//\n"):
                output.write(line[10:])
                if (lines[i+x+1] != "//\n"):
                    output.write("\n")
                x += 1
                line = lines[i+x]
            output.write("//\n\n")
            i = i + x

    input.close()
    output.close()

def make_database():
    files = make_file_list()
    for f in files:
        format_block(f)

def format_block_gibbs(input_file):
    output = open("Gibbs_blocks_output.txt", "a")
    input = open("Gibbs_Output/"+input_file)
    lines = input.readlines()
    lines = [l.strip() for l in lines]
    print (input_file)
    identifier = input_file[input_file.index("2023"):input_file.index(".txt")]

    output.write("ID   "+identifier+"; BLOCK\nAC   AAA; distance from previous block=(1)\nDE   AAA\n")
    output.write("BL   MOTIF "+lines[0])
    output.write(" width="+str(len(lines[0]))+" seq="+str(len(lines))+"\n")
    for l in lines[:-1]:
        output.write("DNAAA_AAA      (    1) ")
        output.write(l+"\n\n")
    output.write("DNAAA_AAA      (    1) ")
    output.write(lines[-1]+"\n")
    output.write("//\n\n")
    input.close()
    output.close()

def make_database_gibbs():
    f = open ("gibbs_files.txt")
    files = f.readlines()
    files = [l.strip() for l in files]
    for f in files:
        format_block_gibbs(f)

make_database_gibbs()