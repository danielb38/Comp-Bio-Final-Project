def process_data():
    f = open ("astral-scopdom-seqres-gd-sel-gs-bib-95-1.55.fa.txt")
    lines = [l.strip() for l in f.readlines()]
    for i in range(len(lines)):
        if "{" in lines[i]:
            print(lines[i])
            family = lines[i][lines[i].index(" ") + 1:lines[i].index(" (")]
            output = open ("Inputs_Astral_95/"+family+".txt","a")
            x = 1
            while (i + x < len(lines) and not "{" in lines[i+x]):
                output.write(lines[i+x])
                x += 1
            output.write("\n")

process_data()