def clean_text(input_file):
    f = open(input_file, "r")
    lines = f.readlines()
    f.close()
    f = open("CLEANED_"+input_file,"w")
    for l in lines:
        line = l.strip().split()
        for c in line[:-1]:
            f.write (c.strip()+" ")
        f.write (line[-1]+"\n")
    f.close()

clean_text("blosum62.txt")