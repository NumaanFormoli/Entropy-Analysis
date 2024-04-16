# This file removes all samples that have any type of - or N's

input_file = open("new_1200_NProtein_Nucleotide.fa", "r")
output_file = open("in-between2.fa","w")

# A file is iterable
# We can read each line with a simple for loop
total_line = list()
keep = True
count = 0
for line in input_file:
    # strips line
    new_line = line.rstrip()

    # resets total_line list every sample 
    if(line[0]!='>'):
        total_line.append(new_line)
    else:
        for line in total_line:
            if(line[0]=='>'):
                continue
            letters = list(line)
            for letter in letters:
                if(letter == 'N' or letter == '-'):
                    if(count > 10):
                        keep = False
                    count += 1
        if(keep):
            count = 0
            for line in total_line:
                output_file.write(line + '\n')
        total_line = list()
        total_line.append(new_line)
        keep = True
        count = 0
