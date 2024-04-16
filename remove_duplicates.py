# Removing the duplicate samples.
import re


input_file = open("new_1200_NProtein_final.fa", "r")
output_file = open("new_1200_NProtein_final_combined1.fa","w")

count = 0
samples_dict = dict()
for line in input_file:
    if re.search('^>', line):
        new_line = line.split('|')
        print(new_line)
        word = new_line[1].rstrip()
        if word not in samples_dict:
            count += 1
            samples_dict[word] = 1
            output_file.write(line)
        else:
            for x in range(16):
                next(input_file)
    else:
        output_file.write(line)

print(count)