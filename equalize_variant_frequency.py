import re



# Opens files of cropped and aligned data
input_file = open("new_1200_NProtein_final_combined1.fa", "r")

# Writes all sample names into file for easier traversing
output_file = open("all_samples_names.fa","w")

# Goes through each sample and writes into file
for line in input_file:
    if re.search('^>', line):
        output_file.write(line)

# Goes through all_samples_names and records the frequency of the variants
samples_dict = dict()
input_file = open("all_samples_names.fa")
for line in input_file:
    new_line = re.split(r'/|>', line)
    word = new_line[len(new_line)-1].rstrip()
    if word in samples_dict:
        samples_dict[word] = samples_dict[word] + 1
    else:
        samples_dict[word] = 1

print(samples_dict)

input_file = open("new_1200_NProtein_final_combined1.fa", "r")
output_file = open("new_1200_NProtein_final_combined_final.fa", "w")

samples_dict = dict()
highest_level = 150
for line in input_file:
    if re.search('^>', line):
        new_line = re.split('/|>', line)
        word = new_line[len(new_line)-1].rstrip()
        if word in samples_dict:
            if(samples_dict[word] >= highest_level):
                for x in range(16):
                    next(input_file)
            else:
                samples_dict[word] = samples_dict[word] + 1
                output_file.write(line)
        else:
                samples_dict[word] = 1
                output_file.write(line)
    else:
        output_file.write(line)