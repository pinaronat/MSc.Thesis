##Generates random sequences in different complexities, writes them to separate files

import random

aa_list = ["A", "M", "W", "V", "T", "G", "Y", "I", "P", "L", "F", "Q", "C", "H", "S", "R", "K", "N", "D", "E"]
seq_dict_1 = {}
seq_dict_2 = {}
seq_dict_3 = {}
seq_dict_4 = {}
seq_dict_5 = {}
seq_dict_6 = {}
#1 complexity
while len(seq_dict_1) < 20:
    aa = random.choice(aa_list)
    
    seq = aa + aa + aa + aa + aa + aa
    
    if not(seq in seq_dict_1.keys()):
        seq_dict_1[seq] = "-"

#2 complexity
while len(seq_dict_2) < 100:
    aa1 = random.choice(aa_list)
    aa2 = random.choice(aa_list)
    aas = [aa1, aa2]
    seq = ""
    while len(seq) < 6:
        seq = seq + random.choice(aas)
    
    if not(seq in seq_dict_2.keys()) and len(set(seq)) == 2:
        seq_dict_2[seq] = "-"   
        
#3 complexity
while len(seq_dict_3) < 100:
    aa1 = random.choice(aa_list)
    aa2 = random.choice(aa_list)
    aa3 = random.choice(aa_list)
    aas = [aa1, aa2, aa3]
    seq = ""
    while len(seq) < 6:
        seq = seq + random.choice(aas)

    if not(seq in seq_dict_3.keys()) and len(set(seq)) == 3:
        seq_dict_3[seq] = "-"

#4 complexity
while len(seq_dict_4) < 100:
    aa1 = random.choice(aa_list)
    aa2 = random.choice(aa_list)
    aa3 = random.choice(aa_list)
    aa4 = random.choice(aa_list)
    aas = [aa1, aa2, aa3, aa4]
    seq = ""
    while len(seq) < 6:
        seq = seq + random.choice(aas)
    
    if not(seq in seq_dict_4.keys()) and len(set(seq)) == 4:
        seq_dict_4[seq] = "-"

#5 complexity
while len(seq_dict_5) < 100:
    aa1 = random.choice(aa_list)
    aa2 = random.choice(aa_list)
    aa3 = random.choice(aa_list)
    aa4 = random.choice(aa_list)
    aa5 = random.choice(aa_list)
    aas = [aa1, aa2, aa3, aa4, aa5]
    seq = ""
    while len(seq) < 6:
        seq = seq + random.choice(aas)
    
    if not(seq in seq_dict_5.keys()) and len(set(seq)) == 5:
        seq_dict_5[seq] = "-"     

#6 complexity
while len(seq_dict_6) < 100:
    aa1 = random.choice(aa_list)
    aa2 = random.choice(aa_list)
    aa3 = random.choice(aa_list)
    aa4 = random.choice(aa_list)
    aa5 = random.choice(aa_list)
    aa6 = random.choice(aa_list)
    
    seq = aa1 + aa2 + aa3 + aa4 + aa5 + aa6
    
    if not(seq in seq_dict_6.keys()) and len(set(seq)) == 6:
        seq_dict_6[seq] = "-"


#write to a file
write_file = open("rand_sequences.txt", "w")
for key in seq_dict_1.keys():
    write_file.write(key + "\n")

for key in seq_dict_2.keys():
    write_file.write(key + "\n")

for key in seq_dict_3.keys():
    write_file.write(key + "\n")

for key in seq_dict_4.keys():
    write_file.write(key + "\n")

for key in seq_dict_5.keys():
    write_file.write(key + "\n")

for key in seq_dict_6.keys():
    write_file.write(key + "\n")
write_file.close()