#calculates the frequencies of each aminoacid in the human proteome using the fasta file from UniProt
freq_dict = {}
prot_file = open("human_proteome.fasta", "r")
lines = prot_file.readlines()
 
count_A = 0
count_M = 0
count_W = 0
count_V = 0
count_T = 0
count_G = 0
count_Y = 0
count_I = 0
count_P = 0
count_L = 0
count_F = 0
count_Q = 0
count_C = 0
count_H = 0
count_S = 0
count_R = 0
count_K = 0
count_N = 0
count_D = 0
count_E = 0
count_total = 0
 
for line in lines:
    if line[0] == ">":
        continue
    for i in range(len(line)):
        if line[i].replace("\n","") == "A":
            count_A = count_A + 1
        elif line[i].replace("\n","") == "M":
            count_M = count_M + 1
        elif line[i].replace("\n","") == "W":
            count_W = count_W + 1
        elif line[i].replace("\n","") == "V":
            count_V = count_V + 1
        elif line[i].replace("\n","") == "T":
            count_T = count_T + 1
        elif line[i].replace("\n","") == "G":
            count_G = count_G + 1
        elif line[i].replace("\n","") == "Y":
            count_Y = count_Y + 1
        elif line[i].replace("\n","") == "I":
            count_I = count_I + 1
        elif line[i].replace("\n","") == "P":
            count_P = count_P + 1
        elif line[i].replace("\n","") == "L":
            count_L = count_L + 1
        elif line[i].replace("\n","") == "F":
            count_F = count_F + 1
        elif line[i].replace("\n","") == "Q":
            count_Q = count_Q + 1
        elif line[i].replace("\n","") == "C":
            count_C = count_C + 1
        elif line[i].replace("\n","") == "H":
            count_H = count_H + 1
        elif line[i].replace("\n","") == "S":
            count_S = count_S + 1
        elif line[i].replace("\n","") == "R":
            count_R = count_R + 1
        elif line[i].replace("\n","") == "K":
            count_K = count_K + 1
        elif line[i].replace("\n","") == "N":
            count_N = count_N + 1
        elif line[i].replace("\n","") == "D":
            count_D = count_D + 1
        elif line[i].replace("\n","") == "E":
            count_E = count_E + 1
        count_total = count_total + 1
freq_dict["A"] = count_A/count_total
freq_dict["M"] = count_M/count_total
freq_dict["W"] = count_W/count_total
freq_dict["V"] = count_V/count_total
freq_dict["T"] = count_T/count_total
freq_dict["G"] = count_M/count_total
freq_dict["Y"] = count_Y/count_total
freq_dict["I"] = count_I/count_total
freq_dict["P"] = count_P/count_total
freq_dict["L"] = count_L/count_total
freq_dict["F"] = count_F/count_total
freq_dict["Q"] = count_Q/count_total
freq_dict["C"] = count_C/count_total
freq_dict["H"] = count_H/count_total
freq_dict["S"] = count_S/count_total
freq_dict["R"] = count_R/count_total
freq_dict["K"] = count_K/count_total
freq_dict["N"] = count_N/count_total
freq_dict["D"] = count_D/count_total
freq_dict["E"] = count_E/count_total
freq_dict["total"] = count_total

#write the dictionary to a file
write_file = open("aa_freq.txt", "w")
for key in freq_dict.keys():
    write_file.write(key + "-" + str(freq_dict[key]) + "\n")

write_file.close()
prot_file.close()