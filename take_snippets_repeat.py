import random
import os

#selects random snippets out of the human proteome, 500 snippets of each complexity (500 three complexity,
#500 four coplexity and 500 five complexity snippets)
def take_snip():
    proteome = open("human_proteome.fasta", "r")
    lines = proteome.readlines()

    #create a proteome dictionary
    count = 0
    prot_dict = {}
    prot_seq = ""

    for line in lines:
        if line[0] == ">":
            if not(prot_seq == ""):
                prot_dict[count] = prot_seq
            count = count + 1
            prot_seq = ""
            continue
        prot_seq = prot_seq + line.replace("\n", "")


    ###use the dictionary to randomly select a protein, then a snippet from that protein
    chosen_hexa3 = []
    chosen_hexa4 = []
    chosen_hexa5 = []
    
    #for 3 complexity
    while len(chosen_hexa3) < 500:
        index = random.randint(1, count-1) #select the protein number
        prot = prot_dict[index]
        lower_bound_max = len(prot) - 7
        if lower_bound_max <= 0:
            continue
        selected_index = random.randint(0, lower_bound_max)
        if not(prot[selected_index] == "B" or prot[selected_index] == "Z" or prot[selected_index] == "X") and not(prot[selected_index+1] == "B" or prot[selected_index+1] == "Z" or prot[selected_index+1] == "X") and not(prot[selected_index+2] == "B" or prot[selected_index+2] == "Z" or prot[selected_index+2] == "X") and not(prot[selected_index+3] == "B" or prot[selected_index+3] == "Z" or prot[selected_index+3] == "X") and not(prot[selected_index+4] == "B" or prot[selected_index+4] == "Z" or prot[selected_index+4] == "X") and not(prot[selected_index+5] == "B" or prot[selected_index+5] == "Z" or prot[selected_index+5] == "X"):
            hexa = prot[selected_index] + prot[selected_index+1] + prot[selected_index+2] + prot[selected_index+3] + prot[selected_index+4] + prot[selected_index+5]
        else:
            continue
        comp = len(list(set(hexa)))
        if comp == 3:
            chosen_hexa3.append(hexa)
        else:
            continue
    
    #for 4 complexity
    while len(chosen_hexa4) < 500:
        index = random.randint(1, count-1) #select the protein number
        prot = prot_dict[index]
        lower_bound_max = len(prot) - 7
        if lower_bound_max <= 0:
            continue
        selected_index = random.randint(0, lower_bound_max)
        if not(prot[selected_index] == "B" or prot[selected_index] == "Z" or prot[selected_index] == "X") and not(prot[selected_index+1] == "B" or prot[selected_index+1] == "Z" or prot[selected_index+1] == "X") and not(prot[selected_index+2] == "B" or prot[selected_index+2] == "Z" or prot[selected_index+2] == "X") and not(prot[selected_index+3] == "B" or prot[selected_index+3] == "Z" or prot[selected_index+3] == "X") and not(prot[selected_index+4] == "B" or prot[selected_index+4] == "Z" or prot[selected_index+4] == "X") and not(prot[selected_index+5] == "B" or prot[selected_index+5] == "Z" or prot[selected_index+5] == "X"):
            hexa = prot[selected_index] + prot[selected_index+1] + prot[selected_index+2] + prot[selected_index+3] + prot[selected_index+4] + prot[selected_index+5]
        else:
            continue
        comp = len(list(set(hexa)))
        if comp == 4:
            chosen_hexa4.append(hexa)
        else:
            continue
 
    #for 5 complexity
    while len(chosen_hexa5) < 500:
        index = random.randint(1, count-1) #select the protein number
        prot = prot_dict[index]
        lower_bound_max = len(prot) - 7
        if lower_bound_max <= 0:
            continue
        selected_index = random.randint(0, lower_bound_max)
        if not(prot[selected_index] == "B" or prot[selected_index] == "Z" or prot[selected_index] == "X") and not(prot[selected_index+1] == "B" or prot[selected_index+1] == "Z" or prot[selected_index+1] == "X") and not(prot[selected_index+2] == "B" or prot[selected_index+2] == "Z" or prot[selected_index+2] == "X") and not(prot[selected_index+3] == "B" or prot[selected_index+3] == "Z" or prot[selected_index+3] == "X") and not(prot[selected_index+4] == "B" or prot[selected_index+4] == "Z" or prot[selected_index+4] == "X") and not(prot[selected_index+5] == "B" or prot[selected_index+5] == "Z" or prot[selected_index+5] == "X"):
            hexa = prot[selected_index] + prot[selected_index+1] + prot[selected_index+2] + prot[selected_index+3] + prot[selected_index+4] + prot[selected_index+5]
        else:
            continue
        comp = len(list(set(hexa)))
        if comp == 5:
            chosen_hexa5.append(hexa)
        else:
            continue   
    
    proteome.close()
    return chosen_hexa3, chosen_hexa4, chosen_hexa5
    
##main
os.chdir("/Users/pinaronat/Desktop/Thesis Work/protein_frq/snippets/repeats")
for i in range(100):
    [peptides3, peptides4, peptides5] = take_snip()
    filename = "random_snippets_{}.txt".format(i)
    write_file = open(filename, "w")
    for i in range(len(peptides3)):
        write_file.write("> random_prot" + "\n" + peptides3[i] + "\n")
    for i in range(len(peptides4)):
        write_file.write("> random_prot" + "\n" + peptides4[i] + "\n")
    for i in range(len(peptides5)):
        write_file.write("> random_prot" + "\n" + peptides5[i] + "\n")
    write_file.close()
