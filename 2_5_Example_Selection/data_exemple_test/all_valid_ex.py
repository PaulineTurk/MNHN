# IMPORTS
from itertools import groupby




# variable de code ...
# accession



alphabet = ["A", "B", "C"]

def pre_treat_seed(accession, seed, alphabet: list): # , seed_treat

    #seed = fasta.read ...(accession)

    for name_seq, seq in seed:
        
        print(f"NAME: {name_seq}")
        print(f"SEQ: {seq}")
        list_in_alphabet_right = [char in alphabet for char in seq]
        print(f"list_in_alphabet_right: {list_in_alphabet_right}")
        list_counter_right = [[k,len(list(v))] for k,v in groupby(list_in_alphabet_right)]
        print(f"list_counter_right: {list_counter_right}")
        
        info_seq = []

        list_counter_left = list_counter_right[::-1]
        print(f"list_counter_left: {list_counter_left}")


        # passage de n**2 à 2n en ordre de grandeur
        # RIGHT
        for elem in list_counter_right:
            count_temp = elem[1]
            if elem[0] == True: # True
                while count_temp > 0:
                    info_seq.append([True, count_temp -1])
                    count_temp -= 1
            else:
                for i in range(count_temp):
                    info_seq.append([False])

        # LEFT
        info_seq_mir = info_seq[::-1]
        position = 0
        for elem in list_counter_left:
            count_temp = elem[1]
            if elem[0] == True: # True
                while count_temp > 0:
                    info_seq_mir[position].insert(0, count_temp -1)
                    count_temp -= 1
                    position += 1
            else:
                for i in range(count_temp):
                    position += 1

        print(info_seq_mir[::-1])

        
        print(f"info_seq: {info_seq}")

                


def ex_valid_in_couple(seed):

    for key, value in seed:




pre_treat_seed(seed, alphabet)



