import numpy as np

def table_3d_count_initialisation(alphabet):
    table_3d_count = {}
    for aa_origine in alphabet:
        table_3d_count[aa_origine] = {}
        for aa_destination in alphabet:
            table_3d_count[aa_origine][aa_destination] = {}
            for aa_context in alphabet:
                table_3d_count[aa_origine][aa_destination][aa_context]  = 0
    return table_3d_count

def table_3d_count_origine(seed, nb_seq,
                           alphabet,
                           index_range,
                           triplet_count,# enlever le delay number
                           start,
                           end,
                           step):    # ajouter argument direction voisinage du range
                                     # et ajouter dans la fonction

    nb_ex_train = 0
    for i in range(nb_seq-1):
        name_origine, seq_origine  = seed[i]
        print(f"NAME ORIGINE: {name_origine}")
        for j in range(i+1, nb_seq):
            name_destination, seq_destination  = seed[j]
            print(f"NAME DESTINATION: {name_destination}")

            for index in index_range:
                aa_origine = seq_origine[index]
                aa_destination = seq_destination[index]
                print(f"aa_o, aa_d: {aa_origine, aa_destination}")

                if all(x in alphabet for x in [aa_origine, aa_destination]):
                    list_aa_context = [seq_origine[position] for position in range(index + start,
                                                                                   index + end,
                                                                                   step)]
                    print(f"list_aa_context: {list_aa_context}")
                    if all(x in alphabet for x in list_aa_context):
                        nb_ex_train += 1
                        aa_context = list_aa_context[-1]
                        print(f"aa_context: {aa_context}")
                        print(nb_ex_train)
                        triplet_count[aa_origine][aa_destination][aa_context] += 1

                    # cas symétrique
                    print("SYM")
                    print(f"aa_o, aa_d: {aa_destination, aa_origine}")
                    list_aa_context_sym = [seq_destination[position] for position in range(index + start,
                                                                                       index + end,
                                                                                       step)]
                    print(f"list_aa_context_sym:: {list_aa_context_sym:}")
                    if all(x in alphabet for x in list_aa_context_sym):
                        nb_ex_train += 1
                        aa_context_sym = list_aa_context_sym[-1]
                        print(f"aa_context_sym: {aa_context_sym}")
                        print(nb_ex_train)
                        triplet_count[aa_destination][aa_origine][aa_context_sym] += 1
    
    print(triplet_count)
    print(nb_ex_train)

    return triplet_count, nb_ex_train





alphabet = ["A", "B", "C"]

seed = (("seq1", "ABBCC-AB"),
        ("seq2", "AABCCAAB"),
        ("seq3", "AABC--AB"))
nb_seq = len(seed)
len_seq = len(seed[0][1])
delay_num = 3

if delay_num < 0:
    index_range = np.arange(- delay_num, len_seq ,1)
    start, end, step = -1, delay_num-1, -1
else:
    index_range = np.arange(0, len_seq - delay_num, 1)
    start, end, step = 1, delay_num +1, 1


triplet_count = table_3d_count_initialisation(alphabet)
    
table_3d_count_origine(seed, nb_seq,
                           alphabet,
                           index_range,
                           triplet_count,
                           start, end, step)