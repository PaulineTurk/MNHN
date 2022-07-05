################################################################################
#                                  Importations                                #    
################################################################################
import numpy as np
from tqdm.auto import tqdm
import time
import os

import sys  
from pathlib import Path 
file = Path(__file__).resolve()
sys.path.append(file.parents[1]) 

import utils.fastaReader as fastaReader
import utils.folder as folder



################################################################################
#                                  Fonctions                                   #    
################################################################################







    

# def table_3d_count(file_fasta, path_folder_pid, triplet_count, nb_ex_train,
#                    delay_num, kp_SeqChoice, 
#                    ALPHABET, accession_num, pid_inf, pid_sup):    
#     """
#     Compte les triplets (aa_origine, aa_destination, aa_voisin) de telles sorte que
#     aa_origine, aa_destination et aa_voisin appartiennent à ALPHABET.
#     Aussi, aa_origine et aa_destination doivent appartenir à un alignement vérifiant:
#     pid_inf <= pid < pid_sup
#     """
#     pid_couple = np.load(f"{path_folder_pid}/{accession_num}.pid.npy", allow_pickle='TRUE').item()    
#     seed = fastaReader.read_multi_fasta(file_fasta)
#     len_seq = len(seed[0][1])

#     if delay_num < 0:   # before                  
#         index_range = - delay_num, len_seq
#     elif delay_num > 0: # after
#         index_range = 0, len_seq - delay_num
                        
#     if seed:
#         for name_k, seq_k in seed:
#             for name_p, seq_p in seed:
#                 if name_k != name_p:
#                     if pid_couple[name_k][name_p] >= pid_inf and pid_couple[name_k][name_p] < pid_sup:
#                         if kp_SeqChoice == "k":
#                             seq_c = seq_k
#                         elif kp_SeqChoice == "p":
#                             seq_c = seq_p
#                         for aa_index in range(index_range[0], index_range[1]):       
#                             aa_k = seq_k[aa_index] 
#                             aa_p = seq_p[aa_index] 
#                             #print("aa_index", aa_index)
#                             if all(x in ALPHABET for x in [aa_k, aa_p]):     
#                                 index_neighbor = aa_index + delay_num
#                                 aa_c = seq_c[index_neighbor] 
#                                 #print("index_neighbor", index_neighbor)
#                                 if aa_c in ALPHABET: 
#                                     #print("aa_k,aa_p,aa_c",aa_k, aa_p,aa_c) 
#                                     nb_ex_train += 1
#                                     triplet_count[aa_k][aa_p][aa_c] += 1
#     else:
#         print(accession_num)

#     return triplet_count, nb_ex_train




def table_3d_count_initialisation(ALPHABET, pseudo_compte):
    """
    Initialisation de la table 3d de taille ALPHABET **3
    un pseudo_compte strictement positif permet d'éviter que l'estimation d'une probabilité soit à 0
    """
    table_3d_count = {}
    for aa_origine in ALPHABET:
        table_3d_count[aa_origine] = {}
        for aa_destination in ALPHABET:
            table_3d_count[aa_origine][aa_destination] = {}
            for aa_voisin in ALPHABET:
                table_3d_count[aa_origine][aa_destination][aa_voisin]  = pseudo_compte
    return table_3d_count




def table_3d_count_origine_gauche(file_fasta, path_folder_pid, triplet_count, nb_ex_train, delay_num,
                                  ALPHABET, accession_num, pid_inf, pid_sup):    
    """
    Compte les triplets (aa_origine, aa_destination, aa_voisin) de telles sorte que
    aa_origine, aa_destination et aa_voisin appartiennent à ALPHABET.
    Aussi, aa_origine et aa_destination doivent appartenir à un alignement vérifiant:
    pid_inf <= pid < pid_sup
    origine_gauche: le voisin est dans la séquence d'origine et delay_num est négatif
    """
    pid_couple = np.load(f"{path_folder_pid}/{accession_num}.pid.npy", allow_pickle='TRUE').item()    
    seed = fastaReader.read_multi_fasta(file_fasta)
    len_seq = len(seed[0][1])
    index_range = np.arange(- delay_num, len_seq ,1)
    nb_seq = len(seed)


    if seed:
        for i in range(nb_seq-1):
            name_origine, seq_origine  = seed[i]
            for j in range(i+1, nb_seq):
                name_destination, seq_destination  = seed[j]

                if pid_inf <= pid_couple[name_origine][name_destination] < pid_sup:
                    for index in index_range:     
                        aa_origine = seq_origine[index] 
                        aa_destination = seq_destination[index] 

                        if all(x in ALPHABET for x in [aa_origine, aa_destination]):  
                            aa_voisin = seq_origine[index + delay_num]   
                            if aa_voisin in ALPHABET:
                                nb_ex_train += 1
                                triplet_count[aa_origine][aa_destination][aa_voisin] += 1
                            # cas symétrique
                            aa_voisin_sym = seq_destination[index + delay_num]
                            if aa_voisin_sym in ALPHABET:
                                nb_ex_train += 1
                                triplet_count[aa_destination][aa_origine][aa_voisin_sym] += 1

    else:
        print(accession_num)

    return triplet_count, nb_ex_train




def table_3d_count_origine_droite(file_fasta, path_folder_pid, triplet_count, nb_ex_train, delay_num,
                                  ALPHABET, accession_num, pid_inf, pid_sup):    
    """
    Compte les triplets (aa_origine, aa_destination, aa_voisin) de telles sorte que
    aa_origine, aa_destination et aa_voisin appartiennent à ALPHABET.
    Aussi, aa_origine et aa_destination doivent appartenir à un alignement vérifiant:
    pid_inf <= pid < pid_sup
    origine_droite: le voisin est dans la séquence d'origine et delay_num est positif
    """
    pid_couple = np.load(f"{path_folder_pid}/{accession_num}.pid.npy", allow_pickle='TRUE').item()    
    seed = fastaReader.read_multi_fasta(file_fasta)
    len_seq = len(seed[0][1])
    index_range = np.arange(0, len_seq - delay_num, 1)
    nb_seq = len(seed)


    if seed:
        for i in range(nb_seq-1):
            name_origine, seq_origine  = seed[i]
            for j in range(i+1, nb_seq):
                name_destination, seq_destination  = seed[j]

                if pid_inf <= pid_couple[name_origine][name_destination] < pid_sup:
                    for index in index_range:     
                        aa_origine = seq_origine[index] 
                        aa_destination = seq_destination[index] 
                       
                        if all(x in ALPHABET for x in [aa_origine, aa_destination]):    
                            aa_voisin = seq_origine[index + delay_num]
                            if aa_voisin in ALPHABET:
                                nb_ex_train += 1
                                triplet_count[aa_origine][aa_destination][aa_voisin] += 1
                            # cas symétrique
                            aa_voisin_sym = seq_destination[index + delay_num]
                            if aa_voisin_sym in ALPHABET:
                                nb_ex_train += 1
                                triplet_count[aa_destination][aa_origine][aa_voisin_sym] += 1

    else:
        print(accession_num)

    return triplet_count, nb_ex_train



def table_3d_count_destination_gauche(file_fasta, path_folder_pid, triplet_count, nb_ex_train, delay_num,
                                  ALPHABET, accession_num, pid_inf, pid_sup):    
    """
    Compte les triplets (aa_origine, aa_destination, aa_voisin) de telles sorte que
    aa_origine, aa_destination et aa_voisin appartiennent à ALPHABET.
    Aussi, aa_origine et aa_destination doivent appartenir à un alignement vérifiant:
    pid_inf <= pid < pid_sup
    destination_gauche: le voisin est dans la séquence d'origine et delay_num est négatif
    """
    pid_couple = np.load(f"{path_folder_pid}/{accession_num}.pid.npy", allow_pickle='TRUE').item()    
    seed = fastaReader.read_multi_fasta(file_fasta)
    len_seq = len(seed[0][1])
    index_range = np.arange(- delay_num, len_seq ,1)
    nb_seq = len(seed)


    if seed:
        for i in range(nb_seq-1):
            name_origine, seq_origine  = seed[i]
            for j in range(i+1, nb_seq):
                name_destination, seq_destination  = seed[j]

                if pid_inf <= pid_couple[name_origine][name_destination] < pid_sup:
                    for index in index_range:     
                        aa_origine = seq_origine[index] 
                        aa_destination = seq_destination[index] 
                      
                        if aa_origine in ALPHABET and aa_destination in ALPHABET:     
                            aa_voisin = seq_destination[index + delay_num]
                            if aa_voisin in ALPHABET:
                                nb_ex_train += 1
                                triplet_count[aa_origine][aa_destination][aa_voisin] += 1
                            # cas symétrique
                            aa_voisin_sym = seq_origine[index + delay_num]
                            if aa_voisin_sym in ALPHABET:
                                nb_ex_train += 1
                                triplet_count[aa_destination][aa_origine][aa_voisin_sym] += 1

    else:
        print(accession_num)

    return triplet_count, nb_ex_train




def table_3d_count_destination_droite(file_fasta, path_folder_pid, triplet_count, nb_ex_train, delay_num,
                                  ALPHABET, accession_num, pid_inf, pid_sup):    
    """
    Compte les triplets (aa_origine, aa_destination, aa_voisin) de telles sorte que
    aa_origine, aa_destination et aa_voisin appartiennent à ALPHABET.
    Aussi, aa_origine et aa_destination doivent appartenir à un alignement vérifiant:
    pid_inf <= pid < pid_sup
    destination_droite: le voisin est dans la séquence d'origine et delay_num est positif
    """
    pid_couple = np.load(f"{path_folder_pid}/{accession_num}.pid.npy", allow_pickle='TRUE').item()    
    seed = fastaReader.read_multi_fasta(file_fasta)
    len_seq = len(seed[0][1])
    index_range = np.arange(0, len_seq - delay_num, 1)
    nb_seq = len(seed)

    if seed:
        for i in range(nb_seq-1):          
            name_origine, seq_origine  = seed[i]
            # compteur_i += 1
            for j in range(i+1, nb_seq): 
                name_destination, seq_destination  = seed[j]
                # compteur_j += 1

                if pid_inf <= pid_couple[name_origine][name_destination] < pid_sup:
                    for index in index_range:     
                        aa_origine = seq_origine[index] 
                        aa_destination = seq_destination[index] 
                        #nb_ex_train += 1
                     
                        if aa_origine in ALPHABET and aa_destination in ALPHABET:  
                            aa_voisin = seq_destination[index + delay_num]
                            if aa_voisin in ALPHABET:
                                nb_ex_train += 1
                                triplet_count[aa_origine][aa_destination][aa_voisin] += 1
                            # cas symétrique
                            aa_voisin_sym = seq_origine[index + delay_num]
                            if aa_voisin_sym in ALPHABET:
                                nb_ex_train += 1
                                triplet_count[aa_destination][aa_origine][aa_voisin_sym] += 1

    else:
        print(accession_num)

    return triplet_count, nb_ex_train




def table_3d_count_destination_droite_version_plus_rapide(file_fasta, path_folder_pid, triplet_count, nb_ex_train, delay_num,
                                  ALPHABET, accession_num, pid_inf, pid_sup):    
    """
    Compte les triplets (aa_origine, aa_destination, aa_voisin) de telles sorte que
    aa_origine, aa_destination et aa_voisin appartiennent à ALPHABET.
    Aussi, aa_origine et aa_destination doivent appartenir à un alignement vérifiant:
    pid_inf <= pid < pid_sup
    destination_droite: le voisin est dans la séquence d'origine et delay_num est positif
    """
    pid_couple = np.load(f"{path_folder_pid}/{accession_num}.pid.npy", allow_pickle='TRUE').item()    
    seed = fastaReader.read_multi_fasta(file_fasta)
    len_seq = len(seed[0][1])
    index_range = np.arange(0, len_seq - delay_num, 1)
    nb_seq = len(seed)

    if seed:
        for i in range(nb_seq-1):          
            name_origine, seq_origine  = seed[i]
            for j in range(i+1, nb_seq): 
                name_destination, seq_destination  = seed[j]

                if pid_inf <= pid_couple[name_origine][name_destination] < pid_sup:
                    for index in index_range:     
                        aa_origine = seq_origine[index] 
                        aa_destination = seq_destination[index] 
                        aa_voisin = seq_destination[index + delay_num]
                     
                        if all(x in ALPHABET for x in [aa_origine, aa_destination,aa_voisin]):  
                                nb_ex_train += 1
                                triplet_count[aa_origine][aa_destination][aa_voisin] += 1
                        # cas symétrique
                        aa_voisin_sym = seq_origine[index + delay_num]
                        if all(x in ALPHABET for x in [aa_origine, aa_destination,aa_voisin_sym]): 
                            nb_ex_train += 1
                            triplet_count[aa_destination][aa_origine][aa_voisin_sym] += 1

    else:
        print(accession_num)

    return triplet_count, nb_ex_train















def multi_table_3d_count(path_folder_fasta, path_folder_pid, 
                         delay_num, RefSeqChoice, 
                         ALPHABET, pid_inf, pid_sup,
                         path_NeighborRes,
                         pseudo_compte):
    """
    Compte le nombre de chauque triplet rencontré dans les alignements multiples (MSA) de path_folder_fasta.
    Seuls les paires d'alignements vérifiant les conditions sur pid et ALPHABET sont comptés.
    Une unique position de voisinage est selectionnée par table de comptage. 
    Cette position est déterminée par les paramètres delay_num et kp_SeqChoice.

    La table3d de comptage qui en résulte est enregistrée dans path_NeighborRes 
    et est en sortie de cette fonction (ainsi que le nom du dossier des MSA)

    RefSeqChoice: origine ou destination
    """

    nb_ex_train = 0

    triplet_count = table_3d_count_initialisation(ALPHABET, pseudo_compte)
    path, dirs, files = next(os.walk(path_folder_fasta))
    nb_files = len(files)


    dico_fonction = {"origine_gauche": table_3d_count_origine_gauche,
                     "origine_droite": table_3d_count_origine_droite,
                     "destination_gauche": table_3d_count_destination_gauche,
                     "destination_droite": table_3d_count_destination_droite}

    if delay_num < 0:
        if RefSeqChoice == "origine":
            fonction_table_3d_count = "origine_gauche"
        if RefSeqChoice == "destination":
            fonction_table_3d_count = "destination_gauche"
    
    if delay_num > 0:
        if RefSeqChoice == "origine":
            fonction_table_3d_count = "origine_droite"
        if RefSeqChoice == "destination":
            fonction_table_3d_count = "destination_droite"



    # liste des PosixPath des alignements d'apprentissage
    files = [x for x in Path(path_folder_fasta).iterdir()]

    start = time.time()

    for file_counter in tqdm(range(nb_files), desc=f'compte des triplets', mininterval=60):  # mininterval: temps refresh (0.1s par défaut)
        file = files[file_counter]
        accession_num = folder.get_accession_number(file)

        triplet_count, nb_ex_train = dico_fonction[fonction_table_3d_count](file, path_folder_pid,
                                                                            triplet_count, nb_ex_train, delay_num,
                                                                            ALPHABET, accession_num, pid_inf, pid_sup)
        #print(file_counter)
        #print("nb ex train:", nb_ex_train)
    end = time.time()
    diff = end - start
    items_per_second = nb_files/diff
    print(f"Compte des triplets : {diff:.2f}s | {items_per_second:.2f}it/s")


    # pourcentage de triplets non évalués
    count_triplet_not_evaluated = 0
    for aa_1 in ALPHABET:
        for aa_2 in ALPHABET:
            for aa_3 in ALPHABET:
                if triplet_count[aa_1][aa_2][aa_3] == pseudo_compte:
                    count_triplet_not_evaluated += 1
    percentage_triplet_not_evaluated = 100*count_triplet_not_evaluated/(len(ALPHABET)**3)
    print(f"Pourcentage de triplets d'acides aminés non estimés: {percentage_triplet_not_evaluated} %")

    print("Nombre d'exemples d'apprentissage :",  '{:_}'.format(nb_ex_train))

    path_triplet_count = f"{path_NeighborRes}/table_3d_count_({str(delay_num)},{RefSeqChoice})"
    np.save(path_triplet_count , triplet_count)

    return triplet_count