# IMPORTS

import csv




def ex_save(seed, info_seq_dico, pid_file, L,
            new_folder_example, pid_inf, pid_sup):
    # if os.path.exists(new_folder_example):
    #     os.remove(new_folder_example)

    with open(new_folder_example, 'w',
            encoding='UTF8', newline='') as f:
        writer = csv.writer(f)
        header_context_ol = [f"aa_ol_{i}" for i in range(1, L+1)]
        header_context_or = [f"aa_or_{i}" for i in range(1, L+1)]
        header_context_dl = [f"aa_dl_{i}" for i in range(1, L+1)]
        header_context_dr = [f"aa_dr_{i}" for i in range(1, L+1)]
        header_info = ["pid", "name_origin", "name_destination", "aa_origin", "aa_destination"]

        header = header_info + header_context_ol + header_context_or + header_context_dl + header_context_dr
        writer.writerow(header)

        for name_seq_1, seq_1 in seed:
            for name_seq_2, seq_2 in seed:
                list_1 = info_seq_dico[name_seq_1]
                list_2 = info_seq_dico[name_seq_2]
                pid = pid_file[name_seq_1][name_seq_2]
                if pid_inf <= pid < pid_sup:
                    for idx, (elem_1, elem_2) in enumerate(zip(list_1, list_2)):
                        # print(f"idx: {idx}")
                        # print(f"elem_1: {elem_1}")
                        # print(f"elem_2: {elem_2}")
                        if all([elem_1 != [False], elem_2 != [False]]):
                            if all([elem_1[0] >= L, elem_2[0] >= L, elem_1[2] >= L, elem_2[2] >= L]):
                                context_ol = [seq_1[i] for i in range(idx-1, idx-1-L, -1)]
                                context_or = [seq_1[i] for i in range(idx+1, idx+1+L, +1)]
                                context_dl = [seq_2[i] for i in range(idx-1, idx-1-L, -1)]
                                context_dr = [seq_2[i] for i in range(idx+1, idx+1+L, +1)]
                                info = [pid, name_seq_1, name_seq_2, seq_1[idx], seq_2[idx]]

                                context = info + context_ol + context_or + context_dl + context_dr

                                data = [elem for elem in context]
                                #print(data)
                                writer.writerow(data)
