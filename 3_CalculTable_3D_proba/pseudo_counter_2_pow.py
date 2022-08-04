def pseudo_counter_generator(pseudo_counter_min, pseudo_counter_max):
    LIST_PSEUDO_COUNTER_2_POW = []

    pseudo_counter = pseudo_counter_min
    while pseudo_counter <= pseudo_counter_max:
        pseudo_counter = pseudo_counter*2
        LIST_PSEUDO_COUNTER_2_POW.append(pseudo_counter)
    return LIST_PSEUDO_COUNTER_2_POW


print(pseudo_counter_generator(0.1, 10))