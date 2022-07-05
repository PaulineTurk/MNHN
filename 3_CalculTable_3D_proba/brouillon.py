def table_3d_remove_1(table_3d_count, alphabet, value):
    """
    Remove a value from each cell of the 3d_table_count
    """
    for aa_origine in alphabet:
        for aa_destination in alphabet:
            # for aa_context in alphabet:
                # table_3d_count[aa_origine][aa_destination][aa_context] -= value
            table_3d_count[aa_origine][aa_destination] -= value
    print(table_3d_count)
    return table_3d_count

dico_test = {"A": {"A":1, "B":3},
             "B": {"A":5, "B":3}}

table_3d_remove_1(dico_test, ["A", "B"], 5)