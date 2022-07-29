# IMPORTS


# FUNCTIONS

def table_3d_count_initialisation(alphabet):
    table_3d_count = {}
    for aa_origine in alphabet:
        table_3d_count[aa_origine] = {}
        for aa_destination in alphabet:
            table_3d_count[aa_origine][aa_destination] = {}
            for aa_context in alphabet:
                table_3d_count[aa_origine][aa_destination][aa_context]  = 0
    return table_3d_count
