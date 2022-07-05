

import numpy as np

list_1 = [np.array([0, 1, 2]), np.array([2, 2, 2])] 
list_2 = [np.array([0, 1, 2]), np.array([4, 4, 4])] 
# print(list_1)
# print(list_1[-1])
total_list_vect = [list_1[-1], list_2[-1]]
print(total_list_vect)
#print(total_list_vect)

# final_vector = np.prod(np.vstack(total_list_vect), axis=0)
# final_vector = final_vector/ np.sum(final_vector)

# print(final_vector[2])