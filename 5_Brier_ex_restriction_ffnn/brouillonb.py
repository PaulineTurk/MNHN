

import numpy as np




list_example = []
ALPHABET = ["A", "B"]
example_selected = ["A", "BB", 'A', "*", ""]

example_selected_str = ''.join(example_selected)
print(example_selected_str)
if all(character in ALPHABET for character in example_selected_str):
    list_example.append(example_selected)
print(list_example)


