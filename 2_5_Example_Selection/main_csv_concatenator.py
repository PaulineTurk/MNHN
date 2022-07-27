"""
https://www.freecodecamp.org/news/how-to-combine-multiple-csv-files-with-8-lines-of-code-265183e0854/

nohup python3 main_csv_concatenator.py > concatenator.out 2>&1 &
"""

import os
import glob
import pandas as pd
import time


import sys
from pathlib import Path
file = Path(__file__).resolve()
sys.path.append(file.parents[0])


DATA = f"{file.parents[2]}/MNHN_RESULT/2_5_EXAMPLES"
L = 6
PID_INF = 40
PID_SUP = 50
NAME_CSV_CONCAT = f"ALL_EXAMPLES_{L}_{PID_INF}_{PID_SUP}"
#NAME_CSV_CONCAT = "test_concat"
os.chdir(f"{DATA}/EXAMPLES_{L}_{PID_INF}_{PID_SUP}")
#os.chdir(f"{DATA}/test_concat")

# get the list of files to concatenate
extension = 'csv'
all_filenames = [i for i in glob.glob('*.{}'.format(extension))]
print(f"NUM_FILE: {len(all_filenames)}")

# combine all files in the list
start = time.time()
combined_csv = pd.concat([pd.read_csv(f) for f in all_filenames])

# export to csv
combined_csv.to_csv(f"{DATA}/{NAME_CSV_CONCAT}.csv", index=False, encoding='utf-8-sig')
end = time.time()
print(f"CONCATENATION TIME: {round(end - start, 2)} s")