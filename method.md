# Method:

## Preprocessing ~ 1h 50min

### Step 1: Data-Preprocessing on Pfam ~ 1h

1. General folder processing + Visualisation ~ 3min OK
2. PID computing ~ 45min OK
3. Clustering > 99% + Visualisation ~ 5min OK
4. Data split (90% train/10% test) ~ 2s OK

### Step 2: Example-Preprocessing ~ 50min

> must be done for Pfam_TRAIN and Pfam_TEST

> up to 6 successive contextual amino-acids in each direction,
and with pid between 40 and 50%

1. Selection of valid example positions: Pfam_train (90%) ~ 5min OK & Pfam_test (10%) ~ 6 s OK
2. Store the valid examples: Pfam_train (90%) ~ 1h 5min OK Pfam_test (10%) ~ 10 min OK
3. Store the fraction of examples per seed ~ Pfam_train (90%) ~ 3 min OK Pfam_test (10%) ~ 26 s OK



## Training Phase ~ 4h

> up to 6 successive contextual amino-acids in one direction,
and with pid between 40 and 50%

### Step 3.1: Table 2D computing ~ 2h

(7 pseudo-counters evaluated)

1. 2D counts: ~ 2h
2. 2D version frequencies, probabilities with pseudo-counters and scores ~ 0.03 s
3. Visualisation 2D tables ~ 12s

### Step 3.2: Table 3D count ~ 3h OK

1. 3D counts: ~ 3h OK

### Step 3.3: Table 3D versions ~ 0s ... OK

1. 3D frequencies ~ 0 s  OK
2. 3D probabilities with pseudo-counters ~ 0 s OK

### Step 4: Example TRAIN/TEST selection

> 1 million examples selected

1. Store random examples TRAIN/TEST (repetition allowed): TRAIN ~ 10min OK TEST ~ 2min OK

### Step 5: Pseudo-counter selection 

> Brier with naive Bayes with each pseudo-counter tested

1. Selection of pseudo-counter 2D + Visualisation  + ~ 1min OK
2. Selection of pseudo-counter 3D with one contextual amino acid (uni) ~ 40min
3. Selection of pseudo-counter 3D with multiple contextual amino acids (multi) ~ 55min OK
4. Visualisation of pseudo-counter 3D uni/multi: /



## Testing Phase

### Step 6: Brier with naive Bayes

> pseudo-counters 2D and 3D fixed according to the results of step 5

1. Brier Score test computing with one contextual amino acid (uni) + Visualisation ~ 5min (one test) + 
1. Brier Score computing with multiple contextual amino acids (multi) + Visualisation ~ 7min (one test) + OK
