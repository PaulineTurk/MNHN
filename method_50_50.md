Method:
Preprocessing ~ 1h 50min
Step 1: Data-Preprocessing on Pfam ~ 1h

    General folder processing + description ~ 3min
    PID computing ~ 45min
    Clustering > 99% + description ~ 5min
    Data split (50% train/50% test) ~ 2s

### Step 2: Example-Preprocessing on half Pfam ~ 50min

    must be done for Pfam_TRAIN and Pfam_TEST

    up to 6 successive contextual amino-acids in one direction, and with pid between 40 and 50%

    Selection of valid example positions ~ 3min
    Store the valid examples ~ 45min
    Store the fraction of examples per seed ~ 2min

Training Phase

    up to 6 successive contextual amino-acids in one direction, and with pid between 40 and 50%

Step 3.1: Table 2D computing ~ 1h

(7 pseudo-counters evaluated)

    2D counts: ~ 1h
    2D version frequencies, probabilities with pseudo-counters and scores ~ 0.01 s
    Visualisation 2D tables ~ 10s

Step 3.2: Table 3D count ~ 1h40min

    3D counts: ~ 1h40min

Step 3.3: Table 3D versions ~ 0.4s

    3D frequencies ~ 0.04 s
    3D probabilities with pseudo-counters ~ 0.34 s

Step 4: Example TRAIN/TEST selection ~ 6min

    1 million examples selected

    Store random examples TRAIN/TEST (repetition allowed) ~ 6min

Step 5: Pseudo-counter selection

    Brier with naive Bayes with each pseudo-counter tested

    Selection of pseudo-counter 2D + Visualisation: +
    Selection of pseudo-counter 3D with one contextual amino acid (uni) ~ 45min
    Selection of pseudo-counter 3D with multiple contextual amino acids (multi) ~ 55min
    Visualisation of pseudo-counter 3D uni/multi: /

Testing Phase
Step 6: Brier with naive Bayes

    pseudo-counters 2D and 3D fixed according to the results of step 5

    Brier Score test computing with one contextual amino acid (uni) + Visualisation ~ 5min (one test) +
    Brier Score computing with multiple contextual amino acids (multi) + Visualisation ~ 7min (one test) +
