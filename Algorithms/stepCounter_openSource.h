#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <stdint.h>

#define SAMPLE_RATE_HZ        52
#define N_SAMPLES             9000   // Max number of samples
#define THRESHOLD_MULTIPLIER  0.2f
#define POST_WINDOW           8     // ~250ms at 52Hz
#define N_SCORE_WINDOW        5
#define FIR_TAPS              101

extern int process_sample_oxford(int Ax, int Ay, int Az);
extern int process_all_samples(int rows);
