#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <stdint.h>

#define SAMPLE_RATE                 52
#define NUM_SAMPLES_PER_ITERATION   48
#define PEAK_MIN_DISTANCE           3 //4 //6   // Allow more frequent steps (~100-120/minute walking)
#define Q                           6                    // Fixed-point: Q6 format
#define ONE_G_FIXED                 (1 << Q)  // 1g in Q6
#define THRESH_SCALE_NUM            1 //2    // 0.2 * std_dev
#define THRESH_SCALE_DEN            10
#define CLAMP_MAX                   3 * ONE_G_FIXED
#define DIP_THRESHOLD               8

#define N 48                        // FIFO batch size
#define WINDOW_SIZE 4               // For scoring window
#define PEAK_WINDOW 48              // Post-processing window
#define THRESH_FACTOR 1             // Threshold multiplier
#define STD_MIN 1                   // Minimum stddev to consider motion

#define FILTER_ORDER                21      // Example:  Adjust as needed
#define Q15_FACTOR                  32768    // 2^15 for Q15 fixed-point representation
#define RUNNING_AVERAGE_WINDOW      30
#define THRESHOLD_MULTIPLIER        2    // Integer multiplier (for Q15)
#define HYSTERESIS_SAMPLES          5

#define FIR_LEN                     33
#define Q_SCALE                     1024
#define DIP_WINDOW                  8
#define MAX_RECENT_STEPS            5
#define MAG_BUF_LEN                 64

#define NORMAL_JOGGING              0
#define BRISK_WALKING               1

#define IDLE_WINDOW_SIZE   			32
#define IDLE_ACCEL_THRESH  			70
#define IDLE_GYRO_THRESH   			3000  // Adjust based on gyro full scale

#define MAG_MIN_THREHSOLD			20

#define GYRO_SUPPRESS_THRESHOLD 	22000

extern int step_count;

extern int16_t process_sample(int Ax, int Ay, int Az, int Gx, int Gy, int Gz);