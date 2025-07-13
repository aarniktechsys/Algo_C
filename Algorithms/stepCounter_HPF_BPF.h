#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <stdint.h>

#define N           48
#define HP_ORDER    4
#define HP_TAPS     5
#define BP_ORDER    8
#define BP_TAPS     9

#define BLOCK_SIZE      48

#define GYRO_THRESH_DEG_PER_SEC 25.0f
#define PITCH_DELTA_THRESH      25.0f
#define MIN_STEP_INTERVAL       16     // 0.3s at 52Hz
#define MAX_STEP_INTERVAL       40     // 0.8s at 52Hz

#define ORIENTATION_THRESHOLD   0.8f  // cos(36 degrees)
#define ORIENTATION_HOLD_COUNT  4   // ~200ms at 52Hz
#define MIN_NORM_THRESHOLD      0.05f    // avoid divide-by-zero or unstable readings

#define GYRO_WINDOW             2  // 5-point smoothing (2 on each side)
#define GYRO_IDLE_THRESHOLD     2000.0f
#define MAX_ACTIVE_GYRO_SAMPLES 3

#define ACCEL_IDLE_THRESH_XY    60
#define ACCEL_IDLE_THRESH_Z     80

#define PITCH_IDLE_TOLERANCE_DEG        7.0f
#define ROLL_IDLE_TOLERANCE_DEG         7.0f
#define ORIENTATION_IDLE_MIN_SAMPLES    40  // at least 40/48 samples must be stable
#define ORIENTATION_DELTA_THRESHOLD     2.0f   // degrees between consecutive samples
#define REF_SAMPLE_COUNT                5  // Number of samples to average for baseline

#define PI 3.14159265f

extern int last_peak_index_global;

int process_step_block(int16_t *Ax_raw, int16_t *Ay_raw, int16_t *Az_raw,
                       int16_t *Gx_raw, int16_t *Gy_raw, int16_t *Gz_raw,
                       int16_t* recent_step_counter);
int process_step_block_SF(int16_t *Ax_raw, int16_t *Ay_raw, int16_t *Az_raw,
                          int16_t *Gx_raw, int16_t *Gy_raw, int16_t *Gz_raw, uint8_t total_samples);