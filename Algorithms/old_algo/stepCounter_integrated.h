#ifndef STEP_COUNTER_INTEGRATED_H
#define STEP_COUNTER_INTEGRATED_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <stdint.h>
#include "activity_profiles.h"

// Type alias for sensor data
typedef struct {
    float x;
    float y;
    float z;
} SENSOR_DATA_F;

#define N           48
#define HP_ORDER    4
#define HP_TAPS     5
#define BP_ORDER    8
#define BP_TAPS     9
#define BLOCK_SIZE  48

//============== PRODUCTION-GRADE CONSTANTS ==============

// Gravity detection (axis-agnostic)
#define GRAVITY_TARGET_MG           1000.0f
#define GRAVITY_TOLERANCE_MG        80.0f
#define GRAVITY_NORM_MIN            0.95f
#define GRAVITY_NORM_MAX            1.05f

// Activity state machine
#define IDLE_ENTRY_BLOCKS           2
#define IDLE_EXIT_BLOCKS            3
#define TRANSITION_THRESHOLD_G      0.25f

// Gyro validation (production-grade)
#define GYRO_IDLE_THRESHOLD_DPS     20.0f
#define MAX_ACTIVE_GYRO_SAMPLES     2
#define GYRO_WINDOW                 2
#define GYRO_RMS_IDLE_DPS           15.0f

// Threshold computation (production-grade)
#define THRESHOLD_BASE_K            1.2f
#define THRESHOLD_FINAL_K_NORMAL    0.90f
#define THRESHOLD_FINAL_K_IDLE      1.30f
#define ABSOLUTE_MIN_THRESHOLD_G    0.10f
#define ABSOLUTE_MIN_STD_G          0.02f

// Peak validation (reduced nuisance steps)
#define PEAK_PROMINENCE_MIN_NORMAL  0.08f
#define PEAK_PROMINENCE_MIN_IDLE    0.12f
#define PEAK_WIDTH_MIN_SAMPLES      3
#define PEAK_ASYMMETRY_RATIO        0.4f

// Cadence enforcement (production-grade)
#define MIN_PEAK_DISTANCE_MIN       8
#define MIN_PEAK_DISTANCE_MAX       35
#define ACTIVE_MODE_REFRACTORY_MS   250

// Idle mode strictness
#define IDLE_HARD_MIN_INTERVAL      20
#define IDLE_PROM_MIN_GAIN          0.60f
#define IDLE_PROM_WIN               5
#define IDLE_PER_WINDOW_CAP         1
#define IDLE_MIN_ENERGY_THRESHOLD   0.030f
#define IDLE_MIN_STD_THRESHOLD      0.025f

// Wrist raise detection (enhanced, production-grade)
#define RAISE_PITCH_MIN_DEG         28.0f
#define RAISE_PITCH_MAX_DEG         75.0f
#define RAISE_ROLL_MAX_DEG          60.0f
#define RAISE_Z_MIN_MG              650
#define RAISE_GYRO_MIN_DPS          12.0f
#define RAISE_GYRO_MAX_DPS          200.0f
#define RAISE_DETECT_MIN_COUNT      4
#define RAISE_ACCEL_MIN_G           0.35f
#define RAISE_SUSTAIN_BLOCKS        2
#define RAISE_COOLDOWN_BLOCKS       5

// Block-to-block continuity
#define BLOCK_CADENCE_MEMORY_DEPTH  5
#define CROSS_BLOCK_DISTANCE_CHECK  1

// Orientation stability (production-grade)
#define ORIENTATION_DELTA_THRESHOLD 3.5f
#define ORIENTATION_IDLE_MIN_SAMPLES 42
#define REF_SAMPLE_COUNT            8

// Debug & telemetry
#define ENABLE_STEP_TELEMETRY       1
#define ENABLE_WRIST_RAISE_DEBUG    1

// Legacy constants
#define GYRO_THRESH_DEG_PER_SEC     GYRO_IDLE_THRESHOLD_DPS
#define PITCH_DELTA_THRESH          25.0f
#define STABLE_STEP_PEAK_MIN        2
#define ORIENTATION_THRESHOLD       0.8f
#define ORIENTATION_HOLD_COUNT      4
#define MIN_NORM_THRESHOLD          0.05f
#define ACCEL_IDLE_THRESH_XY        60
#define ACCEL_IDLE_THRESH_Z         80
#define PITCH_IDLE_TOLERANCE_DEG    7.0f
#define ROLL_IDLE_TOLERANCE_DEG     7.0f
#define CADENCE_MIN_SPM             42
#define CADENCE_MAX_SPM             200
#define MIN_STEPS_STEPTRACK_IDLE    6
#define MIN_CONSEC_NONZERO_IDLE     2

#define RAISE_PITCH_MIN         30.0f
#define RAISE_PITCH_MAX         70.0f
#define RAISE_ROLL_ABS_MAX      50.0f
#define RAISE_Z_MG_MIN          700
#define RAISE_GYRO_THRESH_DPS   15.0f

#define GYRO_WINDOW             2
#define GYRO_IDLE_THRESHOLD     2000.0f

#ifndef GYRO_RMS_BLOCKS
#define GYRO_RMS_BLOCKS 5
#endif

//============== PRODUCTION-GRADE STRUCTURES ==============

typedef struct {
    uint8_t state;
    uint8_t entry_cnt;
    uint8_t exit_cnt;
} ActivityState;

typedef struct {
    float mean;
    float std;
    float energy;
    float max_val;
    float min_val;
    int max_idx;
    int min_idx;
} SignalStats;

typedef struct {
    int sample_idx_global;
    int last_peak_global;
    int last_step_sample;
    uint8_t consecutive_peaks;
    uint16_t cadence_history[BLOCK_CADENCE_MEMORY_DEPTH];
    uint8_t cadence_idx;
    float gyro_rms_history[GYRO_RMS_BLOCKS];
    uint8_t gyro_rms_idx;
    ActivityType_t activity_type;
    uint8_t interval_idx;
    uint8_t energy_idx;
    uint8_t activity_stable_count;
    uint16_t peak_intervals[10];
    float energy_history[5];
} StepCounterState;

typedef struct {
    uint8_t is_active;
    uint8_t sustain_count;
    uint8_t cooldown_count;
    float last_pitch;
    float last_roll;
    uint8_t motion_count;
} WristRaiseState;

typedef struct {
    uint8_t display_went_to_sleep;
} CustomerSettings_t;

extern uint8_t is_wrist_raise;
extern ActivityState g_activity_state;
extern StepCounterState g_step_state;
extern WristRaiseState g_wrist_state;
extern CustomerSettings_t g_st_Customer_Settings;

typedef struct {
	int32_t mean;
	int32_t std;
	int32_t energy;
	int32_t threshold;
	uint8_t peak_count;
	uint8_t consecutive_peaks;
	uint8_t block_steps;
	int16_t last_distance;
	uint8_t is_idle;
} StepCounterDebug;

extern StepCounterDebug g_step_debug;

//============== FUNCTION DECLARATIONS ==============

int process_step_block_integrated(SENSOR_DATA_F *accel_data, SENSOR_DATA_F *gyro_data, uint8_t total_samples);
extern void init_step_counter_state(void);
extern void update_activity_state(float accel_mag, float gyro_rms);
extern void set_selected_activity(ActivityType_t activity);
extern ActivityType_t get_selected_activity(void);
extern const ActivityProfile_t* get_current_activity_profile(void);

extern uint8_t detect_gravity_axis(float ax, float ay, float az);
extern uint8_t is_idle_state_v2(const float *ax, const float *ay, const float *az,
                                const float *gx, const float *gy, const float *gz,
                                uint8_t total_samples, uint8_t gravity_axis);
extern void compute_signal_stats(const float *data, uint8_t len, SignalStats *stats);
extern uint8_t validate_peak(const float *data, int idx, uint8_t len, float min_prom, uint8_t is_idle);
extern void detect_wrist_raise_v2(const float *ax, const float *ay, const float *az,
                                  const float *gx, const float *gy, const float *gz,
                                  uint8_t total_samples);
extern void highpass_filtfilt(const float *input, float *output);
extern void bandpass_filtfilt(const float *input, float *output);

#endif // STEP_COUNTER_INTEGRATED_H
