#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <stdint.h>
#include "activity_profiles.h"

#if(1)

// Type alias for sensor data arrays
typedef struct {
    float x;
    float y;
    float z;
} SENSOR_DATA_F;
extern SENSOR_DATA_F accel_data[48], gyro_data[48];;

#define N           48
#define HP_ORDER    4
#define HP_TAPS     5
#define BP_ORDER    8
#define BP_TAPS     9

#define BLOCK_SIZE 48
//============== PRODUCTION-GRADE CONSTANTS ==============

// Gravity detection (axis-agnostic)
#define GRAVITY_TARGET_MG           1000.0f
#define GRAVITY_TOLERANCE_MG        80.0f   // ±80 mg tolerance
#define GRAVITY_NORM_MIN            0.95f   // norm must be 0.95–1.05g
#define GRAVITY_NORM_MAX            1.05f

// Activity state machine
#define IDLE_ENTRY_BLOCKS           2       // 2 blocks (~1.8s) to enter idle
#define IDLE_EXIT_BLOCKS            3       // 3 blocks (~2.8s) to exit idle
#define TRANSITION_THRESHOLD_G      0.25f   // 250 mg motion to trigger exit

// Gyro validation (production-grade, fixed)
#define GYRO_IDLE_THRESHOLD_DPS     20.0f   // ±20 dps rotation
#define MAX_ACTIVE_GYRO_SAMPLES     2       // Allow only 2 high-gyro samples
#define GYRO_WINDOW                 2       // 5-point moving average
#define GYRO_RMS_IDLE_DPS           15.0f   // RMS threshold for fidget suppression

// Threshold computation (production-grade)
#define THRESHOLD_BASE_K            1.2f    // base_threshold = mean + 1.2·σ
#define THRESHOLD_FINAL_K_NORMAL    0.90f   // final threshold multiplier (active)
#define THRESHOLD_FINAL_K_IDLE      1.30f   // final threshold multiplier (idle)
#define ABSOLUTE_MIN_THRESHOLD_G    0.10f   // Absolute floor: 100 mg
#define ABSOLUTE_MIN_STD_G          0.02f   // Enforce minimum std for stability

// Peak validation (reduced nuisance steps)
#define PEAK_PROMINENCE_MIN_NORMAL  0.08f   // 80 mg minimum prominence
#define PEAK_PROMINENCE_MIN_IDLE    0.12f   // 120 mg (stricter in idle)
#define PEAK_WIDTH_MIN_SAMPLES      3       // Peak must be at least 3 samples wide
#define PEAK_ASYMMETRY_RATIO        0.4f    // Peak must have asymmetric sides

// Cadence enforcement (production-grade)
#define MIN_PEAK_DISTANCE_MIN       8       // Fastest step: ~385 ms (~2.6 Hz)
#define MIN_PEAK_DISTANCE_MAX       35      // Slowest step: ~670 ms (~1.5 Hz)
#define ACTIVE_MODE_REFRACTORY_MS   250     // 250 ms refractory in active mode

// Idle mode strictness (reduce false positives)
#define IDLE_HARD_MIN_INTERVAL      20      // ~385 ms refractory in idle
#define IDLE_PROM_MIN_GAIN          0.60f   // Stricter prominence (60% of std)
#define IDLE_PROM_WIN               5       // Wider prominence window
#define IDLE_PER_WINDOW_CAP         1       // Max 1 step per block in idle
#define IDLE_MIN_ENERGY_THRESHOLD   0.030f  // Minimum energy floor (30 mg²)
#define IDLE_MIN_STD_THRESHOLD      0.025f  // Minimum std (25 mg)

// Wrist raise detection (enhanced, production-grade)
#define RAISE_PITCH_MIN_DEG         28.0f
#define RAISE_PITCH_MAX_DEG         75.0f
#define RAISE_ROLL_MAX_DEG          60.0f
#define RAISE_Z_MIN_MG              650     // Lower threshold
#define RAISE_GYRO_MIN_DPS          12.0f   // More sensitive rotation
#define RAISE_GYRO_MAX_DPS          200.0f  // Reject unrealistic speed
#define RAISE_DETECT_MIN_COUNT      4       // 4 samples minimum
#define RAISE_ACCEL_MIN_G           0.35f   // Minimum upward acceleration
#define RAISE_SUSTAIN_BLOCKS        2       // Must sustain for 2 blocks
#define RAISE_COOLDOWN_BLOCKS       5       // 5-block cooldown

// Block-to-block continuity
#define BLOCK_CADENCE_MEMORY_DEPTH  5       // Track last 5 blocks
#define CROSS_BLOCK_DISTANCE_CHECK  1       // Check cross-block peak spacing

// Orientation stability (production-grade)
#define ORIENTATION_DELTA_THRESHOLD 3.5f    // ±3.5 degrees tolerance
#define ORIENTATION_IDLE_MIN_SAMPLES 42     // 42/48 samples for stable idle
#define REF_SAMPLE_COUNT            8       // Average first 8 samples for baseline

// Debug & telemetry
#define ENABLE_STEP_TELEMETRY       1       // Track step metrics
#define ENABLE_WRIST_RAISE_DEBUG    1       // Debug wrist raise

// Legacy constants (kept for compatibility)
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

#define GYRO_WINDOW             2  // 5-point smoothing (2 on each side)
#define GYRO_IDLE_THRESHOLD     2000.0f


#ifndef GYRO_RMS_BLOCKS
#define GYRO_RMS_BLOCKS 5
#endif

#define GYRO_Z_FIDGET_THRESHOLD     30000.0f   // 30 DPS = 30000 mdps Z-axis (pure wrist rotation = fidgeting)

//============== PRODUCTION-GRADE STRUCTURES ==============

typedef struct {
    uint8_t state;      // 0=ACTIVE, 1=TRANSITION, 2=IDLE
    uint8_t entry_cnt;  // Blocks in current state
    uint8_t exit_cnt;   // Blocks since last motion
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
    int last_step_sample;              // For 250ms hard refractory (13 samples @ 52Hz)
    uint8_t consecutive_peaks;         // Count of consecutive valid peaks
    uint16_t cadence_history[BLOCK_CADENCE_MEMORY_DEPTH];
    uint8_t cadence_idx;
    float gyro_rms_history[GYRO_RMS_BLOCKS];
    uint8_t gyro_rms_idx;
    ActivityType_t activity_type;      // Current selected activity (NORMAL_WALK, BRISK_WALK, etc)
} StepCounterState;

typedef struct {
    uint8_t is_active;
    uint8_t sustain_count;
    uint8_t cooldown_count;
    float last_pitch;
    float last_roll;
    uint8_t motion_count;
} WristRaiseState;

extern uint8_t is_wrist_raise;
extern ActivityState g_activity_state;
extern StepCounterState g_step_state;
extern WristRaiseState g_wrist_state;

// Debug structure for SD card logging (all multiplied by 10000 for int32_t storage)
typedef struct {
	int32_t mean;        // Divide by 10000 to get g
	int32_t std;         // Divide by 10000 to get g
	int32_t energy;      // Divide by 10000 to get g²
	int32_t threshold;   // Divide by 10000 to get g
	uint8_t peak_count;
	uint8_t consecutive_peaks;
	uint8_t block_steps;
	int16_t last_distance;
	uint8_t is_idle;
} StepCounterDebug;

extern StepCounterDebug g_step_debug;

int process_step_block_remaster(SENSOR_DATA_F *accel_data, SENSOR_DATA_F *gyro_data, uint8_t total_samples);
extern void init_step_counter_state(void);
extern void update_activity_state(float accel_mag, float gyro_rms);
extern void set_selected_activity(ActivityType_t activity);
extern ActivityType_t get_selected_activity(void);
extern const ActivityProfile_t* get_current_activity_profile(void);

// Helper function declarations for activity algorithms
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
#endif