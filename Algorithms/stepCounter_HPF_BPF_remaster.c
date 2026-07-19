#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include "stepCounter_HPF_BPF_remaster.h"

#define DEBUG_SAMPLES  1

#if (DEBUG_SAMPLES)
typedef struct {
	float std_d;
	float energy_d;
	float mean_z_d;
	float threshold_d;
	float avg_interval_d;
	int min_peak_distance_d;
	uint8_t is_accel_idle_d;
	uint8_t is_gyro_idle_d;
	uint8_t is_orientation_idle_d;
	uint16_t steps_d;
	uint8_t peak_indices[48];
	uint8_t peak_cnt_d;
	uint8_t gravity_axis;
	float gravity_magnitude;
	uint8_t activity_state;
} STEP_COUNT_PARAMS_DEBUG;

STEP_COUNT_PARAMS_DEBUG g_st_StepCntParams_debug = {0};
#endif

#define FILTER_STAGES 1

// Customer settings structure
typedef struct {
	uint8_t display_went_to_sleep;
} CustomerSettings_t;

CustomerSettings_t g_st_Customer_Settings = {0};

float Ax[BLOCK_SIZE] = {0};
float Ay[BLOCK_SIZE] = {0};
float Az[BLOCK_SIZE] = {0};
float Gx[BLOCK_SIZE] = {0};
float Gy[BLOCK_SIZE] = {0};
float Gz[BLOCK_SIZE] = {0};

uint8_t is_wrist_raise = false;
ActivityState g_activity_state = {0};
StepCounterState g_step_state = {0};
WristRaiseState g_wrist_state = {0};
StepCounterDebug g_step_debug = {0};  // Debug data for SD card logging

// High-pass filter (0.3Hz @ 52Hz)
const float b_hp[HP_TAPS] = { 0.953739f, -3.814958f, 5.722436f, -3.814958f, 0.953739f };
const float a_hp[HP_TAPS] = { 1.0f, -3.905279f, 5.720297f, -3.724635f, 0.909619f };

// Band-pass filter (0.5–4.5Hz @ 52Hz)
const float b_bp[9] = {
    0.00194832f, 0.0f, -0.00779327f, 0.0f, 0.01168990f,
    0.0f, -0.00779327f, 0.0f, 0.00194832f
};

const float a_bp[9] = {
    1.0f, -6.62984475f, 19.36053485f, -32.55730178f,
    34.50841663f, -23.61716055f, 10.19337717f, -2.53669815f, 0.27867724f
};

//===============================================================================
// @brief Initialize step counter state machine
//===============================================================================
void init_step_counter_state(void) {
	memset(&g_activity_state, 0, sizeof(ActivityState));
	memset(&g_step_state, 0, sizeof(StepCounterState));
	memset(&g_wrist_state, 0, sizeof(WristRaiseState));
	g_activity_state.state = 0;  // Start in ACTIVE
	g_activity_state.entry_cnt = 0;
	g_step_state.sample_idx_global = 0;
	g_step_state.last_peak_global = -20;  // Initialize so first peak in block 0 passes cadence (distance will be 0-20)
	g_step_state.activity_type = ACTIVITY_IDLE;  // Default activity
}

//===============================================================================
// @brief Set the selected activity for precise step counting
// @param activity Activity type (NORMAL_WALK, BRISK_WALK, JOGGING, RUNNING)
//===============================================================================
void set_selected_activity(ActivityType_t activity) {
	if (activity >= ACTIVITY_IDLE && activity <= ACTIVITY_RUNNING) {
		g_step_state.activity_type = activity;
#if STEP_DIAG
		printf("[Activity] Set to: %s\n", get_activity_name(activity));
#endif
	}
}

//===============================================================================
// @brief Get the currently selected activity
// @return Current activity type
//===============================================================================
ActivityType_t get_selected_activity(void) {
	return g_step_state.activity_type;
}

//===============================================================================
// @brief Get the activity profile for current selected activity
// @return Pointer to activity-specific thresholds
//===============================================================================
const ActivityProfile_t* get_current_activity_profile(void) {
	return get_activity_profile(g_step_state.activity_type);
}

//===============================================================================
// @brief Forward-backward high-pass filter
//===============================================================================
void highpass_filtfilt(const float *input, float *output){
    float x[HP_TAPS] = {0}, y[HP_TAPS] = {0}, temp[N];

    for (int n = 0; n < N; ++n) {
        for (int i = HP_ORDER; i > 0; --i) { x[i] = x[i - 1]; y[i] = y[i - 1]; }
        x[0] = input[n];
        y[0] = b_hp[0]*x[0];
        for (int i = 1; i <= HP_ORDER; ++i)
            y[0] += b_hp[i]*x[i] - a_hp[i]*y[i];
        temp[n] = y[0];
    }
    for (int i = 0; i < HP_TAPS; ++i) x[i] = y[i] = 0;
    for (int n = N - 1; n >= 0; --n) {
        for (int i = HP_ORDER; i > 0; --i) { x[i] = x[i - 1]; y[i] = y[i - 1]; }
        x[0] = temp[n];
        y[0] = b_hp[0]*x[0];
        for (int i = 1; i <= HP_ORDER; ++i)
            y[0] += b_hp[i]*x[i] - a_hp[i]*y[i];
        output[n] = y[0];
    }
}

//===============================================================================
// @brief Forward-backward band-pass filter
//===============================================================================
void bandpass_filtfilt(const float *input, float *output) {
    float x[BP_TAPS] = {0}, y[BP_TAPS] = {0}, temp[N];

    for (int n = 0; n < N; ++n) {
        for (int i = BP_ORDER; i > 0; --i) { x[i] = x[i-1]; y[i] = y[i-1]; }
        x[0] = input[n];
        y[0] = b_bp[0]*x[0];
        for (int i = 1; i <= BP_ORDER; ++i)
            y[0] += b_bp[i]*x[i] - a_bp[i]*y[i];
        temp[n] = y[0];
    }
    for (int i = 0; i < BP_TAPS; ++i) x[i] = y[i] = 0;
    for (int n = N - 1; n >= 0; --n) {
        for (int i = BP_ORDER; i > 0; --i) { x[i] = x[i-1]; y[i] = y[i-1]; }
        x[0] = temp[n];
        y[0] = b_bp[0]*x[0];
        for (int i = 1; i <= BP_ORDER; ++i)
            y[0] += b_bp[i]*x[i] - a_bp[i]*y[i];
        output[n] = y[0];
    }
}

//===============================================================================
// @brief Detect gravity axis (production-grade, axis-agnostic)
// Returns: 0=X, 1=Y, 2=Z, 3=not detected
//===============================================================================
uint8_t detect_gravity_axis(float ax_g, float ay_g, float az_g) {
	float ax_abs = fabsf(ax_g);
	float ay_abs = fabsf(ay_g);
	float az_abs = fabsf(az_g);

	// Gravity should be strongest on one axis (near 1.0g)
	if (ax_abs > ay_abs && ax_abs > az_abs) {
		if (ax_abs >= GRAVITY_NORM_MIN && ax_abs <= GRAVITY_NORM_MAX)
			return 0;  // Gravity on X
	} else if (ay_abs > ax_abs && ay_abs > az_abs) {
		if (ay_abs >= GRAVITY_NORM_MIN && ay_abs <= GRAVITY_NORM_MAX)
			return 1;  // Gravity on Y
	} else {
		if (az_abs >= GRAVITY_NORM_MIN && az_abs <= GRAVITY_NORM_MAX)
			return 2;  // Gravity on Z
	}
	return 3;  // Not detected
}

//===============================================================================
// @brief Check if device is in idle state (production-grade)
//===============================================================================
uint8_t is_idle_state_v2(
		const float *ax, const float *ay, const float *az,
		const float *gx, const float *gy, const float *gz,
		uint8_t total_samples, uint8_t gravity_axis) {

	int idle_accel_count = 0;
	int idle_orientation_count = 0;
	float gyro_rms = 0.0f;
	float ref_ax = 0.0f, ref_ay = 0.0f, ref_az = 0.0f;

	// Reference orientation: average first 8 samples
	for (int i = 0; i < (REF_SAMPLE_COUNT < total_samples ? REF_SAMPLE_COUNT : total_samples); i++) {
		ref_ax += ax[i];
		ref_ay += ay[i];
		ref_az += az[i];
	}
	ref_ax /= REF_SAMPLE_COUNT;
	ref_ay /= REF_SAMPLE_COUNT;
	ref_az /= REF_SAMPLE_COUNT;

	// Check each sample
	for (int i = 0; i < total_samples; i++) {
		// Accel idle: low motion magnitude
		float accel_mag = sqrtf(ax[i]*ax[i] + ay[i]*ay[i] + az[i]*az[i]);
		if (fabsf(accel_mag - 1.0f) < 0.12f)  // Within ±120 mg of gravity
			idle_accel_count++;

		// Orientation idle: stable pitch/roll
		float pitch = atan2f(ax[i], sqrtf(ay[i]*ay[i] + az[i]*az[i])) * 180.0f / M_PI;
		float roll = atan2f(ay[i], sqrtf(ax[i]*ax[i] + az[i]*az[i])) * 180.0f / M_PI;

		float ref_pitch = atan2f(ref_ax, sqrtf(ref_ay*ref_ay + ref_az*ref_az)) * 180.0f / M_PI;
		float ref_roll = atan2f(ref_ay, sqrtf(ref_ax*ref_ax + ref_az*ref_az)) * 180.0f / M_PI;

		if (fabsf(pitch - ref_pitch) < ORIENTATION_DELTA_THRESHOLD &&
			fabsf(roll - ref_roll) < ORIENTATION_DELTA_THRESHOLD)
			idle_orientation_count++;

		// Gyro RMS
		gyro_rms += gx[i]*gx[i] + gy[i]*gy[i] + gz[i]*gz[i];
	}

	gyro_rms = sqrtf(gyro_rms / total_samples);

	uint8_t is_accel_idle = (idle_accel_count >= (total_samples - 6));
	uint8_t is_orientation_idle = (idle_orientation_count >= ORIENTATION_IDLE_MIN_SAMPLES);
	uint8_t is_gyro_idle = (gyro_rms < GYRO_IDLE_THRESHOLD_DPS);

	return ((is_accel_idle || is_orientation_idle) && is_gyro_idle);
}

//===============================================================================
// @brief Compute signal statistics (production-grade)
//===============================================================================
void compute_signal_stats(
		const float *signal, uint8_t len,
		SignalStats *stats) {
	memset(stats, 0, sizeof(SignalStats));

	// Mean
	for (int i = 0; i < len; i++)
		stats->mean += signal[i];
	stats->mean /= len;

	// Std, energy, min/max
	for (int i = 0; i < len; i++) {
		float diff = signal[i] - stats->mean;
		stats->std += diff * diff;
		stats->energy += signal[i] * signal[i];

		if (signal[i] > stats->max_val || i == 0) {
			stats->max_val = signal[i];
			stats->max_idx = i;
		}
		if (signal[i] < stats->min_val || i == 0) {
			stats->min_val = signal[i];
			stats->min_idx = i;
		}
	}

	stats->std = sqrtf(stats->std / len);
	stats->energy /= len;

	// Enforce minimum std for stability
	if (stats->std < ABSOLUTE_MIN_STD_G)
		stats->std = ABSOLUTE_MIN_STD_G;
}

//===============================================================================
// @brief Detect activity type based on average step rate
// Maps step interval (samples) to activity classification
//===============================================================================
ActivityType_Auto detect_activity_from_rate(uint16_t avg_step_interval_samples) {
	// Convert samples to BPM: BPM = 60 Hz / (samples / 52 Hz) = 3120 / samples
	float bpm = (avg_step_interval_samples > 0) ? (3120.0f / avg_step_interval_samples) : 0;

	if (bpm < 80) return AUTO_SLOW_WALK;       // <80 BPM
	if (bpm < 110) return AUTO_NORMAL_WALK;    // 80-110 BPM
	if (bpm < 140) return AUTO_BRISK_WALK;     // 110-140 BPM
	if (bpm < 180) return AUTO_JOGGING;        // 140-180 BPM
	return AUTO_RUNNING;                       // >180 BPM
}

//===============================================================================
// @brief Get adaptive cadence range based on detected activity
//===============================================================================
void get_adaptive_cadence_range(ActivityType_Auto activity, uint8_t *min, uint8_t *max) {
	// At 52 Hz, X samples = (X/52)*1000 ms
	switch (activity) {
		case AUTO_SLOW_WALK:    // <80 BPM, step ~750ms = 39 samples
			*min = 22;  // 423ms
			*max = 75;  // 1442ms
			break;

		case AUTO_NORMAL_WALK:  // 80-110 BPM, step ~545ms = 28 samples
			*min = 18;  // 346ms (173 BPM max)
			*max = 70;  // 1346ms (44 BPM min)
			break;

		case AUTO_BRISK_WALK:   // 110-140 BPM, step ~430ms = 22 samples
			*min = 15;  // 288ms (207 BPM max)
			*max = 65;  // 1250ms (48 BPM min)
			break;

		case AUTO_JOGGING:      // 140-180 BPM, step ~345ms = 18 samples
			*min = 12;  // 231ms (259 BPM max) - catches faster jogging
			*max = 55;  // 1058ms (57 BPM min)
			break;

		case AUTO_RUNNING:      // >180 BPM, step ~333ms = 17 samples
			*min = 10;  // 192ms (323 BPM max) - catches fast running
			*max = 50;  // 962ms (62 BPM min)
			break;

		default:
			*min = 18;  // Default to normal walk
			*max = 70;
	}
}

//===============================================================================
// @brief Update step interval tracking for activity detection
//===============================================================================
void update_step_interval(int distance_since_last_peak) {
	g_step_state.step_intervals[g_step_state.interval_idx] = (uint16_t)distance_since_last_peak;
	g_step_state.interval_idx = (g_step_state.interval_idx + 1) % 10;
	if (g_step_state.interval_count < 10) {
		g_step_state.interval_count++;
	}
}

//===============================================================================
// @brief Calculate average step interval from tracked intervals
//===============================================================================
uint16_t get_average_step_interval(void) {
	if (g_step_state.interval_count == 0) return 0;

	uint32_t sum = 0;
	for (uint8_t i = 0; i < g_step_state.interval_count; i++) {
		sum += g_step_state.step_intervals[i];
	}
	return (uint16_t)(sum / g_step_state.interval_count);
}

//===============================================================================
// @brief Compute peak "peakiness" score using Mean Difference (Oxford paper)
// Higher score = more peak-like
//===============================================================================
float compute_peak_score(const float *signal, int idx, int window_size, uint8_t total_samples) {
	if (idx < window_size || idx >= (total_samples - window_size))
		return 0.0f;

	float score = 0.0f;
	int count = 0;

	// Mean difference to left and right
	for (int k = 1; k <= window_size; k++) {
		if (idx - k >= 0) {
			score += (signal[idx] - signal[idx - k]);
			count++;
		}
		if (idx + k < total_samples) {
			score += (signal[idx] - signal[idx + k]);
			count++;
		}
	}

	return (count > 0) ? (score / count) : 0.0f;
}

//===============================================================================
// @brief Validate peak prominence and width (production-grade)
//===============================================================================
uint8_t validate_peak(
		const float *signal, int idx, uint8_t total_samples,
		float min_prominence, uint8_t is_idle) {

	if (idx < 1 || idx >= (total_samples - 1))
		return 0;

	float current = signal[idx];
	float left = signal[idx - 1];
	float right = signal[idx + 1];

	// Basic peak shape
	if (!(current > left && current > right))
		return 0;

	// Peak width (must be at least 3 samples wide)
	int width = 1;
	for (int i = 1; (idx - i) >= 0 && signal[idx - i] > signal[idx] * 0.5f; i++)
		width++;
	for (int i = 1; (idx + i) < total_samples && signal[idx + i] > signal[idx] * 0.5f; i++)
		width++;

	if (width < PEAK_WIDTH_MIN_SAMPLES)
		return 0;

	// Prominence (peak height above surrounding valley)
	float left_valley = signal[idx - 1];
	float right_valley = signal[idx + 1];
	for (int j = -4; j <= -1; j++) {
		int jj = idx + j;
		if (jj >= 0 && signal[jj] < left_valley)
			left_valley = signal[jj];
	}
	for (int j = 1; j <= 4; j++) {
		int jj = idx + j;
		if (jj < total_samples && signal[jj] < right_valley)
			right_valley = signal[jj];
	}

	float prominence = current - (left_valley > right_valley ? left_valley : right_valley);

	if (prominence < min_prominence)
		return 0;

	// Asymmetry check (peak shouldn't be perfectly symmetric = noise)
	float left_slope = current - left;
	float right_slope = current - right;
	float ratio = (left_slope > right_slope) ? (right_slope / left_slope) : (left_slope / right_slope);

	if (ratio < PEAK_ASYMMETRY_RATIO)
		return 0;

	return 1;
}

//===============================================================================
// @brief Update activity state machine (production-grade hysteresis)
//===============================================================================
void update_activity_state(float accel_mag, float gyro_rms) {
	uint8_t is_motion = (accel_mag > TRANSITION_THRESHOLD_G || gyro_rms > GYRO_IDLE_THRESHOLD_DPS);

	switch (g_activity_state.state) {
		case 0:  // ACTIVE
			if (!is_motion) {
				g_activity_state.exit_cnt++;
				if (g_activity_state.exit_cnt >= IDLE_ENTRY_BLOCKS)
					g_activity_state.state = 2;  // → IDLE
			} else {
				g_activity_state.exit_cnt = 0;
			}
			break;

		case 2:  // IDLE
			if (is_motion) {
				g_activity_state.entry_cnt++;
				if (g_activity_state.entry_cnt >= IDLE_EXIT_BLOCKS)
					g_activity_state.state = 0;  // → ACTIVE
			} else {
				g_activity_state.entry_cnt = 0;
			}
			break;
	}
}

//===============================================================================
// @brief Detect wrist raise (production-grade, enhanced)
//===============================================================================
void detect_wrist_raise_v2(
		const float *ax, const float *ay, const float *az,
		const float *gx, const float *gy, const float *gz,
		uint8_t total_samples) {

	if (!g_st_Customer_Settings.display_went_to_sleep)
		return;  // Only when display is asleep

	int raise_count = 0;
	float avg_pitch = 0.0f, avg_roll = 0.0f, avg_gy = 0.0f;

	for (int i = 0; i < total_samples; i++) {
		float pitch = atan2f(ax[i], sqrtf(ay[i]*ay[i] + az[i]*az[i])) * 180.0f / M_PI;
		float roll = atan2f(ay[i], sqrtf(ax[i]*ax[i] + az[i]*az[i])) * 180.0f / M_PI;
		float accel_z_norm = az[i];

		// Check all criteria
		if (pitch >= RAISE_PITCH_MIN_DEG && pitch <= RAISE_PITCH_MAX_DEG &&
			fabsf(roll) <= RAISE_ROLL_MAX_DEG &&
			accel_z_norm > (RAISE_Z_MIN_MG / 1000.0f) &&
			gy[i] > RAISE_GYRO_MIN_DPS && gy[i] < RAISE_GYRO_MAX_DPS) {

			raise_count++;
			avg_pitch += pitch;
			avg_roll += roll;
			avg_gy += gy[i];
		}
	}

	if (raise_count >= RAISE_DETECT_MIN_COUNT) {
		avg_pitch /= raise_count;
		avg_roll /= raise_count;
		avg_gy /= raise_count;

		g_wrist_state.motion_count++;
		g_wrist_state.last_pitch = avg_pitch;
		g_wrist_state.last_roll = avg_roll;

		if (g_wrist_state.motion_count >= RAISE_SUSTAIN_BLOCKS) {
			g_wrist_state.is_active = 1;
			g_wrist_state.sustain_count = 1;
			is_wrist_raise = 1;
		}
	} else {
		g_wrist_state.motion_count = 0;
	}

	// Cooldown handling
	if (g_wrist_state.is_active) {
		g_wrist_state.cooldown_count++;
		if (g_wrist_state.cooldown_count >= RAISE_COOLDOWN_BLOCKS) {
			g_wrist_state.is_active = 0;
			g_wrist_state.cooldown_count = 0;
			is_wrist_raise = 0;
		}
	}
}

//===============================================================================
// @brief Main step counter algorithm (production-grade)
// Processes 48 samples for >95% accuracy with reduced nuisance steps
//===============================================================================
//===============================================================================
// @brief REFACTORED STEP COUNTER - SIMPLIFIED PIPELINE
// Architecture: Raw XYZ → Magnitude → Bandpass → Stats → Threshold → Peaks →
//              Global State → Cadence Validation → Steps (@ 52 Hz)
//===============================================================================
int process_step_block_remaster(SENSOR_DATA_F* accel_data, SENSOR_DATA_F* gyro_data, uint8_t total_samples) {
	static float magnitude[BLOCK_SIZE], filtered_mag[BLOCK_SIZE];
	static float gyro_x[BLOCK_SIZE], gyro_y[BLOCK_SIZE], gyro_z[BLOCK_SIZE];
	static float gyro_mag[BLOCK_SIZE];
	static uint16_t idle_reject_count = 0;
	static uint16_t energy_reject_count = 0;
	static uint16_t total_blocks = 0;

	int block_step_count = 0;
	total_blocks++;

	// ===== PHASE 1: Copy sensor data and compute gyro Z-axis RMS =====
	float gyro_z_rms = 0.0f;
	for (int i = 0; i < total_samples; ++i) {
		gyro_x[i] = gyro_data[i].x;
		gyro_y[i] = gyro_data[i].y;
		gyro_z[i] = gyro_data[i].z;
		gyro_z_rms += gyro_z[i] * gyro_z[i];
	}
	gyro_z_rms = sqrtf(gyro_z_rms / total_samples);

	// ===== PHASE 1b: WRIST FIDGET REJECTION - DISABLED =====
	// NOTE: Cannot reliably distinguish walking arm swing from fidgeting using Z-gyro alone
	// Both produce high Z-axis rotation (50,000-100,000+ mdps)
	// Walking arm swing is NECESSARY for wrist-worn step detection
	// Idle detection + energy checks provide sufficient filtering
	// TODO: Implement gyro+accel correlation for future improvement

	// Disabled for wrist-worn accuracy:
	// if (gyro_z_rms > GYRO_Z_FIDGET_THRESHOLD) { return 0; }

	// ===== PHASE 2: Compute magnitude directly from raw accel =====
	// Do NOT filter individual axes - compute magnitude first!
	for (int i = 0; i < total_samples; ++i) {
		magnitude[i] = sqrtf(accel_data[i].x * accel_data[i].x +
		                      accel_data[i].y * accel_data[i].y +
		                      accel_data[i].z * accel_data[i].z);
		gyro_mag[i] = sqrtf(gyro_x[i]*gyro_x[i] + gyro_y[i]*gyro_y[i] + gyro_z[i]*gyro_z[i]);
	}

	// ===== PHASE 3: Bandpass filter (0.5-4 Hz) removes DC component =====
	// No separate HP filter needed - BP removes DC at 0.5 Hz
	bandpass_filtfilt(magnitude, filtered_mag);

	// ===== PHASE 4: Compute statistics (accel magnitude only) =====
	float mean = 0.0f, std = 0.0f, energy = 0.0f;

	// Mean
	for (int i = 0; i < total_samples; ++i) {
		mean += filtered_mag[i];
	}
	mean /= total_samples;

	// Std and Energy (for step detection threshold and idle check)
	for (int i = 0; i < total_samples; ++i) {
		float diff = filtered_mag[i] - mean;
		std += diff * diff;
		energy += filtered_mag[i] * filtered_mag[i];
	}
	std = sqrtf(std / total_samples);
	energy /= total_samples;

	// ===== PHASE 10: Idle detection (accel magnitude only) =====
	// If acceleration signal has very low variation, device is idle/stationary
	// Walking requires minimum energy in the filtered signal
	// Energy thresholds in mg² (after conversion to mg units)
	// Idle: std < 8 mg or energy < 10 (mg)²
	uint8_t is_idle = (std < 8.0f || energy < 10.0f);

	if (is_idle) {
#if (DEBUG_SAMPLES)
		printf("[BLOCK IDLE] std=%.1f mg, energy=%.0f mg² (too low for walking) → returning 0 steps\n", std, energy);
#endif
		idle_reject_count++;
		return 0;
	}

	// ===== PHASE 10b: Reject extreme fidgeting (too much energy) =====
	// Energy is in (mg)² = 1,000,000 × (g)²
	// Normal walking: ~30,000-160,000 mg² (0.03-0.16 g²)
	// Vigorous walking/running: ~100,000-300,000 mg² (normal for active movement)
	// Device drop or extreme shock: >500,000 mg² (unrealistic)
	// Threshold raised to 1.5M to allow jogging (600k-1M mg²) and vigorous running (up to 1.2M mg²)
	// Only reject if energy exceeds physical device limits or true motion artifacts
	if (energy > 1500000.0f) {
#if (DEBUG_SAMPLES)
		printf("[EXTREME ENERGY] energy=%.0f mg² > 1.5M → returning 0 steps\n", energy);
#endif
		// Reset consecutive peaks when rejecting extreme-energy block
		energy_reject_count++;
		g_step_state.consecutive_peaks = 0;
		return 0;
	}

	// ===== PHASE 5: Single adaptive threshold =====
	// threshold = mean + 0.8 * std  (commercial algorithms use 0.7-1.0)
	float threshold = mean + 0.8f * std;

#if (DEBUG_SAMPLES)
	printf("[Block] std=%.4f energy=%.4f mean=%.3f threshold=%.3f\n",
		std, energy, mean, threshold);
#endif

	// ===== PHASE 5b: Detect peaks - find all local maxima with scoring =====
	// Find ALL local maxima and compute peakiness score (Oxford paper method)
	// Use Mean Difference scoring to enhance peak quality assessment
	int peak_indices[BLOCK_SIZE];
	float peak_scores[BLOCK_SIZE];
	int peak_count_in_block = 0;

	int score_window = 4;  // Window size for peakiness calculation

	for (int i = 1; i < total_samples - 1; ++i) {
		// Simple local maximum
		if (filtered_mag[i] > filtered_mag[i - 1] &&
			filtered_mag[i] > filtered_mag[i + 1]) {

			// Compute peakiness score (Oxford Mean Difference method)
			float score = compute_peak_score(filtered_mag, i, score_window, total_samples);

			// Be permissive: Accept all local maxima (score filtering done in cadence check)
			// This catches small-amplitude valid steps that have good cadence
			peak_indices[peak_count_in_block] = i;
			peak_scores[peak_count_in_block] = score;
			peak_count_in_block++;
		}
	}

#if (DEBUG_SAMPLES)
	if (peak_count_in_block > 0) {
		printf("[Peaks] %d detected in block (scored): ", peak_count_in_block);
		for (int i = 0; i < peak_count_in_block; i++) {
			printf("%d(%.3f) ", peak_indices[i], peak_scores[i]);
		}
		printf("\n");
	}
#endif

	// ===== PHASE 7: Global state processing =====
	// This is the CORE of the algorithm: maintain state across blocks
	// Do NOT reset g_step_state between blocks!

	// Process each peak found in this block
	for (int i = 0; i < peak_count_in_block; ++i) {
		int local_idx = peak_indices[i];
		int global_idx = g_step_state.sample_idx_global + local_idx;
		int distance_since_last_peak = global_idx - g_step_state.last_peak_global;
		float peak_score = peak_scores[i];

#if (DEBUG_SAMPLES)
		printf("[PeakAnalysis] idx=%d global=%d score=%.3f\n",
			local_idx, global_idx, peak_score);
#endif

		// ===== PHASE 10b: Peak quality filter =====
		// Reject low-score peaks (noise, weak artifacts from filtering)
		// High-quality steps have strong peakiness signatures
		// Threshold: 5.0 optimized for broad cadence range
		// Scores 0-5: obvious noise, 5+: accept if cadence valid
		// Cadence validation (17-75 samples) provides primary filtering
		float min_peak_score = 5.0f;  // Very lenient: cadence validation is primary filter

		if (peak_score < min_peak_score) {
#if (DEBUG_SAMPLES)
			printf("  ↳ [REJECTED] Low peak quality: score=%.1f < min=%.1f (weak peak/noise)\n",
				peak_score, min_peak_score);
#endif
			continue;
		}

		// ===== PHASE 11: Refractory period (block boundary protection) =====
		// Skip peaks that occur within 154ms (~8 samples @ 52 Hz) of the last accepted peak
		// This prevents double-counting the same step while allowing multi-peak jogging signal
		// At 52 Hz: 154ms = 8 samples (allows faster cadences like running/jogging)
		// Note: Cadence range (10-75) provides primary validation, refractory is secondary
		if (distance_since_last_peak < 8) {
#if (DEBUG_SAMPLES)
			printf("  ↳ [REJECTED] Refractory period: dist=%d < 8 samples (double-count protection)\n",
				distance_since_last_peak);
#endif
			continue;
		}

		// ===== PHASE 8: Global Cadence Validation =====
		// Single permissive range works across all activities:
		// - Slow walking: ~39 samples (75 BPM)
		// - Normal walking: ~28 samples (111 BPM)
		// - Brisk walking: ~22 samples (142 BPM)
		// - Jogging: ~18 samples (173 BPM)
		// - Running: ~17 samples (183 BPM)
		//
		// At 52 Hz per sample interval:
		//   10 samples = 192ms (312 BPM) - fast running
		//   20 samples = 385ms (156 BPM) - brisk walk
		//   30 samples = 577ms (104 BPM) - normal walk
		//   40 samples = 769ms (78 BPM) - slow walk
		//   75 samples = 1442ms (42 BPM) - very slow walk

		// Track step intervals for activity detection (debug info only)
		update_step_interval(distance_since_last_peak);

		// Every 10 valid peaks, show detected activity for diagnostics
		if ((g_step_state.consecutive_peaks % 10) == 0 && g_step_state.interval_count >= 5) {
			uint16_t avg_interval = get_average_step_interval();
			g_step_state.detected_activity = detect_activity_from_rate(avg_interval);

#if (DEBUG_SAMPLES)
			const char* activity_names[] = {"UNKNOWN", "SLOW_WALK", "NORMAL_WALK", "BRISK_WALK", "JOGGING", "RUNNING"};
			float bpm = (avg_interval > 0) ? (3120.0f / avg_interval) : 0;
			printf("[ACTIVITY_INFO] %s (avg_interval=%d samples, ~%.0f BPM)\n",
				activity_names[g_step_state.detected_activity], avg_interval, bpm);
#endif
		}

		// Single global cadence range (permissive to catch all activities)
		uint8_t cadence_min = 10;   // 192ms - allows fast running
		uint8_t cadence_max = 75;   // 1442ms - allows very slow walking

		uint8_t cadence_valid = (distance_since_last_peak >= cadence_min &&
								 distance_since_last_peak <= cadence_max);

		if (!cadence_valid) {
#if (DEBUG_SAMPLES)
			if (distance_since_last_peak < cadence_min) {
				printf("  ↳ [REJECTED] TOO FAST: dist=%d < min=%d (noise or arm swing)\n",
					distance_since_last_peak, cadence_min);
			} else {
				printf("  ↳ [REJECTED] TOO SLOW: dist=%d > max=%d (probably noise or missed peak)\n",
					distance_since_last_peak, cadence_max);
			}
#endif
			// Reset counter only for unconfirmed patterns
			if (g_step_state.consecutive_peaks < 5) {
				g_step_state.consecutive_peaks = 0;
#if (DEBUG_SAMPLES)
				printf("  ↳ [RESET] consecutive_peaks → 0\n");
#endif
			}
			g_step_state.last_peak_global = global_idx;
			continue;
		}

		// ===== Valid peak with valid cadence =====
		g_step_state.consecutive_peaks++;
		g_step_state.last_peak_global = global_idx;
		g_step_state.last_step_sample = global_idx;

#if (DEBUG_SAMPLES)
		printf("  ↳ [ACCEPTED] Valid cadence: dist=%d samples, consecutive=%d, score=%.1f\n",
			distance_since_last_peak, g_step_state.consecutive_peaks, peak_score);
#endif

		// ===== PHASE 8b: Count step immediately after 1st valid cadence peak =====
		// CHANGED: Count every peak that passes cadence check
		// This prevents losing valid steps waiting for "confirmation"
		// Cadence validation itself acts as the filter
		if (g_step_state.consecutive_peaks >= 1) {
			block_step_count++;
#if (DEBUG_SAMPLES)
			printf("  ↳ [STEP_COUNTED] peak_score=%.1f, dist=%d samples\n",
				peak_score, distance_since_last_peak);
#endif
		}
	}

	// ===== Update global sample index for next block =====
	// IMPORTANT: This advances the global time counter
	// At 52 Hz with 48 samples per block: 48 samples = 923 ms (0.923 seconds)
	g_step_state.sample_idx_global += total_samples;

	// ===== Populate debug structure for SD card logging (multiply by 10000 for int32_t) =====
	g_step_debug.mean = (int32_t)(mean);
	g_step_debug.std = (int32_t)(std);
	g_step_debug.energy = (int32_t)(energy);
	g_step_debug.threshold = (int32_t)(threshold);
	g_step_debug.peak_count = peak_count_in_block;
	g_step_debug.consecutive_peaks = g_step_state.consecutive_peaks;
	g_step_debug.block_steps = block_step_count;
	g_step_debug.last_distance = (g_step_state.last_peak_global > 0) ?
		(g_step_state.sample_idx_global - g_step_state.last_peak_global) : 0;
	g_step_debug.is_idle = is_idle;

#if (DEBUG_SAMPLES)
	printf("[BlockEnd] sample_idx_global now=%d steps_this_block=%d\n",
		g_step_state.sample_idx_global, block_step_count);

	// Print diagnostic summary every 100 blocks
	if ((total_blocks % 100) == 0) {
		float idle_percent = (idle_reject_count * 100.0f) / total_blocks;
		float energy_percent = (energy_reject_count * 100.0f) / total_blocks;
		printf("\n[DIAGNOSTIC] After %d blocks:\n", total_blocks);
		printf("  - Idle rejection: %d blocks (%.1f%%)\n", idle_reject_count, idle_percent);
		printf("  - Extreme energy rejection: %d blocks (%.1f%%)\n", energy_reject_count, energy_percent);
		printf("  - Active blocks: %d (%.1f%%)\n",
			total_blocks - idle_reject_count - energy_reject_count,
			((total_blocks - idle_reject_count - energy_reject_count) * 100.0f) / total_blocks);
		printf("\n");
	}
#endif

	return block_step_count;
}
