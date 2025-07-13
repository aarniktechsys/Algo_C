#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <stdint.h>
#include "stepCounter.h"


int PEAK_MIN_DIST = 24;
int PEAK_THRESHOLD = 225;
int DIP_DROP = 50;

// Bandpass FIR filter coefficients for 1–4 Hz range
const int16_t fir_bandpass_1_4hz[FIR_LEN] = {
    -5, -6, -6, -4, 0, 7, 16, 27, 37, 44,
     45, 38, 23, 0, -27, -54, -73, -77, -63, -33,
      9, 54, 93, 117, 118, 95, 53, 0, -52, -95,
    -117, -113, -84
};

int16_t accel_fir_buf[FIR_LEN] = {0};
int16_t gyro_fir_buf[FIR_LEN] = {0};
int fir_index = 0;

//int last_step_sample = -PEAK_MIN_DIST;
int step_count = 0;
int16_t fused_mag_buf[MAG_BUF_LEN];
int mag_buf_index = 0;

int sample_counter = 0;
int last_step_sample = -30;

int recent_step_deltas[MAX_RECENT_STEPS] = {14, 14, 14, 14, 14};
int step_delta_index = 0;

static float filtered_interval = 20.0f;

static int32_t a_avg = 0, g_avg = 0;


int16_t process_sample(int Ax, int Ay, int Az, int Gx, int Gy, int Gz);
static int16_t apply_fir_filter(int16_t *buffer, int16_t sample);
static int16_t compute_magnitude(int16_t x, int16_t y, int16_t z);
static void update_adaptive_thresholds(void);

//===============================================================================
// @brief  FIR filter
//
// @param[in]    : Accelerometer/Gyroscope X,Y,Z
// @param[out]   :
// @param[inout] :
// @retval       :
//
//
//===============================================================================
static int16_t apply_fir_filter(int16_t *buffer, int16_t sample) {
    static const int16_t fir_coeffs[FIR_LEN] = {83, 135, 178, 198, 178, 135, 83 };
    buffer[fir_index] = sample;
    int32_t acc = 0;
    for (int i = 0; i < FIR_LEN; i++) {
        int idx = (fir_index + i) % FIR_LEN;
        acc += fir_coeffs[i] * buffer[idx];
    }
    fir_index = (fir_index + 1) % FIR_LEN;
    return (int16_t)(acc / Q_SCALE);
}

//===============================================================================
// @brief  compute_magnitude
//
// @param[in]    : Accelerometer/Gyroscope X,Y,Z
// @param[out]   :
// @param[inout] :
// @retval       :
//
//
//===============================================================================
static int16_t compute_magnitude(int16_t x, int16_t y, int16_t z) {
    int64_t mag_sq = (int64_t)x * x + (int64_t)y * y + (int64_t)z * z;
	int64_t low = 0, high = mag_sq, mid, sqrt_val = 0;

	while (low <= high) {
		mid = (low + high) / 2;
		if (mid * mid <= mag_sq) {
			sqrt_val = mid;
			low = mid + 1;
		} else {
			high = mid - 1;
		}
	}
	return (int16_t)sqrt_val;
}

//===============================================================================
// @brief  Compute threshold value
//
// @param[in]    : 
// @param[out]   :
// @param[inout] :
// @retval       :
//
//
//===============================================================================
static void update_adaptive_thresholds(void) {
	int total = 0;
	int valid_steps = 0;

	for (int i = 0; i < MAX_RECENT_STEPS; i++) {
		if (recent_step_deltas[i] > 0) {
			total += recent_step_deltas[i];
			valid_steps++;
		}
	}

	if (valid_steps < MAX_RECENT_STEPS) return;

	int avg_interval = total / valid_steps;

	if (avg_interval < 14) avg_interval = 14;
	if (avg_interval > 40) avg_interval = 40;

	//printf("avg_interval = %d\n", avg_interval);

	if (avg_interval < 16) {
		PEAK_MIN_DIST = 20;
		PEAK_THRESHOLD = 475;
		DIP_DROP = 100;
	} else if (avg_interval < 26) {
		PEAK_MIN_DIST = 24;
		PEAK_THRESHOLD = 440;
		DIP_DROP = 90;
	} else {
		PEAK_MIN_DIST = 16;
		PEAK_THRESHOLD = 400;
		DIP_DROP = 60;
	}
	if (PEAK_THRESHOLD < 400) PEAK_THRESHOLD = 400;	
}

//===============================================================================
// @brief  Step counter algorithm
//
// @param[in]    : Accelerometer/Gyroscope X,Y,Z
// @param[out]   :
// @param[inout] :
// @retval       :
//
//
//===============================================================================
int16_t process_sample(int Ax, int Ay, int Az, int Gx, int Gy, int Gz) {
    static int32_t idle_accel_sum = 0, idle_gyro_sum = 0;
    static int idle_counter = 0;  
    int steps = 0; 

    // Compute magnitude of accelerometer and gyroscope
    int a_mag = compute_magnitude(Ax, Ay, Az);
    int g_mag = compute_magnitude(Gx, Gy, Gz);

    // Idle suppression logic
	idle_accel_sum += a_mag;
	idle_gyro_sum += g_mag;
	idle_counter++;

	if (idle_counter >= IDLE_WINDOW_SIZE) {
		int32_t a_mean = idle_accel_sum / IDLE_WINDOW_SIZE;
		int32_t g_mean = idle_gyro_sum / IDLE_WINDOW_SIZE;
		idle_accel_sum = 0;
		idle_gyro_sum = 0;
		idle_counter = 0;
		if (a_mean < IDLE_ACCEL_THRESH && g_mean < IDLE_GYRO_THRESH) {
			return 0;
		}
	}

    a_avg = (a_avg * 15 + a_mag) / 16;
	g_avg = (g_avg * 15 + g_mag) / 16;
	int16_t a_hp = a_mag - a_avg;
	int16_t g_hp = g_mag - g_avg;

    // FIR filters
	int16_t a_filtered = apply_fir_filter(accel_fir_buf, a_hp);
	int16_t g_filtered = apply_fir_filter(gyro_fir_buf, g_hp);

    // Fusion and clamping
	int16_t fused = (3 * a_filtered + g_filtered) / 4;
	int16_t abs_fused = (fused > 0) ? fused : -fused;
	if (abs_fused > 2000) abs_fused = 2000;

	fused_mag_buf[mag_buf_index % MAG_BUF_LEN] = abs_fused;

	update_adaptive_thresholds();

    if (sample_counter > FIR_LEN + DIP_WINDOW) {
		int i = mag_buf_index % MAG_BUF_LEN;
		int i_prev = (mag_buf_index - 1 + MAG_BUF_LEN) % MAG_BUF_LEN;
		int i_next = (mag_buf_index + 1) % MAG_BUF_LEN;

		if (fused_mag_buf[i] > PEAK_THRESHOLD &&
			fused_mag_buf[i_prev] < fused_mag_buf[i] &&
			fused_mag_buf[i] > fused_mag_buf[i_next] &&
			(sample_counter - last_step_sample) >= PEAK_MIN_DIST) {

			int peak_val = fused_mag_buf[i];
			bool dip_found = false;
			for (int j = 1; j <= DIP_WINDOW; j++) {
				int dip_idx = (mag_buf_index + j) % MAG_BUF_LEN;
				if (fused_mag_buf[dip_idx] < (peak_val - DIP_DROP)) {
					dip_found = true;
					break;
				}
			}

			int delta = sample_counter - last_step_sample;
			recent_step_deltas[step_delta_index % MAX_RECENT_STEPS] = delta;
			//printf("Recent step deta = %d\n", recent_step_deltas[step_delta_index % MAX_RECENT_STEPS]);
			step_delta_index++;

			int total = 0, valid = 0;
			for (int k = 0; k < MAX_RECENT_STEPS; k++) {
				if (recent_step_deltas[k] > 0 && recent_step_deltas[k] < 64) {
					total += recent_step_deltas[k];
					valid++;
				}
			}

			int mean = valid ? total / valid : 0;
			int variance = 0;
			for (int k = 0; k < MAX_RECENT_STEPS; k++) {
				if (recent_step_deltas[k] > 0 && recent_step_deltas[k] < 64) {
					int diff = recent_step_deltas[k] - mean;
					variance += diff * diff;
				}
			}

			if (dip_found) {
				(steps)++;
				last_step_sample = sample_counter;
			}
		}
	}

	mag_buf_index++;
	sample_counter++;

	return steps;    
}