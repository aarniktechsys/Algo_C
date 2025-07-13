#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include "stepCounter_HPF_BPF.h"

float Ax[BLOCK_SIZE] = {0};
float Ay[BLOCK_SIZE] = {0};
float Az[BLOCK_SIZE] = {0};
float Gx[BLOCK_SIZE] = {0};
float Gy[BLOCK_SIZE] = {0};
float Gz[BLOCK_SIZE] = {0};

float last_orientation_deg = 0.0f;
uint8_t is_wrist_raise = false;

void highpass_filtfilt(const float *input, float *output);
void bandpass_filtfilt(const float *input, float *output);
int process_step_block_SF(int16_t *Ax_raw, int16_t *Ay_raw, int16_t *Az_raw,
    int16_t *Gx_raw, int16_t *Gy_raw, int16_t *Gz_raw, uint8_t total_samples);

int process_step_block_SF(int16_t *Ax_raw, int16_t *Ay_raw, int16_t *Az_raw,
    int16_t *Gx_raw, int16_t *Gy_raw, int16_t *Gz_raw, uint8_t total_samples) {

    static float fx[BLOCK_SIZE], fy[BLOCK_SIZE], fz[BLOCK_SIZE];
    static float magnitude[BLOCK_SIZE], filtered_mag[BLOCK_SIZE], magnitude_gyro[BLOCK_SIZE];
    static float temp_Ax[BLOCK_SIZE] = {0}, temp_Ay[BLOCK_SIZE] = {0}, temp_Az[BLOCK_SIZE] = {0};
    static int16_t temp_Ax_mg[BLOCK_SIZE] = {0}, temp_Ay_mg[BLOCK_SIZE] = {0}, temp_Az_mg[BLOCK_SIZE] = {0};
    static float gyro_x[BLOCK_SIZE] = {0}, gyro_y[BLOCK_SIZE] = {0}, gyro_z[BLOCK_SIZE] = {0};
    float gyro_mag_sq[BLOCK_SIZE] = {0}, gyro_mag_smooth[BLOCK_SIZE] = {0};
    float pitch_arr[BLOCK_SIZE] = {0}, roll_arr[BLOCK_SIZE] = {0};

    int active_gyro_samples = 0, idle_sample_count = 0, orientation_idle_count = 0;
    int16_t steps = 0;
    float mean_z = 0;
    int orientation_stable_count = 0;

    // Load and convert raw samples
    for (int i = 0; i < total_samples; ++i) {
        temp_Ax[i] = ((float)Ax_raw[i]) / 1000.0f;
        temp_Ay[i] = ((float)Ay_raw[i]) / 1000.0f;
        temp_Az[i] = ((float)Az_raw[i]) / 1000.0f;

        gyro_x[i] = ((float)Gx_raw[i]) / 1000.0f;
        gyro_y[i] = ((float)Gy_raw[i]) / 1000.0f;
        gyro_z[i] = ((float)Gz_raw[i]) / 1000.0f;

        temp_Ax_mg[i] = (int16_t)(temp_Ax[i] * 1000.0f + 0.5f);
        temp_Ay_mg[i] = (int16_t)(temp_Ay[i] * 1000.0f + 0.5f);
        temp_Az_mg[i] = (int16_t)(temp_Az[i] * 1000.0f + 0.5f);

        mean_z += temp_Az[i];
    }
    mean_z /= total_samples;

    // Apply high-pass filter
    highpass_filtfilt(temp_Ax, temp_Ax);
    highpass_filtfilt(temp_Ay, temp_Ay);
    highpass_filtfilt(temp_Az, temp_Az);

    float ref_pitch = 0.0f, ref_roll = 0.0f;

    // Compute magnitude, orientation, idle
    for (int i = 0; i < total_samples; ++i) {
        magnitude[i] = sqrtf(temp_Ax[i]*temp_Ax[i] + temp_Ay[i]*temp_Ay[i] + temp_Az[i]*temp_Az[i]);
        magnitude_gyro[i] = (gyro_x[i]*gyro_x[i] + gyro_y[i]*gyro_y[i] + gyro_z[i]*gyro_z[i]);

        // Accel idle based on tilt thresholds
        bool tilt_only = (abs(temp_Az_mg[i] - 1000) <= 100 &&
                            abs(temp_Ax_mg[i]) < 200 &&
                            abs(temp_Ay_mg[i]) < 200);
        if (tilt_only) idle_sample_count++;

        // Compute pitch/roll per sample
        pitch_arr[i] = atan2f(temp_Ax[i], sqrtf(temp_Ay[i]*temp_Ay[i] + temp_Az[i]*temp_Az[i])) * 180.0f / M_PI;
        roll_arr[i]  = atan2f(temp_Ay[i], sqrtf(temp_Ax[i]*temp_Ax[i] + temp_Az[i]*temp_Az[i])) * 180.0f / M_PI;

        if (i < REF_SAMPLE_COUNT) {
            ref_pitch += pitch_arr[i];
            ref_roll  += roll_arr[i];
        }
    }

    ref_pitch /= REF_SAMPLE_COUNT;
    ref_roll  /= REF_SAMPLE_COUNT;

    // Orientation idle check
    for (int i = 0; i < total_samples; ++i) {
        float pitch_delta = fabsf(pitch_arr[i] - pitch_arr[i - 1]);
        float roll_delta  = fabsf(roll_arr[i] - roll_arr[i - 1]);

        if (pitch_delta < ORIENTATION_DELTA_THRESHOLD  &&
            roll_delta  < ORIENTATION_DELTA_THRESHOLD) {
                orientation_stable_count++;
        } else {
            /*printf("pitch[%d]=%.2f pitch[%d]=%.2f pitch_delta=%.2f | roll[%d]=%.2f roll[%d]=%.2f delta_roll=%.2f\n",
                i - 1, pitch_arr[i - 1], i, pitch_arr[i], pitch_delta,
                i - 1, roll_arr[i - 1], i, roll_arr[i], roll_delta);*/
        }
    }

    // Smoothed gyro magnitude
    for (int i = GYRO_WINDOW; i < BLOCK_SIZE - GYRO_WINDOW; ++i) {
        float sum = 0;
        for (int j = -GYRO_WINDOW; j <= GYRO_WINDOW; ++j)
            sum += magnitude_gyro[i + j];
        gyro_mag_smooth[i] = sum / (2 * GYRO_WINDOW + 1);
    }

    for (int i = GYRO_WINDOW; i < BLOCK_SIZE - GYRO_WINDOW; ++i) {
        if (gyro_mag_smooth[i] > GYRO_IDLE_THRESHOLD)
            active_gyro_samples++;
    }

    // Idle flags
    uint8_t is_accel_idle = (idle_sample_count >= (total_samples - 8));
    uint8_t is_gyro_idle  = (active_gyro_samples <= MAX_ACTIVE_GYRO_SAMPLES);
    uint8_t is_orientation_idle = (orientation_stable_count >= ORIENTATION_IDLE_MIN_SAMPLES);

    printf("is_accel_idle= %d, is_gyro_idle= %d, is_orientation_idle= %d\n", is_accel_idle, is_gyro_idle, is_orientation_idle);
    printf("idle_sample_count = %d, orientation_idle_count = %d\n", idle_sample_count, orientation_idle_count);

    if ((is_orientation_idle || is_accel_idle) && is_gyro_idle) {
        printf("Watch is idle\n");
        return 0;
    }

    // Band-pass filter for step detection
    bandpass_filtfilt(magnitude, filtered_mag);

    float mean = 0.0f, std = 0.0f;
    for (int i = 0; i < BLOCK_SIZE; ++i)
        mean += filtered_mag[i];
    mean /= BLOCK_SIZE;

    for (int i = 0; i < BLOCK_SIZE; ++i)
        std += (filtered_mag[i] - mean) * (filtered_mag[i] - mean);
    std = sqrtf(std / BLOCK_SIZE);

    float energy = 0.0f;
    for (int i = 0; i < BLOCK_SIZE; ++i)
        energy += filtered_mag[i] * filtered_mag[i];
    energy /= BLOCK_SIZE;

   /* if (std < 0.015f || energy < 0.02f || (mean_z > 950.0f && mean_z < 1050.0f && std < 0.02f))
        return 0;*/

    //float threshold = mean + 0.15f * std;
    //float threshold = mean + 1.5f * std;
    //if (threshold < 0.08f) threshold = 0.08f;
    float base_threshold = mean + 1.2f * std;    

    /*int last_peak_index = -1000;
    const int min_peak_distance = 3;
    int peak_count = 0;*/

    int peak_indices[BLOCK_SIZE];
    int peak_count = 0;

    for (int i = 1; i < BLOCK_SIZE - 1; i++) {
        if (filtered_mag[i] > base_threshold &&
            filtered_mag[i] > filtered_mag[i - 1] &&
            filtered_mag[i] > filtered_mag[i + 1]) {

            // Save peak index
            peak_indices[peak_count++] = i;
        }
    }

    for (int i = 0; i < peak_count; i++){
        printf("peak_indices=%d", peak_indices[i]);
    }
    printf("\n");

    // Estimate average step interval from detected peaks
    int total_interval = 0;
    int valid_intervals = 0;
    for (int i = 1; i < peak_count; i++) {
        int interval = peak_indices[i] - peak_indices[i - 1];
        if (interval >= 3 && interval <= 40) {
            total_interval += interval;
            valid_intervals++;
        }
    }

    float avg_interval = (valid_intervals > 0) ? ((float)total_interval / valid_intervals) : 0.0f;

    float k = (avg_interval > 18.0f) ? 0.9f : 1.5f;
    float threshold = mean + k * std;
    float threshold_floor = 0.3f * std; 

    // Dynamically select min_peak_distance
    int min_peak_distance = 20;
    if (avg_interval > 35.0f) {  // slow walking (~1.5 Hz)
        min_peak_distance = 30;  
    } else if (avg_interval > 25.0f) {
        min_peak_distance = 24;   // normal walking (~2 Hz)
    } else if (avg_interval > 18.0f) {
        min_peak_distance = 18;   // Brisk walking (~2.8 Hz)
    } else if (avg_interval > 12.0f) {
        min_peak_distance = 14;   // jogging (~4 Hz)
    }else{
        min_peak_distance = 10;  // fast running (~5.2 Hz max)
    }

    if (min_peak_distance < 12) min_peak_distance = 12;

    printf("avg_interval= %f\n",avg_interval);
    printf("min_peak_distance= %d\n",min_peak_distance);

    if (avg_interval > 14.0f && 
        (std < 0.015f || energy < 0.02f || 
        (mean_z > 950.0f && mean_z < 1050.0f && std < 0.02f))){
            return 0;
     }        

    int last_peak_index = -min_peak_distance;   

    for (int i = 1; i < BLOCK_SIZE - 1; ++i) {
        if (filtered_mag[i] > threshold &&
            filtered_mag[i] > filtered_mag[i - 1] &&
            filtered_mag[i] > filtered_mag[i + 1] &&
            filtered_mag[i] > threshold_floor) {

            float current = filtered_mag[i];
            float prev = filtered_mag[i - 1];
            float next = filtered_mag[i + 1];            

            /*if ((i - last_peak_index) >= min_peak_distance){ // && filtered_mag[i] > 0.2f) {
                peak_count++;
                last_peak_index = i;
            }*/
            if ((i - last_peak_index) >= min_peak_distance ||
                ((i - last_peak_index) == min_peak_distance - 1 && current > 1.2f * threshold)) {
                peak_count++;
                last_peak_index = i;
            }
        }
    }       

    return peak_count;
}
    


