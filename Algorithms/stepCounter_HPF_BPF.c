#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include "stepCounter_HPF_BPF.h"

// High-pass filter (0.3Hz @ 52Hz)
const float b_hp[HP_TAPS] = { 0.953739f, -3.814958f, 5.722436f, -3.814958f, 0.953739f };
const float a_hp[HP_TAPS] = { 1.0f, -3.905279f, 5.720297f, -3.724635f, 0.909619f };

#if(0)
// Band-pass filter (0.3–8Hz @ 52Hz)
const float b_bp[BP_TAPS] = {
    0.01780133f, 0.0f, -0.07120533f, 0.0f, 0.10680800f,
    0.0f, -0.07120533f, 0.0f, 0.01780133f
};
const float a_bp[BP_TAPS] = {
    1.0f, -5.495348f, 13.258485f, -18.532786f, 16.566255f,
    -9.734108f, 3.659269f, -0.800733f, 0.078967f
};
#endif

#if(1)
// Band-pass filter (0.5–4.5Hz @ 52Hz)
const float b_bp[9] = {
    0.00194832f, 0.0f, -0.00779327f, 0.0f, 0.01168990f,
    0.0f, -0.00779327f, 0.0f, 0.00194832f
};

const float a_bp[9] = {
    1.0f, -6.62984475f, 19.36053485f, -32.55730178f,
    34.50841663f, -23.61716055f, 10.19337717f, -2.53669815f, 0.27867724f
};
#endif

typedef enum {
    ORIENT_UNKNOWN = 0,
    ORIENT_X_POS,
    ORIENT_X_NEG,
    ORIENT_Y_POS,
    ORIENT_Y_NEG,
    ORIENT_Z_POS,
    ORIENT_Z_NEG
} ORIENTATION6D;

void highpass_filtfilt(const float *input, float *output);
void bandpass_filtfilt(const float *input, float *output);
int detect_steps(const float *data, float threshold, int min_distance, float min_prominence);
int process_step_block(int16_t *Ax_raw, int16_t *Ay_raw, int16_t *Az_raw,
                        int16_t *Gx_raw, int16_t *Gy_raw, int16_t *Gz_raw,
                        int16_t* recent_step_counter);
float compute_pitch(float ax, float ay, float az);

// Forward-backward high-pass filter
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

// Forward-backward band-pass filter
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

// Peak detection (threshold + distance + prominence)
int detect_steps(const float *data, float threshold, int min_distance, float min_prominence){
    int steps = 0, last_peak = -min_distance;

    for (int i = 1; i < N - 1; ++i) {
        if (data[i] > data[i-1] && data[i] > data[i+1] && data[i] >= threshold) {
            float lbase = data[i], rbase = data[i];
            for (int j = i-1; j >= 0 && data[j] < lbase; --j) lbase = data[j];
            for (int j = i+1; j < N && data[j] < rbase; ++j) rbase = data[j];
            float prom = data[i] - fmaxf(lbase, rbase);
            if (prom < min_prominence) continue;

            if (i - last_peak >= min_distance) {
                steps++;
                last_peak = i;
            }
        }
    }
    return steps;
}

float compute_pitch(float ax, float ay, float az) {
    return atan2f(-ax, sqrtf(ay * ay + az * az)) * (180.0f / PI);
}

ORIENTATION6D get_6D_orientation(float ax, float ay, float az) {
    float norm = sqrtf(ax * ax + ay * ay + az * az);
    if (norm < MIN_NORM_THRESHOLD) norm = 1.0f;

    float nx = ax / norm;
    float ny = ay / norm;
    float nz = az / norm;

    float abs_x = fabsf(nx);
    float abs_y = fabsf(ny);
    float abs_z = fabsf(nz);

    if (abs_x > abs_y && abs_x > abs_z) {
        return (nx > 0) ? ORIENT_X_POS : ORIENT_X_NEG;
    } else if (abs_y > abs_x && abs_y > abs_z) {
        return (ny > 0) ? ORIENT_Y_POS : ORIENT_Y_NEG;
    } else if (abs_z > abs_x && abs_z > abs_y) {
        return (nz > 0) ? ORIENT_Z_POS : ORIENT_Z_NEG;
    } else {
        return ORIENT_UNKNOWN;
    }
}

// === MAIN ENTRY POINT FOR EACH 48-SAMPLE BLOCK ===
int process_step_block(int16_t *Ax_raw, int16_t *Ay_raw, int16_t *Az_raw,
                        int16_t *Gx_raw, int16_t *Gy_raw, int16_t *Gz_raw,
                        int16_t* recent_step_counter) {
    float Ax[N], Ay[N], Az[N], mag[N], filtered_mag[N],pitch[N],gyro_x[N],gyro_y[N];
    float mean = 0.0f, std = 0.0f;
    bool is_wrist_raise = false;
    int steps = 0;
    int hold_count = 0;
    float norm = 0;

    // Scale input (from mg to g)
    for (int i = 0; i < N; ++i) {
        Ax[i] = ((float)Ax_raw[i]) / 1000.0f;
        Ay[i] = ((float)Ay_raw[i]) / 1000.0f;
        Az[i] = ((float)Az_raw[i]) / 1000.0f;
        gyro_x[i] = (float)Gx_raw[i] / 1000.0f;
		gyro_y[i] = (float)Gy_raw[i] / 1000.0f;
    }

    // High-pass filter each axis
    highpass_filtfilt(Ax, Ax);
    highpass_filtfilt(Ay, Ay);
    highpass_filtfilt(Az, Az);

    // Compute magnitude
    for (int i = 0; i < N; ++i){
        mag[i] = sqrtf(Ax[i]*Ax[i] + Ay[i]*Ay[i] + Az[i]*Az[i]);                
    }

    // Band-pass filter magnitude
    bandpass_filtfilt(mag, filtered_mag);

    // Compute mean and std of filtered magnitude
    for (int i = 0; i < N; ++i) mean += filtered_mag[i];
    mean /= N;

    for (int i = 0; i < N; ++i) std += (filtered_mag[i] - mean) * (filtered_mag[i] - mean);
    std = sqrtf(std / N);

    float threshold = mean + 0.15f * std;

    // Detect steps
    //return detect_steps(filtered_mag, threshold, 9 /*distance*/, 0.01f /*prominence*/);
    int last_peak_index = -1000;
    const int min_peak_distance = 3;
    int peak_count = 0;
    int valid_peak_count = 0;
    int step_count = 0;
    int last_valid_peak_index = -1000;

    for (int i = 1; i < N - 1; i++) {
        if (filtered_mag[i] > threshold &&
            filtered_mag[i] > filtered_mag[i - 1] &&
            filtered_mag[i] > filtered_mag[i + 1]) {
#if(0)
            ORIENTATION6D ori = get_6D_orientation(Ax[i], Ay[i], Az[i]);

            printf("Ori=%d\n", ori);

            if (ori == ORIENT_Z_NEG) {
                hold_count++;
            } 
    
            if (hold_count >= ORIENTATION_HOLD_COUNT) {
                printf("Display UP orientation confirmed\n");
                is_wrist_raise = true;
            }

            if (is_wrist_raise) {
                printf("Wrist raise detected, skipping peak\n");
                continue;
            }  
#endif          
            
            if ((i - last_peak_index) >= min_peak_distance && filtered_mag[i] > 0.2f) {
                peak_count++;
                last_peak_index = i;
            }            
            
        }
    }   

    steps = peak_count;
    return steps;
    
}
