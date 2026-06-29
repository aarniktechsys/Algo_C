#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <stdint.h>
#include <math.h>
#include "../main.h"
#include "stepCounter.h"
#include "PedestrianDeadRecon.h"

// Madgwick quaternion
volatile float q0 = 1.0f, q1 = 0.0f, q2 = 0.0f, q3 = 0.0f;
float beta = BETA_DEF;

// Step detection state
int pdr_step_count = 0;
float last_mag = 0.0f;
int step_timer = 0;

// Position tracking
float pos_x = 0.0f, pos_y = 0.0f;

void process_sensor_data(float ax, float ay, float az,
    float gx, float gy, float gz,
    float mx, float my, float mz);
void MadgwickAHRSupdate(float gx, float gy, float gz, float ax, float ay, float az,
                        float mx, float my, float mz);
float get_heading_deg();
bool detect_step(float ax, float ay, float az);
void update_position(float step_length, int peak_count);
void pdr_to_latlon(double lat0_deg, double lon0_deg, float heading_deg, float step_len_m,
                   double *new_lat_deg, double *new_lon_deg);
//===============================================================================
// @brief  Madgwick Filter Update
//
// @param[in]    : Sensors - Accel, Gyro, Mag
// @param[out]   :
// @param[inout] :
// @retval       :
//===============================================================================
void MadgwickAHRSupdate(float gx, float gy, float gz, float ax, float ay, float az,
                        float mx, float my, float mz)
{
    float recipNorm;
    float s0, s1, s2, s3;
    float qDot1, qDot2, qDot3, qDot4;
    float hx, hy;
    float _2q0mx, _2q0my, _2q0mz, _2q1mx;
    float _2q0 = 2.0f * q0;
    float _2q1 = 2.0f * q1;
    float _2q2 = 2.0f * q2;
    float _2q3 = 2.0f * q3;
    float _2q0q2 = 2.0f * q0 * q2;
    float _2q2q3 = 2.0f * q2 * q3;
    float q0q0 = q0 * q0;
    float q0q1 = q0 * q1;
    float q0q2 = q0 * q2;
    float q0q3 = q0 * q3;
    float q1q1 = q1 * q1;
    float q1q2 = q1 * q2;
    float q1q3 = q1 * q3;
    float q2q2 = q2 * q2;
    float q2q3 = q2 * q3;
    float q3q3 = q3 * q3;

    // Normalize accelerometer
    recipNorm = 1.0f / sqrtf(ax * ax + ay * ay + az * az);
    ax *= recipNorm;
    ay *= recipNorm;
    az *= recipNorm;

    // Normalize magnetometer
    recipNorm = 1.0f / sqrtf(mx * mx + my * my + mz * mz);
    mx *= recipNorm;
    my *= recipNorm;
    mz *= recipNorm;

    // Reference direction of Earth's magnetic field
    _2q0mx = 2.0f * q0 * mx;
    _2q0my = 2.0f * q0 * my;
    _2q0mz = 2.0f * q0 * mz;
    _2q1mx = 2.0f * q1 * mx;
    hx = mx * q0q0 - _2q0my * q3 + _2q0mz * q2 + mx * q1q1 + _2q1 * my * q2 + _2q1 * mz * q3 - mx * q2q2 - mx * q3q3;
    hy = _2q0mx * q3 + my * q0q0 - _2q0mz * q1 + _2q1mx * q2 - my * q1q1 + my * q2q2 + _2q2 * mz * q3 - my * q3q3;

    // Gradient descent algorithm corrective step (simplified)
    float f1 = 2*(q1*q3 - q0*q2) - ax;
    float f2 = 2*(q0*q1 + q2*q3) - ay;
    float f3 = 2*(0.5f - q1*q1 - q2*q2) - az;
    s0 = f1 * (-2*q2) + f2 * (2*q1);
    s1 = f1 * (2*q3) + f2 * (2*q0);
    s2 = f1 * (2*q0) + f2 * (2*q3);
    s3 = f1 * (-2*q1) + f2 * (-2*q2);
    recipNorm = 1.0f / sqrtf(s0*s0 + s1*s1 + s2*s2 + s3*s3);
    s0 *= recipNorm;
    s1 *= recipNorm;
    s2 *= recipNorm;
    s3 *= recipNorm;

    // Integrate rate of change of quaternion
    qDot1 = 0.5f * (-q1 * gx - q2 * gy - q3 * gz) - beta * s0;
    qDot2 = 0.5f * (q0 * gx + q2 * gz - q3 * gy) - beta * s1;
    qDot3 = 0.5f * (q0 * gy - q1 * gz + q3 * gx) - beta * s2;
    qDot4 = 0.5f * (q0 * gz + q1 * gy - q2 * gx) - beta * s3;

    // Integrate to yield quaternion
    q0 += qDot1 / SAMPLE_FREQ;
    q1 += qDot2 / SAMPLE_FREQ;
    q2 += qDot3 / SAMPLE_FREQ;
    q3 += qDot4 / SAMPLE_FREQ;

    // Normalize quaternion
    recipNorm = 1.0f / sqrtf(q0*q0 + q1*q1 + q2*q2 + q3*q3);
    q0 *= recipNorm;
    q1 *= recipNorm;
    q2 *= recipNorm;
    q3 *= recipNorm;

    /* printf("q0 = %f", q0);
    printf("q1 = %f", q1);
    printf("q2 = %f", q2);
    printf("q3 = %f", q3); */
}
//===============================================================================
// @brief  Interpolate with known longitude and latitude
//
// @param[in]    : 
// @param[out]   :
// @param[inout] :
// @retval       :
//===============================================================================
void pdr_to_latlon(double lat0_deg, double lon0_deg, float heading_deg, float step_len_m,
                   double* new_lat_deg, double* new_lon_deg) {
    // Convert inputs to radians
    double lat0_rad = DEG2RAD(lat0_deg);
    double lon0_rad = DEG2RAD(lon0_deg);
    double heading_rad = DEG2RAD(heading_deg);

    double delta = step_len_m / EARTH_RADIUS_M;

    // Compute new latitude
    double lat1_rad = asin(sin(lat0_rad) * cos(delta) +
                           cos(lat0_rad) * sin(delta) * cos(heading_rad));

    // Compute new longitude
    double lon1_rad = lon0_rad + atan2(sin(heading_rad) * sin(delta) * cos(lat0_rad),
                                       cos(delta) - sin(lat0_rad) * sin(lat1_rad));

    // Convert back to degrees
    *new_lat_deg = RAD2DEG(lat1_rad);
    *new_lon_deg = RAD2DEG(lon1_rad);
}

// Convert total (x, y) meters to lat/lon offset
void displacement_to_latlon(double lat0_deg, double lon0_deg, float dx_m, float dy_m,
                            double* new_lat_deg, double* new_lon_deg) {
    double lat0_rad = DEG2RAD(lat0_deg);
    double lon0_rad = DEG2RAD(lon0_deg);

    // Latitude offset (north-south)
    double dlat_rad = dy_m / EARTH_RADIUS_M;
    double dlon_rad = dx_m / (EARTH_RADIUS_M * cos(lat0_rad));

    *new_lat_deg = lat0_deg + RAD2DEG(dlat_rad);
    *new_lon_deg = lon0_deg + RAD2DEG(dlon_rad);
}
//===============================================================================
// @brief  Compute heading from quaternion
//
// @param[in]    : 
// @param[out]   :
// @param[inout] :
// @retval       :
//===============================================================================
float get_heading_deg() {
    double lat0 = 17.5039247;  // Start Latitude (Hyderabad)
    double lon0 = 78.3584966;  // Start Longitude
    float heading = atan2f(2.0f * (q0 * q3 + q1 * q2),
                           1.0f - 2.0f * (q2 * q2 + q3 * q3));
    heading *= (180.0f / 3.14159265f);
    heading  = (heading) - 176.24;
    if (heading < 0) heading += 360.0f;
    
    printf("new_lat =  %f,new_lon = %f\n", heading);
    return heading;
}
//===============================================================================
// @brief   Step Detection
//
// @param[in]    : 
// @param[out]   :
// @param[inout] :
// @retval       :
//===============================================================================
bool detect_step(float ax, float ay, float az) {
    float mag = sqrtf(ax*ax + ay*ay + az*az);
    float diff = mag - last_mag;
    last_mag = mag;

    // Simple peak detection with debounce
    if (diff > 0.7f && step_timer <= 0) {
        step_timer = (int)(0.4f * SAMPLE_FREQ); // 400ms debounce
        pdr_step_count++;
        printf("detect_step_pdr_step_count =   %d ", pdr_step_count);
        return true;
    }
    if (step_timer > 0) step_timer--;
    return false;
}
//===============================================================================
// @brief   Update Position using PDR (Pedestrian Dead Reckoning)
//
// @param[in]    : 
// @param[out]   :
// @param[inout] :
// @retval       :
//===============================================================================
void update_position(float step_length, int peak_count) {    
    double new_lat, new_lon;
    double lat0 = 17.5019634;  // Start Latitude (Hyderabad)
    double lon0 = 78.3601644;  // Start Longitude

    float heading_deg = get_heading_deg();                // Keep heading in degrees for logging
    float heading_rad = DEG2RAD(heading_deg);             // Convert only for math
    pos_x += step_length * cosf(heading_rad);
    pos_y += step_length * sinf(heading_rad);

    peak_count = steps;

    //pdr_to_latlon(lat0, lon0, heading_deg, step_length, &new_lat, &new_lon);

    displacement_to_latlon(lat0, lon0, pos_x, pos_y, &new_lat, &new_lon);

    fprintf(pdr_log, "%d,%f,%f,%f,%f,%f\n", peak_count, pos_x, pos_y, heading_deg,new_lat,new_lon);

    printf("Step %d, Heading=%f, PosX=%f, PosY=%f, Lattitude=%f, Longitude=%f\n", peak_count, heading_deg, pos_x, pos_y,
                        new_lat,new_lon);
}
//===============================================================================
// @brief  PDR Algorithm
//
// @param[in]    : Data logger file
// @param[out]   :
// @param[inout] :
// @retval       :
//===============================================================================
void process_sensor_data(float ax, float ay, float az,
                         float gx, float gy, float gz,
                         float mx, float my, float mz){

    static uint32_t pres_step_cnt = 0;
    static uint32_t prev_step_cnt = 0;

    MadgwickAHRSupdate(gx, gy, gz, ax, ay, az, mx, my, mz);
    #if(0)
    pres_step_cnt = process_sample(ax, ay, az, gx, gy, gz);    
    if(prev_step_cnt != pres_step_cnt){        
        prev_step_cnt = pres_step_cnt;
        pdr_step_count += pres_step_cnt;
        printf("pdr_step_count =   %d ", pdr_step_count);
        update_position(0.7f); // Assume 0.7m per step
    }
    #endif
#if(1)
    if(detect_step(ax,ay,az)){
        update_position(0.7f,1); // Assume 0.7m per step
    }
#endif
    
}