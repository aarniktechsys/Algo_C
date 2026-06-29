#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <stdint.h>


#define SAMPLE_FREQ    52  //100.0f   // Hz
#define BETA_DEF       0.1f     // 2 * proportional gain
#define DEG2RAD(x)     ((x) * 0.01745329251f)
#define RAD2DEG(x)     ((x) * 57.32484076433f)

#define EARTH_RADIUS_M 6378137.0f

extern volatile float q0, q1, q2, q3;
extern void process_sensor_data(float ax, float ay, float az,
                                float gx, float gy, float gz,
                                float mx, float my, float mz);
extern void MadgwickAHRSupdate(float gx, float gy, float gz, float ax, float ay, float az,
                               float mx, float my, float mz);
extern void update_position(float step_length, int peak_count);