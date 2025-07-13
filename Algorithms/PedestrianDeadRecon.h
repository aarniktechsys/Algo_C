#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <stdint.h>


#define SAMPLE_FREQ    52  //100.0f   // Hz
#define BETA_DEF       0.1f     // 2 * proportional gain
#define DEG2RAD(x)     ((x) * 0.01745329251f)

extern void process_sensor_data(float ax, float ay, float az,
                                float gx, float gy, float gz,
                                float mx, float my, float mz);