/* Real-Time + Batch Oxford Step Detection Pipeline in C
 * Adapted for CSV Test Harness
 * Author: Aarnik Tech Sys
 * Description: Integrates Oxford StepCountingAlgo into memory-based C test pipeline
 */

 #include <math.h>
 #include <stdint.h>
 #include <string.h>
 #include "stepCounter_openSource.h"
 
 // Oxford StepCounter state
 static float prevAccMag = 0;
 static float prevAccMag2 = 0;
 static float stepThreshold = 1.0f; // Adjustable threshold
 static int lastStepIndex = -999;
 static int sampleIndex = 0;
 static int oxfordStepCount = 0;

 int process_sample_oxford(int Ax, int Ay, int Az);
 
 // Reset Oxford algorithm
 void StepCountingAlgo_Init(void) {
     prevAccMag = prevAccMag2 = 0.0f;
     stepThreshold = 1.0f;
     lastStepIndex = -999;
     sampleIndex = 0;
     oxfordStepCount = 0;
 }
 
 // Process one accelerometer sample
 int StepCountingAlgo_Process(float accX, float accY, float accZ) {
     float mag = sqrtf(accX * accX + accY * accY + accZ * accZ);
 
     // Simple peak detection
     if (prevAccMag2 < prevAccMag && prevAccMag > mag && prevAccMag > stepThreshold) {
         if (sampleIndex - lastStepIndex > 10) { // ~200ms min distance
             oxfordStepCount++;
             lastStepIndex = sampleIndex;
             prevAccMag2 = prevAccMag;
             prevAccMag = mag;
             sampleIndex++;
             return 1; // step detected
         }
     }
 
     prevAccMag2 = prevAccMag;
     prevAccMag = mag;
     sampleIndex++;
     return 0;
 }
 
 // CSV integration test wrapper
 int process_sample_oxford(int Ax, int Ay, int Az) {
     float ax_g = Ax / 1000.0f;
     float ay_g = Ay / 1000.0f;
     float az_g = Az / 1000.0f;
 
     return StepCountingAlgo_Process(ax_g, ay_g, az_g);
 } 

 