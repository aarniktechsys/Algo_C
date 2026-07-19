#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <stdint.h>
#include <math.h>
#include "main.h"
//#include "Algorithms/stepCounter.h"
//#include "Algorithms/PedestrianDeadRecon.h"
//#include "Algorithms/stepCounter_openSource.h"
//#include "Algorithms/stepCounter_HPF_BPF.h"

//===============================================================================
// @brief Select which step counter algorithm to use
// Set to 0 for stepCounter_integrated
// Set to 1 for stepCounter_HPF_BPF_remaster
//===============================================================================
#define USE_REMASTER_VERSION 1

#if USE_REMASTER_VERSION
#include "Algorithms/stepCounter_HPF_BPF_remaster.h"
#define PROCESS_STEP_BLOCK process_step_block_remaster
#define ALGO_NAME "HPF_BPF_Remaster"
#else
#include "Algorithms/stepCounter_integrated.h"
#define PROCESS_STEP_BLOCK process_step_block_integrated
#define ALGO_NAME "Integrated"
#endif

int steps = 0;
FILE *pdr_log;
FILE *step_detection_log;
//===============================================================================
// @brief  Count number of rows in a data logger file
//
// @param[in]    : Data logger file
// @param[out]   :
// @param[inout] :
// @retval       :
//
//
//===============================================================================
int count_rows(const char *filename){
    char line[MAX_LENGTH];
    int count = 0;

    FILE *file = fopen(filename, "r");

    if (!file){
        printf("ERROR: Failed to open file: %s\n", filename);
        return 0;
    }
    printf("File opened successfully\n");

    //Skip header
    if (fgets(line, sizeof(line), file) == NULL) {
        printf("WARNING: File is empty or could not read header\n");
        fclose(file);
        return 0;
    }

    //Count the rows
    while (fgets(line,sizeof(line),file)){
        count++;
    }

    printf("Number of rows in file = %d\n", count);

    fclose(file);
    return count;
}

//===============================================================================
// @brief  Approximate square root function for Q8.24 format, scaled to Q4.12
//
// @param[in]    : Accelerometer/Gyroscope X,Y,Z
// @param[out]   :
// @param[inout] :
// @retval       :
//
//
//===============================================================================
int sqrt_approx(int value) {
	if (value <= 0) return 0; // Return 0 for non-positive values

    int res = value;
    int bit;

    for (bit = 1 << 30; bit > value;bit >>= 2);        

    while (bit != 0) {
        if (value >= res + bit) {
        	value -= res + bit;
            res = (res >> 1) + bit;
        } else {
            res >>= 1;
        }
        bit >>= 2;
    }

    return res;
}

//===============================================================================
// @brief  main
//
// @param[in]    : Data logger file
// @param[out]   :
// @param[inout] :
// @retval       :
//
//
//===============================================================================
int main(){
    printf("hello world\n");
    const char *filename = "D:\\AarnikTechSys\\Firmware\\Repository\\Design\\Algorithms\\People_Data\\19072026\\Sekhar_42_Male_171Cm_82Kg_Slowwalking_Outdoor_cool_1784462389803_1784462979979_1010Steps_15_18_min_km.csv";
     
    int rows = count_rows(filename);

    if (rows <=0){
        printf("No data found to read in file\n");
        return 1;
    }

    long *unique_id = (long *) malloc(rows * sizeof(long));  
    int  *record_count  = (int *)malloc(rows * sizeof(int));
    int  *ax            = (int *)malloc(rows * sizeof(int));
    int  *ay            = (int *)malloc(rows * sizeof(int));
    int  *az            = (int *)malloc(rows * sizeof(int));
    int  *gx            = (int *)malloc(rows * sizeof(int));
    int  *gy            = (int *)malloc(rows * sizeof(int));
    int  *gz            = (int *)malloc(rows * sizeof(int));
    int  *mx            = (int *)malloc(rows * sizeof(int));
    int  *my            = (int *)malloc(rows * sizeof(int));
    int  *mz            = (int *)malloc(rows * sizeof(int));
    int  *direction     = (int *)malloc(rows * sizeof(int));
    int  *alt_x         = (int *)malloc(rows * sizeof(int));
    int  *alt_y         = (int *)malloc(rows * sizeof(int));
    int  *hx            = (int *)malloc(rows * sizeof(int));
    int  *spo2          = (int *)malloc(rows * sizeof(int));

    int16_t Accel_x[48] = {0};
    int16_t Accel_y[48] = {0};
    int16_t Accel_z[48] = {0};    

    int16_t Gyro_x[48] = {0};
    int16_t Gyro_y[48] = {0};
    int16_t Gyro_z[48] = {0};   

    int16_t Mag_x[48] = {0};
    int16_t Mag_y[48] = {0};
    int16_t Mag_z[48] = {0};

    int16_t Ax[48] = {0};
    int16_t Ay[48] = {0};
    int16_t Az[48] = {0};
    int16_t Gx[48] = {0};
    int16_t Gy[48] = {0};
    int16_t Gz[48] = {0};
    int16_t Mx[48] = {0};
    int16_t My[48] = {0};
    int16_t Mz[48] = {0};

    char line[256];
    int row = 0;
    int total_step_count = 0;
    int total_step_count_OF = 0;    
	static int32_t idle_accel_sum = 0,idle_gyro_sum = 0;
	static int idle_counter = 0;
    int16_t recent_step_counter = 0;
    int block_offset  = 0;    
    const float GYRO_SENSITIVITY_MDPS_PER_COUNT = 500.0f * 1000.0f / 32768.0f;  // ~15.26 mdps/count

    FILE *file = fopen(filename, "r");

    if (file == NULL) {
        printf("ERROR: Failed to open file: %s\n", filename);
        printf("Possible causes: File not found, path too long, or invalid characters\n");
        return 1;
    }

    //skip header
    fgets(line, sizeof(line), file);

    while (fgets(line,sizeof(line),file) && row < rows){
        sscanf(line, "%ld,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d",
               &unique_id[row], &record_count[row],
               &ax[row], &ay[row], &az[row],
               &gx[row], &gy[row], &gz[row],
               &mx[row], &my[row], &mz[row],
               &direction[row], &alt_x[row], &alt_y[row],
               &hx[row], &spo2[row]);
        row++;
    }

    fclose(file);

    pdr_log = fopen("pdr_log.csv", "w");
    if (pdr_log != NULL) {
        // Write the CSV header
        fprintf(pdr_log, "steps,x,y,heading,latitude,longitude\n");
    } else {
        // Handle file open error
        printf("Failed to open pdr_log.csv\n");
    }

    step_detection_log = fopen("step_detection_log.csv", "w");
    if (step_detection_log != NULL) {
        // Write the CSV header
        fprintf(step_detection_log, "avg_interval,min_peak_distance,std,energy,mean_z\n");
    } else {
        // Handle file open error
        printf("Failed to open step_detection_log.csv\n");
    }


    int total_iterations = rows / N;
    //total_iterations = 20;
    steps = 0;

    printf("================== STEP COUNTER DEBUG SESSION ==================\n");
    printf("Algorithm: %s\n", ALGO_NAME);
    printf("Rows=%d, Total Iterations=%d, N=%d\n", rows, total_iterations, N);
    printf("Data file: %s\n", filename);
    printf("================================================================\n\n");

    // Initialize step counter state before processing
    init_step_counter_state();

    // === VALIDATION: Check first block data ===
    printf(">>> VALIDATION: Checking first block (iteration 0)...\n");
    for (int i = 0; i < N; i++){
        int idx = i + 0 * N;
        if (i < 5) {  // Print first 5 samples
            printf("  Sample[%d]: ax=%d ay=%d az=%d | gx=%d gy=%d gz=%d\n",
                   idx, ax[idx], ay[idx], az[idx], gx[idx], gy[idx], gz[idx]);
        }
    }
    printf(">>> Accel data range: min/max check\n");
    int ax_min = ax[0], ax_max = ax[0];
    for (int i = 0; i < N; i++) {
        if (ax[i] < ax_min) ax_min = ax[i];
        if (ax[i] > ax_max) ax_max = ax[i];
    }
    printf("    Accel X range: [%d, %d] (expect ±2000-4000 for normal motion)\n", ax_min, ax_max);
    printf("\n");

    // === MAIN PROCESSING LOOP ===
    printf(">>> PROCESSING %d blocks of %d samples each...\n\n", total_iterations, N);

    // Create float sensor data structures
    SENSOR_DATA_F accel_data[N], gyro_data[N];

    init_step_counter_state();

    for (int iteration_cnt = 0; iteration_cnt < total_iterations; iteration_cnt++){

        // Convert sensor data to float and calculate dynamic acceleration
        float accel_mean_mag = 0.0f;
        float gyro_mean_mag = 0.0f;

        for (int i = 0; i < N; i++){
            int idx = i + iteration_cnt * N;
            
            accel_data[i].x = (float)ax[idx];
            accel_data[i].y = (float)ay[idx];
            accel_data[i].z = (float)az[idx];

            // Gyro data is in DPS (degrees per second)
            gyro_data[i].x = (float)gx[idx];
            gyro_data[i].y = (float)gy[idx];
            gyro_data[i].z = (float)gz[idx];
        }        

        // Process block with selected step counter algorithm
        printf("[Block %2d] Processing...\n", iteration_cnt);
        printf("          Dynamic accel mean=%.2f mg, Gyro mean=%.2f DPS\n",
               accel_mean_mag, gyro_mean_mag);
        int block_steps = PROCESS_STEP_BLOCK(accel_data, gyro_data, N);
        steps += block_steps;

        // Log results
        printf("[Block %2d] Result: +%d steps, Total: %d\n",
               iteration_cnt, block_steps, steps);
        printf("          Debug: mean=%d std=%d energy=%d threshold=%d peaks=%d\n",
               g_step_debug.mean, g_step_debug.std, g_step_debug.energy,
               g_step_debug.threshold, g_step_debug.peak_count);
        printf("\n");

        block_offset += N;
    }

    printf("\n================== FINAL RESULTS ==================\n");
    printf("Total step count: %d steps over %d blocks\n", steps, total_iterations);
    printf("Average steps per block: %.2f\n", (float)steps / total_iterations);
    printf("===================================================\n");
    printf("Press any key to exit...\n");
    getchar();

    return 0;
}



