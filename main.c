#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <stdint.h>
#include "main.h"
#include "Algorithms/stepCounter.h"
#include "Algorithms/PedestrianDeadRecon.h"
#include "Algorithms/stepCounter_openSource.h"
#include "Algorithms/stepCounter_HPF_BPF.h"

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
        printf("error opening a file\n");
    }
    printf("file opened successfully\n");

    //Skip header
    fgets(line, sizeof(line), file);

    //Count the rows
    while (fgets(line,sizeof(line),file)){
        count++;
    };

    printf("number of rows in file = %ld", count);

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
     //const char *filename = "D:\\AarnikTechSys\\Firmware\\Repository\\Design\\Algorithms\\People_Data\\14062025\\Chakri_42_5_7_96_Male_Walking_1749920228_1749920957_1749901177768_1749901422859.csv";
    const char *filename = "C:\\Users\\Admin\\Downloads\\J Chandra Sekhar_41_Male_171(Cm)_85(Kg)_Brisk walking_Outdoor_cool_1752065182312_1752066144350.csv";
     
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
    int steps = 0;
	static int32_t idle_accel_sum = 0,idle_gyro_sum = 0;
	static int idle_counter = 0;
    int16_t recent_step_counter = 0;
    int block_offset  = 0;    

    FILE *file = fopen(filename, "r");

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

    int total_iterations = rows / N;

    printf("Rows=%d\n", rows);   
# if(0)
    for (int iteration_cnt = 0; iteration_cnt < rows;iteration_cnt++){

        int Ax = (int)ax[iteration_cnt];
        int Ay = (int)ay[iteration_cnt];
        int Az = (int)az[iteration_cnt];
        int Gx = (int)gx[iteration_cnt];
        int Gy = (int)gy[iteration_cnt];
        int Gz = (int)gz[iteration_cnt];
        int Mx = (int)mx[iteration_cnt];
        int My = (int)my[iteration_cnt];
        int Mz = (int)mz[iteration_cnt];

    #if (EXECUTE_STEP_COUNTER_ALGO)
        steps += process_sample(Ax, Ay, Az, Gx, Gy, Gz);        
    #endif
    #if (EXECUTE_PDR_ALGO)
        process_sensor_data(Ax, Ay, Ax, Gx, Gy, Gz, Mx, My, Mz);
    #endif
    #if (EXECUTE_OPENSOURCE_STEP_CNT)
        if(process_sample_oxford(Ax, Ay, Az) == 1){
            steps++;
        }
    #endif
    }
#endif       

    for (int iteration_cnt = 0; iteration_cnt < total_iterations;iteration_cnt++){

        for (int i = 0; i < N; i++){
            Ax[i] = (int)ax[i + iteration_cnt * N];
            Ay[i] = (int)ay[i + iteration_cnt * N];
            Az[i] = (int)az[i + iteration_cnt * N];
            Gx[i] = (int)gx[i + iteration_cnt * N];
            Gy[i] = (int)gy[i + iteration_cnt * N];
            Gz[i] = (int)gz[i + iteration_cnt * N];
            Mx[i] = (int)mx[i + iteration_cnt * N];
            My[i] = (int)my[i + iteration_cnt * N];
            Mz[i] = (int)mz[i + iteration_cnt * N];
        }
        //steps += process_step_block(&Ax[0], &Ay[0], &Az[0],&Gx[0],&Gy[0],&Gz[0],&recent_step_counter);
        steps += process_step_block_SF(&Ax[0], &Ay[0], &Az[0],&Gx[0],&Gy[0],&Gz[0],N);
        block_offset += N; 
    }

    printf("Total step count= %d\n",steps);
    getchar();

    return 0;
}
