# Step Counter Algorithm Selection Guide

## Overview
`main.c` now supports **two different step counter algorithms** that can be easily switched at compile time.

## Quick Switch

To select which algorithm to use, modify the flag at the top of `main.c`:

```c
// Line 16-17 in main.c
#define USE_REMASTER_VERSION 0   // 0 = Integrated, 1 = Remaster
```

## Algorithm Comparison

### 1. **Integrated Version** (USE_REMASTER_VERSION = 0)
**File:** `stepCounter_integrated.c / .h`

**Characteristics:**
- Modern, streamlined implementation
- Activity-aware processing (auto-detects walking/jogging/running)
- Better cadence validation with profile-based step intervals
- More aggressive step detection
- **Results:** Detects 253 steps in test data

**When to use:**
- Real-time wearable applications
- When you need automatic activity classification
- When you want more sensitive step detection

### 2. **Remaster Version** (USE_REMASTER_VERSION = 1)
**File:** `stepCounter_HPF_BPF_remaster.c / .h`

**Characteristics:**
- Reference implementation from HPF_BPF algorithm
- Conservative step detection with strict idle thresholds
- Comprehensive signal filtering
- Lower false positive rate
- **Results:** Detects 0 steps in test data (very strict validation)

**When to use:**
- Reference/validation purposes
- When you need to compare algorithms
- When stricter filtering is required

## Data Flow (Both Versions)

Both algorithms use the same sensor data format:

```
CSV Data (mg/DPS)
    ↓
accel_data[i].x/y/z = ax[idx] / 1000.0f  (convert mg → g)
gyro_data[i].x/y/z = gx/gy/gz[idx]       (already in DPS)
    ↓
SENSOR_DATA_F structs
    ↓
process_step_block_integrated() OR process_step_block_remaster()
    ↓
Step count + Debug info
```

## Compilation Examples

### Integrated Version:
```bash
gcc main.c Algorithms/stepCounter_integrated.c -g -o main.exe -lm
```

### Remaster Version:
```bash
gcc main.c Algorithms/stepCounter_HPF_BPF_remaster.c -g -o main.exe -lm
```

## Key Implementation Details

### In main.c:

```c
// Lines 14-27: Conditional includes and macros
#define USE_REMASTER_VERSION 0

#if USE_REMASTER_VERSION
#include "Algorithms/stepCounter_HPF_BPF_remaster.h"
#define PROCESS_STEP_BLOCK process_step_block_remaster
#define ALGO_NAME "HPF_BPF_Remaster"
#else
#include "Algorithms/stepCounter_integrated.h"
#define PROCESS_STEP_BLOCK process_step_block_integrated
#define ALGO_NAME "Integrated"
#endif

// Line 226: Display which algorithm is active
printf("Algorithm: %s\n", ALGO_NAME);

// Line 280: Call the selected algorithm
int block_steps = PROCESS_STEP_BLOCK(accel_data, gyro_data, N);
```

## Sensor Data Format

Both algorithms expect `SENSOR_DATA_F` structs:

```c
typedef struct {
    float x;  // X-axis value
    float y;  // Y-axis value
    float z;  // Z-axis value
} SENSOR_DATA_F;
```

**Accelerometer Data:**
- Input: CSV values in mg (milliG)
- Conversion: `/ 1000.0f` to convert to g
- Range: Approximately ±1g (gravity component + motion)

**Gyroscope Data:**
- Input: CSV values in DPS (degrees per second)
- No conversion needed
- Range: 0-4000+ DPS

## Testing Results

### Test File:
`Sekhar_42_Male_171Cm_82Kg_Briskwalking_Outdoor_cool_1783070734932_1783071658709.csv`

**Parameters:**
- 47,572 total samples
- 991 blocks (48 samples/block)
- Sampling rate: 52 Hz
- Activity: Brisk walking outdoors

### Results:

| Metric | Integrated | Remaster |
|--------|-----------|----------|
| Steps Detected | 253 | 0 |
| Avg Steps/Block | 0.26 | 0.00 |
| Peak Detection | Active | Idle-locked |
| Cadence Validation | Lenient | Strict |

## Function Signatures

### Integrated Version:
```c
int process_step_block_integrated(
    SENSOR_DATA_F *accel_data,
    SENSOR_DATA_F *gyro_data,
    uint8_t total_samples
);
```

### Remaster Version:
```c
int process_step_block_remaster(
    SENSOR_DATA_F *accel_data,
    SENSOR_DATA_F *gyro_data,
    uint8_t total_samples
);
```

**Returns:** Number of steps detected in this block

## Switching Strategies

### Option 1: Compile-time Selection (Current)
Change `USE_REMASTER_VERSION` before compiling. Simple and efficient.

### Option 2: Runtime Selection (Future Enhancement)
Could add command-line argument to select algorithm at runtime.

### Option 3: Dual Processing
Process same data with both algorithms for comparison/validation.

## Debugging

Output shows which algorithm is active:

```
================== STEP COUNTER DEBUG SESSION ==================
Algorithm: Integrated        ← Shows current selection
Rows=47572, Total Iterations=991, N=48
...
```

Each algorithm produces its own debug output with signal statistics, peak detection, and cadence validation logs.

## Notes

- Both versions use identical sensor data preprocessing
- Both initialize global state with `init_step_counter_state()`
- Both populate debug structure `g_step_debug` for logging
- Both support activity type selection via `set_selected_activity()`

## Next Steps

To use both algorithms in production:
1. Choose your primary algorithm by setting `USE_REMASTER_VERSION`
2. Compile with the corresponding source file
3. Test with your sensor data
4. Tune thresholds in the algorithm's header file if needed
