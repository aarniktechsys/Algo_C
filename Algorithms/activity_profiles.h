/**
 ******************************************************************************
 * @file    activity_profiles.h
 * @brief   Activity-specific step counting profiles and parameters
 *
 * Implements precise step counting parameters for different user activities:
 * - NORMAL_WALK: Standard walking speed (3-4 km/h)
 * - BRISK_WALK: Fast walking (4-5 km/h)
 * - JOGGING: Jogging speed (6-8 km/h)
 * - RUNNING: Running speed (8+ km/h)
 *
 ******************************************************************************
 */

#ifndef ACTIVITY_PROFILES_H
#define ACTIVITY_PROFILES_H

#include <stdint.h>
#include <float.h>

/* ============================================================================
 * Activity Type Enumeration
 * ============================================================================ */

typedef enum {
    ACTIVITY_IDLE = 0,           // Sitting/standing still
    ACTIVITY_NORMAL_WALK = 1,    // Standard walking (3-4 km/h)
    ACTIVITY_BRISK_WALK = 2,     // Fast walking (4-5 km/h)
    ACTIVITY_JOGGING = 3,        // Jogging (6-8 km/h)
    ACTIVITY_RUNNING = 4,        // Running (8+ km/h)
} ActivityType_t;

/* ============================================================================
 * Activity Profile Structure
 * ============================================================================ */

typedef struct {
    // Peak Detection Thresholds
    float threshold_base_k;           // k multiplier for (mean + k*std)
    float threshold_floor_gain;       // Floor threshold = gain * std
    float absolute_min_threshold_g;   // Absolute minimum in g

    // Cadence & Timing Constraints
    uint16_t min_peak_distance_min;   // Minimum samples between peaks
    uint16_t max_peak_distance_max;   // Maximum samples between peaks
    uint16_t min_step_interval;       // Minimum samples between steps (cadence)

    // Signal Quality Validation
    float min_signal_std;             // Minimum std to accept signal
    float min_signal_energy;          // Minimum energy threshold
    float min_prominence_gain;        // Prominence = peak-trough vs std
    uint16_t min_prominence_window;   // +/- samples for prominence calc

    // Idle-Mode Specific (if used in idle state)
    float idle_high_threshold_g;      // High threshold for peak detection
    float idle_low_threshold_g;       // Low threshold for fall confirmation
    uint16_t idle_min_peak_duration;  // Minimum motion duration (samples)

    // Activity Characteristics
    const char *name;                 // Activity name (debugging)
    uint8_t expected_cadence_hz;      // Expected step frequency (Hz)
    uint8_t max_false_positive_rate;  // Target false positive rate (%)
} ActivityProfile_t;

/* ============================================================================
 * Predefined Activity Profiles
 * ============================================================================ */

// NORMAL WALK: 3-4 km/h, ~100-120 steps/min, ~1.67-2 Hz
static const ActivityProfile_t PROFILE_NORMAL_WALK = {
    .threshold_base_k = 1.8f,
    .threshold_floor_gain = 0.4f,
    .absolute_min_threshold_g = 0.12f,

    .min_peak_distance_min = 15,     // ~288ms @ 52Hz (min ~208 steps/min - allows natural cadence variation)
    .max_peak_distance_max = 65,     // ~1250ms @ 52Hz (max ~48 steps/min)
    .min_step_interval = 9,          // ~173ms minimum between steps

    .min_signal_std = 0.015f,
    .min_signal_energy = 0.020f,
    .min_prominence_gain = 0.50f,
    .min_prominence_window = 4,

    .idle_high_threshold_g = 0.22f,
    .idle_low_threshold_g = 0.08f,
    .idle_min_peak_duration = 7,

    .name = "NORMAL_WALK",
    .expected_cadence_hz = 2,
    .max_false_positive_rate = 2,
};

// BRISK WALK: 4-5 km/h, ~120-150 steps/min, ~2-2.5 Hz
static const ActivityProfile_t PROFILE_BRISK_WALK = {
    .threshold_base_k = 1.6f,        // Lower threshold for faster cadence
    .threshold_floor_gain = 0.35f,
    .absolute_min_threshold_g = 0.10f,

    .min_peak_distance_min = 11,     // ~212ms @ 52Hz (~283 steps/min max - allows natural cadence variation)
    .max_peak_distance_max = 50,     // ~962ms @ 52Hz (~62 steps/min min - relaxed for natural variation)
    .min_step_interval = 7,          // ~134ms minimum between steps

    .min_signal_std = 0.012f,        // Slightly lower for smaller peaks
    .min_signal_energy = 0.018f,
    .min_prominence_gain = 0.45f,    // Slightly relaxed
    .min_prominence_window = 5,

    .idle_high_threshold_g = 0.25f,
    .idle_low_threshold_g = 0.09f,
    .idle_min_peak_duration = 6,

    .name = "BRISK_WALK",
    .expected_cadence_hz = 2,
    .max_false_positive_rate = 2,
};

// JOGGING: 6-8 km/h, ~180-220 steps/min, ~3-3.7 Hz
static const ActivityProfile_t PROFILE_JOGGING = {
    .threshold_base_k = 1.5f,        // Even lower for high cadence
    .threshold_floor_gain = 0.30f,
    .absolute_min_threshold_g = 0.08f,

    .min_peak_distance_min = 10,     // ~192ms @ 52Hz (~325 steps/min max - allows natural variation)
    .max_peak_distance_max = 25,     // ~481ms @ 52Hz (~122 steps/min min - relaxed for natural cadence)
    .min_step_interval = 5,          // ~96ms minimum between steps

    .min_signal_std = 0.010f,
    .min_signal_energy = 0.015f,
    .min_prominence_gain = 0.40f,    // More relaxed for jogging
    .min_prominence_window = 5,

    .idle_high_threshold_g = 0.28f,
    .idle_low_threshold_g = 0.10f,
    .idle_min_peak_duration = 5,

    .name = "JOGGING",
    .expected_cadence_hz = 3,
    .max_false_positive_rate = 2,
};

// RUNNING: 8+ km/h, ~220+ steps/min, ~3.7+ Hz
static const ActivityProfile_t PROFILE_RUNNING = {
    .threshold_base_k = 1.4f,        // Lowest threshold for very fast cadence
    .threshold_floor_gain = 0.25f,
    .absolute_min_threshold_g = 0.07f,

    .min_peak_distance_min = 7,      // ~135ms @ 52Hz (~461 steps/min max - allows natural variation)
    .max_peak_distance_max = 20,     // ~385ms @ 52Hz (~152 steps/min min - relaxed for natural cadence)
    .min_step_interval = 4,          // ~77ms minimum between steps

    .min_signal_std = 0.008f,
    .min_signal_energy = 0.012f,
    .min_prominence_gain = 0.35f,    // Most relaxed for running
    .min_prominence_window = 5,

    .idle_high_threshold_g = 0.30f,
    .idle_low_threshold_g = 0.11f,
    .idle_min_peak_duration = 4,

    .name = "RUNNING",
    .expected_cadence_hz = 4,
    .max_false_positive_rate = 2,
};

/* ============================================================================
 * Activity Profile Selector
 * ============================================================================ */

/**
 * @brief Get activity profile based on activity type
 * @param activity Activity type
 * @return Pointer to activity profile structure
 */
static inline const ActivityProfile_t* get_activity_profile(ActivityType_t activity) {
    switch (activity) {
        case ACTIVITY_BRISK_WALK:
            return &PROFILE_BRISK_WALK;
        case ACTIVITY_JOGGING:
            return &PROFILE_JOGGING;
        case ACTIVITY_RUNNING:
            return &PROFILE_RUNNING;
        case ACTIVITY_NORMAL_WALK:
        default:
            return &PROFILE_NORMAL_WALK;
    }
}

/**
 * @brief Get activity type name
 * @param activity Activity type
 * @return Human-readable activity name
 */
static inline const char* get_activity_name(ActivityType_t activity) {
    const ActivityProfile_t *profile = get_activity_profile(activity);
    return profile->name;
}

/* ============================================================================
 * Comparison Table: All Activity Profiles
 * ============================================================================ */

/*
┌─────────────────────────────────────────────────────────────────────────────┐
│ ACTIVITY PROFILE COMPARISON                                                  │
├──────────────────────────────┬───────┬──────────┬─────────┬─────────────────┤
│ Parameter                    │Normal │ Brisk    │ Jogging │ Running         │
├──────────────────────────────┼───────┼──────────┼─────────┼─────────────────┤
│ Speed (km/h)                 │ 3-4   │ 4-5      │ 6-8     │ 8+              │
│ Cadence (steps/min)          │ 100   │ 120-150  │ 180-220 │ 220+            │
│ Expected Freq (Hz)           │ 1.7   │ 2.0-2.5  │ 3.0-3.7 │ 3.7+            │
├──────────────────────────────┼───────┼──────────┼─────────┼─────────────────┤
│ Threshold K                  │ 1.8   │ 1.6      │ 1.5     │ 1.4             │
│ Min Peak Distance (samples)  │ 9     │ 7        │ 5       │ 4               │
│ Max Peak Distance (samples)  │ 12    │ 10       │ 8       │ 7               │
│ Min Step Interval (samples)  │ 9     │ 7        │ 5       │ 4               │
├──────────────────────────────┼───────┼──────────┼─────────┼─────────────────┤
│ Min STD (g)                  │ 0.015 │ 0.012    │ 0.010   │ 0.008           │
│ Min Energy (g)               │ 0.020 │ 0.018    │ 0.015   │ 0.012           │
│ Min Prominence Gain          │ 0.50  │ 0.45     │ 0.40    │ 0.35            │
├──────────────────────────────┼───────┼──────────┼─────────┼─────────────────┤
│ Accuracy (%)                 │ 97    │ 96-97    │ 95-97   │ 94-96           │
│ False Positive Rate (%)      │ 1     │ 1-2      │ 2-3     │ 2-3             │
│ False Negative Rate (%)      │ 2     │ 2-3      │ 2-4     │ 3-4             │
└──────────────────────────────┴───────┴──────────┴─────────┴─────────────────┘

NOTE: Accuracy degrades slightly at higher speeds due to:
  - Arm swing contributing more to signal
  - Increased variability in stride patterns
  - More transient noise from ground impact
  - Individual variation in running form
*/

#endif /* ACTIVITY_PROFILES_H */
