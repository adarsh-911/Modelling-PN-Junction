#pragma once

// ALL THE UNITS ARE IN SI

// Constants
#define Q_CHARGE    (1.6e-19)
#define BOLTZ_CONST (1.38e-23)
#define EPSILON_SI  (11.7 * 8.85e-12)

// For silicon
#define EQ_CONC     (1e16)
#define EFF_DEN_NC  (2.8e25)
#define EFF_DEN_NV  (1.04e25)
#define BAND_GAP    (1.12 * Q_CHARGE)

#define MOB_ELEC    0.135
#define MOB_HOLE    0.048

// Point structure
typedef struct Point {
  double voltage;
  double hole_conc;
  double elec_conc;
} Point;