#include "../inc/equillibrium.h"
#include <math.h>

#define MAX_IT 10000

void run_to_equillibrium() {
  int id = 0;
  double prev = 0;
  while (id != MAX_IT) {

    // Run the poisson solver
    poisson_update_main();
    
    // Boltzmann Approximation (at equillibrium)
    for (int i = 0 ; i < MESH_SIZE ; i++) {
      points[i].elec_conc = EFF_DEN_NC * exp((Q_CHARGE * points[i].voltage) /(BOLTZ_CONST * TEMP));
      points[i].hole_conc = EFF_DEN_NV * exp(-(Q_CHARGE * points[i].voltage + BAND_GAP) /(BOLTZ_CONST * TEMP));
    }

    // Convergence check
    if (fabs((points[MESH_SIZE/2].voltage - prev) / prev) < TOL) break;

    prev = points[MESH_SIZE/2].voltage;
    id++;

  }
  return ;
}