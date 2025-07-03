#include "../inc/equillibrium.h"
#include <math.h>

#define MAX_IT 10000

void run_to_equillibrium() {
  int id = 0;
  while (id != MAX_IT) {

    // Run the poisson solver
    poisson_update_main();
    
    // Boltzmann Approximation (at equillibrium)
    for (int i = 1 ; i < MESH_SIZE - 1 ; i++) {
      points[i].elec_conc = EFF_DEN_NC * exp((Q_CHARGE * points[i].voltage) /(BOLTZ_CONST * TEMP));
      points[i].hole_conc = EFF_DEN_NV * exp(-(Q_CHARGE * points[i].voltage + BAND_GAP) /(BOLTZ_CONST * TEMP));
    }

    id++;
  }
  return ;
}