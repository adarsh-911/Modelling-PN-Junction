#include "../../inc/steady_state/steady_state.h"
#include <math.h>

#define MAX_IT 2000

void run_to_steady_state() {
  int id = 0;
  double prev = 0;

  // Initialize the concentration array to zero
  init_conc_arr();
  
  while (id != MAX_IT) {

    // Poisson update
    poisson_update_main();

    // Update concentrations
    scharfetter_gummel_main();

    // Convergence check
    if (fabs((points[MESH_SIZE/2].voltage - prev) / prev) < TOL) break;

    prev = points[MESH_SIZE/2].voltage;
    id++;
  }
}