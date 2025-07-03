#include "inc/init.h"
#include "inc/equillibrium.h"
#include "inc/steady_state/steady_state.h"
#include <stdio.h>

int main() {

  // Initialize the space
  init_points("equillibrium");

  // Run the device to equillibrium
  run_to_equillibrium();

  // Print converged values
  print_elements("Equillibrium condition");

  // Copy to new array
  write_eq_points();

  // Inititalize again with applied voltage
  init_points("steady_state");

  // Run the device to steady state
  run_to_steady_state();

  // Print converged values
  print_elements("Steady state condition");

  // Export to json
  export_to_json();

  return 0;
}