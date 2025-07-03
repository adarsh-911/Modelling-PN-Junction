#include "inc/init.h"
#include "inc/equillibrium.h"
#include "inc/steady_state/steady_state.h"
#include <stdio.h>
#include <string.h>

int main() {

  // Initialize the space
  init_points("equillibrium");

  // Run the device to equillibrium
  run_to_equillibrium();

  // Copy to new array
  write_eq_points();

  // Inititalize again with applied voltage
  init_points("steady_state");

  // Run the device to steady state
  run_to_steady_state();

  // Export to json
  export_to_json();
  
  // Optional print
  char op[2];
  printf("---------------------------------------------------------------------------------\n");
  printf("Print converged values? (y/n): ");
  if (scanf("%s", op) != 1) return 1;

  if (!strcmp(op, "y")) {
    print_elements("Equillibrium condition");
    print_elements("Steady state condition");
  }

  return 0;
}