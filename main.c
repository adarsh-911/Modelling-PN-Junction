#include "inc/init.h"
#include "inc/equillibrium.h"
#include <stdio.h>

int main() {

  // Initialize the space
  init_points();

  // Run the device to equillibrium
  run_to_equillibrium();

  // Export to json
  export_to_json();

  // Print converged values
  printf("---------------------------------------------------------------------------------\n");
  printf("Equillibrium condition\n");
  printf("---------------------------------------------------------------------------------\n");

  print_elements();

  return 0;
}