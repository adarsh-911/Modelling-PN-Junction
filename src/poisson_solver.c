#include "poisson_solver.h"
#include <string.h>
#include <math.h>

double jacobian [MESH_SIZE][MESH_SIZE];
double poisson_func [MESH_SIZE];
double potential_change [MESH_SIZE];

void generate_jacobian() {
  for (int i = 0 ; i < MESH_SIZE ; i++) {
    for (int j = 0 ; j < MESH_SIZE ; j++) {
      if (((i - j) == 1 && i < MESH_SIZE - 1) || ((j - i) == 1 && i > 0)) jacobian[i][j] = 1/pow(DEL, 2);
      else if (i == j) 
        jacobian[i][j] = - (2 / pow(DEL, 2)) - (1 / (EPSILON_SI * BOLTZ_CONST * TEMP)) * (pow(Q_CHARGE, 2) * (points[i].hole_conc + points[i].elec_conc));
      else jacobian[i][j] = 0;
    }
  }
  jacobian[0][0] = 1;
  jacobian[MESH_SIZE - 1][MESH_SIZE - 1] = 1;

  return ;
}

double poisson_func_val (double prev_v, double curr_v, double next_v, int idx, char* side) {
  double rho = 0.0;
  if (!strcmp(side, "n-side")) rho = Q_CHARGE * (DONOR_CONC + points[idx].hole_conc - points[idx].elec_conc);
  else if (!strcmp(side, "p-side")) rho = Q_CHARGE * (-ACC_CONC + points[idx].hole_conc - points[idx].elec_conc);

  return (prev_v - 2*curr_v + next_v) / pow(DEL, 2) + (rho / EPSILON_SI);
}

void update_func_array() {
  for (int i = 1 ; i < MESH_SIZE - 1 ; i++) {
    if (i <= floor(MESH_SIZE/2)) 
      poisson_func[i] = -1 * poisson_func_val(points[i-1].voltage, points[i].voltage, points[i+1].voltage, i, "p-side");
    else
      poisson_func[i] = -1 * poisson_func_val(points[i-1].voltage, points[i].voltage, points[i+1].voltage, i, "n-side");
  }
  poisson_func[0] = 0.0;
  poisson_func[MESH_SIZE - 1] = 0.0;

  return ;
}

void row_reduce_matrix() {
  for (int i = 0 ; i < MESH_SIZE ; i++) {
    for (int j = 0 ; j < MESH_SIZE ; j++) {
      if ((i - j) == 1) {
        double factor = jacobian[i][j] / jacobian[i - 1][j];
        for (int k = 0 ; k < MESH_SIZE ; k++) jacobian[i][k] -= jacobian[i - 1][k] * factor;
        jacobian[i][j] = 0.0;
        if (i > 0) poisson_func[i] -= poisson_func[i - 1] * factor;
      }
    }
    poisson_func[MESH_SIZE - 1] = 0.0;
  }

  return ;
}

void find_potential_update() {
  potential_change[0] = 0.0;
  potential_change[MESH_SIZE - 1] = poisson_func[MESH_SIZE - 1] / jacobian[MESH_SIZE - 1][MESH_SIZE - 1];

  for (int i = MESH_SIZE - 1 ; i > 0 ; i--) 
    potential_change[i - 1] = (poisson_func[i - 1] - (jacobian[i - 1][i] * potential_change[i]))/ jacobian[i - 1][i - 1];

  return ;
}

void poisson_update_main() {

  // Generate jacobian matrix
  generate_jacobian();

  // Update the poisson function array
  update_func_array();

  // Simplify the equations to two variables by reducing rows
  row_reduce_matrix();

  // Find the change in potential
  find_potential_update();

  // Update the change
  for (int i = 0 ; i < MESH_SIZE ; i++) points[i].voltage += potential_change[i];
}