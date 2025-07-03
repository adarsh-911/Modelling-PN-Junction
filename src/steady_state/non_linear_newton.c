#include "../../inc/steady_state/non_linear_newton.h"
#include <stdio.h>
#include <string.h>
#include <math.h>

double jacobian_holes [MESH_SIZE][MESH_SIZE];
double jacobian_elec [MESH_SIZE][MESH_SIZE];

double sharfetter_func_holes [MESH_SIZE];
double sharfetter_func_elec [MESH_SIZE];

double hole_conc_inst [MESH_SIZE];
double elec_conc_inst [MESH_SIZE];

double B_func (double vol) {
  double z = vol / (BOLTZ_CONST * TEMP / Q_CHARGE);

  // Taylor series approximation
  if (fabs(z) < 1e-6) {
    double sum = 1 + z / 2;
    return 1 / sum;
  }
  else if (z == 0) return 1;
  else return z / (exp(z) - 1);
}

void init_conc_arr() {
  for (int i = 0 ; i < MESH_SIZE ; i++) {
    hole_conc_inst[i] = 0.0;
    elec_conc_inst[i] = 0.0;
  }

  return ;
}

void init_func() {
  for (int i = 0 ; i < MESH_SIZE ; i++) {
    sharfetter_func_holes[i] = points_eq[i].hole_conc;
    sharfetter_func_elec[i] = points_eq[i].elec_conc;
  }

  return ;
}

void generate_jacobian_holes() {
  const double Kp = - (BOLTZ_CONST * TEMP * MOB_HOLE * LYF_TM_HOLE / (Q_CHARGE * pow(DEL, 2)));

  for (int i = 0 ; i < MESH_SIZE ; i++) {
    for (int j = 0 ; j < MESH_SIZE ; j++) {
      if ((i - j) == 1) {
        if (i == 0 || i == MESH_SIZE - 1) jacobian_holes[i][j] = 0.0;
        else jacobian_holes[i][j] = Kp * B_func(points[i].voltage - points[i - 1].voltage);
      }
      else if ((j - i) == 1) {
        if (i == 0 || i == MESH_SIZE - 1) jacobian_holes[i][j] = 0.0;
        else jacobian_holes[i][j] = Kp * B_func(points[i].voltage - points[i + 1].voltage);
      }
      else if (i == j && (i != 0 && i != MESH_SIZE - 1))
        jacobian_holes[i][j] = 1.0 - (Kp * (B_func(points[i + 1].voltage - points[i].voltage) + B_func(points[i - 1].voltage - points[i].voltage)));
      else jacobian_holes[i][j] = 0.0;
    }
  }
  jacobian_holes[0][0] = 1.0;
  jacobian_holes[MESH_SIZE - 1][MESH_SIZE - 1] = 1.0;

  return ;
}

void generate_jacobian_elec() {
  const double Kn = (BOLTZ_CONST * TEMP * MOB_ELEC * LYF_TM_ELEC / (Q_CHARGE * pow(DEL, 2)));

  for (int i = 0 ; i < MESH_SIZE ; i++) {
    for (int j = 0 ; j < MESH_SIZE ; j++) {
      if ((i - j) == 1) {
        if (i == 0 || i == MESH_SIZE - 1) jacobian_elec[i][j] = 0.0;
        else jacobian_elec[i][j] = Kn * B_func(points[i - 1].voltage - points[i].voltage);
      }
      else if ((j - i) == 1) {
        if (i == 0 || i == MESH_SIZE - 1) jacobian_elec[i][j] = 0.0;
        else jacobian_elec[i][j] = Kn * B_func(points[i + 1].voltage - points[i].voltage);
      }
      else if (i == j && (i != 0 && i != MESH_SIZE - 1))
        jacobian_elec[i][j] = - 1.0 - (Kn * (B_func(points[i].voltage - points[i + 1].voltage) + B_func(points[i].voltage - points[i - 1].voltage)));
      else jacobian_elec[i][j] = 0.0;
    }
  }
  jacobian_elec[0][0] = 1.0;
  jacobian_elec[MESH_SIZE - 1][MESH_SIZE - 1] = 1.0;

  return ;
}

void reduce_rows_matrix (char* name) {
  double (*jacobian)[MESH_SIZE];
  double (*func);

  if (!strcmp(name, "holes")) {
    jacobian = jacobian_holes;
    func = sharfetter_func_holes;
  }
  else if (!strcmp(name, "electrons")) {
    jacobian = jacobian_elec;
    func = sharfetter_func_elec;
  }
  else return ;

  // Reduce rows
  for (int i = 0 ; i < MESH_SIZE ; i++) {
    for (int j = 0 ; j < MESH_SIZE ; j++) {
      if ((i > j) && fabs(jacobian[i][j]) > 1e-6) {
        double factor = jacobian[j][j] / jacobian[i][j];
        for (int k = 0 ; k < MESH_SIZE ; k++) jacobian[i][k] = jacobian[i][k] * factor - jacobian[j][k];

        func[i] = func[i] * factor - func[j];
      }
    }
  }

  return ;
}

void find_conc_update (char* name) {
  double (*jacobian)[MESH_SIZE], (*func), (*update_arr);

  if (!strcmp(name, "holes")) {
    jacobian = jacobian_holes;
    func = sharfetter_func_holes;
    update_arr = hole_conc_inst;
  }
  else if (!strcmp(name, "electrons")) {
    jacobian = jacobian_elec;
    func = sharfetter_func_elec;
    update_arr = elec_conc_inst;
  }
  else return ;

  for (int i = MESH_SIZE - 1 ; i >= 0 ; i--) {
    double sum = func[i];
    for (int j = i + 1 ; j < MESH_SIZE ; j++) sum -= jacobian[i][j] * update_arr[j];

    update_arr[i] = sum / jacobian[i][i];
  }

  return ;
}

void scharfetter_gummel_main() {

  // Initialize the functions
  init_func();

  // Generate jacobian matrices
  generate_jacobian_holes();
  generate_jacobian_elec();

  // Reduce rows
  reduce_rows_matrix("holes");
  reduce_rows_matrix("electrons");

  // Update the concentrations
  find_conc_update("holes");
  find_conc_update("electrons");

  // Update the changes
  for (int i = 0 ; i < MESH_SIZE ; i++) {
    points[i].hole_conc = hole_conc_inst[i];
    points[i].elec_conc = elec_conc_inst[i];
  }

  return ;
}
