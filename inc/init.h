#pragma once
#include "def.h"

// Simulation and device Parameters
#define VAPP       0.0
#define LENGTH     (11e-6)
#define MESH_SIZE  311
#define TEMP       300
#define ACC_CONC   (1e21)
#define DONOR_CONC (1e20)

#define DEL        (LENGTH/(MESH_SIZE - 1))

void init_points();
void print_elements();
int export_to_json();
extern Point points[MESH_SIZE];