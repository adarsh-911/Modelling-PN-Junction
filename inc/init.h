#pragma once
#include "def.h"

// Simulation and device Parameters
#define VAPP       0.3
#define LENGTH     (11e-6)
#define MESH_SIZE  300
#define TEMP       300
#define ACC_CONC   (1e21)
#define DONOR_CONC (1e20)

#define DEL        (LENGTH/(MESH_SIZE - 1))

void init_points(char* state);
void print_elements(char* name);
int export_to_json();
void write_eq_points();

extern Point points[MESH_SIZE];
extern Point points_eq [MESH_SIZE];