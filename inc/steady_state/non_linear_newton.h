#pragma once

#include "../equillibrium.h"

void init_conc_arr();
void scharfetter_gummel_main();

extern double hole_conc_inst [MESH_SIZE];
extern double elec_conc_inst [MESH_SIZE];