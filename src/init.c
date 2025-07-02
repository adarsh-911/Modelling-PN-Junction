#include "../inc/init.h"
#include <stdio.h>
#include <math.h>

Point points [MESH_SIZE];

void init_points() {
  for (int i = 1 ; i < MESH_SIZE - 1 ; i++) {
    points[i].voltage = 0.0;
    points[i].hole_conc = 0.0;
    points[i].elec_conc = 0.0;
  }
  double pot_start = VAPP - (BAND_GAP/Q_CHARGE) - (BOLTZ_CONST * TEMP / Q_CHARGE) * log(ACC_CONC/EFF_DEN_NV);
  double pot_end = (BOLTZ_CONST * TEMP / Q_CHARGE) * log(DONOR_CONC/EFF_DEN_NC);

  points[0] = (Point){pot_start, ACC_CONC, (EQ_CONC * EQ_CONC)/ACC_CONC};
  points[MESH_SIZE - 1] = (Point){pot_end, (EQ_CONC * EQ_CONC)/DONOR_CONC, DONOR_CONC}; 

  return ;
}

void print_elements() {
  for (int i = 0 ; i < MESH_SIZE ; i++)
    printf("V = %e\tn = %e\tp = %e\n", points[i].voltage, points[i].elec_conc, points[i].hole_conc);

  return ;
}

int export_to_json() {
  FILE *fp = fopen("out.json", "w");

  if (!fp) {
    perror("fopen");
    return 1;
  }

  fprintf(fp, "[\n");

  for (int i = 0 ; i < MESH_SIZE ; i++)
    fprintf(fp,"  {\"V\": %.8e, \"n\": %.8e, \"p\": %.8e}%s\n", points[i].voltage, points[i].elec_conc, points[i].hole_conc, (i == MESH_SIZE - 1) ? "" : ",");

  fprintf(fp, "]\n");

  fclose(fp);
}