#ifndef CORE_H
#define CORE_H

#include "parameters.h"

typedef struct {
    double x[N];
    double y[N];
    double z[N];
} Vector_SOA;

void init_pos(Vector_SOA* v_positions, const double rho);
void init_vel(Vector_SOA* v_velocities, double* temp, double* ekin);
void forces(Vector_SOA* v_positions, Vector_SOA* v_forces, double* epot, double* pres,
            const double* temp, const double rho, const double V, const double L);
void velocity_verlet(Vector_SOA* v_positions, Vector_SOA* v_velocities, Vector_SOA* v_forces, double* epot,
                     double* ekin, double* pres, double* temp, const double rho,
                     const double V, const double L);

#endif
