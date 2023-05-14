#ifndef CORE_H
#define CORE_H

#include "parameters.h"
#include <immintrin.h>
#include <stdint.h> // uint32_t
#include <stdlib.h> // rand(), RAND_MAX

typedef struct {
    float x[N];
    float y[N];
    float z[N];
} Vector_SOA;

void init_pos(Vector_SOA* v_positions, const float rho);
void init_vel(Vector_SOA* v_velocities, float* temp, float* ekin);
void forces(Vector_SOA* v_positions, Vector_SOA* v_forces, float* epot, float* pres,
            const float* temp, const float rho, const float V, const float L);
void velocity_verlet(Vector_SOA* v_positions, Vector_SOA* v_velocities, Vector_SOA* v_forces, float* epot,
                     float* ekin, float* pres, float* temp, const float rho,
                     const float V, const float L);

#endif
