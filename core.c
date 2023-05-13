#include "core.h"
#include "parameters.h"

#include <math.h>
#include <stdlib.h> // rand()

#define ECUT (4.0 * (pow(RCUT, -12) - pow(RCUT, -6)))

/*
void init_pos(Vector_SOA* restrict v_positions, const float rho)
{
    // inicialización de las posiciones de los átomos en un cristal FCC

    float a = cbrt(4.0 / rho);
    int nucells = ceil(cbrt((float)N / 4.0));
    int idx = 0;

    for (int i = 0; i < nucells; i++) {
        for (int j = 0; j < nucells; j++) {
            for (int k = 0; k < nucells; k++) {
                v_positions->x[idx] = i * a; // x
                v_positions->y[idx] = j * a; // y
                v_positions->z[idx] = k * a; // z
                v_positions->w[idx] = k * a; // z
                    // del mismo átomo
                v_positions->x[idx + 1] = (i + 0.5) * a;
                v_positions->y[idx + 1] = (j + 0.5) * a;
                v_positions->z[idx + 1] = k * a;
                v_positions->w[idx + 1] = k * a;

                v_positions->x[idx + 2] = (i + 0.5) * a;
                v_positions->y[idx + 2] = j * a;
                v_positions->z[idx + 2] = (k + 0.5) * a;
                v_positions->w[idx + 2] = (k + 0.5) * a;

                v_positions->x[idx + 3] = i * a;
                v_positions->y[idx + 3] = (j + 0.5) * a;
                v_positions->z[idx + 3] = (k + 0.5) * a;
                v_positions->w[idx + 3] = (k + 0.5) * a;

                idx += 4;
            }
        }
    }
}*/


void init_pos(Vector_SOA* restrict v_positions, const float rho)
{
    // Initialization of atom positions in an FCC crystal

    float a = cbrt(4.0 / rho);
    int nucells = ceil(cbrt((float)N / 4.0));
    int idx = 0;

    for (int i = 0; i < nucells * nucells * nucells; i++) {
        int k = i % nucells;
        int j = (i / nucells) % nucells;
        int i_ = i / (nucells * nucells);

        v_positions->x[idx] = i_ * a; // x
        v_positions->y[idx] = j * a; // y
        v_positions->z[idx] = k * a; // z

        v_positions->x[idx + 1] = (i_ + 0.5) * a;
        v_positions->y[idx + 1] = (j + 0.5) * a;
        v_positions->z[idx + 1] = k * a;

        v_positions->x[idx + 2] = (i_ + 0.5) * a;
        v_positions->y[idx + 2] = j * a;
        v_positions->z[idx + 2] = (k + 0.5) * a;

        v_positions->x[idx + 3] = i_ * a;
        v_positions->y[idx + 3] = (j + 0.5) * a;
        v_positions->z[idx + 3] = (k + 0.5) * a;

        idx += 4;
    }
}

/*
void init_pos(Vector_SOA* restrict v_positions, const float rho)
{
    // Inicialización de las posiciones de los átomos en un cristal FCC

    float a = cbrt(4.0 / rho);
    int nucells = ceil(cbrt((float)N / 4.0));
    int idx = 0;

    for (int i = 0; i < nucells * nucells * nucells; i++) {
        int x_idx = i % nucells;
        int y_idx = (i / nucells) % nucells;
        int z_idx = (i / (nucells * nucells)) % nucells;

        float x_offset = (i % 2 == 0 || (i / (nucells * nucells)) % 2 == 0) ? 0.0 : 0.5;
        float y_offset = (i % 4 < 2) ? 0.0 : 0.5;
        float z_offset = (i % 8 < 4) ? 0.0 : 0.5;

        v_positions->x[idx] = (x_idx + x_offset) * a;
        v_positions->y[idx] = (y_idx + y_offset) * a;
        v_positions->z[idx] = (z_idx + z_offset) * a;
        v_positions->w[idx] = (z_idx + z_offset) * a;

        idx++;
    }
}*/


void init_vel(Vector_SOA* restrict v_velocities, float* restrict temp, float* restrict ekin)
{
    // inicialización de velocidades aleatorias

    float sf, sumvx = 0.0, sumvy = 0.0, sumvz = 0.0, sumv2 = 0.0, sumvw = 0.0;

    for (int i = 0; i < N; i++) {
        v_velocities->x[i] = rand() / (float)RAND_MAX - 0.5;
        v_velocities->y[i] = rand() / (float)RAND_MAX - 0.5;
        v_velocities->z[i] = rand() / (float)RAND_MAX - 0.5;
        v_velocities->w[i] = rand() / (float)RAND_MAX - 0.5;

        sumvx += v_velocities->x[i];
        sumvy += v_velocities->y[i];
        sumvz += v_velocities->z[i];
        sumvw += v_velocities->w[i];
        sumv2 += v_velocities->x[i] * v_velocities->x[i] + v_velocities->y[i] * v_velocities->y[i]
               + v_velocities->z[i] * v_velocities->z[i] + v_velocities->w[i] * v_velocities->w[i] * 0; //podria no estar w?
    }

    sumvx /= (float)N;
    sumvy /= (float)N;
    sumvz /= (float)N;
    sumvw /= (float)N;
    *temp = sumv2 / (3.0 * N);
    *ekin = 0.5 * sumv2;
    sf = sqrt(T0 / *temp);

    for (int i = 0; i < N; i++) { // elimina la velocidad del centro de masa
        // y ajusta la temperatura
        v_velocities->x[i] = (v_velocities->x[i] - sumvx) * sf;
        v_velocities->y[i] = (v_velocities->y[i] - sumvy) * sf;
        v_velocities->z[i] = (v_velocities->z[i] - sumvz) * sf;
        v_velocities->w[i] = (v_velocities->w[i] - sumvw) * sf;
    }
}


static float minimum_image(float cordi, const float cell_length)
{
    // imagen más cercana

    if (cordi <= -0.5 * cell_length) {
        cordi += cell_length;
    } else if (cordi > 0.5 * cell_length) {
        cordi -= cell_length;
    }
    return cordi;// + cell_length * (cordi <= -0.5 * cell_length) - cell_length * (cordi > 0.5 * cell_length)
}


void forces(Vector_SOA* restrict v_positions, Vector_SOA* restrict v_forces, float* restrict epot, float* restrict pres,
            const float* restrict temp, const float rho, const float V, const float L)
{
    // calcula las fuerzas LJ (12-6)

    for (int i = 0; i < N; i++) {
        v_forces->x[i] = 0.0;
        v_forces->y[i] = 0.0;
        v_forces->z[i] = 0.0;
        v_forces->w[i] = 0.0;
    }
    float pres_vir = 0.0;
    float rcut2 = RCUT * RCUT;
    *epot = 0.0;

    for (int i = 0; i < N - 1; i++) {

        float xi = v_positions->x[i];
        float yi = v_positions->y[i];
        float zi = v_positions->z[i];
        float wi = v_positions->w[i];

        for (int j = i + 1; j < N; j++) {

            float xj = v_positions->x[j];
            float yj = v_positions->y[j];
            float zj = v_positions->z[j];
            float wj = v_positions->w[j];

            // distancia mínima entre r_i y r_j
            float rx = xi - xj;
            float ry = yi - yj;
            float rz = zi - zj;
            float rw = wi - wj;
            rx = minimum_image(rx, L);
            ry = minimum_image(ry, L);
            rz = minimum_image(rz, L);
            rw = minimum_image(rw, L);

            float rij2 = rx * rx + ry * ry + rz * rz;

            float condition = (rij2 <= rcut2) ? 1.0 : 0.0;
            // if (rij2 <= rcut2) {
                float r2inv = condition * (1.0 / rij2);
                float r6inv = condition * r2inv * r2inv * r2inv;

                float fr = condition * 24.0 * r2inv * r6inv * (2.0 * r6inv - 1.0);

                v_forces->x[i] += condition * fr * rx;
                v_forces->y[i] += condition * fr * ry;
                v_forces->z[i] += condition * fr * rz;
                v_forces->w[i] += condition * fr * rw;

                v_forces->x[j] -= condition * fr * rx;
                v_forces->y[j] -= condition * fr * ry;
                v_forces->z[j] -= condition * fr * rz;
                v_forces->w[j] -= condition * fr * rw;

                *epot += condition * (4.0 * r6inv * (r6inv - 1.0) - ECUT);
                pres_vir += condition * fr * rij2;
            //}
        }
    }
    pres_vir /= (V * 3.0);
    *pres = *temp * rho + pres_vir;
}


static float pbc(float cordi, const float cell_length)
{
    // condiciones periodicas de contorno coordenadas entre [0,L)
    if (cordi <= 0) {
        cordi += cell_length;
    } else if (cordi > cell_length) {
        cordi -= cell_length;
    }
    return cordi;
}


void velocity_verlet(Vector_SOA* restrict v_positions, Vector_SOA* restrict v_velocities, Vector_SOA* restrict v_forces, float* restrict epot,
                     float* restrict ekin, float* restrict pres, float* restrict temp, const float rho,
                     const float V, const float L)
{

    for (int i = 0; i < N; i++) { // actualizo posiciones
        v_positions->x[i] += v_velocities->x[i] * DT + 0.5 * v_forces->x[i] * DT * DT;
        v_positions->y[i] += v_velocities->y[i] * DT + 0.5 * v_forces->y[i] * DT * DT;
        v_positions->z[i] += v_velocities->z[i] * DT + 0.5 * v_forces->z[i] * DT * DT;
        v_positions->w[i] += v_velocities->w[i] * DT + 0.5 * v_forces->w[i] * DT * DT;

        v_positions->x[i] = pbc(v_positions->x[i], L);
        v_positions->y[i] = pbc(v_positions->y[i], L);
        v_positions->z[i] = pbc(v_positions->z[i], L);
        v_positions->w[i] = pbc(v_positions->w[i], L);

        v_velocities->x[i] += 0.5 * v_forces->x[i] * DT;
        v_velocities->y[i] += 0.5 * v_forces->y[i] * DT;
        v_velocities->z[i] += 0.5 * v_forces->z[i] * DT;
        v_velocities->w[i] += 0.5 * v_forces->w[i] * DT;
    }

    forces(v_positions, v_forces, epot, pres, temp, rho, V, L); // actualizo fuerzas

    float sumv2 = 0.0;
    for (int i = 0; i < N; i++) { // actualizo velocidades
        v_velocities->x[i] += 0.5 * v_forces->x[i] * DT;
        v_velocities->y[i] += 0.5 * v_forces->y[i] * DT;
        v_velocities->z[i] += 0.5 * v_forces->z[i] * DT;
        v_velocities->w[i] += 0.5 * v_forces->w[i] * DT;
        
        sumv2 += v_velocities->x[i] * v_velocities->x[i] + v_velocities->y[i] * v_velocities->y[i]
               + v_velocities->z[i] * v_velocities->z[i] + v_velocities->w[i] * v_velocities->w[i] * 0;
    }

    *ekin = 0.5 * sumv2;
    *temp = sumv2 / (3.0 * N);
}
