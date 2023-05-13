#include "core.h"
#include "parameters.h"

#include <math.h>
#include <stdlib.h> // rand()

#define ECUT (4.0 * (pow(RCUT, -12) - pow(RCUT, -6)))


void init_pos(Vector_SOA* v_positions, const double rho)
{
    // inicialización de las posiciones de los átomos en un cristal FCC

    double a = cbrt(4.0 / rho);
    int nucells = ceil(cbrt((double)N / 4.0));
    int idx = 0;

    for (int i = 0; i < nucells; i++) {
        for (int j = 0; j < nucells; j++) {
            for (int k = 0; k < nucells; k++) {
                v_positions->x[idx] = i * a; // x
                v_positions->y[idx] = j * a; // y
                v_positions->z[idx] = k * a; // z
                    // del mismo átomo
                v_positions->x[idx + 1] = (i + 0.5) * a;
                v_positions->y[idx + 1] = (j + 0.5) * a;
                v_positions->z[idx + 1] = k * a;

                v_positions->x[idx + 2] = (i + 0.5) * a;
                v_positions->y[idx + 2] = j * a;
                v_positions->z[idx + 2] = (k + 0.5) * a;

                v_positions->x[idx + 3] = i * a;
                v_positions->y[idx + 3] = (j + 0.5) * a;
                v_positions->z[idx + 3] = (k + 0.5) * a;

                idx += 4;
            }
        }
    }
}


void init_vel(Vector_SOA* v_velocities, double* temp, double* ekin)
{
    // inicialización de velocidades aleatorias

    double sf, sumvx = 0.0, sumvy = 0.0, sumvz = 0.0, sumv2 = 0.0;

    for (int i = 0; i < N; i++) {
        v_velocities->x[i] = rand() / (double)RAND_MAX - 0.5;
        v_velocities->y[i] = rand() / (double)RAND_MAX - 0.5;
        v_velocities->z[i] = rand() / (double)RAND_MAX - 0.5;

        sumvx += v_velocities->x[i];
        sumvy += v_velocities->y[i];
        sumvz += v_velocities->z[i];
        sumv2 += v_velocities->x[i] * v_velocities->x[i] + v_velocities->y[i] * v_velocities->y[i]
            + v_velocities->z[i] * v_velocities->z[i];
    }

    sumvx /= (double)N;
    sumvy /= (double)N;
    sumvz /= (double)N;
    *temp = sumv2 / (3.0 * N);
    *ekin = 0.5 * sumv2;
    sf = sqrt(T0 / *temp);

    for (int i = 0; i < N; i++) { // elimina la velocidad del centro de masa
        // y ajusta la temperatura
        v_velocities->x[i] = (v_velocities->x[i] - sumvx) * sf;
        v_velocities->y[i] = (v_velocities->y[i] - sumvy) * sf;
        v_velocities->z[i] = (v_velocities->z[i] - sumvz) * sf;
    }
}


static double minimum_image(double cordi, const double cell_length)
{
    // imagen más cercana

    if (cordi <= -0.5 * cell_length) {
        cordi += cell_length;
    } else if (cordi > 0.5 * cell_length) {
        cordi -= cell_length;
    }
    return cordi;
}


void forces(Vector_SOA* v_positions, Vector_SOA* v_forces, double* epot, double* pres,
            const double* temp, const double rho, const double V, const double L)
{
    // calcula las fuerzas LJ (12-6)

    for (int i = 0; i < N; i++) {
        v_forces->x[i] = 0.0;
        v_forces->y[i] = 0.0;
        v_forces->z[i] = 0.0;
    }
    double pres_vir = 0.0;
    double rcut2 = RCUT * RCUT;
    *epot = 0.0;

    for (int i = 0; i < N - 1; i++) {

        double xi = v_positions->x[i];
        double yi = v_positions->y[i];
        double zi = v_positions->z[i];

        for (int j = i + 1; j < N; j++) {

            double xj = v_positions->x[j];
            double yj = v_positions->y[j];
            double zj = v_positions->z[j];

            // distancia mínima entre r_i y r_j
            double rx = xi - xj;
            double ry = yi - yj;
            double rz = zi - zj;
            rx = minimum_image(rx, L);
            ry = minimum_image(ry, L);
            rz = minimum_image(rz, L);

            double rij2 = rx * rx + ry * ry + rz * rz;

            if (rij2 <= rcut2) {
                double r2inv = 1.0 / rij2;
                double r6inv = r2inv * r2inv * r2inv;

                double fr = 24.0 * r2inv * r6inv * (2.0 * r6inv - 1.0);

                v_forces->x[i] += fr * rx;
                v_forces->y[i] += fr * ry;
                v_forces->z[i] += fr * rz;

                v_forces->x[j] -= fr * rx;
                v_forces->y[j] -= fr * ry;
                v_forces->z[j] -= fr * rz;

                *epot += 4.0 * r6inv * (r6inv - 1.0) - ECUT;
                pres_vir += fr * rij2;
            }
        }
    }
    pres_vir /= (V * 3.0);
    *pres = *temp * rho + pres_vir;
}


static double pbc(double cordi, const double cell_length)
{
    // condiciones periodicas de contorno coordenadas entre [0,L)
    if (cordi <= 0) {
        cordi += cell_length;
    } else if (cordi > cell_length) {
        cordi -= cell_length;
    }
    return cordi;
}


void velocity_verlet(Vector_SOA* v_positions, Vector_SOA* v_velocities, Vector_SOA* v_forces, double* epot,
                     double* ekin, double* pres, double* temp, const double rho,
                     const double V, const double L)
{

    for (int i = 0; i < N; i++) { // actualizo posiciones
        v_positions->x[i] += v_velocities->x[i] * DT + 0.5 * v_forces->x[i] * DT * DT;
        v_positions->y[i] += v_velocities->y[i] * DT + 0.5 * v_forces->y[i] * DT * DT;
        v_positions->z[i] += v_velocities->z[i] * DT + 0.5 * v_forces->z[i] * DT * DT;

        v_positions->x[i] = pbc(v_positions->x[i], L);
        v_positions->y[i] = pbc(v_positions->y[i], L);
        v_positions->z[i] = pbc(v_positions->z[i], L);

        v_velocities->x[i] += 0.5 * v_forces->x[i] * DT;
        v_velocities->y[i] += 0.5 * v_forces->y[i] * DT;
        v_velocities->z[i] += 0.5 * v_forces->z[i] * DT;
    }

    forces(v_positions, v_forces, epot, pres, temp, rho, V, L); // actualizo fuerzas

    double sumv2 = 0.0;
    for (int i = 0; i < N; i++) { // actualizo velocidades
        v_velocities->x[i] += 0.5 * v_forces->x[i] * DT;
        v_velocities->y[i] += 0.5 * v_forces->y[i] * DT;
        v_velocities->z[i] += 0.5 * v_forces->z[i] * DT;

        sumv2 += v_velocities->x[i] * v_velocities->x[i] + v_velocities->y[i] * v_velocities->y[i]
            + v_velocities->z[i] * v_velocities->z[i];
    }

    *ekin = 0.5 * sumv2;
    *temp = sumv2 / (3.0 * N);
}
