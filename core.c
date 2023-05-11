#include "core.h"
#include "parameters.h"

#include <math.h>
#include <stdlib.h> // rand()

#define ECUT (4.0 * (pow(RCUT, -12) - pow(RCUT, -6)))


void init_pos(struct Atoms* positions, const double rho)
{
    // inicialización de las posiciones de los átomos en un cristal FCC

    double a = cbrt(4.0 / rho);
    int nucells = ceil(cbrt((double)N / 4.0));
    int idx = 0;

    for (int i = 0; i < nucells; i++) {
        for (int j = 0; j < nucells; j++) {
            for (int k = 0; k < nucells; k++) {
                positions->x[idx] = i * a;
                positions->y[idx] = j * a;
                positions->z[idx] = k * a;

                positions->x[idx+1] = (i + 0.5) * a;
                positions->y[idx+1] = (j + 0.5) * a;
                positions->z[idx+1] = k * a;

                positions->x[idx+2] = (i + 0.5) * a;
                positions->y[idx+2] = j * a;
                positions->z[idx+2] = (k + 0.5) * a;

                positions->x[idx+3] = i * a;
                positions->y[idx+3] = (j + 0.5) * a;
                positions->z[idx+3] = (k + 0.5) * a;

                idx += 4;
            }
        }
    }
}


void init_vel(struct Atoms* velocities, double* temp, double* ekin)
{
    // inicialización de velocidades aleatorias

    double sf, sumvx = 0.0, sumvy = 0.0, sumvz = 0.0, sumv2 = 0.0;

    for (int i = 0; i < N; i++) {
        velocities->x[i] = rand() / (double)RAND_MAX - 0.5;
        velocities->y[i] = rand() / (double)RAND_MAX - 0.5;
        velocities->z[i] = rand() / (double)RAND_MAX - 0.5;

        sumvx += velocities->x[i];
        sumvy += velocities->y[i];
        sumvz += velocities->z[i];
        sumv2 += velocities->x[i] * velocities->x[i] + velocities->y[i] * velocities->y[i] + velocities->z[i] * velocities->z[i];
    }

    sumvx /= (double)N;
    sumvy /= (double)N;
    sumvz /= (double)N;
    *temp = sumv2 / (3.0 * N);
    *ekin = 0.5 * sumv2;
    sf = sqrt(T0 / *temp);

    for (int i = 0; i < N; i++) { // elimina la velocidad del centro de masa
        // y ajusta la temperatura
        velocities->x[i] = (velocities->x[i] - sumvx) * sf;
        velocities->y[i] = (velocities->y[i] - sumvy) * sf;
        velocities->z[i] = (velocities->y[i] - sumvz) * sf;
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


void forces(const struct Atoms* positions, struct Atoms* forcess, double* epot, double* pres,
            const double* temp, const double rho, const double V, const double L)
{
    // calcula las fuerzas LJ (12-6)

    for (int i = 0; i < N; i++) {
        forcess->x[i] = 0.0;
        forcess->y[i] = 0.0;
        forcess->z[i] = 0.0;
    }

    double pres_vir = 0.0;
    double rcut2 = RCUT * RCUT;
    *epot = 0.0;

    for (int i = 0; i < N-1; i++) {

        double xi = positions->x[i];
        double yi = positions->y[i];
        double zi = positions->z[i];

        for (int j = i + 1; j < N; j++) {

            double xj = positions->x[j];
            double yj = positions->y[j];
            double zj = positions->z[j];

            // distancia mínima entre r_i y r_j
            double rx = xi - xj;
            rx = minimum_image(rx, L);
            double ry = yi - yj;
            ry = minimum_image(ry, L);
            double rz = zi - zj;
            rz = minimum_image(rz, L);

            double rij2 = rx * rx + ry * ry + rz * rz;

            if (rij2 <= rcut2) {
                double r2inv = 1.0 / rij2;
                double r6inv = r2inv * r2inv * r2inv;

                double fr = 24.0 * r2inv * r6inv * (2.0 * r6inv - 1.0);

                forcess->x[i] += fr * rx;
                forcess->y[i] += fr * ry;
                forcess->z[i] += fr * rz;

                forcess->x[j] -= fr * rx;
                forcess->y[j] -= fr * ry;
                forcess->z[j] -= fr * rz;

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


void velocity_verlet(struct Atoms* positions, struct Atoms* velocities, struct Atoms* forcess, double* epot,
                     double* ekin, double* pres, double* temp, const double rho,
                     const double V, const double L)
{

    for (int i = 0; i < N; i++) { // actualizo posiciones
        positions->x[i] += velocities->x[i] * DT + 0.5 * forcess->x[i] * DT * DT;
        positions->y[i] += velocities->y[i] * DT + 0.5 * forcess->y[i] * DT * DT;
        positions->z[i] += velocities->z[i] * DT + 0.5 * forcess->z[i] * DT * DT;

        positions->x[i] = pbc(positions->x[i], L);
        positions->y[i] = pbc(positions->y[i], L);
        positions->z[i] = pbc(positions->z[i], L);

        velocities->x[i] += 0.5 * forcess->x[i] * DT;
        velocities->y[i] += 0.5 * forcess->y[i] * DT;
        velocities->z[i] += 0.5 * forcess->z[i] * DT;
    }

    forces(positions, forcess, epot, pres, temp, rho, V, L); // actualizo fuerzas

    double sumv2 = 0.0;
    for (int i = 0; i < N; i++) { // actualizo velocidades
        velocities->x[i] += 0.5 * forcess->x[i] * DT;
        velocities->y[i] += 0.5 * forcess->y[i] * DT;
        velocities->z[i] += 0.5 * forcess->z[i] * DT;

        sumv2 += velocities->x[i] * velocities->x[i] + velocities->y[i] * velocities->y[i]
            + velocities->z[i] * velocities->z[i];
    }

    *ekin = 0.5 * sumv2;
    *temp = sumv2 / (3.0 * N);
}
