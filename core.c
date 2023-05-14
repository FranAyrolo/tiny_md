#include "core.h"
#include "parameters.h"

#include <math.h>
#include <stdlib.h> // rand()

#define ECUT (4.0 * (pow(RCUT, -12) - pow(RCUT, -6)))

void init_pos(Vector_SOA* restrict v_positions, const double rho)
{
    // Initialization of atom positions in an FCC crystal

    double a = cbrt(4.0 / rho);
    int nucells = ceil(cbrt((double)N / 4.0));
    int idx = 0;

    for (int i = 0; i < nucells * nucells * nucells; i++) {
        int k = i % nucells;
        int j = (i / nucells) % nucells;
        int i_ = i / (nucells * nucells);

        v_positions->x[idx] = i_ * a; // x
        v_positions->x[idx + 1] = (i_ + 0.5) * a;
        v_positions->x[idx + 2] = (i_ + 0.5) * a;
        v_positions->x[idx + 3] = i_ * a;

        v_positions->y[idx] = j * a; // y
        v_positions->y[idx + 1] = (j + 0.5) * a;
        v_positions->y[idx + 2] = j * a;
        v_positions->y[idx + 3] = (j + 0.5) * a;

        v_positions->z[idx] = k * a; // z
        v_positions->z[idx + 1] = k * a;
        v_positions->z[idx + 2] = (k + 0.5) * a;
        v_positions->z[idx + 3] = (k + 0.5) * a;

        idx += 4;
    }
}


void init_vel(Vector_SOA* restrict v_velocities, double* restrict temp, double* restrict ekin)
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

inline static double minimum_image(double cordi, const double cell_length)
{
    const double half_cell_length = 0.5 * cell_length;
    cordi -= cell_length * (cordi > half_cell_length);
    cordi += cell_length * (cordi <= -half_cell_length);
    return cordi;
}


void forces(Vector_SOA* restrict v_positions, Vector_SOA* restrict v_forces, double* restrict epot, double* restrict pres,
            const double* restrict temp, const double rho, const double V, const double L)
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
            rx = minimum_image(rx, L);
            double ry = yi - yj;
            ry = minimum_image(ry, L);
            double rz = zi - zj;
            rz = minimum_image(rz, L);        

            double rij2 = rx * rx + ry * ry + rz * rz;

            double condition = (rij2 <= rcut2) ? 1.0 : 0.0;
            double r2inv = condition * (1.0 / rij2);
            double r6inv = condition * r2inv * r2inv * r2inv;

            double fr = condition * 24.0 * r2inv * r6inv * (2.0 * r6inv - 1.0);

            v_forces->x[i] += condition * fr * rx;
            v_forces->y[i] += condition * fr * ry;
            v_forces->z[i] += condition * fr * rz;

            v_forces->x[j] -= condition * fr * rx;
            v_forces->y[j] -= condition * fr * ry;
            v_forces->z[j] -= condition * fr * rz;

            *epot += condition * (4.0 * r6inv * (r6inv - 1.0) - ECUT);
            pres_vir += condition * fr * rij2;
        }
    }

    pres_vir /= (V * 3.0);
    *pres = *temp * rho + pres_vir;
}


inline static double pbc(double cordi, const double cell_length)
{
    cordi += cell_length * (cordi <= 0);
    cordi -= cell_length * (cordi > cell_length && cordi > 0);
    return cordi;
}



void velocity_verlet(Vector_SOA* restrict v_positions, Vector_SOA* restrict v_velocities, Vector_SOA* restrict v_forces, double* restrict epot,
                     double* restrict ekin, double* restrict pres, double* restrict temp, const double rho,
                     const double V, const double L)
{

    for (int i = 0; i < N; i++) { // actualizo posiciones
        v_positions->x[i] += v_velocities->x[i] * DT + 0.5 * v_forces->x[i] * DT * DT;
        v_positions->x[i] = pbc(v_positions->x[i], L);
        v_velocities->x[i] += 0.5 * v_forces->x[i] * DT;

        v_positions->y[i] += v_velocities->y[i] * DT + 0.5 * v_forces->y[i] * DT * DT;
        v_positions->y[i] = pbc(v_positions->y[i], L);
        v_velocities->y[i] += 0.5 * v_forces->y[i] * DT;

        v_positions->z[i] += v_velocities->z[i] * DT + 0.5 * v_forces->z[i] * DT * DT;
        v_positions->z[i] = pbc(v_positions->z[i], L);
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
