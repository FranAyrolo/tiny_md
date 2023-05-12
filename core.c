#include "core.h"
#include "parameters.h"

#include <math.h>
#include <stdlib.h> // rand()

#define ECUT (4.0 * (pow(RCUT, -12) - pow(RCUT, -6)))

AtomSystem *createAtomSystem(int numAtoms) {
    AtomSystem *system = malloc(sizeof(AtomSystem));
    system->positions = malloc(numAtoms * sizeof(Vector3D));
    system->velocities = malloc(numAtoms * sizeof(Vector3D));
    system->forces = malloc(numAtoms * sizeof(Vector3D));
    return system;
}

void destroyAtomSystem(AtomSystem* system) {
    free(system->positions);
    free(system->velocities);
    free(system->forces);
    free(system);
}

void init_pos(AtomSystem* system, const double rho)
{
    double a = cbrt(4.0 / rho);
    int nucells = ceil(cbrt((double)N / 4.0));
    int numAtoms = nucells * nucells * nucells;
    int idx = idx = 0;
    for(int linearIdx = 0; linearIdx < numAtoms; linearIdx++) {
        int i = linearIdx / (nucells * nucells);
        int j = (linearIdx / nucells) % nucells;
        int k = linearIdx % nucells;

        system->positions[idx] = (Vector3D) { i * a, j * a, k * a };
        system->positions[idx + 1] = (Vector3D) { (i + 0.5) * a, (j + 0.5) * a, k * a };
        system->positions[idx + 2] = (Vector3D) { (i + 0.5) * a, j * a, (k + 0.5) * a };
        system->positions[idx + 3] = (Vector3D) { i * a, (j + 0.5) * a, (k + 0.5) * a };
        idx += 4;
    }
}

void init_vel(AtomSystem* system, double* temp, double* ekin)
{
    // inicialización de velocidades aleatorias

    double sf, sumvx = 0.0, sumvy = 0.0, sumvz = 0.0, sumv2 = 0.0;

    for (int i = 0; i < N; i++) {
        system->velocities[i] = (Vector3D) {
            rand() / (double)RAND_MAX - 0.5,
            rand() / (double)RAND_MAX - 0.5,
            rand() / (double)RAND_MAX - 0.5
        };

        sumvx += system->velocities[i].x;
        sumvy += system->velocities[i].y;
        sumvz += system->velocities[i].z;
        sumv2 += system->velocities[i].x * system->velocities[i].x + system->velocities[i].y *
            system->velocities[i].y + system->velocities[i].z * system->velocities[i].z;
    }

    sumvx /= (double)N;
    sumvy /= (double)N;
    sumvz /= (double)N;
    *temp = sumv2 / (3.0 * N);
    *ekin = 0.5 * sumv2;
    sf = sqrt(T0 / *temp);

    for (int i = 0; i < N; i++) { // elimina la velocidad del centro de masa
        // y ajusta la temperatura
        system->velocities[i].x = (system->velocities[i].x - sumvx) * sf;
        system->velocities[i].y = (system->velocities[i].y - sumvy) * sf;
        system->velocities[i].z = (system->velocities[i].z - sumvz) * sf;
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


void forces(const AtomSystem* system, double* epot, double* pres,
            const double* temp, const double rho, const double V, const double L)
{
    // calcula las fuerzas LJ (12-6)

    for (int i = 0; i < N; i++) {
        system->forces[i] = (Vector3D) { 0.0, 0.0, 0.0 };
    }
    double pres_vir = 0.0;
    double rcut2 = RCUT * RCUT;
    *epot = 0.0;

    for (int i = 0; i < N - 1; i++) {

        double xi = system->positions[i].x;
        double yi = system->positions[i].y;
        double zi = system->positions[i].z;

        for (int j = i + 1; j < N; j++) {

            double xj = system->positions[j].x;
            double yj = system->positions[j].y;
            double zj = system->positions[j].z;

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

                system->forces[i].x += fr * rx;
                system->forces[i].y += fr * ry;
                system->forces[i].z += fr * rz;

                system->forces[j].x -= fr * rx;
                system->forces[j].y -= fr * ry;
                system->forces[j].z -= fr * rz;

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


void velocity_verlet(AtomSystem* system, double* epot,
                     double* ekin, double* pres, double* temp, const double rho,
                     const double V, const double L)
{

    for (int i = 0; i < N; i++) { // actualizo posiciones
        system->positions[i].x += system->velocities[i].x * DT + 0.5 * system->forces[i].x * DT * DT;
        system->positions[i].y += system->velocities[i].y * DT + 0.5 * system->forces[i].y * DT * DT;
        system->positions[i].z += system->velocities[i].z * DT + 0.5 * system->forces[i].z * DT * DT;

        system->positions[i].x = pbc(system->positions[i].x, L);
        system->positions[i].y = pbc(system->positions[i].y, L);
        system->positions[i].z = pbc(system->positions[i].z, L);

        system->velocities[i].x += 0.5 * system->forces[i].x * DT;
        system->velocities[i].y += 0.5 * system->forces[i].y * DT;
        system->velocities[i].z += 0.5 * system->forces[i].z * DT;
    }

    forces(system, epot, pres, temp, rho, V, L); // actualizo fuerzas

    double sumv2 = 0.0;
    for (int i = 0; i < N; i++) { // actualizo velocidades
        system->velocities[i].x += 0.5 * system->forces[i].x * DT;
        system->velocities[i].y += 0.5 * system->forces[i].y * DT;
        system->velocities[i].z += 0.5 * system->forces[i].z * DT;

        sumv2 += system->velocities[i].x * system->velocities[i].x + system->velocities[i].y * system->velocities[i].y
            + system->velocities[i].z * system->velocities[i].z;
    }

    *ekin = 0.5 * sumv2;
    *temp = sumv2 / (3.0 * N);
}
