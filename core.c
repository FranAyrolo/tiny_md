#include "core.h"
#include "parameters.h"

#include <math.h>
#include <stdlib.h> // rand()
#include <immintrin.h>

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

void init_pos(Vector3D* positions, const double rho)
{
    double a = cbrt(4.0 / rho);
    int nucells = ceil(cbrt((double)N / 4.0));
    int numAtoms = nucells * nucells * nucells;
    int idx = idx = 0;
    for(int linearIdx = 0; linearIdx < numAtoms; linearIdx++) {
        int i = linearIdx / (nucells * nucells);
        int j = (linearIdx / nucells) % nucells;
        int k = linearIdx % nucells;

        positions[idx] = (Vector3D) { i * a, j * a, k * a };
        positions[idx + 1] = (Vector3D) { (i + 0.5) * a, (j + 0.5) * a, k * a };
        positions[idx + 2] = (Vector3D) { (i + 0.5) * a, j * a, (k + 0.5) * a };
        positions[idx + 3] = (Vector3D) { i * a, (j + 0.5) * a, (k + 0.5) * a };
        idx += 4;
    }
}

void init_vel(Vector3D* velocities, double* temp, double* ekin)
{
    // inicialización de velocidades aleatorias

    double sf, sumvx = 0.0, sumvy = 0.0, sumvz = 0.0, sumv2 = 0.0;

    for (int i = 0; i < N; i++) {
        velocities[i] = (Vector3D) {
            rand()/(double)RAND_MAX-0.5,
            rand()/(double)RAND_MAX-0.5,
            rand()/(double)RAND_MAX-0.5
        };

        sumvx += velocities[i].x;
        sumvy += velocities[i].y;
        sumvz += velocities[i].z;
        sumv2 += velocities[i].x * velocities[i].x + velocities[i].y *
            velocities[i].y + velocities[i].z * velocities[i].z;
    }

    sumvx /= (double)N;
    sumvy /= (double)N;
    sumvz /= (double)N;
    *temp = sumv2 / (3.0 * N);
    *ekin = 0.5 * sumv2;
    sf = sqrt(T0 / *temp);

    for (int i = 0; i < N; i++) { // elimina la velocidad del centro de masa
        // y ajusta la temperatura
        velocities[i].x = (velocities[i].x - sumvx) * sf;
        velocities[i].y = (velocities[i].y - sumvy) * sf;
        velocities[i].z = (velocities[i].z - sumvz) * sf;
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


void forces(const Vector3D* positions, Vector3D* forces3D, double* epot, double* pres,
            const double* temp, const double rho, const double V, const double L)
{
    // calcula las fuerzas LJ (12-6)

    for (int i = 0; i < N; i++) {
        forces3D[i] = (Vector3D) { 0.0, 0.0, 0.0 };
    }
    double pres_vir = 0.0;
    double rcut2 = RCUT * RCUT;
    *epot = 0.0;

    for (int i = 0; i < N - 1; i++) {

        double xi = positions[i].x;
        double yi = positions[i].y;
        double zi = positions[i].z;

        double forces_ix = 0.0;
        double forces_iy = 0.0;
        double forces_iz = 0.0;

        for (int j = i + 1; j < N; j++) {

            double xj = positions[j].x;
            double yj = positions[j].y;
            double zj = positions[j].z;

            // distancia mínima entre r_i y r_j
            double rx = xi - xj;
            rx = minimum_image(rx, L);
            double ry = yi - yj;
            ry = minimum_image(ry, L);
            double rz = zi - zj;
            rz = minimum_image(rz, L);

            double rij2 = rx * rx + ry * ry + rz * rz;

            double r2inv = 1.0 / rij2;
            double r6inv = r2inv * r2inv * r2inv;

            double fr = 24.0 * r2inv * r6inv * (2.0 * r6inv - 1.0);

            double condition = (rij2 <= rcut2) ? 1.0 : 0.0;

            forces_ix += condition * fr * rx;
            forces_iy += condition * fr * ry;
            forces_iz += condition * fr * rz;

            forces3D[j].x -= condition * fr * rx;
            forces3D[j].y -= condition * fr * ry;
            forces3D[j].z -= condition * fr * rz;

            *epot += condition * (4.0 * r6inv * (r6inv - 1.0) - ECUT);
            pres_vir += condition * fr * rij2;
        }
        forces3D[i].x += forces_ix;
        forces3D[i].y += forces_iy;
        forces3D[i].z += forces_iz;
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


void velocity_verlet(Vector3D* restrict positions, Vector3D* restrict velocities, Vector3D* restrict forces3D,
        double* restrict epot, double* restrict ekin, double* restrict pres, double* restrict temp, const double rho, const double V, const double L)
{

    for (int i = 0; i < N; i++) { // actualizo posiciones
        positions[i].x += velocities[i].x * DT + 0.5 * forces3D[i].x * DT * DT;
        positions[i].y += velocities[i].y * DT + 0.5 * forces3D[i].y * DT * DT;
        positions[i].z += velocities[i].z * DT + 0.5 * forces3D[i].z * DT * DT;

        positions[i].x = pbc(positions[i].x, L);
        positions[i].y = pbc(positions[i].y, L);
        positions[i].z = pbc(positions[i].z, L);

        velocities[i].x += 0.5 * forces3D[i].x * DT;
        velocities[i].y += 0.5 * forces3D[i].y * DT;
        velocities[i].z += 0.5 * forces3D[i].z * DT;
    }

    forces(positions, forces3D, epot, pres, temp, rho, V, L); // actualizo fuerzas

    double sumv2 = 0.0;
    for (int i = 0; i < N; i++) { // actualizo velocidades
        velocities[i].x += 0.5 * forces3D[i].x * DT;
        velocities[i].y += 0.5 * forces3D[i].y * DT;
        velocities[i].z += 0.5 * forces3D[i].z * DT;

        sumv2 += velocities[i].x * velocities[i].x + velocities[i].y * velocities[i].y
            + velocities[i].z * velocities[i].z;
    }

    *ekin = 0.5 * sumv2;
    *temp = sumv2 / (3.0 * N);
}
