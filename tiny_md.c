#define _XOPEN_SOURCE 500  // M_PI
#include "core.h"
#include "parameters.h"
#include "wtime.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

int main()
{
    FILE *file_xyz, *file_thermo;
    file_xyz = fopen("trajectory.xyz", "w");
    file_thermo = fopen("thermo.log", "w");
    double Ekin, Epot, Temp, Pres; // variables macroscopicas
    double Rho, cell_V, cell_L, tail, Etail, Ptail;

    struct Atoms* positions = malloc(sizeof(struct Atoms));
    struct Atoms* velocities = malloc(sizeof(struct Atoms));
    struct Atoms* forcess = malloc(sizeof(struct Atoms));

    positions->x = (double*)malloc(sizeof(double) * N);
    positions->y = (double*)malloc(sizeof(double) * N);
    positions->z = (double*)malloc(sizeof(double) * N);
    velocities->x = (double*)malloc(sizeof(double) * N);
    velocities->y = (double*)malloc(sizeof(double) * N);
    velocities->z = (double*)malloc(sizeof(double) * N);
    forcess->x = (double*)malloc(sizeof(double) * N);
    forcess->y = (double*)malloc(sizeof(double) * N);
    forcess->z = (double*)malloc(sizeof(double) * N);

    printf("# Número de partículas:      %d\n", N);
    printf("# Temperatura de referencia: %.2f\n", T0);
    printf("# Pasos de equilibración:    %d\n", TEQ);
    printf("# Pasos de medición:         %d\n", TRUN - TEQ);
    printf("# (mediciones cada %d pasos)\n", TMES);
    printf("# densidad, volumen, energía potencial media, presión media\n");
    fprintf(file_thermo, "# t Temp Pres Epot Etot\n");

    srand(SEED);
    double t = 0.0, sf = 0.0;
    double Rhob;
    Rho = RHOI;
    init_pos(positions, Rho);
    double start = wtime();
    for (int m = 0; m < 9; m++) {
        Rhob = Rho;
        Rho = RHOI - 0.1 * (double)m;
        cell_V = (double)N / Rho;
        cell_L = cbrt(cell_V);
        tail = 16.0 * M_PI * Rho * ((2.0 / 3.0) * pow(RCUT, -9) - pow(RCUT, -3)) / 3.0;
        Etail = tail * (double)N;
        Ptail = tail * Rho;

        int i = 0;
        sf = cbrt(Rhob / Rho);
        for (int k = 0; k < N; k++) { // reescaleo posiciones a nueva densidad
            positions->x[k] *= sf;
            positions->y[k] *= sf;
            positions->z[k] *= sf;
        }
        init_vel(velocities, &Temp, &Ekin);
        forces(positions, forcess, &Epot, &Pres, &Temp, Rho, cell_V, cell_L);

        for (i = 1; i < TEQ; i++) { // loop de equilibracion

            velocity_verlet(positions, velocities, forcess, &Epot, &Ekin, &Pres, &Temp, Rho, cell_V, cell_L);

            sf = sqrt(T0 / Temp);
            for (int k = 0; k < N; k++) { // reescaleo de velocidades
                velocities->x[k] *= sf;
                velocities->y[k] *= sf;
                velocities->z[k] *= sf;
            }
        }

        int mes = 0;
        double epotm = 0.0, presm = 0.0;
        for (i = TEQ; i < TRUN; i++) { // loop de medicion

            velocity_verlet(positions, velocities, forcess, &Epot, &Ekin, &Pres, &Temp, Rho, cell_V, cell_L);

            sf = sqrt(T0 / Temp);
            for (int k = 0; k < N; k++) { // reescaleo de velocidades
                velocities->x[k] *= sf;
                velocities->y[k] *= sf;
                velocities->z[k] *= sf;
            }

            if (i % TMES == 0) {
                Epot += Etail;
                Pres += Ptail;

                epotm += Epot;
                presm += Pres;
                mes++;

                fprintf(file_thermo, "%f %f %f %f %f\n", t, Temp, Pres, Epot, Epot + Ekin);
                fprintf(file_xyz, "%d\n\n", N);
                for (int k = 0; k < N; k++) {
                    fprintf(file_xyz, "Ar %e %e %e\n", positions->x[k], positions->y[k], positions->z[k]);
                }
            }

            t += DT;
        }
        printf("%f\t%f\t%f\t%f\n", Rho, cell_V, epotm / (double)mes, presm / (double)mes);
    }

    double elapsed = wtime() - start;
    printf("# Tiempo total de simulación = %f segundos\n", elapsed);
    printf("# Tiempo simulado = %f [fs]\n", t * 1.6);
    printf("# ns/day = %f\n", (1.6e-6 * t) / elapsed * 86400);
    //                       ^1.6 fs -> ns       ^sec -> day
    free(positions->x);
    free(positions->y);
    free(positions->z);
    free(positions);
    free(velocities->x);
    free(velocities->y);
    free(velocities->z);
    free(velocities);
    free(forcess->x);
    free(forcess->y);
    free(forcess->z);
    free(forcess);
    return 0;
}
