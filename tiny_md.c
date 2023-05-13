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
    float Ekin, Epot, Temp, Pres; // variables macroscopicas
    float Rho, cell_V, cell_L, tail, Etail, Ptail;

    Vector_SOA* v_positions = (Vector_SOA*)malloc(sizeof(Vector_SOA));
    /*v_positions->x = (float*)malloc(N * sizeof(float));
    v_positions->y = (float*)malloc(N * sizeof(float));
    v_positions->z = (float*)malloc(N * sizeof(float));
    v_positions->w = (float*)malloc(N * sizeof(float));*/
    Vector_SOA* v_velocities = (Vector_SOA*)malloc(sizeof(Vector_SOA));
    /*v_velocities->x = (float*)malloc(N * sizeof(float));
    v_velocities->y = (float*)malloc(N * sizeof(float));
    v_velocities->z = (float*)malloc(N * sizeof(float));
    v_velocities->w = (float*)malloc(N * sizeof(float));*/
    Vector_SOA* v_forces = (Vector_SOA*)malloc(sizeof(Vector_SOA));
    /*v_forces->x = (float*)malloc(N * sizeof(float));
    v_forces->y = (float*)malloc(N * sizeof(float));
    v_forces->z = (float*)malloc(N * sizeof(float));
    v_forces->w = (float*)malloc(N * sizeof(float));*/
    

    printf("# Número de partículas:      %d\n", N);
    printf("# Temperatura de referencia: %.2f\n", T0);
    printf("# Pasos de equilibración:    %d\n", TEQ);
    printf("# Pasos de medición:         %d\n", TRUN - TEQ);
    printf("# (mediciones cada %d pasos)\n", TMES);
    printf("# densidad, volumen, energía potencial media, presión media\n");
    fprintf(file_thermo, "# t Temp Pres Epot Etot\n");

    srand(SEED);
    float t = 0.0, sf;
    float Rhob;
    Rho = RHOI;
    init_pos(v_positions, Rho);
    double start = wtime();
    for (int m = 0; m < 9; m++) {
        Rhob = Rho;
        Rho = RHOI - 0.1 * (float)m;
        cell_V = (float)N / Rho;
        cell_L = cbrt(cell_V);
        tail = 16.0 * M_PI * Rho * ((2.0 / 3.0) * pow(RCUT, -9) - pow(RCUT, -3)) / 3.0;
        Etail = tail * (float)N;
        Ptail = tail * Rho;

        int i = 0;
        sf = cbrt(Rhob / Rho);
        for (int k = 0; k < N; k++) { // reescaleo posiciones a nueva densidad
            v_positions->x[k] *= sf;
            v_positions->y[k] *= sf;
            v_positions->z[k] *= sf;
            v_positions->w[k] *= sf;
        }
        init_vel(v_velocities, &Temp, &Ekin);
        forces(v_positions, v_forces, &Epot, &Pres, &Temp, Rho, cell_V, cell_L);

        for (i = 1; i < TEQ; i++) { // loop de equilibracion

            velocity_verlet(v_positions, v_velocities, v_forces, &Epot, &Ekin, &Pres, &Temp, Rho, cell_V, cell_L);

            sf = sqrt(T0 / Temp);
            for (int k = 0; k < N; k++) { // reescaleo de velocidades
                v_velocities->x[k] *= sf;
                v_velocities->y[k] *= sf;
                v_velocities->z[k] *= sf;
                v_velocities->w[k] *= sf;
            }
        }

        int mes = 0;
        float epotm = 0.0, presm = 0.0;
        for (i = TEQ; i < TRUN; i++) { // loop de medicion

            velocity_verlet(v_positions, v_velocities, v_forces, &Epot, &Ekin, &Pres, &Temp, Rho, cell_V, cell_L);

            sf = sqrt(T0 / Temp);
            for (int k = 0; k < N; k++) { // reescaleo de velocidades
                v_velocities->x[k] *= sf;
                v_velocities->y[k] *= sf;
                v_velocities->z[k] *= sf;
                v_velocities->w[k] *= sf;
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
                    //fprintf(file_xyz, "Ar %e %e %e\n", v_positions->x[k + 0], v_positions->y[k + 1], v_positions->z[k + 2]);
                }
            }

            t += DT;
        }
        printf("%f\t%f\t%f\t%f\n", Rho, cell_V, epotm / (float)mes, presm / (float)mes);
    }

    double elapsed = wtime() - start;
    printf("# Tiempo total de simulación = %f segundos\n", elapsed);
    printf("# Tiempo simulado = %f [fs]\n", t * 1.6);
    printf("# ns/day = %f\n", (1.6e-6 * t) / elapsed * 86400);
    //                       ^1.6 fs -> ns       ^sec -> day
    return 0;
}
