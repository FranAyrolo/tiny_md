/*
 * Tiny Molecular Dynamics
 *
 * Unidades: Lennard-Jones
 *
 */

#ifndef PARAMETERS_H
#define PARAMETERS_H

#include <time.h> // time(NULL)

// Valores que podemos usar de N = 256, 500, 864, 1372, 2048, 2916, 4000, 5324, 6912
// 8788 10976 13500 16384 19652 23328 27436 32000
#ifndef N
#define N 256
#endif

#ifndef SEED // rand SEED para las velocidades
#define SEED (time(NULL))
#endif

#ifndef T0 // isoterma
#define T0 2.0
#endif

#ifndef RHOI // densidad inicial
#define RHOI 1.2
#endif

#ifndef RCUT // radio de corte
#define RCUT 2.5
#endif

#ifndef DT // paso temporal ~ 1.6 fs para el Ar
#define DT 0.005
#endif

#ifndef TEQ // pasos de equilibracion
#define TEQ 500
#endif

#ifndef TRUN // TRUN - TEQ: pasos de medicion
#define TRUN 2000
#endif

#ifndef TMES // cada cuantos pasos se mide
#define TMES 10
#endif

#endif
