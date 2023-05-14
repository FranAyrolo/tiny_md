#include "core.h"
#include "parameters.h"

#include <math.h>
#include <stdlib.h> // rand()

#define ECUT (4.0 * (pow(RCUT, -12) - pow(RCUT, -6)))

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

/*
static uint32_t s[4];

static void xoshiro_seed(uint32_t seed) {
    s[0] = seed;
    s[1] = seed + 1;
    s[2] = seed + 2;
    s[3] = seed + 3;
}

static __m128i xoshiro_next() {
    const __m128i s0 = _mm_set_epi32(s[3], s[2], s[1], s[0]);

    const __m128i s1 = _mm_xor_si128(s0, _mm_slli_epi32(s0, 5));
    const __m128i s2 = _mm_xor_si128(s1, _mm_slli_epi32(s1, 15));
    const __m128i s3 = _mm_xor_si128(s2, _mm_slli_epi32(s2, 27));

    const __m128i result = _mm_add_epi32(s0, s3);

    // Update the state
    s[0] = _mm_extract_epi32(s1, 3);
    s[1] = _mm_extract_epi32(s2, 3);
    s[2] = _mm_extract_epi32(s3, 3);
    s[3] = _mm_extract_epi32(result, 3);

    return result;
}

static void generate_random_numbers(int* array, size_t len) {
    const size_t num_iterations = len / 4;
    for(size_t i = 0; i < num_iterations; i++) {
        __m128i random = xoshiro_next();
        __m128i scaled_random = _mm_cvtps_epi32(_mm_mul_ps(_mm_cvtepi32_ps(random), _mm_set1_ps(RAND_MAX)));
        _mm_storeu_si128((__m128i*)(array + 4 * i), scaled_random);
    }

    // If len is not divisible by 4, generate the remaining random numbers
    const size_t remaining = len % 4;
    if(remaining > 0) {
        __m128i random = xoshiro_next();
        __m128i scaled_random = _mm_cvtps_epi32(_mm_mul_ps(_mm_cvtepi32_ps(random), _mm_set1_ps(RAND_MAX)));
        int* last_values = (int*)&scaled_random;
        for(size_t i = 0; i < remaining; i++) {
            array[len - remaining + i] = last_values[i];
        }
    }
}
*/


void init_vel(Vector_SOA* restrict v_velocities, float* restrict temp, float* restrict ekin)
{
    // inicialización de velocidades aleatorias

    float sf, sumvx = 0.0, sumvy = 0.0, sumvz = 0.0, sumv2 = 0.0;

    // int random_array[3*N];
    // xoshiro_seed(1234);
    //generate_random_numbers(random_array, 3*N);
    //int k = 0;

    for (int i = 0; i < N; i++) {
        v_velocities->x[i] = rand() / (float)RAND_MAX - 0.5;
        v_velocities->y[i] = rand() / (float)RAND_MAX - 0.5;
        v_velocities->z[i] = rand() / (float)RAND_MAX - 0.5;
        //k += 3;

        sumvx += v_velocities->x[i];
        sumvy += v_velocities->y[i];
        sumvz += v_velocities->z[i];
        sumv2 += v_velocities->x[i] * v_velocities->x[i] + v_velocities->y[i] * v_velocities->y[i]
               + v_velocities->z[i] * v_velocities->z[i];
    }

    sumvx /= (float)N;
    sumvy /= (float)N;
    sumvz /= (float)N;

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

inline static float minimum_image(float cordi, const float cell_length)
{
    const float half_cell_length = 0.5 * cell_length;
    cordi -= cell_length * (cordi > half_cell_length);
    cordi += cell_length * (cordi <= -half_cell_length);
    return cordi;
}


void forces(Vector_SOA* restrict v_positions, Vector_SOA* restrict v_forces, float* restrict epot, float* restrict pres,
            const float* restrict temp, const float rho, const float V, const float L)
{
    // calcula las fuerzas LJ (12-6)
    for (int i = 0; i < N; i++) {
        v_forces->x[i] = 0.0;
        v_forces->y[i] = 0.0;
        v_forces->z[i] = 0.0;
    }
    float pres_vir = 0.0;
    float rcut2 = RCUT * RCUT;
    *epot = 0.0;

    for (int i = 0; i < N - 1; i++) {

        float xi = v_positions->x[i];
        float yi = v_positions->y[i];
        float zi = v_positions->z[i];

        for (int j = i + 1; j < N; j++) {

            float xj = v_positions->x[j];
            float yj = v_positions->y[j];
            float zj = v_positions->z[j];

            // distancia mínima entre r_i y r_j
            float rx = xi - xj;
            rx = minimum_image(rx, L);
            float ry = yi - yj;
            ry = minimum_image(ry, L);
            float rz = zi - zj;
            rz = minimum_image(rz, L);        

            float rij2 = rx * rx + ry * ry + rz * rz;

            float condition = (rij2 <= rcut2) ? 1.0 : 0.0;
            float r2inv = condition * (1.0 / rij2);
            float r6inv = condition * r2inv * r2inv * r2inv;

            float fr = condition * 24.0 * r2inv * r6inv * (2.0 * r6inv - 1.0);

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


inline static float pbc(float cordi, const float cell_length)
{
    cordi += cell_length * (cordi <= 0);
    cordi -= cell_length * (cordi > cell_length && cordi > 0);
    return cordi;
}



void velocity_verlet(Vector_SOA* restrict v_positions, Vector_SOA* restrict v_velocities, Vector_SOA* restrict v_forces, float* restrict epot,
                     float* restrict ekin, float* restrict pres, float* restrict temp, const float rho,
                     const float V, const float L)
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

    float sumv2 = 0.0;
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
