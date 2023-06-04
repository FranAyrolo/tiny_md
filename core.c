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

// Arreglo para el algoritmo de Xoshiro
static uint32_t s[8];

static void xoshiro_seed(uint32_t seed) {
    s[0] = seed;
    s[1] = seed + 1;
    s[2] = seed + 2;
    s[3] = seed + 3;
    s[4] = seed + 4;
    s[5] = seed + 5;
    s[6] = seed + 6;
    s[7] = seed + 7;
}

// Genera flotantes entre -0.5 y 0.5 en batchs
static __m256 xoshiro_next() {
    const __m256i s0 = _mm256_set_epi32(s[7], s[6], s[5], s[4], s[3], s[2], s[1], s[0]);

    const __m256i s1 = _mm256_xor_si256(s0, _mm256_slli_epi32(s0, 5));
    const __m256i s2 = _mm256_xor_si256(s1, _mm256_slli_epi32(s1, 15));
    const __m256i s3 = _mm256_xor_si256(s2, _mm256_slli_epi32(s2, 27));

    const __m256i result = _mm256_add_epi32(s0, s3);

    // Actualizamos el estado
    s[0] = _mm256_extract_epi32(s1, 7);
    s[1] = _mm256_extract_epi32(s2, 7);
    s[2] = _mm256_extract_epi32(s3, 7);
    s[3] = _mm256_extract_epi32(result, 7);
    s[4] = _mm256_extract_epi32(s1, 6);
    s[5] = _mm256_extract_epi32(s2, 6);
    s[6] = _mm256_extract_epi32(s3, 6);
    s[7] = _mm256_extract_epi32(result, 6);

    const __m256 scale = _mm256_set1_ps(1.0f / (float)UINT32_MAX);
    const __m256 scaled_result = _mm256_mul_ps(_mm256_cvtepi32_ps(result), scale);
    const __m256 scaled_minus_half = _mm256_set1_ps(-0.5f);
    const __m256 final_result = _mm256_add_ps(scaled_result, scaled_minus_half);

    return final_result;
}

static void generate_random_floats(float* array, size_t len) {
    const size_t num_iterations = len / 8;
    for(size_t i = 0; i < num_iterations; i++) {
        __m256 random = xoshiro_next();
        _mm256_storeu_ps(array + 8 * i, random);
    }

    // Si len no es divisible por 8, generamos los restantes numeros aleatorios
    const size_t remaining = len % 8;
    if(remaining > 0) {
        __m256 random = xoshiro_next();
        float* last_values = (float*)&random;
        for(size_t i = 0; i < remaining; i++) {
            array[len - remaining + i] = last_values[i];
        }
    }
}


void init_vel(Vector_SOA* restrict v_velocities, float* restrict temp, float* restrict ekin)
{
    // inicialización de velocidades aleatorias

    float sf, sumvx = 0.0, sumvy = 0.0, sumvz = 0.0, sumv2 = 0.0;
    float random_array[N*3];
    xoshiro_seed(SEED);
    generate_random_floats(random_array, 3*N);
    int k = 0;

    for (int i = 0; i < N; i++) {
        v_velocities->x[i] = random_array[k++];
        v_velocities->y[i] = random_array[k++];
        v_velocities->z[i] = random_array[k++];

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

void forces(float* restrict pos_x, float* restrict pos_y, float* restrict pos_z, 
            float* restrict forces_x, float* restrict forces_y, float* restrict forces_z, 
            float* epot, float* pres, const float* temp, const float rho, const float V, const float L)
{
    // calcula las fuerzas LJ (12-6)
    for (int i = 0; i < N; i++) {
        forces_x[i] = 0.0;
        forces_y[i] = 0.0;
        forces_z[i] = 0.0;
    }
    
    float pres_vir = 0.0;
    float rcut2 = RCUT * RCUT;
    *epot = 0.0;
    float temp_epot = *epot;
    
    int i, j, k;
    float priv_forces_x[N] = {0.0}, priv_forces_y[N] = {0.0}, priv_forces_z[N] = {0.0};
    
    #pragma omp parallel default(private) reduction(+:temp_epot, pres_vir) \
    firstprivate(priv_forces_x, priv_forces_y, priv_forces_z, L, rcut2) private(i, j, k) \
    shared(pos_x, pos_y, pos_z, forces_x, forces_y, forces_z)
    {
    #pragma omp for nowait
    for (k = 0; k < N * (N - 1) / 2; k++) {
        j = (int)((sqrt(8 * k + 1) + 1) / 2);
        i = k - j * (j - 1) / 2;
        
        float xi = pos_x[i];
        float yi = pos_y[i];
        float zi = pos_z[i];
        
        float xj = pos_x[j];
        float yj = pos_y[j];
        float zj = pos_z[j];
        
        // distancia mínima entre r_i y r_j
        float rx = minimum_image(xi - xj, L);
        float ry = minimum_image(yi - yj, L);
        float rz = minimum_image(zi - zj, L);        
        float rij2 = rx * rx + ry * ry + rz * rz;
        float condition = (rij2 <= rcut2) ? 1.0 : 0.0;
        float r2inv = condition * (1.0 / rij2);
        float r6inv = condition * r2inv * r2inv * r2inv;
        float fr = condition * 24.0 * r2inv * r6inv * (2.0 * r6inv - 1.0);
        
        priv_forces_x[i] += condition * fr * rx;
        priv_forces_y[i] += condition * fr * ry;
        priv_forces_z[i] += condition * fr * rz;
        
        priv_forces_x[j] -= condition * fr * rx;
        priv_forces_y[j] -= condition * fr * ry;
        priv_forces_z[j] -= condition * fr * rz;

        temp_epot += condition * (4.0 * r6inv * (r6inv - 1.0) - ECUT);
        pres_vir += condition * fr * rij2;
    }
    for(i = 0; i < N; i++) {
    //#pragma omp atomic
        forces_x[i] += priv_forces_x[i];
        forces_y[i] += priv_forces_y[i];
        forces_z[i] += priv_forces_z[i];
    }
        
    }
    *epot += temp_epot;
    pres_vir /= (V * 3.0);
    *pres = *temp * rho + pres_vir;
}


inline static float pbc(float cordi, const float cell_length)
{
    cordi += cell_length * (cordi <= 0);
    cordi -= cell_length * (cordi > cell_length && cordi > 0);
    return cordi;
}


void velocity_verlet(Vector_SOA* restrict v_positions, Vector_SOA* restrict v_velocities, 
                        Vector_SOA* restrict v_forces, float* restrict epot, float* restrict ekin, 
                        float* restrict pres, float* restrict temp, const float rho, const float V, const float L)
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

    forces(v_positions->x, v_positions->y, v_positions->z, v_forces->x, v_forces->y, v_forces->z, 
            epot, pres, temp, rho, V, L); // actualizo fuerzas

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
