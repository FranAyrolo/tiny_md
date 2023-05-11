#ifndef CORE_H
#define CORE_H

struct Atoms {
    double *x;
    double *y;
    double *z;
};

void init_pos(struct Atoms* positions, const double rho);
void init_vel(struct Atoms* velocities, double* temp, double* ekin);
void forces(const struct Atoms* positions, struct Atoms* forcess, double* epot, double* pres,
            const double* temp, const double rho, const double V, const double L);
void velocity_verlet(struct Atoms* positions, struct Atoms* velocities, struct Atoms* forcess, double* epot,
                     double* ekin, double* pres, double* temp, const double rho,
                     const double V, const double L);

#endif
