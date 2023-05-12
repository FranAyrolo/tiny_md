#ifndef CORE_H
#define CORE_H

typedef struct {
    double x;
    double y;
    double z;
} Vector3D;

typedef struct {
    Vector3D* positions;
    Vector3D* velocities;
    Vector3D* forces;
} AtomSystem;

AtomSystem* createAtomSystem(int numAtoms);

void destroyAtomSystem(AtomSystem* system);

void init_pos(AtomSystem* system, const double rho);
void init_vel(AtomSystem* system, double* temp, double* ekin);
void forces(const AtomSystem* system, double* epot, double* pres,
            const double* temp, const double rho, const double V, const double L);
void velocity_verlet(AtomSystem* system, double* epot,
                     double* ekin, double* pres, double* temp, const double rho,
                     const double V, const double L);

#endif
