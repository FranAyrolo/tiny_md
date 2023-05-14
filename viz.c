#define _XOPEN_SOURCE 500  // M_PI
#include "core.h"
#include "parameters.h"

#include <GL/glut.h> // OpenGL
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

// variables globales
static double Ekin, Epot, Temp, Pres; // variables macroscopicas
static double Rho, V, box_size, tail, Etail, Ptail;
static double Rhob, sf, epotm, presm;
static int switcher = 0, frames = 0, mes;
static Vector_SOA* v_positions;
static Vector_SOA* v_velocities;
static Vector_SOA* v_forces;


// OpenGL specific drawing routines
static int win_id;
static int win_x = 900, win_y = 900;


static void pre_display(void)
{ // 3D
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();

    gluPerspective(45.0, (double)win_x / win_y, 1.0, 0.0);
    gluLookAt(1.0, 1.0, 1.0, 0.5, 0.5, 0.5, 1.0, 0.0, 0.0);

    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();

    glClearColor(0.0f, 0.0f, 0.0f, 1.0f);
    glClear(GL_COLOR_BUFFER_BIT);
}


static void post_display(void)
{
    glutSwapBuffers();
}


static void draw_atoms(void)
{
    double glL = cbrt((double)N / (RHOI - 0.8));

    double resize = 0.5;

    // grafico las lineas que delimitan la caja de simulación
    glBegin(GL_LINES);

    double box_line = resize * (box_size / glL);
    glColor3d(0.0, 0.0, 1.0);

    glVertex3d(0.0, 0.0, 0.0);
    glVertex3d(0.0, 0.0, box_line);

    glVertex3d(0.0, 0.0, 0.0);
    glVertex3d(0.0, box_line, 0.0);

    glVertex3d(0.0, 0.0, 0.0);
    glVertex3d(box_line, 0.0, 0.0);

    glVertex3d(box_line, box_line, box_line);
    glVertex3d(box_line, box_line, 0.0);

    glVertex3d(box_line, box_line, box_line);
    glVertex3d(box_line, 0.0, box_line);

    glVertex3d(box_line, box_line, box_line);
    glVertex3d(0.0, box_line, box_line);

    glVertex3d(0.0, box_line, 0.0);
    glVertex3d(box_line, box_line, 0.0);

    glVertex3d(0.0, box_line, box_line);
    glVertex3d(0.0, 0.0, box_line);

    glVertex3d(box_line, 0.0, box_line);
    glVertex3d(box_line, 0.0, 0.0);

    glVertex3d(box_line, 0.0, box_line);
    glVertex3d(0.0, 0.0, box_line);

    glVertex3d(0.0, box_line, box_line);
    glVertex3d(0.0, box_line, 0.0);

    glVertex3d(box_line, box_line, 0.0);
    glVertex3d(box_line, 0.0, 0.0);

    glEnd();

    // grafico las particulas (x, y, z) en el punto (dx, dy, dx), son reescaleadas
    // a [0, 1] y luego multiplicadas con un factor que las achica para poder
    // apreciar mejor el cambio en el volumen
    glBegin(GL_POINTS);

    int di;

    double dx;
    double dy;
    double dz;

    for (di = 0; di < N; di++) {
        dx = (v_positions->x[di] / glL) * resize;
        dy = (v_positions->y[di] / glL) * resize;
        dz = (v_positions->z[di] / glL) * resize;

        glColor3d(0.0, 1.0, 0.0);
        glVertex3d(dx, dy, dz);
    }

    glEnd();
}


/*static void reshape_func ( int width, int height )
{
	glutSetWindow ( win_id );
	glutReshapeWindow ( width, height );

	win_x = width;
	win_y = height;
}*/


static void idle_func(void)
{

    if (switcher == 3) {

        Rho = RHOI;
        V = (double)N / Rho;
        box_size = cbrt(V);
        tail = 16.0 * M_PI * Rho * ((2.0 / 3.0) * pow(RCUT, -9) - pow(RCUT, -3)) / 3.0;
        Etail = tail * (double)N;
        Ptail = tail * Rho;

        init_pos(v_positions, Rho);
        init_vel(v_velocities, &Temp, &Ekin);
        forces(v_positions, v_forces, &Epot, &Pres, &Temp, Rho, V, box_size);

        switcher = 0;

    } else if (switcher == 2) { // imprimo propiedades en la terminal y cambio la densidad

        printf("%f\t%f\t%f\t%f\n", Rho, V, epotm / (double)mes,
               presm / (double)mes);

        Rhob = Rho;
        Rho = Rho - 0.1;


        V = (double)N / Rho;
        box_size = cbrt(V);
        tail = 16.0 * M_PI * Rho * ((2.0 / 3.0) * pow(RCUT, -9) - pow(RCUT, -3)) / 3.0;
        Etail = tail * (double)N;
        Ptail = tail * Rho;

        sf = cbrt(Rhob / Rho);
        for (int k = 0; k < N; k++) { // reescaleo posiciones a nueva densidad
            v_positions->x[k] *= sf;
            v_positions->y[k] *= sf;
            v_positions->z[k] *= sf;
        }
        init_vel(v_velocities, &Temp, &Ekin);
        forces(v_positions, v_forces, &Epot, &Pres, &Temp, Rho, V, box_size);

        switcher = 0;
        if (fabs(Rho - (RHOI - 0.9f)) < 1e-6) {
            printf("\n");
            switcher = 3;
        }

    } else if (switcher == 1) { // loop de medición

        for (int i = frames; i < frames + TMES; i++) {

            velocity_verlet(v_positions, v_velocities, v_forces, &Epot, &Ekin, &Pres, &Temp, Rho,
                            V, box_size);

            sf = sqrt(T0 / Temp);
            for (int k = 0; k < N; k++) { // reescaleo de velocidades
                v_velocities->x[k] *= sf;
                v_velocities->y[k] *= sf;
                v_velocities->z[k] *= sf;
            }
        }

        Epot += Etail;
        Pres += Ptail;

        epotm += Epot;
        presm += Pres;
        mes++;

        frames += TMES;
        if (frames % TRUN == 0) {
            switcher = 2;
        }

    } else if (switcher == 0) { // loop de equilibración

        while (frames % TEQ != 0) {

            velocity_verlet(v_positions, v_velocities, v_forces, &Epot, &Ekin, &Pres, &Temp, Rho,
                            V, box_size);

            sf = sqrt(T0 / Temp);
            for (int k = 0; k < N; k++) { // reescaleo de velocidades
                v_velocities->x[k] *= sf;
                v_velocities->y[k] *= sf;
                v_velocities->z[k] *= sf;
            }

            frames++;
        }

        mes = 0;
        epotm = 0.0;
        presm = 0.0;

        switcher = 1;
    }
    glutSetWindow(win_id);
    glutPostRedisplay();
}


static void display_func(void)
{
    pre_display();
    draw_atoms();
    post_display();
}


static void open_glut_window(void)
{
    glutInitDisplayMode(GLUT_RGBA | GLUT_DOUBLE);

    glutInitWindowPosition(0, 0);
    glutInitWindowSize(win_x, win_y);
    win_id = glutCreateWindow("tiny molecular dynamics | visualization");

    glClearColor(0.0f, 0.0f, 0.0f, 1.0f);
    glClear(GL_COLOR_BUFFER_BIT);
    glutSwapBuffers();
    glClear(GL_COLOR_BUFFER_BIT);
    glutSwapBuffers();

    pre_display();

    // glutKeyboardFunc ( key_func );
    // glutMouseFunc ( mouse_func );
    // glutMotionFunc ( motion_func );
    //glutReshapeFunc ( reshape_func );

    glutIdleFunc(idle_func);
    glutDisplayFunc(display_func);
}


// viz main

int main(int argc, char** argv)
{
    v_positions = (Vector_SOA*)malloc(sizeof(Vector_SOA));
    v_velocities = (Vector_SOA*)malloc(sizeof(Vector_SOA));
    v_forces = (Vector_SOA*)malloc(sizeof(Vector_SOA));

    glutInit(&argc, argv);

    // parametros iniciales para que los pueda usar (antes de modificar)
    // `idle_func`
    srand(SEED);
    Rho = RHOI;
    Rhob = Rho;
    V = (double)N / Rho;
    box_size = cbrt(V);
    tail = 16.0 * M_PI * Rho * ((2.0 / 3.0) * pow(RCUT, -9) - pow(RCUT, -3)) / 3.0;
    Etail = tail * (double)N;
    Ptail = tail * Rho;

    init_pos(v_positions, Rho);
    init_vel(v_velocities, &Temp, &Ekin);
    forces(v_positions, v_forces, &Epot, &Pres, &Temp, Rho, V, box_size);
    //
    //

    printf("# Número de partículas:      %d\n", N);
    printf("# Temperatura de referencia: %.2f\n", T0);
    printf("# Pasos de equilibración:    %d\n", TEQ);
    printf("# Pasos de medición:         %d\n", TRUN - TEQ);
    printf("# (mediciones cada %d pasos)\n", TMES);

    open_glut_window();

    glutMainLoop();

    exit(0);
}
