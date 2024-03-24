#ifndef MISC_HPP
#define MISC_HPP

#include <iostream>

// auxiliary macros to access the fields
#define U(i, j) u[(i) * prm.NY + (j)]
#define V(i, j) v[(i) * prm.NY + (j)]
#define W(i, j) w[(i) * prm.NY + (j)]
#define ADV_U(i, j) adv_u[(i) * prm.NY + (j)]
#define ADV_V(i, j) adv_v[(i) * prm.NY + (j)]
#define Ustar(i, j) ustar[(i) * prm.NY + (j)]
#define Vstar(i, j) vstar[(i) * prm.NY + (j)]
#define P(i, j) p[(i) * prm.NY + (j)]
#define T(i, j) T[(i) * prm.NY + (j)]
#define S(i, j) S[(i) * prm.NY + (j)]
#define DIV(i, j) div((i) * prm.ny + (j))
#define x(i) (-0.5 * prm.dx + (i) * prm.dx)
#define y(j) (-0.5 * prm.dy + (j) * prm.dy)

typedef struct prm {
  int NX;      // number of points in the x direction (including the ghost points)
  int NY;      // number of points in the y direction (including the ghost points)
  int nx;      // number of points in the x direction (excluding the ghost points)
  int ny;      // number of points in the y direction (excluding the ghost points)
  double L;    // length of the domain in the x direction
  double H;    // length of the domain in the y direction
  double dt;   // time step
  double dx;   // grid spacing in the x direction
  double dy;   // grid spacing in the y direction
  double T;    // final time
  double A_T;  // amplitude of the flow for temperature
  double A_S;  // amplitude of the flow for salinity
  // dimensionless numbers
  double Pr;     // Prandtl number
  double Sc;     // Schmidt number
  double Le;     // Lewis number
  double Ra_T;   // Rayleigh number for temperature
  double Ra_S;   // Rayleigh number for salinity
  double R_rho;  // density stability ratio
  size_t NXNY;   // total number of points in the domain (including the ghost points)
  size_t nxny;   // total number of points in the domain (excluding the ghost points)
} Prm;

#endif  // MISC_HPP