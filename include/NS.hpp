#ifndef NS_HPP
#define NS_HPP

#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/Sparse>

#include "misc.hpp"

// auxiliary macros to access the fields
#define p(jj) (obstacle.InterpolatingPoints[6 * k + (jj)])
#define m_x (obstacle.MirrorPoints[k].x)
#define m_y (obstacle.MirrorPoints[k].y)

typedef Eigen::SparseMatrix<double> SpMat;  // declares a column-major sparse matrix type of double
typedef Eigen::Triplet<double> Trip;        // declares a triplet type of double

using namespace std;

// @brief 1st order semi-Lagrangian for advection
// @param u velocity in x direction
// @param v velocity in y direction
// @param q scalar field to be advected
// @param prm parameters of the simulation (dx, dy, dt, etc.)
// @param sign sign of the direction of the velocity fields (to be easily inverted). If unsure, use 1.
// @param obstacle object representing the obstacle in the domain
void Semilag(double* u, double* v, double* q, Prm prm, int sign);

// @brief 2nd order semi-Lagrangian for advection
// @param u velocity in x direction
// @param v velocity in y direction
// @param q0 initial scalar field to be advected
// @param q1 final scalar field to be advected
// @param prm parameters of the simulation (dx, dy, dt, etc.)
// @param obstacle object representing the obstacle in the domain
void Semilag2(double* u, double* v, double* q0, double* q1, Prm prm);

// @brief Set the boundary conditions for the velocity fields, but on the ghost points (i.e. fictitious points outside the domain) by linear interpolation
// @param u velocity in x direction
// @param v velocity in y direction
// @param prm parameters of the simulation (dx, dy, dt, etc.)
// @param obstacle object representing the obstacle in the domain
void BC_velocity(double* u, double* v, Prm prm);

// @brief Set the boundary conditions for the pressure field, but on the ghost points (i.e. fictitious points outside the domain) by linear interpolation
// @param p pressure
// @param prm parameters of the simulation (dx, dy, dt, etc.)
// @param obstacle object representing the obstacle in the domain
void BC_pressure(double* p, Prm prm);

// @brief Set the boundary conditions for the temperature field, but on the ghost points (i.e. fictitious points outside the domain) by linear interpolation
// @param T temperature
// @param prm parameters of the simulation (dx, dy, dt, etc.)
// @param obstacle object representing the obstacle in the domain
void BC_temperature(double* T, Prm prm);

// @brief Set the boundary conditions for the salinity field, but on the ghost points (i.e. fictitious points outside the domain) by linear interpolation
// @param S salinity
// @param prm parameters of the simulation (dx, dy, dt, etc.)
// @param obstacle object representing the obstacle in the domain
void BC_salinity(double* S, Prm prm);

// @brief Returns the value of the income or outflow of temperature at the upper boundary
// @param x x-coordinate of the point in the upper boundary
// @param prm parameters of the simulation (dx, dy, dt, etc.)
// @return value of the income or outflow of temperature at the upper boundary
double flux_T(double x, Prm prm);

// @brief Returns the value of the income or outflow of salinity at the upper boundary
// @param x x-coordinate of the point in the upper boundary
// @param prm parameters of the simulation (dx, dy, dt, etc.)
// @return value of the income or outflow of salinity at the upper boundary
double flux_S(double x, Prm prm);

// @brief Build the Poisson matrix for the pressure Poisson equation poisson(p) = div(u)
// @param coeffs vector of triplets to store the coefficients of the Possion matrix. Triplets consist of (i, j, value), where i and j are the indices of the matrix and value is the value of the coefficient
// @param prm parameters of the simulation (dx, dy, dt, etc.)
void buildPoissonMatrix(vector<Trip>& coeffs, Prm prm);
void buildPoissonMatrix_V(vector<Trip>& coeffs, double lambda, double mu, Prm prm);
void buildPoissonMatrix_UTS(vector<Trip>& coeffs, double lambda, double mu, Prm prm);
void crankNicholson(double* T, double* adv_T, double* S, double* adv_S, Prm prm);
#endif  // NS_HPP