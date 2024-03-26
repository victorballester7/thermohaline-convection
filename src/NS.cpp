#include "../include/NS.hpp"

#include <cstring>

void Semilag(double* u, double* v, double* q, Prm prm, int sign) {
  int sign_u, sign_v;
  double a, b;
  double* aux = (double*)calloc(prm.NXNY, sizeof(double));
  for (int i = 1; i < prm.NX - 1; i++) {
    for (int j = 1; j < prm.NY - 1; j++) {
      if (sign * U(i, j) > 0) {
        a = 1 - sign * U(i, j) * prm.dt / prm.dx;
        sign_u = 1;
      } else {
        a = 1 + sign * U(i, j) * prm.dt / prm.dx;
        sign_u = -1;
      }
      if (sign * V(i, j) > 0) {
        b = 1 - sign * V(i, j) * prm.dt / prm.dy;
        sign_v = 1;
      } else {
        b = 1 + sign * V(i, j) * prm.dt / prm.dy;
        sign_v = -1;
      }

      aux[i * prm.NY + j] = a * b * q[i * prm.NY + j] +
                            (1 - a) * b * q[(i - sign_u) * prm.NY + j] +
                            a * (1 - b) * q[i * prm.NY + j - sign_v] +
                            (1 - a) * (1 - b) * q[(i - sign_u) * prm.NY + j - sign_v];
    }
  }
  memcpy(q, aux, prm.NXNY * sizeof(double));
  free(aux);
}

void Semilag2(double* u, double* v, double* q0, double* q1, Prm prm) {
  memcpy(q1, q0, prm.NXNY * sizeof(double));
  Semilag(u, v, q1, prm, 1);
  Semilag(u, v, q1, prm, -1);
  for (int i = 0; i < prm.NX * prm.NY; i++)
    q1[i] = q0[i] + (q0[i] - q1[i]) / 2;
  Semilag(u, v, q1, prm, 1);
}

void BC_velocity(double* u, double* v, Prm prm) {
  for (int i = 0; i < prm.NX; i++) {
    // bottom boundary: u = 0, v = 0
    U(i, 0) = -U(i, 1);
    V(i, 0) = -V(i, 1);
    // top boundary: \partial_y u = 0, v = 0
    U(i, prm.NY - 1) = U(i, prm.NY - 2);
    V(i, prm.NY - 1) = -V(i, prm.NY - 2);
  }

  for (int j = 0; j < prm.NY; j++) {
    // left boundary: \partial_x u = 0, \partial_x v = 0
    U(0, j) = U(1, j);
    V(0, j) = V(1, j);
    // right boundary: \partial_x u = 0, \partial_x v = 0
    U(prm.NX - 1, j) = U(prm.NX - 2, j);
    V(prm.NX - 1, j) = V(prm.NX - 2, j);
  }
}

void BC_pressure(double* p, Prm prm) {
  for (int i = 0; i < prm.NX; i++) {
    // bottom boundary: \partial_y p = 0
    P(i, 0) = P(i, 1);
    // top boundary: \partial_y p = 0
    P(i, prm.NY - 1) = P(i, prm.NY - 2);
  }
  for (int j = 0; j < prm.NY; j++) {
    // left boundary: \partial_x p = 0
    P(0, j) = P(1, j);
    // right boundary: \partial_x p = 0
    P(prm.NX - 1, j) = P(prm.NX - 2, j);
  }
  // CHANGE ON POINT TO DIRICHLET
  // P(1, 0) = -P(1, 1);
}

void BC_temperature(double* T, Prm prm) {
  for (int i = 0; i < prm.NX; i++) {
    // bottom boundary: T = 0 (once normalized)
    T(i, 0) = -T(i, 1);
    // top boundary: \partial_y T = flow_T
    T(i, prm.NY - 1) = T(i, prm.NY - 2) + prm.dy * flux_T(x(i), prm);
    // T(i, prm.NY - 1) = -T(i, prm.NY - 2);
  }
  for (int j = 0; j < prm.NY; j++) {
    // left boundary: \partial_x T = 0
    T(0, j) = T(1, j);
    // right boundary: \partial_x T = 0
    T(prm.NX - 1, j) = T(prm.NX - 2, j);
  }
}

void BC_salinity(double* S, Prm prm) {
  for (int i = 0; i < prm.NX; i++) {
    // bottom boundary: S = 0 (once normalized)
    S(i, 0) = -S(i, 1);
    // top boundary: \partial_y S = flow_S
    S(i, prm.NY - 1) = S(i, prm.NY - 2) + prm.dy * flux_S(x(i), prm);
  }
  for (int j = 0; j < prm.NY; j++) {
    // left boundary: \partial_x S = 0
    S(0, j) = S(1, j);
    // right boundary: \partial_x S = 0
    S(prm.NX - 1, j) = S(prm.NX - 2, j);
  }
}

double flux_T(double x, Prm prm) {
  return -prm.A_T * cos(2 * M_PI * x / prm.L);
}

double flux_S(double x, Prm prm) {
  return -prm.A_S * cos(2 * M_PI * x / prm.L);
}

void buildPoissonMatrix(vector<Trip>& coeffs, Prm prm) {
  // Remeber we are taking the following BC for the pressure:
  // on the top, bottom and left boundaries, the normal derivative of the pressure is 0, i.e.:  P(i, NY - 1) = P(i, NY - 2)
  //                                                                                            P(i, 0) = P(i, 1);
  //                                                                                            P(0, j) = P(1, j);
  // respectively. On the right boundary, the pressure is 0, i.e.: P(NX - 1, j) = P(NX - 2, j).
  // Also, note that the coefficients are ordered by column, i.e. the first NY elements are the first column, the next NY elements are the second column, and so on.
  // The matrix has size (nx * ny) x (nx * ny), where nx = NX - 2 and ny = NY - 2 (i.e. the number of points in the domain minus the ghost points).
  // Matrix for the 2nd derivative in x (except for the division by dx^2)
  // NOTE: in order for the system to have a unique solution, we need to impose a Dirichlet BC on the pressure on at least one point. We will impose it on the point (0, 0) in the y-direction.
  //  ------ ny ------
  // ||-----------------------------------------------------------------------------------------||
  // || -1           |  1           |                                                           || |
  // ||   -1         |    1         |                                                           || |
  // ||      .       |      .       |                                                           || ny
  // ||        .     |        .     |                                                           || |
  // ||          .   |          .   |                                                           || |
  // ||           -1 |            1 |                                                           || |
  // ||-----------------------------------------------------------------------------------------||
  // || 1            | -2           | 1            |                                            ||
  // ||   1          |   -2         |   1          |                                            ||
  // ||      .       |      .       |      .       |                                            ||
  // ||        .     |        .     |        .     |                                            ||
  // ||          .   |          .   |          .   |                                            ||
  // ||            1 |           -2 |            1 |                                            ||
  // ||-----------------------------------------------------------------------------------------||
  // ||              | 1            | -2           | 1            |                             ||
  // ||              |   1          |   -2         |   1          |                             ||
  // ||              |      .       |      .       |      .       |                             ||
  // ||              |        .     |        .     |        .     |                             ||
  // ||              |          .   |          .   |          .   |                             ||
  // ||              |            1 |           -2 |            1 |                             ||
  // ||-----------------------------------------------------------------------------------------||
  // ||                             |  .           |  .           |  .           |              ||
  // ||                             |    .         |    .         |    .         |              ||
  // ||                             |      .       |      .       |      .       |              ||
  // ||                             |        .     |        .     |        .     |              ||
  // ||                             |          .   |          .   |          .   |              ||
  // ||                             |            . |            . |            . |              ||
  // ||-----------------------------------------------------------------------------------------||
  // ||                                            | 1            | -2           | 1            ||
  // ||                                            |   1          |   -2         |   1          ||
  // ||                                            |      .       |      .       |      .       ||
  // ||                                            |        .     |        .     |        .     ||
  // ||                                            |          .   |          .   |          .   ||
  // ||                                            |            1 |           -2 |            1 ||
  // ||-----------------------------------------------------------------------------------------||
  // ||                                                           | 1            | -1           ||
  // ||                                                           |   1          |   -1         ||
  // ||                                                           |      .       |      .       ||
  // ||                                                           |        .     |        .     ||
  // ||                                                           |          .   |          .   ||
  // ||                                                           |            1 |           -1 ||
  // ||-----------------------------------------------------------------------------------------||

  // Matrix for the 2nd derivative in y (except for the division by dy^2)
  //  --------- ny --------
  // ||------------------------------------------------------------------------------------||
  // || -1  1              |                                                               || |
  // ||  1 -2  1           |                                                               || |
  // ||    1 -2  1         |                                                               || |
  // ||      .  .  .       |                                                               || ny
  // ||        .  .  .     |                                                               || |
  // ||          .  .  .   |                                                               || |
  // ||            1 -2  1 |                                                               || |
  // ||               1 -1 |                                                               || |
  // ||------------------------------------------------------------------------------------||
  // ||                    | -1  1              |                                          ||
  // ||                    |  1 -2  1           |                                          ||
  // ||                    |    1 -2  1         |                                          ||
  // ||                    |      .  .  .       |                                          ||
  // ||                    |        .  .  .     |                                          ||
  // ||                    |          .  .  .   |                                          ||
  // ||                    |            1 -2  1 |                                          ||
  // ||                    |               1 -1 |                                          ||
  // ||------------------------------------------------------------------------------------||
  // ||                                          |   .                |                    ||
  // ||                                          |     .              |                    ||
  // ||                                          |       .            |                    ||
  // ||                                          |         .          |                    ||
  // ||                                          |           .        |                    ||
  // ||                                          |             .      |                    ||
  // ||                                          |               .    |                    ||
  // ||                                          |                 .  |                    ||
  // ||------------------------------------------------------------------------------------||
  // ||                                                               | -1  1              ||
  // ||                                                               |  1 -2  1           ||
  // ||                                                               |    1 -2  1         ||
  // ||                                                               |      .  .  .       ||
  // ||                                                               |        .  .  .     ||
  // ||                                                               |          .  .  .   ||
  // ||                                                               |            1 -2  1 ||
  // ||                                                               |               1 -1 ||
  // ||------------------------------------------------------------------------------------||

  // number of non-zero coefficients in the Poi matrix
  // dx -> 3 * ny * (nx - 2) + 2 * (2 * ny) = ny * (3 * nx - 2)
  // dy -> (3 * ny - 2) * nx
  // common (diag) -> nx * ny
  // total = dx + dy - common = 5 * nx * ny - 2 * ny - 2 * nx
  // we add one for the Dirichlet BC

  coeffs.reserve((uint)(5 * prm.nx * prm.ny - 2 * (prm.ny + prm.nx) + 2));
  int dim = prm.nx * prm.ny;
  double dx_2 = 1. / (prm.dx * prm.dx);
  double dy_2 = 1. / (prm.dy * prm.dy);
  double diagX, diagY;

  for (int i = 0; i < dim; i++) {
    // diagonal
    diagX = -2. * dx_2;
    diagY = -2. * dy_2;
    if (i % prm.ny == 0 || i % prm.ny == prm.ny - 1) diagY = -dy_2;
    if (i < prm.ny || i >= dim - prm.ny) diagX = -dx_2;
    coeffs.push_back(Trip(i, i, -diagX - diagY));

    // Dxx part (secondary diagonals)
    if (i < dim - prm.ny) {
      coeffs.push_back(Trip(i, i + prm.ny, -dx_2));  // upper diagonal
      coeffs.push_back(Trip(i + prm.ny, i, -dx_2));  // lower diagonal
    }

    // Dyy part (secondary diagonals)
    if ((i + 1) % prm.ny != 0) {
      coeffs.push_back(Trip(i, i + 1, -dy_2));  // upper diagonal
      coeffs.push_back(Trip(i + 1, i, -dy_2));  // lower diagonal
    }
  }
  // coeffs.push_back(Trip(0, 0, -2 * dy_2));
  coeffs.push_back(Trip(0, 1, dy_2));
  coeffs.push_back(Trip(1, 0, dy_2));
}