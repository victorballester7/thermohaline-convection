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
    if (prm.periodic_x) {  // periodic BC
      // left boundary
      U(0, j) = U(prm.NX - 2, j);
      V(0, j) = V(prm.NX - 2, j);
      // right boundary
      U(prm.NX - 1, j) = U(1, j);
      V(prm.NX - 1, j) = V(1, j);
    } else {
      // left boundary: \partial_x u = 0, \partial_x v = 0
      U(0, j) = U(1, j);
      V(0, j) = V(1, j);

      // right boundary: \partial_x u = 0, \partial_x v = 0
      U(prm.NX - 1, j) = U(prm.NX - 2, j);
      V(prm.NX - 1, j) = V(prm.NX - 2, j);
    }
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
    if (prm.periodic_x) {  // periodic BC
      // left boundary
      P(0, j) = P(prm.NX - 2, j);
      // right boundary
      P(prm.NX - 1, j) = P(1, j);
    } else {
      // left boundary: \partial_x p = 0
      P(0, j) = P(1, j);
      // right boundary: \partial_x p = 0
      P(prm.NX - 1, j) = P(prm.NX - 2, j);
      // CHANGE ON POINT TO DIRICHLET
      P(1, 0) = -P(1, 1);
    }
  }
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
    if (prm.periodic_x) {  // periodic BC
      // left boundary
      T(0, j) = T(prm.NX - 2, j);
      // right boundary
      T(prm.NX - 1, j) = T(1, j);
    } else {
      // left boundary: \partial_x T = 0
      T(0, j) = T(1, j);
      // right boundary: \partial_x T = 0
      T(prm.NX - 1, j) = T(prm.NX - 2, j);
    }
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
    if (prm.periodic_x) {  // periodic BC
      // left boundary
      S(0, j) = S(prm.NX - 2, j);
      // right boundary
      S(prm.NX - 1, j) = S(1, j);
    } else {
      // left boundary: \partial_x T = 0
      S(0, j) = S(1, j);
      // right boundary: \partial_x T = 0
      S(prm.NX - 1, j) = S(prm.NX - 2, j);
    }
  }
}

double flux_T(double x, Prm prm) {
  double aux = 2 * M_PI / prm.L;
  return -aux * prm.A_T * cos(aux * x) / tanh(aux * prm.H);
}

double flux_S(double x, Prm prm) {
  double aux = 2 * M_PI / prm.L;
  return -aux * prm.A_S * cos(aux * x) / tanh(aux * prm.H);
}

void buildPoissonMatrix(vector<Trip>& coeffs, Prm prm) {
  // Remeber we are taking the following BC for the pressure:
  // on the top, bottom the normal derivative of the pressure is 0, i.e.:  P(i, NY - 1) = P(i, NY - 2)
  //                                                                       P(i, 0) = P(i, 1);
  // On the right and left boundary, the pressure is peridic, i.e.: P(0, j) = P(NX - 2, j)
  //                                                                P(NX - 1, j) = P(1, j)
  // Also, note that the coefficients are ordered by column, i.e. the first NY elements are the first column, the next NY elements are the second column, and so on.
  // The matrix has size (nx * ny) x (nx * ny), where nx = NX - 2 and ny = NY - 2 (i.e. the number of points in the domain minus the ghost points).
  // Matrix for the 2nd derivative in x (except for the division by dx^2)
  //  ------ ny ------
  // ||-----------------------------------------------------------------------------------------||
  // || -2           |  1           |                                            | 1            || |
  // ||   -2         |    1         |                                            |   1          || |
  // ||      .       |      .       |                                            |     .        || ny
  // ||        .     |        .     |                                            |       .      || |
  // ||          .   |          .   |                                            |         .    || |
  // ||           -2 |            1 |                                            |           1  || |
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
  // || 1            |                                            | 1            | -2           ||
  // ||   1          |                                            |   1          |   -2         ||
  // ||     .        |                                            |      .       |      .       ||
  // ||       .      |                                            |        .     |        .     ||
  // ||         .    |                                            |          .   |          .   ||
  // ||           1  |                                            |            1 |           -2 ||
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

  coeffs.reserve((uint)(5 * prm.nx * prm.ny - 2 * (prm.ny + prm.nx) + 2 * prm.ny));
  int dim = prm.nx * prm.ny;
  double dx_2 = 1. / (prm.dx * prm.dx);
  double dy_2 = 1. / (prm.dy * prm.dy);
  double diagX, diagY;

  for (int i = 0; i < dim; i++) {
    // diagonal
    diagX = -2. * dx_2;
    diagY = -2. * dy_2;
    if (i % prm.ny == 0 || i % prm.ny == prm.ny - 1) diagY = -dy_2;
    if (!prm.periodic_x) {
      if (i < prm.ny || i >= dim - prm.ny) diagX = -dx_2;
    }

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
  if (prm.periodic_x) {  // periodic BC
    for (int i = 0; i < prm.ny; i++) {
      coeffs.push_back(Trip(i + prm.ny * (prm.nx - 1), i, -dx_2));  // top row
      coeffs.push_back(Trip(i, i + prm.ny * (prm.nx - 1), -dx_2));  // bottom row
    }
  } else {
    coeffs.push_back(Trip(0, 0, 2 * dy_2));  // we impose Dirichlet BC on the first element (bottom left corner) Thus, we change -1 | 2 -1 = 0 | 1 -1 (in Neumann) to -1 | 2 -1 = 0 | 2 -1 (in Dirichlet)
  }
}
void buildPoissonMatrix_V(vector<Trip>& coeffs, double lambda, double mu, Prm prm) {
  coeffs.reserve((uint)(5 * prm.nx * prm.ny - 2 * (prm.ny + prm.nx) + 2 * prm.ny));
  int dim = prm.nx * prm.ny;
  double diagX, diagY;

  cout << "lambda: " << lambda << " mu: " << mu << endl;

  for (int i = 0; i < dim; i++) {
    // diagonal
    diagX = 2 * lambda;
    diagY = 2 * mu;
    if (i % prm.ny == 0 || i % prm.ny == prm.ny - 1) diagY = 3 * mu;
    if (!prm.periodic_x) {
      if (i < prm.ny || i >= dim - prm.ny) diagX = lambda;
    }

    coeffs.push_back(Trip(i, i, 1 + diagX + diagY));

    // Dxx part (secondary diagonals)
    if (i < dim - prm.ny) {
      coeffs.push_back(Trip(i, i + prm.ny, -lambda));  // upper diagonal
      coeffs.push_back(Trip(i + prm.ny, i, -lambda));  // lower diagonal
    }

    // Dyy part (secondary diagonals)
    if ((i + 1) % prm.ny != 0) {
      coeffs.push_back(Trip(i, i + 1, -mu));  // upper diagonal
      coeffs.push_back(Trip(i + 1, i, -mu));  // lower diagonal
    }
  }
  if (prm.periodic_x) {  // periodic BC
    for (int i = 0; i < prm.ny; i++) {
      coeffs.push_back(Trip(i + prm.ny * (prm.nx - 1), i, -lambda));  // top row
      coeffs.push_back(Trip(i, i + prm.ny * (prm.nx - 1), -lambda));  // bottom row
    }
  }
}

void buildPoissonMatrix_UTS(vector<Trip>& coeffs, double lambda, double mu, Prm prm) {
  // Remeber we are taking the following BC for the pressure:
  // on the top, bottom the normal derivative of the pressure is 0, i.e.:  P(i, NY - 1) = P(i, NY - 2)
  //                                                                   m    P(i, 0) = P(i, 1);
  // On the right and left boundary, the pressure is peridic, i.e.: P(0, j) = P(NX - 2, j)
  //                                                                P(NX - 1, j) = P(1, j)
  // Also, note that the coefficients are ordered by column, i.e. the first NY elements are the first column, the next NY elements are the second column, and so on.
  // The matrix has size (nx * ny) x (nx * ny), where nx = NX - 2 and ny = NY - 2 (i.e. the number of points in the domain minus the ghost points).
  // Matrix for the 2nd derivative in x (except for the division by dx^2)
  //  ------ ny ------
  // ||-----------------------------------------------------------------------------------------||
  // || -2           |  1           |                                            | 1            || |
  // ||   -2         |    1         |                                            |   1          || |
  // ||      .       |      .       |                                            |     .        || ny
  // ||        .     |        .     |                                            |       .      || |
  // ||          .   |          .   |                                            |         .    || |
  // ||           -2 |            1 |                                            |           1  || |
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
  // || 1            |                                            | 1            | -2           ||
  // ||   1          |                                            |   1          |   -2         ||
  // ||     .        |                                            |      .       |      .       ||
  // ||       .      |                                            |        .     |        .     ||
  // ||         .    |                                            |          .   |          .   ||
  // ||           1  |                                            |            1 |           -2 ||
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

  coeffs.reserve((uint)(5 * prm.nx * prm.ny - 2 * (prm.ny + prm.nx) + 2 * prm.ny));
  int dim = prm.nx * prm.ny;
  double diagX, diagY;

  cout << "lambda: " << lambda << " mu: " << mu << endl;

  for (int i = 0; i < dim; i++) {
    // diagonal
    diagX = 2 * lambda;
    diagY = 2 * mu;
    if (i % prm.ny == 0) diagY = 3 * mu;
    if (i % prm.ny == prm.ny - 1) diagY = mu;
    if (!prm.periodic_x) {
      if (i < prm.ny || i >= dim - prm.ny) diagX = lambda;
    }

    coeffs.push_back(Trip(i, i, 1 + diagX + diagY));

    // Dxx part (secondary diagonals)
    if (i < dim - prm.ny) {
      coeffs.push_back(Trip(i, i + prm.ny, -lambda));  // upper diagonal
      coeffs.push_back(Trip(i + prm.ny, i, -lambda));  // lower diagonal
    }

    // Dyy part (secondary diagonals)
    if ((i + 1) % prm.ny != 0) {
      coeffs.push_back(Trip(i, i + 1, -mu));  // upper diagonal
      coeffs.push_back(Trip(i + 1, i, -mu));  // lower diagonal
    }
  }
  if (prm.periodic_x) {  // periodic BC
    for (int i = 0; i < prm.ny; i++) {
      coeffs.push_back(Trip(i + prm.ny * (prm.nx - 1), i, -lambda));  // top row
      coeffs.push_back(Trip(i, i + prm.ny * (prm.nx - 1), -lambda));  // bottom row
    }
  }
}

void crankNicholson(double* T, double* adv_T, double* S, double* adv_S, Prm prm) {
  const double w = 1.18;  //
  const double max_iter = 20;
  double errorLoo = TOL + 1;
  double mu_x = prm.dt / (prm.dx * prm.dx);
  double mu_y = prm.dt / (prm.dy * prm.dy);
  double tmp;

  // we will store the each new iterate in T, so we start with initial guess the just advected field
  mempcpy(T, adv_T, prm.NXNY * sizeof(double));
  BC_temperature(T, prm);

  mempcpy(S, adv_S, prm.NXNY * sizeof(double));
  BC_salinity(S, prm);

  // SOR method
  for (int m = 1; m < max_iter && errorLoo > TOL; m++) {
    errorLoo = 0;
    for (int i = 1; i < prm.NX - 1; i++) {
      for (int j = 1; j < prm.NY - 1; j++) {
        tmp = T(i, j);
        T(i, j) = (1 - w) * T(i, j) + w / (1 + mu_x + mu_y) * ((1 - mu_x - mu_y) * ADV_T(i, j) + mu_x / 2 * (T(i + 1, j) + T(i - 1, j) + ADV_T(i + 1, j) + ADV_T(i - 1, j)) + mu_y / 2 * (T(i, j + 1) + T(i, j - 1) + ADV_T(i, j + 1) + ADV_T(i, j - 1)));
        errorLoo = fmax(errorLoo, fabs(tmp - T(i, j)));
      }
    }
    // printf("Error: %f\n", errorLoo);
    BC_temperature(T, prm);
  }

  if (errorLoo > TOL) {
    std::cerr << "Temperature did not converge in " << max_iter << " iterations. Error: " << errorLoo << std::endl;
    exit(1);
  }

  errorLoo = TOL + 1;
  for (int m = 1; m < max_iter && errorLoo > TOL; m++) {
    errorLoo = 0;
    for (int i = 1; i < prm.NX - 1; i++) {
      for (int j = 1; j < prm.NY - 1; j++) {
        tmp = (1 - w) * S(i, j) + w / (1 + mu_x + mu_y) * ((1 - mu_x - mu_y) * ADV_S(i, j) + mu_x / 2 * (S(i + 1, j) + S(i - 1, j) + ADV_S(i + 1, j) + ADV_S(i - 1, j)) + mu_y / 2 * (S(i, j + 1) + S(i, j - 1) + ADV_S(i, j + 1) + ADV_S(i, j - 1)));
        errorLoo = fmax(errorLoo, fabs(tmp - S(i, j)));
        S(i, j) = tmp;
      }
    }
    BC_salinity(S, prm);
  }

  if (errorLoo > TOL) {
    std::cerr << "Salinity did not converge in " << max_iter << " iterations. Error: " << errorLoo << std::endl;
    exit(1);
  }
}