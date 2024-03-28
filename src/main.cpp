#include <chrono>
#include <cstring>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <vector>

// hdf5 library
#include <H5Cpp.h>

#include "../include/NS.hpp"
#include "../include/misc.hpp"
#include "../include/write.hpp"

using namespace Eigen;
using namespace std;

#define EPS 1e-8
// #define write 1

#define START_TIMER() begin = chrono::steady_clock::now()
#define END_TIMER() end = chrono::steady_clock::now()
#define ADD_TIME_TO(x) x += chrono::duration_cast<chrono::microseconds>(end - begin).count()

#define WRITE_ANIM() ((t - EPS < plot_count * plot_dt) && (t + EPS > plot_count * plot_dt))

int main(void) {
  const string filename_input = "config/input.txt";  // name of the input file to read the parameters
  const string space = "    ";                       // space to print
  const uint per = 10;                               // progress percentage interval to print (each count%)
  uint count = per;
  Prm prm;
  bool animation;       // flag to enable or disable the animation
  double plot_dt;       // frequency to write the animation file
  uint plot_count = 0;  // counter to write the animation file
  string plot_var;      // variable to plot in the animation
  int64_t total_files = 0, total_advection = 0, total_diffusion = 0,
          total_others = 0, total_build_poisson = 0, total_solve_laplace = 0,
          total_constr_obstacle = 0, total_boundary = 0, total_buoyancy = 0;  // variables to measure the time
  chrono::steady_clock::time_point begin, end;                                // variables to measure the time

  // ------------- File input setup ----------------
  START_TIMER();
  ifstream file_input;
  file_input.open(filename_input);
  string tmp;
  if (file_input.is_open()) {
    file_input >> tmp >> prm.L;
    file_input >> tmp >> prm.H;
    file_input >> tmp >> prm.nx;
    file_input >> tmp >> prm.ny;
    file_input >> tmp >> prm.dt;
    file_input >> tmp >> prm.T;
    file_input >> tmp >> prm.Pr;
    file_input >> tmp >> prm.Le;
    file_input >> tmp >> prm.Ra_T;
    file_input >> tmp >> prm.R_rho;
    file_input >> tmp >> prm.A_T;
    file_input >> tmp >> prm.A_S;
    file_input >> tmp >> plot_var;
    file_input >> tmp >> plot_dt;
    file_input >> tmp >> animation;
  }
  file_input.close();

  if (plot_var != "u" && plot_var != "v" && plot_var != "uv" && plot_var != "T" && plot_var != "S" && plot_var != "p") {
    cout << "Invalid variable to plot. Exiting." << endl;
    return 1;
  }

  prm.NX = prm.nx + 2;
  prm.NY = prm.ny + 2;
  prm.dx = prm.L / prm.nx;
  prm.dy = prm.H / prm.ny;
  prm.NXNY = (uint)(prm.NX * prm.NY);
  prm.nxny = (uint)(prm.nx * prm.ny);
  prm.Ra_S = prm.Ra_T / prm.R_rho * prm.Le;
  prm.Sc = prm.Pr * prm.Le;
  // -------------------------------------------

  // ------------- Print plot setup -----------------
  cout << "dx:            " << space << prm.dx << endl;
  cout << "dy:            " << space << prm.dy << endl;
  cout << "dt:            " << space << prm.dt << endl;
  cout << "T:             " << space << prm.T << endl;
  cout << "Pr:            " << space << prm.Pr << endl;
  cout << "Sc:            " << space << prm.Sc << endl;
  cout << "Le:            " << space << prm.Le << endl;
  cout << "Ra_T:          " << space << prm.Ra_T << endl;
  cout << "Ra_S:          " << space << prm.Ra_S << endl;
  cout << "R_rho:         " << space << prm.R_rho << endl;
  cout << "Plot animation?" << space << (animation ? "yes" : "no") << ((plot_var == "u" || plot_var == "v" || plot_var == "uv") ? " (velocity)" : (plot_var == "T" ? " (temperature)" : (plot_var == "S" ? " (salinity)" : " (pressure)"))) << endl;
  // ------------------------------------------------

  double t = 0.0, mean = 0;
  uint numSteps = 0;
  double fraction_completed = prm.T / 100.;  // fraction of the integration time to print

  // allocate memory
  double* u = new double[prm.NXNY];      // x component of velocity
  double* ustar = new double[prm.NXNY];  // x component of velocity
  double* adv_u = new double[prm.NXNY];  // advection term of u
  double* v = new double[prm.NXNY];      // y component of velocity
  double* vstar = new double[prm.NXNY];  // y component of velocity
  double* adv_v = new double[prm.NXNY];  // advection term of v
  double* p = new double[prm.NXNY];      // pressure
  double* T = new double[prm.NXNY];      // temperature
  // double* T_eq = new double[prm.NXNY];     // equilibrium temperature
  // double* Delta_T = new double[prm.NXNY];  // Delta_T
  double* S = new double[prm.NXNY];  // salinity
  // double* S_eq = new double[prm.NXNY];     // equilibrium salinity
  // double* Delta_S = new double[prm.NXNY];  // Delta_S

  // initialize to 0
  for (int i = 0; i < prm.NXNY; i++) {
    u[i] = 0;
    ustar[i] = 0;
    adv_u[i] = 0;
    v[i] = 0;
    vstar[i] = 0;
    adv_v[i] = 0;
    p[i] = 0;
    T[i] = 0;
    // T_eq[i] = 0;
    // Delta_T[i] = 0;
    S[i] = 0;
    // S_eq[i] = 0;
    // Delta_S[i] = 0;
  }
  // Set T(x, y) = y * A_T * cos(pi * x / L) * sinh(pi * y / L)
  // Set S(x, y) = y * A_S * cos(pi * x / L) * sinh(pi * y / L)
  double aux = 2 * M_PI / prm.L;
  for (int i = 0; i < prm.NX; i++) {
    for (int j = 0; j < prm.NY; j++) {
      // the sinh(pi * y / L) term is a normalization factor for the derivative (it)
      // T_eq(i, j) = y(j) * prm.A_T * cos(M_PI * x(i) / prm.L) * sinh(M_PI * y(j) / prm.L) / sinh(M_PI * prm.H / prm.L);
      // T(i, j) = T_eq(i, j);
      // S_eq(i, j) = y(j) * prm.A_S * cos(M_PI * x(i) / prm.L) * sinh(M_PI * y(j) / prm.L) / sinh(M_PI * prm.H / prm.L);
      // S(i, j) = S_eq(i, j);

      // U(i, j) = cos(M_PI * x(i) / prm.L) * (cos(M_PI * y(j) / prm.H) - 1);

      // for Rayleigh-Benard convection
      // T_eq(i, j) = y(j) * prm.A_T / prm.H;
      // T_eq(i, j) = (prm.H - y(j)) * prm.A_T / prm.H;
      // T_eq(i, j) = prm.A_T * cos(2 * M_PI * x(i) / prm.L) * sinh(2 * M_PI * y(j) / prm.L);
      // T(i, j) = T_eq(i, j);
      // S_eq(i, j) = prm.A_S * cos(2 * M_PI * x(i) / prm.L) * sinh(2 * M_PI * y(j) / prm.L);
      // S(i, j) = S_eq(i, j);
      // T_eq(i, j) = prm.A_T * cos(aux * x(i)) * sinh(aux * y(j)) / sinh(aux * prm.H);
      // T(i, j) = T_eq(i, j);
      // S_eq(i, j) = prm.A_S * cos(aux * x(i)) * sinh(aux * y(j)) / sinh(aux * prm.H);
      // S(i, j) = S_eq(i, j);

      // T(i, j) = -prm.A_T * cos(aux * x(i)) * sinh(aux * y(j)) / sinh(aux * prm.H);
      // S(i, j) = -prm.A_S * cos(aux * x(i)) * sinh(aux * y(j)) / sinh(aux * prm.H);

      T(i, j) = sin(M_PI * y(j) / prm.H);
      S(i, j) = sin(M_PI * y(j) / prm.H);
    }
  }
  // hot spot in the middle
  // for (int i = prm.NX / 4; i < 3 * prm.NX / 4; i++) {
  //   for (int j = prm.NY / 4; j < 3 * prm.NY / 4; j++) {
  //     T(i, j) += prm.A_T;
  //   }
  // }

  string filename_out = "output/prova.txt";
  ofstream file_out;
  file_out.open(filename_out);

  saveSetupToHDF5(prm, plot_var, animation);
  saveDataToHDF5(plot_count, u, v, T, S, p, prm.NX, prm.NY, t);
#ifdef write
  write_sol(file_out, T, t, prm, true);
#endif
  plot_count++;
  END_TIMER();
  ADD_TIME_TO(total_files);
  // -------------------------------------------

  // ------------- Build Poission matrix -----------
  START_TIMER();
  vector<Trip> coeffs;
  buildPoissonMatrix(coeffs, prm);

  SpMat A(prm.nx * prm.ny, prm.nx * prm.ny);
  VectorXd div(prm.nx * prm.ny);
  A.setFromTriplets(coeffs.begin(), coeffs.end());

  SimplicialLDLT<SpMat> chol;
  chol.compute(A);  // performs a Cholesky factorization of A. The matrix has to be symmetric and positive definite for this to work as expected
  if (chol.info() != Success) {
    cout << "Cholesky decomposition failed" << endl;
    return 1;
  }

  VectorXd p_solved;
  END_TIMER();
  ADD_TIME_TO(total_build_poisson);
  // -------------------------------------------

  numSteps++;
  double max_u_v = 0, lap, buoy, rhs_u, rhs_v;

  while (t < prm.T - EPS) {
    // adaptative time step
    // cout << "dt: " << prm.dt << endl;
    // for (int i = 0; i < prm.NXNY; i++) {
    //   max_u_v = max(max_u_v, max(abs(u[i]), abs(v[i])));
    // }
    // if (max_u_v > 1e-3) {
    //   // update CFL number
    //   prm.dt = prm.dx * prm.dy / (2 * max_u_v * (prm.dx + prm.dy));  // the factor 2 is to be conservative
    // }

    // ---------- advection term ------------
    START_TIMER();
    // Semilag(u, v, u, prm, 1);
    // Semilag(u, v, v, prm, 1);
    // memcpy(adv_u, u, prm.NXNY * sizeof(double));
    // memcpy(adv_v, v, prm.NXNY * sizeof(double));
    Semilag2(u, v, u, adv_u, prm);
    Semilag2(u, v, v, adv_v, prm);
    END_TIMER();
    ADD_TIME_TO(total_advection);
    // --------------------------------------
#ifdef write
    file_out << "adevected u" << endl;
    write_sol(file_out, adv_u, t, prm, true);

    file_out << "adevected v" << endl;
    write_sol(file_out, adv_v, t, prm, true);
#endif

    // ---------- diffusion term ------------
    START_TIMER();
    for (int i = 1; i < prm.NX - 1; i++) {
      for (int j = 1; j < prm.NY - 1; j++) {
        lap = (U(i + 1, j) - 2 * U(i, j) + U(i - 1, j)) / (prm.dx * prm.dx) +
              (U(i, j + 1) - 2 * U(i, j) + U(i, j - 1)) / (prm.dy * prm.dy);
        rhs_u = ADV_U(i, j) + prm.dt * prm.Pr * lap;
        lap = (V(i + 1, j) - 2 * V(i, j) + V(i - 1, j)) / (prm.dx * prm.dx) +
              (V(i, j + 1) - 2 * V(i, j) + V(i, j - 1)) / (prm.dy * prm.dy);
        buoy = prm.Ra_T * (T(i, j) - S(i, j) / prm.R_rho);
        // only temperature
        // buoy = prm.Ra_T * (T(i, j) - T_eq(i, j));
        // buoy = prm.Ra_T * (T(i, j));
        // buoy = 0;
        rhs_v = ADV_V(i, j) + prm.dt * prm.Pr * (lap + buoy);
        Ustar(i, j) = rhs_u;
        Vstar(i, j) = rhs_v;
      }
    }
    END_TIMER();
    ADD_TIME_TO(total_diffusion);
    // --------------------------------------
#ifdef write
    file_out << "diffused u" << endl;
    write_sol(file_out, ustar, t, prm, true);

    file_out << "diffused v" << endl;
    write_sol(file_out, vstar, t, prm, true);
#endif

    // update velocity ghost points
    START_TIMER();
    BC_velocity(ustar, vstar, prm);
    END_TIMER();
    ADD_TIME_TO(total_boundary);

    // ---------- pressure term ------------
    START_TIMER();
    // compute the (minus) divergence of the velocity field (we omit multiplication by dt, because later on we will divide by dt) and store it in div
    mean = 0;
    for (int i = 1; i < prm.NX - 1; i++) {
      for (int j = 1; j < prm.NY - 1; j++) {
        DIV(i - 1, j - 1) = -(Ustar(i + 1, j) - Ustar(i - 1, j)) / (2 * prm.dx) -
                            (Vstar(i, j + 1) - Vstar(i, j - 1)) / (2 * prm.dy);
        mean += DIV(i - 1, j - 1);
      }
    }
    // subtract the mean value of the divergence in order to ensure that the system is solvable
    mean /= prm.nxny;
    for (int i = 0; i < prm.nxny; i++) div[i] -= mean;
    END_TIMER();
    ADD_TIME_TO(total_others);

    START_TIMER();
    // solve the pressure Poisson equation
    p_solved = chol.solve(div);
    if (chol.info() != Success) {
      cout << "Cholesky solve failed" << endl;
      return 1;
    }
    // convert p_solved back to pointer
    for (int i = 1; i < prm.NX - 1; i++) {
      for (int j = 1; j < prm.NY - 1; j++) {
        P(i, j) = p_solved((i - 1) * prm.ny + (j - 1));
      }
    }

    END_TIMER();
    ADD_TIME_TO(total_solve_laplace);

    // update pressure ghost points
    START_TIMER();
    BC_pressure(p, prm);
    END_TIMER();
    ADD_TIME_TO(total_boundary);

    START_TIMER();
    // ---------- projection step ------------
    for (int i = 1; i < prm.NX - 1; i++) {
      for (int j = 1; j < prm.NY - 1; j++) {
        U(i, j) = Ustar(i, j) - (P(i + 1, j) - P(i - 1, j)) / (2 * prm.dx);
        V(i, j) = Vstar(i, j) - (P(i, j + 1) - P(i, j - 1)) / (2 * prm.dy);
      }
    }
    BC_velocity(u, v, prm);
    END_TIMER();
    ADD_TIME_TO(total_others);
#ifdef write
    file_out << "final u" << endl;
    write_sol(file_out, u, t, prm, true);

    file_out << "final v" << endl;
    write_sol(file_out, v, t, prm, true);
#endif
    // --------------------------------------

    // ---------- temperature Navier-Stokes ------------
    // advection term
    START_TIMER();
    // Semilag(u, v, T, prm, 1);
    // Semilag(u, v, S, prm, 1);
    // memcpy(adv_u, T, prm.NXNY * sizeof(double));
    // memcpy(adv_v, S, prm.NXNY * sizeof(double));
    Semilag2(u, v, T, adv_u, prm);
    Semilag2(u, v, S, adv_v, prm);
#ifdef write
    file_out << "adevected T" << endl;
    write_sol(file_out, adv_u, t, prm, true);

    file_out << "adevected S" << endl;
    write_sol(file_out, adv_v, t, prm, true);
#endif

    END_TIMER();
    ADD_TIME_TO(total_advection);

    // diffusion term
    START_TIMER();

    // crankNicholson(T, adv_u, S, adv_v, prm);

    for (int i = 1; i < prm.NX - 1; i++) {
      for (int j = 1; j < prm.NY - 1; j++) {
        lap = (T(i + 1, j) - 2 * T(i, j) + T(i - 1, j)) / (prm.dx * prm.dx) +
              (T(i, j + 1) - 2 * T(i, j) + T(i, j - 1)) / (prm.dy * prm.dy);
        ADV_U(i, j) += prm.dt * lap;
        lap = (S(i + 1, j) - 2 * S(i, j) + S(i - 1, j)) / (prm.dx * prm.dx) +
              (S(i, j + 1) - 2 * S(i, j) + S(i, j - 1)) / (prm.dy * prm.dy);
        ADV_V(i, j) += prm.dt * lap / prm.Le;
      }
    }
    memcpy(T, adv_u, prm.NXNY * sizeof(double));
    memcpy(S, adv_v, prm.NXNY * sizeof(double));

    END_TIMER();
    ADD_TIME_TO(total_diffusion);
#ifdef write
    file_out << "diffused T" << endl;
    write_sol(file_out, T, t, prm, true);

    file_out << "diffused S" << endl;
    write_sol(file_out, S, t, prm, true);
#endif

    // boundary conditions
    START_TIMER();
    BC_temperature(T, prm);
    BC_salinity(S, prm);
    // for (int i = 0; i < prm.NX; i++) {
    //   for (int j = 0; j < prm.NY; j++) {
    //     Delta_T(i, j) = T(i, j) - T_eq(i, j);
    //     Delta_S(i, j) = S(i, j) - S_eq(i, j);
    //   }
    // }
    END_TIMER();
    ADD_TIME_TO(total_boundary);

    t += prm.dt;

    // ---------- printing to file ------------
    START_TIMER();
    if (WRITE_ANIM()) {
      saveDataToHDF5(plot_count, u, v, T, S, p, prm.NX, prm.NY, t);
#ifdef write
      file_out << "final T" << endl;
      write_sol(file_out, T, t, prm, true);

      file_out << "final S" << endl;
      write_sol(file_out, S, t, prm, true);
#endif
      plot_count++;
    }

    if (t > fraction_completed * count - EPS) {
      cout << count << "%"
           << " dt: " << prm.dt << endl;
      count += per;
    }
    END_TIMER();
    ADD_TIME_TO(total_files);
    // --------------------------------------

    numSteps++;
  }

  print("Total time for files:                   " + space, total_files);
  print("Total time for construction of obstacle:" + space, total_constr_obstacle);
  print("Total time for build laplace:           " + space, total_build_poisson);
  print("Total time for advection:               " + space, total_advection);
  print("Total time for diffusion:               " + space, total_diffusion);
  print("Total time for buoyancy:                " + space, total_buoyancy);
  print("Total time for solve laplace:           " + space, total_solve_laplace);
  print("Total time for boundary:                " + space, total_boundary);
  print("Total time for rest:                    " + space, total_others);
  printf("-----------------------------------------------------------\n");
  print("Total time:                             " + space, total_files + total_advection + total_diffusion + total_others + total_build_poisson + total_solve_laplace + total_constr_obstacle + total_boundary + total_buoyancy);

  // free memory
  delete[] u;
  delete[] ustar;
  delete[] v;
  delete[] vstar;
  delete[] T;
  delete[] S;
  delete[] p;
  delete[] adv_u;
  delete[] adv_v;

  return 0;
}
