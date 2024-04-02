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

#define START_TIMER() begin = chrono::steady_clock::now()
#define END_TIMER() end = chrono::steady_clock::now()
#define ADD_TIME_TO(x) x += chrono::duration_cast<chrono::microseconds>(end - begin).count()

#define WRITE_ANIM() (animation && (t - EPS < plot_count * plot_dt) && (t + EPS > plot_count * plot_dt))

int main(int argc, char* argv[]) {
  const string filename_input = "config/input.txt";  // name of the input file to read the parameters
  const string filename_output_Le_Rrho = "output/Le_Rrho.txt";
  const string space = "    ";  // space to print
  const uint per = 10;          // progress percentage interval to print (each count%)
  uint count = per;
  Prm prm;
  bool animation;       // flag to enable or disable the animation
  double plot_dt;       // frequency to write the animation file
  uint plot_count = 0;  // counter to write the animation file
  string plot_var;      // variable to plot in the animation
  int64_t total_files = 0, total_advection = 0, total_diffusion = 0,
          total_others = 0, total_build_poisson = 0, total_solve_laplace = 0,
          total_boundary = 0;                   // variables to measure the time
  chrono::steady_clock::time_point begin, end;  // variables to measure the time

  bool remove_cout = false;
  // if one argument is passed
  if (argc == 2) {
    remove_cout = true;  // we are running the full simulation
  }

  // ------------- File input setup ----------------
  START_TIMER();
  ifstream file_input;
  ofstream file_output_Le_Rrho;
  file_input.open(filename_input);
  file_output_Le_Rrho.open(filename_output_Le_Rrho, ios::app);  // append to the file
  string tmp;
  if (file_input.is_open()) {
    file_input >> tmp >> prm.L;
    file_input >> tmp >> prm.H;
    file_input >> tmp >> prm.nx;
    file_input >> tmp >> prm.ny;
    file_input >> tmp >> prm.dt;
    file_input >> tmp >> prm.Tfinal;
    file_input >> tmp >> prm.Pr;
    file_input >> tmp >> prm.Le;
    file_input >> tmp >> prm.Ra_T;
    file_input >> tmp >> prm.R_rho;
    file_input >> tmp >> prm.A_T;
    file_input >> tmp >> prm.A_S;
    file_input >> tmp >> prm.implicit;
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
  // -------------------------------------------

  if (!remove_cout) {
    // ------------- Print plot setup -----------------
    cout << "dx:            " << space << prm.dx << endl;
    cout << "dy:            " << space << prm.dy << endl;
    cout << "dt:            " << space << prm.dt << endl;
    cout << "Tfinal:        " << space << prm.Tfinal << endl;
    cout << "Pr:            " << space << prm.Pr << endl;
    cout << "Le:            " << space << prm.Le << endl;
    cout << "Ra_T:          " << space << prm.Ra_T << endl;
    cout << "R_rho:         " << space << prm.R_rho << endl;
    cout << "Plot animation?" << space << (animation ? "yes" : "no") << ((plot_var == "u" || plot_var == "v" || plot_var == "uv") ? " (velocity)" : (plot_var == "T" ? " (temperature)" : (plot_var == "S" ? " (salinity)" : " (pressure)"))) << endl;
    // ------------------------------------------------
  }
  double t = 0.0, mean = 0;
  double meanFluxV_0 = 0, meanFluxV_1 = 0;
  bool changed = false;
  uint numSteps = 0;
  double fraction_completed = prm.Tfinal / 100.;  // fraction of the integration time to print

  // allocate memory
  double* u = new double[prm.NXNY];      // x component of velocity
  double* ustar = new double[prm.NXNY];  // x component of velocity
  double* adv_u = new double[prm.NXNY];  // advection term of u
  double* v = new double[prm.NXNY];      // y component of velocity
  double* vstar = new double[prm.NXNY];  // y component of velocity
  double* adv_v = new double[prm.NXNY];  // advection term of v
  double* p = new double[prm.NXNY];      // pressure
  double* T = new double[prm.NXNY];      // temperature
  double* S = new double[prm.NXNY];      // salinity

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
    S[i] = 0;
  }
  double aux = 2 * M_PI / prm.L;
  for (int i = 0; i < prm.NX; i++) {
    for (int j = 0; j < prm.NY; j++) {
      T(i, j) = -prm.A_T * cos(aux * x(i)) * sinh(aux * y(j)) / sinh(aux * prm.H);
      S(i, j) = -prm.A_S * cos(aux * x(i)) * sinh(aux * y(j)) / sinh(aux * prm.H);
    }
  }

  saveSetupToHDF5(prm, plot_var, animation, plot_dt);
  saveDataToHDF5(plot_count, u, v, T, S, p, prm.NX, prm.NY, t);
  plot_count++;
  END_TIMER();
  ADD_TIME_TO(total_files);
  // -------------------------------------------

  // ------------- Build Poission matrix -----------
  START_TIMER();
  vector<Trip> coeffsP, coeffsT, coeffsS, coeffsU, coeffsV;
  double lambdaU, lambdaV, muU, muV;
  double lambdaT, lambdaS, muT, muS;
  lambdaU = prm.dt / (prm.dx * prm.dx) * prm.Pr;
  lambdaV = lambdaU;
  muU = prm.dt / (prm.dx * prm.dx) * prm.Pr;
  muV = muU;
  lambdaT = prm.dt / (prm.dx * prm.dx);
  lambdaS = prm.dt / (prm.dx * prm.dx) / prm.Le;
  muT = prm.dt / (prm.dy * prm.dy);
  muS = prm.dt / (prm.dy * prm.dy) / prm.Le;
  buildPoissonMatrixP(coeffsP, prm);
  buildHeatMatrix_UTS(coeffsU, lambdaU, muU, prm);
  buildHeatMatrix_UTS(coeffsT, lambdaT, muT, prm);
  buildHeatMatrix_UTS(coeffsS, lambdaS, muS, prm);
  buildHeatMatrix_V(coeffsV, lambdaV, muV, prm);

  SpMat A(prm.nx * prm.ny, prm.nx * prm.ny);
  VectorXd div(prm.nx * prm.ny);
  A.setFromTriplets(coeffsP.begin(), coeffsP.end());

  SimplicialLDLT<SpMat> chol;
  chol.compute(A);  // performs a Cholesky factorization of A. The matrix has to be symmetric and positive definite for this to work as expected
  if (chol.info() != Success) {
    cout << "Cholesky decomposition failed" << endl;
    return 1;
  }

  A.setFromTriplets(coeffsU.begin(), coeffsU.end());
  SimplicialLDLT<SpMat> cholU;
  cholU.compute(A);  // performs a Cholesky factorization of A. The matrix has to be symmetric and positive definite for this to work as expected
  if (cholU.info() != Success) {
    cout << "Cholesky decomposition failed" << endl;
    return 1;
  }

  A.setFromTriplets(coeffsV.begin(), coeffsV.end());
  SimplicialLDLT<SpMat> cholV;
  cholV.compute(A);  // performs a Cholesky factorization of A. The matrix has to be symmetric and positive definite for this to work as expected
  if (cholV.info() != Success) {
    cout << "Cholesky decomposition failed" << endl;
    return 1;
  }

  A.setFromTriplets(coeffsT.begin(), coeffsT.end());
  SimplicialLDLT<SpMat> cholT;
  cholT.compute(A);  // performs a Cholesky factorization of A. The matrix has to be symmetric and positive definite for this to work as expected
  if (cholT.info() != Success) {
    cout << "Cholesky decomposition failed" << endl;
    return 1;
  }

  A.setFromTriplets(coeffsS.begin(), coeffsS.end());
  SimplicialLDLT<SpMat> cholS;
  cholS.compute(A);  // performs a Cholesky factorization of A. The matrix has to be symmetric and positive definite for this to work as expected
  if (cholS.info() != Success) {
    cout << "Cholesky decomposition failed" << endl;
    return 1;
  }

  VectorXd p_solved;
  END_TIMER();
  ADD_TIME_TO(total_build_poisson);
  // -------------------------------------------

  numSteps++;
  double max_u_v = 0, lap, buoy, rhs_u, rhs_v;

  while (t < prm.Tfinal - EPS) {
    // ---------- advection term ------------
    START_TIMER();
    Semilag2(u, v, u, adv_u, prm);
    Semilag2(u, v, v, adv_v, prm);
    END_TIMER();
    ADD_TIME_TO(total_advection);
    // --------------------------------------

    // ---------- diffusion term ------------
    START_TIMER();
    if (prm.implicit) {
      // solve the u heat equation
      for (int i = 1; i < prm.NX - 1; i++) {
        for (int j = 1; j < prm.NY - 1; j++) {
          DIV(i - 1, j - 1) = ADV_U(i, j);
        }
      }
      p_solved = cholU.solve(div);
      if (cholU.info() != Success) {
        cout << "Cholesky solve failed" << endl;
        return 1;
      }
      // convert p_solved back to pointer
      for (int i = 1; i < prm.NX - 1; i++) {
        for (int j = 1; j < prm.NY - 1; j++) {
          Ustar(i, j) = p_solved((i - 1) * prm.ny + (j - 1));
        }
      }

      // solve the v heat equation with external forcing (buoyancy)
      for (int i = 1; i < prm.NX - 1; i++) {
        for (int j = 1; j < prm.NY - 1; j++) {
          DIV(i - 1, j - 1) = ADV_V(i, j) + prm.dt * prm.Pr * prm.Ra_T * (T(i, j) - S(i, j) / prm.R_rho);
        }
      }
      p_solved = cholV.solve(div);
      if (cholV.info() != Success) {
        cout << "Cholesky solve failed" << endl;
        return 1;
      }
      // convert p_solved back to pointer
      for (int i = 1; i < prm.NX - 1; i++) {
        for (int j = 1; j < prm.NY - 1; j++) {
          Vstar(i, j) = p_solved((i - 1) * prm.ny + (j - 1));
        }
      }
    } else {
      for (int i = 1; i < prm.NX - 1; i++) {
        for (int j = 1; j < prm.NY - 1; j++) {
          lap = (U(i + 1, j) - 2 * U(i, j) + U(i - 1, j)) / (prm.dx * prm.dx) +
                (U(i, j + 1) - 2 * U(i, j) + U(i, j - 1)) / (prm.dy * prm.dy);
          rhs_u = ADV_U(i, j) + prm.dt * prm.Pr * lap;
          lap = (V(i + 1, j) - 2 * V(i, j) + V(i - 1, j)) / (prm.dx * prm.dx) +
                (V(i, j + 1) - 2 * V(i, j) + V(i, j - 1)) / (prm.dy * prm.dy);
          buoy = prm.Ra_T * (T(i, j) - S(i, j) / prm.R_rho);
          rhs_v = ADV_V(i, j) + prm.dt * prm.Pr * (lap + buoy);
          Ustar(i, j) = rhs_u;
          Vstar(i, j) = rhs_v;
        }
      }
    }
    END_TIMER();
    ADD_TIME_TO(total_diffusion);
    // --------------------------------------

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
    // --------------------------------------

    // ---------- temperature & salinity Navier-Stokes ------------
    // advection term
    START_TIMER();
    Semilag2(u, v, T, adv_u, prm);
    Semilag2(u, v, S, adv_v, prm);
    END_TIMER();
    ADD_TIME_TO(total_advection);

    // diffusion term
    START_TIMER();
    if (prm.implicit) {  // implicit time stepping
      // solve the Temperature heat equation
      for (int i = 1; i < prm.NX - 1; i++) {
        for (int j = 1; j < prm.NY - 1; j++) {
          DIV(i - 1, j - 1) = ADV_U(i, j);
          if (j == prm.NY - 2) {
            DIV(i - 1, j - 1) += flux_T(x(i), prm) * prm.dt / prm.dy;
          }
        }
      }
      p_solved = cholT.solve(div);
      if (cholT.info() != Success) {
        cout << "Cholesky solve failed" << endl;
        return 1;
      }
      // convert p_solved back to pointer
      for (int i = 1; i < prm.NX - 1; i++) {
        for (int j = 1; j < prm.NY - 1; j++) {
          T(i, j) = p_solved((i - 1) * prm.ny + (j - 1));
        }
      }

      // solve the Salinity heat equation
      for (int i = 1; i < prm.NX - 1; i++) {
        for (int j = 1; j < prm.NY - 1; j++) {
          DIV(i - 1, j - 1) = ADV_V(i, j);
          if (j == prm.NY - 2) {
            DIV(i - 1, j - 1) += flux_S(x(i), prm) * prm.dt / prm.dy / prm.Le;
          }
        }
      }
      p_solved = cholS.solve(div);
      if (cholS.info() != Success) {
        cout << "Cholesky solve failed" << endl;
        return 1;
      }
      // convert p_solved back to pointer
      for (int i = 1; i < prm.NX - 1; i++) {
        for (int j = 1; j < prm.NY - 1; j++) {
          S(i, j) = p_solved((i - 1) * prm.ny + (j - 1));
        }
      }
    } else {  // explicit time stepping
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
    }
    END_TIMER();
    ADD_TIME_TO(total_diffusion);

    // boundary conditions
    START_TIMER();
    BC_temperature(T, prm);
    BC_salinity(S, prm);
    END_TIMER();
    ADD_TIME_TO(total_boundary);

    if (meanFluxV_0 * meanFluxV_1 < 0) {
      file_output_Le_Rrho << prm.R_rho << " " << prm.Le << (meanFluxV_0 < 0 ? " ST " : " TS ") << t << endl;
      changed = true;
    }
    meanFluxV_0 = meanFluxV_1;
    meanFluxV_1 = meanFluxV(v, prm);

    t += prm.dt;

    // ---------- printing to file ------------
    START_TIMER();
    if (WRITE_ANIM()) {
      saveDataToHDF5(plot_count, u, v, T, S, p, prm.NX, prm.NY, t);
      plot_count++;
    }

    if (t > fraction_completed * count - EPS && !remove_cout) {
      cout << count << "%"
           << " dt: " << prm.dt << endl;
      count += per;
    }
    END_TIMER();
    ADD_TIME_TO(total_files);
    // --------------------------------------

    numSteps++;
  }

  if (!remove_cout) {
    print("Total time for files:                   " + space, total_files);
    print("Total time for build laplace:           " + space, total_build_poisson);
    print("Total time for advection:               " + space, total_advection);
    print("Total time for diffusion:               " + space, total_diffusion);
    print("Total time for solve laplace:           " + space, total_solve_laplace);
    print("Total time for boundary:                " + space, total_boundary);
    print("Total time for rest:                    " + space, total_others);
    printf("-----------------------------------------------------------\n");
    print("Total time:                             " + space, total_files + total_advection + total_diffusion + total_others + total_build_poisson + total_solve_laplace + total_boundary);
  } else {
    if (!changed) {
      file_output_Le_Rrho << prm.R_rho << " " << prm.Le << (meanFluxV_0 < 0 ? " S " : " T ") << nan("") << endl;
    }
  }

  // free memory
  delete[] u;
  delete[] ustar;
  delete[] adv_u;
  delete[] v;
  delete[] vstar;
  delete[] adv_v;
  delete[] T;
  delete[] S;
  delete[] p;

  return 0;
}
