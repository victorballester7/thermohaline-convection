#ifndef WRITE_HPP
#define WRITE_HPP

#include <H5Cpp.h>

#include <cmath>
#include <fstream>
#include <iostream>
#include <string>

#include "misc.hpp"

using namespace H5;
using namespace std;

// @brief Print a string and a time in a human-readable format (in the form: "string xx.xx units")
// @param str string to be printed
// @param time time to be printed
void print(string str, int64_t time) {
  // default time is in microseconds
  string microsec = "\u03BCs";
  string millisec = "ms";
  string sec = "s";
  string min = "min";
  string hour = "h";
  double time_d;
  if (time > 10000) {
    time /= 1000;  // now time is in milliseconds
    if (time > 10000) {
      time /= 1000;  // now time is in seconds
      if (time > 100) {
        time_d = (double)time;
        time_d /= 60;  // now time is in minutes
        if (time_d > 100) {
          time_d /= 60;  // now time is in hours
          cout << str << time_d << " " << hour << endl;
          return;
        }
        cout << str << time_d << " " << min << endl;
        return;
      }
      cout << str << time << " " << sec << endl;
      return;
    }
    cout << str << time << " " << millisec << endl;
    return;
  }
  cout << str << time << " " << microsec << endl;
  return;
}

// @brief Write the setup parameters to a file
// @param prm parameters of the simulation
// @param object_type type of the object
// @param vorticity_on boolean to check whether we want to plot the vorticity on the animation by default or not
// @param animation_on boolean to check whether we want to plot the animation or not
void saveSetupToHDF5(Prm prm, string plot_var, bool animation_on) {
  string filename = "output/setup.h5";
  const H5std_string DATASET_NAMES[] = {"NX", "NY", "L", "H", "Pr", "Le", "Ra_T", "R_rho", "dx", "dy", "dt", "Tfinal", "animation_on"};
  int num_datasets = 13;
  try {
    H5File file(filename, H5F_ACC_TRUNC);
    hsize_t dims_scalar[1] = {1};
    DataSpace scalar_dataspace(1, dims_scalar);

    // Create and write datasets for the variables
    double datasets[] = {(double)prm.NX, (double)prm.NY, prm.L, prm.H, prm.Pr, prm.Le, prm.Ra_T, prm.R_rho, prm.dx, prm.dy, prm.dt, prm.Tfinal, (double)animation_on};
    for (int i = 0; i < num_datasets; i++) {
      DataSet dataset = file.createDataSet(DATASET_NAMES[i], PredType::NATIVE_DOUBLE, scalar_dataspace);
      dataset.write(&datasets[i], PredType::NATIVE_DOUBLE);
    }
    // Saving the object type
    // Saving the object type
    H5std_string object_type_name("plot_var");
    H5::StrType strdatatype(H5::PredType::C_S1, plot_var.size());
    DataSet obstacle_dataset = file.createDataSet(object_type_name, strdatatype, scalar_dataspace);
    obstacle_dataset.write(plot_var.c_str(), strdatatype);
  } catch (FileIException& error) {
    error.printErrorStack();
    throw;  // Re-throw the exception to notify the caller of the failure.
  } catch (DataSetIException& error) {
    error.printErrorStack();
    throw;
  } catch (DataSpaceIException& error) {
    error.printErrorStack();
    throw;
  }
}

// @brief Write the data to a file
// @param plot_count index of the plot
// @param u velocity in the x direction
// @param v velocity in the y direction
// @param w voriticity
// @param p pressure
// @param Nx number of cells in the x direction (including ghost cells)
// @param Ny number of cells in the y direction (including ghost cells)
// @param t time
void saveDataToHDF5(uint plot_count, double* u, double* v, double* T, double* S, double* p, int Nx, int Ny, double t) {
  string filename = "output/results/sol_";
  filename = filename + to_string(plot_count) + ".h5";
  const H5std_string DATASET_NAMES[] = {"u", "v", "T", "S", "p"};
  const H5std_string TIMESTAMP_NAME("t");

  try {
    H5File file(filename, H5F_ACC_TRUNC);
    hsize_t dims[2] = {static_cast<hsize_t>(Nx), static_cast<hsize_t>(Ny)};
    DataSpace dataspace(2, dims);

    // Create and write datasets for u, v, w, p
    double* datasets[] = {u, v, T, S, p};
    for (int i = 0; i < 5; i++) {
      DataSet dataset = file.createDataSet(DATASET_NAMES[i], PredType::NATIVE_DOUBLE, dataspace);
      dataset.write(datasets[i], PredType::NATIVE_DOUBLE);
    }

    // Saving the timestamp 't'
    hsize_t dims_scalar[1] = {1};
    DataSpace scalar_dataspace(1, dims_scalar);
    DataSet dataset_t = file.createDataSet(TIMESTAMP_NAME, PredType::NATIVE_DOUBLE, scalar_dataspace);
    dataset_t.write(&t, PredType::NATIVE_DOUBLE);
  } catch (FileIException& error) {
    error.printErrorStack();
    throw;  // Re-throw the exception to notify the caller of the failure.
  } catch (DataSetIException& error) {
    error.printErrorStack();
    throw;
  } catch (DataSpaceIException& error) {
    error.printErrorStack();
    throw;
  }
}

void write_sol(ofstream& file, double* phi, double t, Prm prm, bool ghost) {
  file << t << endl;
  if (ghost) {
    for (int i = 0; i < prm.NX; i++) {
      for (int j = 0; j < prm.NY; j++) {
        file << phi[i * prm.NY + j] << " ";
      }
      file << endl;
    }
  } else {
    for (int i = 1; i < prm.NX - 1; i++) {
      for (int j = 1; j < prm.NY - 1; j++) {
        file << phi[i * prm.NY + j] << " ";
      }
      file << endl;
    }
  }
}

#endif  // WRITE_HPP