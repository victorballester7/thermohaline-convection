# Numerical simulation of the thermohaline circulation in the ocean

This is a simple numerical simulation of the thermohaline circulation in the ocean. The simulation assumes that the density of the ocean is a linear function of temperature and salinity only in order to use then the Bousinessq approximation for the Navier-Stokes equations. You can find more information in the report of the project, inside the folder `latex/report`.

## Compilation and execution of the code

The code is written in `C++` and the animation is done in `Python`. You will need the library `Eigen` to compile it, as it is used for the linear algebra operations (with sparse matrices). You can check its webpage from [here](http://eigen.tuxfamily.org/index.php?title=Main_Page).

In order to properly compile the code, you will to ensure that the `eigen3` folder is correctly linked in the `Makefile` file. Either you have the `eigen3` folder in the path `/usr/include/`, being this folder in the path of the system, or you will need to specify the library path in the `Makefile` in the `INCLUDES` variable. Keep in mind that you may also need to change `python` to `python3` in the `run.sh` file if you are using `python3` instead of `python`.

To compile the code, first download it:

```
git clone git@github.com:victorballester7/thermohaline-convection.git
cd thermohaline-convection
```

Then, you can compile and run the code (having previously check your input parameters in the file `config/input.txt`) with the following command:

```
./run.sh
```

If you get an error, you may need to activate the execution permissions of the file `run.sh` with the following command:

```
chmod +x run.sh
```

If activated in the file `config/input.txt` (check the readme of the input file in `docs/input_README.txt`), the code will also animate the solution.

## Other information

The code is mainly in the `src` folder, and the headers are in the `include` folder. In the `config` folder you can find the input file and the output files are found in the `output` folder. Finally, in the `docs` folder you can find the documentation of the input file and other papers with information about the thermohaline circulation.

## Results

Here you can see some results of the simulation.
