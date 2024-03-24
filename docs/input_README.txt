L                       Length of the domain in the x-direction. Positive real number.
H                       Length of the domain in the y-direction. Positive real number.
nx                      Number of grid points in the x-direction (excluding ghost cells). Positive integer.
ny                      Number of grid points in the y-direction (excluding ghost cells). Positive integer.
dt                      Initial time step (if all goes well, it should be adaptative). Positive real number.
Tfinal                  Time at which the simulation should stop. Positive real number.
Pr                      Prandtl number of the fluid. Positive real number. For salt water, Pr = 13.4 at 0°C and Pr = 7.2 at 20°C.
Le                      Lewis number of the fluid. Positive real number.
Ra_T                    Rayleigh number for temperature. Positive real number.
Ra_S                    Rayleigh number for salinity. Positive real number.
A_T                     Amplitude of the temperature flow at the top boundary. Positive real number.
A_S                     Amplitude of the salinity flow at the top boundary. Positive real number.
plot                    Variable to be plotted. String: 'T' (temperature), 'S' (salinity), 'u' (modulus of the velocity), 'p' (pressure)
plot_dt                 Time step for the animation (print solution each animation_dt time-steps). Positive integer.
animation               Whether to animate the solution or not. Boolean: 1 (yes) or 0 (no).