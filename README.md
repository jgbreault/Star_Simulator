Here I summarize the theoretical background behind equations.py and Star_Simulator.py. equations.py, is only used to define all equations used throughout Star_Simulator.py.

There are 5 coupled differential equations (DEs) that can be used describe a star. These equations describe how the star's density (rho), tempurature (T), mass within radius r (M), luminosty (L), and optical depth (tau) change with distance from the star's centre. These 5 coupled DEs can be solved using numerical integration starting from the star's centre if we have an inital condition corresponding to each DE. Tau is 0.0 at the star's centre. The initial conditions for M and L can be determined once we have the initial conditions for T and rho. To produce a star using Star_Simulator.py, The inital condition for T (T0) is simply set to be within the possible range of values (10.0^6.6 Kelvin to 10.0^7.4 Kelvin). Now the only unkown initial condition is rho0, so in the get_rho0 function the value for rho0 is geussed, checked, then refined using a bisection meathod for a given number of iterations. With all 5 initial conditions obtained, 4th order Runge-Kutta integration is used to solve the 5 coupled DEs. Producing a star takes 10-20 seconds. Larger stars take longer to integrate.

T0 can be varied over the possible range of values to reproduce the Main Sequence. The main sequence is a continuous and distinctive band of stars that appears on plots of stellar surface tempurature versus luminosity. Main_Sequence_from_Experiment.gif shows approximate values for surface tempurature and luminosity for real observed Main Sequence stars. The Main Sequence reproduced using this code (with Newtonian gravity) very closely resembles the shape and range of the real Main Sequence.

I then take things a step farther and toy with gravity. Newtonian gravity is given by G*M/r^2. What if in reality gravity is given by G*M/r^2.0*(1.0 + lam/r) instead. In this case, below some scale lam the gravitational force has a component that scales as r^−3. What constraints can we put on lam before we see noticible changes in the reproduced main sequence. In Main_Sequence_Small-Scale.png it can be seen that lam must be constrained within about -10^8 and 10^8 meters. Beyond this there is significant deviation from the Main Sequence reproduced using Newtonian Gravity.

What if in reality gravity is given by G*M/r^2.0*(1.0 + r/LAM) instead. In this case, above some scale LAM the gravitational force has a component that scales as r^−1. What constraints can we put on LAM before we see noticible changes in the reproduced main sequence. In Main_Sequence_Large-Scale.png it can be seen that LAM must be larger than 10^8 meters. Smaller values of LAM leads to significant deviation from the Main Sequence reproduced using Newtonian Gravity.
