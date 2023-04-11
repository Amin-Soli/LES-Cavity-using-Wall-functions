# An Object-Oriented program in C++ is written to simulate a 2D Cavity flow using LES model with wall functions.

# This code was written in a way that We can run it in both Windows and Linux environment.

# Discretization of the governing equations is based on finite volume method.

# SIMPLE algorithm is used to solve velocity-pressure coupling.

# Due to the higher stability, the convection term in Momentum equation, is discretized based on first order upwind.

# In order to solve u, v, and Pressure correction equations, Line by Line method is used in every iteration.

# Proper choice of under-relaxation factors is needed for convergence.

# Eddy viscosity model “Smagorinsky” is used to calculate Subgrid-scale viscosity. 

# The damping effect of walls on sub-grid viscosity is taken into account using Van Driest damping function.

# For computing walls shear stress, two wall functions called “Standard wall function” and “Enhanced wall function” is used to calculate effective wall viscosity for each wall boundary face.

# This program consists of four folders: 
  * A folder which is named "libs", includes all required classes to solve the problem
  * A folder which is named "solver", includes your main code file ("main.cpp"). There are also some “ .H” files included in “main.cpp”, which are used to somehow make the code clean.
  * A folder which is named "inputData", includes "data.txt" file to enter required data for solving the problem 
  * A folder which is named "post-processing", includes a python script named “plotResults.py” which is responsible for plotting all calculated results. (It is important to note that you should choose a time that exists in the “results” directory in “solver” folder)

# If you want to run this code, you should:
  1. Go to the "inputData" directiory and open the "data.txt" file
  2. Set input variables according to your desire
  3. Go to the "solver" directory
  4. Compile "main.cpp" using a C++ compiler
  5. Run the executable file that will be produced after compiling

# After running the code, you can do post-processing by:
  1. Going to the "post-processing" directory and open the "plotResults.py" 
  2. Setting the time you want to plot the results
  3. Run "plotResults.py" file using a python compiler
  
# You can run this code in different modes such as:
  * Laminar model
  * Smagorinsky LES model without considering wall damping effects
  * Smagorinsky LES model considering wall damping effects using Van Driest damping function
  * Calculate walls shear stress without wall functions
  * Calculate walls shear stress using Standard wall function
  * Calculate walls shear stress using Enhanced wall function
  
# It is important to note that when we have laminar flow, the solver automatically calculates walls shear stress without wall functions, as in laminar flow, 𝜈_w is the same as 𝜈, and there is no need to use wall function to calculate 𝜈_w . However, when we have turbulent flow, you can choose how to calculate walls shear stress. Without wall functions, Standard wall function, and Enhanced wall function are methods that you can use to calculate 𝜈_w.

# Another important functionality of this code is you can resume running the code from what time that you want. For example, imagine that first, you want to solve the problem until t=5s. At t=5s, the code stops running. Then, you realize that you want to solve the problem until t=10s. If you want to do this, you do not need to start solving the problem from t=0. You can resume your run from t=5s. Only thing that you need to do is just set a new time for Startime (in this example Startime=5) in “data.txt” file in “inputData” folder. It is important to note that you should choose a time whose results exist in “results” directory in “solver” folder, as u, v, and pressure are initialized from there.
