———————————————————————————————————————

Single_Double_Pendulum

———————————————————————————————————————

The motion of a single or double pendulum are simulated by solving equations of motion through numerical methods.

———————————————————————————————————————

Most of the code for all of the files below are similar. The ones made for the report
are slightly modified to produce better images for the report.
For tests and more experiments, 'Project A.py' and 'Project A animated.py' are used.
'Project A animated.py' contains almost the same code as 'Project A.py',
but it is an animated version.

———————————————————————————————————————

How to use the files (they have been tested using Spyder in Blackett Computer Suite and Canopy in the Library):
- The variables which can be changed are under ###### CONFIGURATIONS ######
- Make approx = 'y' or 'n' to approximate sin(theta) = theta for single pendulum.
- Make method = method_list[n] where n is an integer from 0 to 3. This chooses an FDM from method_list.
- Type = type_list[m] where m is 0 or 1. This chooses whether the system is double or single pendulum.

- find_critical_step_length = 'y' or 'n' to choose whether to find critical step length.
- critical_step_start is the initial h value to start at.
- critical_step_end is the final h value to finish at.
- number_of_steps is the number of intervals to find the critical step length.
- energy_percentage is a percentage of initial energy. Explained in code and in report.
- theta_difference is a small value. Explained in code and report.

- animation_speed is the speed of the animation. The number is the time in milliseconds between each frame.
- h is step size
- tmax is maximum time value

-Rest should be self explanatory.

———————————————————————————————————————

The Python code that is used to create results for the report are:
Single Pendulum FDMs undamped.py
Single Pendulum FDMs critical step-size.py
Double Pendulum FDMs undamped.py

These codes are create results that are briefly talked about:
Single Pendulum FDMs damped.py
Double Pendulum FDMs damped.py

Some parts of code for the report are commented out to show which part of original code is not needed to produce the result.

———————————————————————————————————————

Contact details:

- kelvin.chan14@imperial.ac.uk