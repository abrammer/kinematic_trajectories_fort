README for trajectory program

A somewhat dumb fortran kinematic trajectory program I wrote in a couple days many years ago. 
Backend is all there to handle multiple WRF files, reading data and tracing trajectories. 
Could/Should be adapted to take input from python or NCL to make the file-io simpler and more generic. 

Runs based on a namelist file for file definitions and trajectory options.  

Utilises bicubic interpolation in the horizontal and polynomial in the vertical, with time-integration following Petterssen (similar to RK). 

