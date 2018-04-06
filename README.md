# rrlsq-matlab
relaxed regularized least squares in MATLAB

## set up and testing
In the top-level folder, run the setup.m and 
runalltests.m scripts. This should build 
and test the mex files used to expose the 
more efficient versions of the LAPACK QR
routines. To force recompilation, create a
variable recompile=1 before calling setup.

## demo
See the demo folder. The results can be 
turned into an HTML page using MATLAB's
publish feature.