#!/usr/bin/env python
# Script for calibrating beta using MPAS 
# Note that the script currenly assumes that the target velocity is an observed surface velocity, as opposed to a balance velocity - but that distinction may not matter too much.
# Matt Hoffman, March 2013

################ SET THESE FOR YOUR SETUP ###########
runcmd = 'time openmpirun -n 4 land_ice_model.exe'   # command to run MPAS on your system
vobsfile='observedSpeed.nc'  # This is the file that has the observed surface velocity field - it can have anything else, too.  
vobsvariable='observedSpeed'  # this is the name of the observed velocity field variable in the above file.  It should have dimensions (Time, nCells)
numiters = 5            # Number of iterations to run the algorithm.  the first time through uses a constant no-sliding beta to get the first guess of the deformation velocity
#####################################################

# These are probably fine.
infile='land_ice_grid.nc'
outfile='output.nc'
modelfiles='namelist.input graph.info.part* dataIce *xml ' + infile   # land_ice_model.exe gets symlinked


# Import modules, get command line input ###########
import os, sys
import subprocess
import numpy as np
#import matplotlib.pyplot as plt
try:
  from Scientific.IO.NetCDF import NetCDFFile
  netCDF_module = 'Scientific.IO.NetCDF'
except ImportError:
  try:
    from netCDF4 import Dataset as NetCDFFile
    netCDF_module = 'netCDF4'
  except ImportError:
      print 'Unable to import any of the following python modules:'
      print '  Scientific.IO.NetCDF \n  netcdf4 '
      print 'One of them must be installed.'
      raise ImportError('No netCDF module found')
from optparse import OptionParser

print "** Gathering information.  (Invoke with --help for more details. All arguments are optional)"
parser = OptionParser(epilog='You should set the obs file and variable names (and run command if you will have the script run MPAS) in the script before running.  Your namelist file should be set to do only a diagnostic solve (no time stepping).  Currently this is setup to tune to observed surface velocities (as opposed to balance velocities).')
parser.add_option("-s", "--start", dest="start", help="start iteration to perform, first iteration is 0.  REQUIRED.", metavar="START_ITERATION")
parser.add_option("-e", "--end", dest="end", help="end iteration to perform.  To perform a single iteration, give START_ITERATION, or don't include this flag and it will default to START_ITERATION", metavar="END_ITERATION")
parser.add_option("-r", "--run", action="store_true", dest="run", help="include this flag to run the model from within the script (useful for non-queue machines).  Default behavior is to NOT run the model.", metavar="RUN")
parser.add_option("-b", "--balance", action="store_true", dest="balancevel", help="include this flag if you are using balance velocities.  Default behavior is to assume observed surface velocities.", metavar="BALANCEVEL")
options, args = parser.parse_args()


print 'You should set the obs file and variable names (and run command if you will have the script run MPAS) in the script before running.  Your namelist file should be set to do only a diagnostic solve (no time stepping).  Currently this is setup to tune to observed surface velocities (as opposed to balance velocities).'

if not options.start:
    sys.exit('Start time is required.  Invoke with -h for more information.')
if not options.end:
    options.end = options.start
if not options.run and options.end > options.start:
   sys.exit('If you are not having the script run the model, then you can only perform one iteration at a time.')
if options.balancevel:
   print 'Assuming balance velocities are being supplied.'
else:
   print 'Assuming surface velocities are being supplied.'


# ==================
# Setup a function
def check_results(outfile, vobsfile, vobsvariable, balancevel):
  ###### check the results
  fo=NetCDFFile(outfile, 'r')
  thickness=fo.variables['thickness'][0,:]
  ux=fo.variables['uReconstructX'][0,:,:]  # use the first time level 
  uy=fo.variables['uReconstructY'][0,:,:]  # use the first time level 
  if balancevel:
    vModel = ((ux[:,0]**2 + uy[:,0]**2) ** 0.5).mean(1)  # get vertically averaged velocity
  else:
    vModel = (ux[:,0]**2 + uy[:,0]**2) ** 0.5  # surface velocity (top vertical layer)

  cellMask=fo.variables['cellMask'][0,:]
  ind = ( (cellMask & 4) == 4 )             # 4 bit = grounded ice  # this gets the indices that are grounded ice.

  fv = NetCDFFile('../' + vobsfile,'r')
  vOb = fv.variables[vobsvariable][0,:]  

  rmse = ( ((vModel[ind] - vOb[ind])**2).sum() / ind.sum() )**0.5   # TODO  check if this is calculating right  
  print "RMSE="+str(rmse)
  fo.close()
  fv.close()

# ===================




# If we aren't running the model, then first calculate the results from the previous iteration (if there is one)
if not options.run and options.start>0:
  # Show results from previous run
  check_results('iter{0:3d}'.format(options.start-1) + '/' + outfile, vobsfile, vobsvariable)


# Loop and iterate these steps - 
# (add relaxation or smarter guess later)
for i in range(options.start, options.end+1):
    dirname = 'iter{0:3d}'.format(i)
    print 'starting iteration ' + str(i) + ' in directory ' + dirname

    # remove previous attempts
    os.system('rm -rf ' + dirname )

    # Setup the directory for the first run
    os.system('mkdir ' + dirname)
    os.system('cp ' + modelfiles + ' ' + dirname)
    os.chdir(dirname)
    if options.run:
       os.system('ln -s ../land_ice_model.exe .')
    print 'Working in ' + os.getcwd()

    ########### Calculate a guess for beta ########
    if i == 0:
        print 'iteration 0 is an initial guess!  Using beta=1e8 as the initial guess.'
        # Make an initial guess for beta
        fi=NetCDFFile(infile,'r+')
        fi.variables['betaTimeSeries'][:]=1.0e8
        fi.close()
        
    else:
        ########## Get Taub ###############
        # Read the previous input, output files
        fo_old=NetCDFFile('../iter{0:3d}'.format(i-1) + '/' + outfile, 'r')
        fi_old=NetCDFFile('../iter{0:3d}'.format(i-1) + '/' + infile,  'r')
        # Get ub and beta
        beta_old=fi_old.variables['betaTimeSeries'][:,0]  # just use the first time level from the input file since beta isn't in the output  (the time level chosen shouldn't matter)
        #cellMask=fi.variables['cellMask'][0,:]  # need to treat floating ice separately?  i guess not since it gets handled independently...
        # use the reconstructed velocities because beta is on cell centers
        ux=fo_old.variables['uReconstructX'][0,:,:]  # use the first time level 
        uy=fo_old.variables['uReconstructY'][0,:,:]  # use the first time level 
        # Calculate previous model's Taub=beta*ub since we aren't outputting Taub
        ub_old = (ux[:,-1]**2 + uy[:,-1]**2) ** 0.5  # sliding velocity - bottom vertical layer
        Taub = beta_old * ub_old
        if options.balancevel:
          ud_old = ((ux[:,:]**2 + uy[:,:]**2) ** 0.5).mean(1) - ub_old  # depth-averaged deformation velocity = depth-averaged velocity - sliding velocity
        else:
          ud_old = (ux[:,0]**2 + uy[:,-0]**2) ** 0.5 - ub_old  # deformation velocity at surface = surface velocity - sliding velocity

        # calculate what the sliding velocity target is now that we have deformation velocity
        fv = NetCDFFile('../' + vobsfile,'r')
        vOb = fv.variables[vobsvariable][0,:]  
        ubTarget = vOb - ud_old  # target sliding = observed velocity - modeled deformation  (see above for how calculation of ud_old)
        # check 0 and negative targets to something small so we always have some sliding
        ubTarget[ubTarget<0.1] = 0.1

        # Calculate new beta guess:   beta = Taub_M / ub_T   M=modeled, T=target
        betaguess = Taub / ubTarget
        betaguess[betaguess>1.0e8] = 1.0e8   # Cap how big beta can be
        fi=NetCDFFile(infile,'r+')
        fi.variables['betaTimeSeries'][:,0] = betaguess[:]

        # use the old velo soln as the guess for the next run
        fi.variables['normalVelocity'][:] = fo_old.variables['normalVelocity'][:]   # TODO make sure this works properly

        fi.close()
        fi_old.close()
        fo_old.close()
        fv.close()


    if options.run:
        ########### Run the model - diagnostic solve
        print 'Starting model run'
        try:
            subprocess.check_call(runcmd, shell=True)
        except:
            sys.exit('model run failed. exiting...')

        ###### check the results now that we have run the model
        check_results(outfile, vobsfile, vobsvariable)

        os.chdir('..')
    else:
        print 'Please run the model now in directory ' + dirname




print 'Completed all requested iterations.'



