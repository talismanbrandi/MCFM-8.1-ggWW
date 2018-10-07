######################################################
# HTCondor Submit Description File.                  #
# COMMON TEMPLATE FOR AFS                            #
# Next commands that can be added to submit files    #
# -- Ayan Paul -- July 2018 -- DESY, Hamburg --      #
######################################################
# Picking up the RUNNAME in a variable
if !defined FNAME
  FNAME = RUNNAME
endif

# Picking up the ID
ID                      = $(Cluster).$(Process)

# Define universe and transfer protocols
# NOTE: Be careful with getenv. 
# $HOME is not defined so we have to get it if we need it.
universe                = vanilla
should_transfer_files   = NO
getenv                  = True
home_dir                = $ENV(HOME)
logs_dir		= $(home_dir)/Github/MCFM-8.1-ggWW/WORK/logs

# Define directories and IO
# initialdir              = $(home_dir)/Github/MCFM-8.1-ggWW/WORK/logs
output                  = $(logs_dir)/$(FNAME).$(ID).out
error                   = $(logs_dir)/$(FNAME).$(ID).err
log                     = $(logs_dir)/$(FNAME).$(Cluster).log

# Set up notifications
notify_user             = ayan.paul@desy.de
notification            = Complete

# Hardware requests
request_cpus            = 12
request_memory          = 8 GB
requirements            = OpSysAndVer == "CentOS7"
+RequestRuntime         = 54000

# Set up the run with exe and arguments. "config" is a directory.
# This is an OMP run.
environment             = OMP_NUM_THREADS=12;OMP_STACKSIZE=16000
executable              = mcfm_omp
arguments               = ../WORK/$(RUNDIR) $(INPUTFILE)

# generate the queue from a list in another file
queue RUNDIR,INPUTFILE from run_list.txt
# FOR MORE DETIALS:
# https://indico.cern.ch/event/611296/contributions/2604376/attachments/1471164/2276521/TannenbaumT_UserTutorial.pdf
