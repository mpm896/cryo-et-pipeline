# !/bin/bash

# Matthew Martinez
# Sanofi - US Dept. of Large Molecules Research, Protein Engineering Group, Structural Biology

# 08/09/2024 - So far this is only functional for IMOD recontruction with serieswatcher and batchruntomo
# TODO: 
#   1. Incorporate each dataset into a SQL database
#   2. Provide functionality with Warp
#   3. For the IMOD pipeline, allow for splitting datasets into half sets for denoising models (i.e. cryo-CARE)


# RUN FROM YOUR PROJECT DIRECTORY! Chances are, if starting a project, this directory should be empty
# Usage:
#   Fill out the parameters in the first section, below, then save
#   Run with: sh submit_imod_pipeline.sh

# Parent script to start the IMOD tomogram reconstruction pipeline 
# For more details on IMOD and the processes here, see: 
#   IMOD: https://bio3d.colorado.edu/imod/
#    Batchruntomo guide: https://bio3d.colorado.edu/imod/doc/batchGuide.html
#    Batchruntomo: https://bio3d.colorado.edu/imod/doc/man/batchruntomo.html
#    Motion correction:
#        framewatcher: https://bio3d.colorado.edu/imod/doc/man/framewatcher.html
#        alignframes: https://bio3d.colorado.edu/imod/doc/man/alignframes.html
#    Alignment and reconstruction:
#        serieswatcher: https://bio3d.colorado.edu/imod/doc/man/serieswatcher.html
#        Batchruntomo: see above


# 1. Script will first copy over all the raw frames from the Glacios PC, as this is a read-only filesystem
# 2. Copy all frames to a subdiretory called Frames (if collected in Tomography 5. If collected in serialEM, frames will already by in a subdirectory). 
# 3. Rename all the mdoc files according the serialEM/IMOD convention if collected in Tomography 5
# 4. Modify the mdoc files (in place) with the correct, user-specified tilt axis if collected in Thermo Tomography 5 instead of serialEM
#     - Thermo Tomography 5 and serialEM have different conventions for recording the tilt axis. 
#     - serialEM is the de facto standard in the field. Here, we are using Tomography 5. The conversion between the two tilt axis numbers are:
#         serialEM tilt axis = -90 - Tomo 5 tilt axis
#         - i.e. if the recorded Tomo 5 tilt axis is -87, then the serialEM tilt axis is -90 - -87 = -3
# 5. Run framewatcher to perform motion correction of frames and output a tiltseries made of motion corrected frames 
#       - This process runs the IMOD alignframes program. Cannot be submitted with nohup
#       - Optional normalized dose weighting of frames
# 6. Run serieswatcher at the same time via a python script db_reconstruct.py, to watch for motion corrected tilt series to run automated reconstruction
# 7. Transfer to the database upon completion of tomogram reconstruction via a python script db_reconstruct.py


# ====================== #
# ENTER PARAMETERS BELOW
# ====================== #

# IMPORTANT NOTE:
#   - WHEN ASSIGNING VARIABLE, DO NOT ADD SPACES AROUND THE = SIGN
#     BASH DOES NOT LIKE THIS!

# This is most likely being run at Sanofi. This parameter lets the script know if the Pipe CLI is configured:
PIPE_CLI=1          # 0: Pipe CLI is not configured, do not use pipe commands
                    # 1: Pipe CLI is configured (if using a Magellan instance, it is most likely configured!)

# Collection software - Important for interpretation of mdocs and tilt axis
SOFTWARE=1          # 0: serialEM
                    # 1: Thermo Scientific Tomography 5

# Tilt axis from mdoc file
TILTAXIS=-87    # Leave as "EMPTY" if data was collected with serialEM


# *** SETUP SETTINGS ***
TRANFER_RAW_DATA=0          # Whether or not to transfer raw data from the microscope PC (i.e. you already did the transfer)
                                # 0: Do not do data transfer
                                # 1: Do the data transfer
RAW_DATA_DIR=/cloud-data/its-cmo-darwin-magellan-workspaces-folders/WS_Cryoem/CX_LMR/Project_directories/cryo-et-test  # Path to raw frames
FRAMES_NAME=Fractions       # Common word in the filenames of ALL THE FRAMES
                                #   - If collected in Tomography 5 with dose fractionation, frames should have "Fractions" in their name
                                #   - If collected in serialEM, not as important since you should have been set to save raw frames in a "Frames" subdirectory
MDOC_DUPLICATE=override     # This is used in the case of Tomography 5 software. For some reason it creates duplicate mdocs containing the word "override"
                                # Do not copy over these mdocs
READ_MDOC=1                 # 0: Do not read mdoc
                                # 1: Read values from the mdoc
GAIN_PATH="EMPTY"           # Path to gain reference file.
                                # If raw frames already gain corrected, set to "EMPTY" (hint: MRC are likely gain corrected, TIF are likely NOT gain corrected)
EXTENSION=*.mrc         # Extension of the frames. Options:
                            #    *.tif
                            #    *.mrc
                            #    *.eer
PIXEL_SIZE=1.19      # Enter pixel size, or leave as "EMPTY" to read from the mdoc
EXPOSURE=2.34        # Enter the exposure, per tilt angle (i.e 2.34). Leave as "EMPTY" to read from mdoc

# *** MOTION CORRECTION WITH FRAMEWATCHER AND IMOD ALIGNFRAMES *** 
RUN_FRAMEWATCHER=1  # 0: Do not do motion correction (i.e. if you already motion corrected your frames)
                        # 1: Do motion correction with framewatcher
DO_MC_DOSEWEIGHT=0	# Do normalized dose weighting during alignframes motion correction. ONLY run this if the dosage in the mdoc looks normal. If the dose is too high (i.e. total dose over 150 or 200), then the output aligned frames will be messed up
			            # 0: Do not do dose weighting
			            # 1: Do dose weighting
DROP_MEAN=100		# Mean signal, below which frames will be dropped and excluded during motion correction. May need to change this to a lower number for cellular targets - not sure

# *** COMMAND FILE SETTINGS - RECONSTRUCTION ***
CPUS=4                  # Number of CPUs to use for parallel processing
GPUS=1                  # Number of GPUs

# *** PREPROCESSING ***
REMOVE_XRAYS=1          # 0: Do not remove X-rays
                        # 1: Remove X-rays/bad pixels (recommended)

# *** COARSE ALIGNMENT ***
PREALIGN_BIN=4          # Bin factor for dataset DURING alignment

# *** TRACKING - GENERAL ***
TRACK_METHOD=1          # 0: Fiducial tracking
                        # 1: Patch tracking
SIZE_GOLD=0             # Size of gold fiducials, in nm. Only need to enter if doing fiducial tracking. Automatically set to 0 if doing patch tracking

# *** PATCH TRACKING ***
PATCH_SIZE_X=200        # Size of patch in X (pixels)
PATCH_SIZE_Y=200        # Size of patch in Y (pixels)
PATCH_OVERLAP_X=0.4     # Fraction overlap of patches in X
PATCH_OVERLAP_Y=0.4     # Fraction overlap of patches in Y

# *** FIDUCIAL TRACKING ***
NUM_BEADS=25            # Number of beads to seed automatically. 25 should be more than enough.
USE_SOBEL=1             # Use the sobel filter to find center of beads. This should be used for cryo data
                        # 0: Do not use sobel filter
                        # 1: Use sobel filter
SOBEL_SIGMA=1.5         # Sigma for sobel filter in binned pixels. Set to 1.5 for cryo data

# *** FINAL ALIGNMENT ***
FINAL_BIN=4             # Bin factor for the output tomogram
DO_CTF=1                # 0: Don't do CTF correction
                        # 1: Do CTF correction
DO_DOSE_WEIGHTING=1     # 0: Don't do dose weighting
                        # 1: Do dose weighting (note: This is important if trying to determine structures to high resolution)
                        #   - For doseweighting, the dose will be taken from the MDOC

# *** FINAL ALIGNMENT - CTF CORRECTION ***
DEFOCUS_RANGE_LOW=1000.0    # Low end of the defocus range to scan, in nm underfocus
DEFOCUS_RANGE_HIGH=10000.0   # High end of the defocus range to scan, in nm underfocus
AUTOFIT_RANGE=12            # Degree range of tilts to average for fitting the CTF (investigate before changing)
AUTOFIT_STEP=1              # Increment for ranges for autofitting
                            #   - Together, autofit range and step mean:
                            #   - If collected tilt series in 3 degree increments, it will average 12 degrees worth of tilts (i.e. 0 degree tilt to 12 degree tilt -> 5 tilt angles)
                            #   - For the next fitting, it will increment by 1 tilt and repeat (i.e. 3 degree tilt to 15 degree tilt -> 5 tilt angles), and so on
TUNE_FITTING_SAMPLING=1     # 0: Do not tune the fitting and sampling during CTF
                            # 1: Tune the fitting and sampling during CTF
                            #   - I don't know if this is actually correct because this parameter is not described on the IMOD directives page. But I know it should be 1                    
VOLTAGE=200                 # Voltage of the microscoe, in keV
CS=2.7                      # Spherical abberation, in mm

# *** FINAL ALIGNMENT - DOSE WEIGHTING ***
DOSE_SYM=0                  # 0: Dose-symmetric data collection
                            # 1: Bidiretional data collection? This info is not on the batchruntomo directives page

# *** RECONSTRUCTION PARAMETERS ***
RECONSTRUCT_METHOD=1        # 0: R-weighted backprojection
                            # 1: R-weighted backprojection with SIRT-like filter
                            # 2: SIRT
FAKE_SIRT_ITERS=10          # SIRT-like radial filter equivalent to given # of SIRT iterations
                            #   - Only necessary of RECONSTRUCT_METHOD=1 (SIRT-like filter)
SIRT_ITERS=10               # Number of SIRT iterations. Don't worry about this if not using RECONSTRUC_METHOD=2 (SIRT)
THICKNESS_UNBINNED=1200     # Thickness of the final tomogram in UNBINNED pixels
                            #   - Reference:
                            #       single particle / nanoparticle sample: ~150 nm is more than enough (with pixel size of 1.19 A, this is 1260 unbinned pixels)
                            #       thick cellular sample: 630 - 850 nm may be necessary (with pixel size of 2.65 A, this is 2400 - 3200 unbinned pixels)
THICKNESS_BINNED="EMPTY"    # Thickness of the final tomogram in BINNED pixels. If a number is input here, it will override THICKNESS_UNBINNED. Otherwise leave as "EMPTY"

# *** POST PROCESSING ***
DO_TRIMVOL=1                # Run trimvol after reconstruction. Should probably keep this as 1
                            #   0: Don't run trimvol
                            #   1: Run trimvol
REORIENT=2                  # Reorientation to do during trimvol
                            #   0: None
                            #   1: Flip
                            #   2: Rotate around X


# *** DENOISING AND MISSING WEDGE CORRECTION ***
DO_DENOISING=1		    # 0: Do not setup files for denoising
			            # 1: Setup files for denoising and missing wedge correction with DeepDeWedge

USER_DB_ID=1

# Scripts for tomogram reconstruction and data transfer
PY_RECONSTRUCT=/cloud-data/its-cmo-darwin-magellan-workspaces-folders/WS_Cryoem/CX_LMR/Project_directories/cryo-et-pipeline/db_reconstruct.py
PY_TRANSFER=/cloud-data/its-cmo-darwin-magellan-workspaces-folders/WS_Cryoem/CX_LMR/Project_directories/cryo-et-pipeline/db_transfer.py
PY_HALFTOMO=/cloud-data/its-cmo-darwin-magellan-workspaces-folders/WS_Cryoem/CX_LMR/Project_directories/cryo-et-pipeline/halftomo_reconstruct.py
PY_DDW=/cloud-data/its-cmo-darwin-magellan-workspaces-folders/WS_Cryoem/CX_LMR/Project_directories/cryo-et-pipeline/pyDDW.py




# ================================================
# ***  DO NOT MODIFY ANYTHING BELOW THIS POINT ***
# ================================================
function typewriter
{
    text="$1"
    delay="$2"

    for i in $(seq 0 $(expr length "${text}")) ; do
        echo -n "${text:$i:1}"
        sleep ${delay}
    done
    echo
}

if [[ $RECONSTRUC_METHOD -eq 2 ]]; then
    DO_SIRT=1
else
    DO_SIRT=0
fi

DATA_DIR=$(basename $RAW_DATA_DIR)
if [ ! -d $DATA_DIR ]; then
    mkdir -p $DATA_DIR/Frames
fi

if [[ $TRANFER_RAW_DATA -eq 1 ]]; then
    echo "==================================="
    typewriter "Copying data from the Glacios PC..." 0.02
    echo "==================================="

    # Setup tmux session
    # TRANSFER_PROC=transfer
    # tmux new-session -d -s $TRANSFER_PROC

    # Move everything from the Glacios PC
    if [[ $PIPE_CLI -eq 1 ]]; then
        MAGELLAN_DIR=$(echo $RAW_DATA_DIR | sed 's/.*\(its\)/\1/g')
        pipe storage cp -r --force --skip-existing cp://$MAGELLAN_DIR $DATA_DIR
        # tmux send-keys "pipe storage cp -r --force --skip-existing cp://$MAGELLAN_DIR $DATA_DIR" C-m
    elif [[ $PIPE_CLI -eq 0 ]]; then
        rsync --progress --ignore-existing -avr $RAW_DATA_DIR/ $DATA_DIR
        # tmux send-keys "rsync --progress --ignore-existing -avr $RAW_DATA_DIR/ $DATA_DIR" C-m
    fi
    # typewriter "See progress of data transfer with: tmux a -t $TRANSFER_PROC" 0.05
fi

# Perform a series of operations if collected in Tomography 5
if [[ $SOFTWARE -eq 0 ]]; then
    typewriter "===== Data collected with serialEM =====" 0.02
elif [[ $SOFTWARE -eq 1 ]]; then
    typewriter "===== Data collected with Thermo Scientific Tomography 5 =====" 0.02
    TILTAXIS=$(echo "scale=2; -90.0 - $TILTAXIS" | bc)
    typewriter " -> New tilt axis is ${TILTAXIS}" 0.02

    # Renamd all the mdocs to take one the extension ".mrc.mdoc"
    typewriter "===== Renaming the mdocs =====" 0.02
    find $DATA_DIR -name "*mdoc" -type f -exec sh -c '
    for f; do
        if [[ $f != *.mrc.mdoc ]]; then
            mv "$f" "${f/.mdoc/.mrc.mdoc}"
        fi
    done' sh {} +

    # Remove potential duplicate MDOCs
    find $DATA_DIR -iname "*$MDOC_DUPLICATE*" -exec rm {} +

    # Modify the mdocs to display the new tilt axis
    typewriter "===== Adjusting the tilt axis in the mdoc =====" 0.02
    for file in "$DATA_DIR/*.mrc.mdoc"; do
        head $file | sed -i -E "/TiltAxisAngle *= */s/-?[0-9]+\.[0-9]+/$TILTAXIS/" $file
    done
    
    # Copy all the frames to a subdirectory called "Frames"
    echo "================================"
    typewriter "Moving frames to subdirectory..." 0.02
    echo "================================"
    mv "$DATA_DIR/*$FRAMES_NAME*" $DATA_DIR/Frames
else
    echo "SOFTWARE=${SOFTWARE} is not a valid option."
fi

# ====================================
# *** SETTING UP THE COMMAND FILES ***
# ====================================
# MASTER FILES CAN BE FOUND AT: /PATH/TO/FILES
# The following commands to produce command files are based off those files

typewriter "========== CONSTRUCTING COMMAND FILES ==========" 0.02

# FIX STUFF WITH DATA_DIR, WATCH_DIR, ETC
WATCH_DIR=$DATA_DIR
MASTER_FRAMEWATCHER_COM=coms/FRAMEWATCHER_MASTER.pcm
OUT_DIR=Aligned
PROC_DIR=$WATCH_DIR/Processed
THUMB_DIR=$OUT_DIR/alignedJPG

if [[ ! -d $OUT_DIR ]]; then
    mkdir -p $THUMB_DIR
fi

if [[ $PIXEL_SIZE == "EMPTY" ]]; then 
	PIX_SIZE=$(header $(ls $WATCH_DIR/$EXTENSION | head -1) | grep Pixel | awk 'NF>1{print $NF}')
fi 

if [ ! -d coms ]; then
    mkdir coms
fi

# Master command file for framewatcher - alignframes
echo "\$alignframes -StandardInput
OutputImageFile
AdjustAndWriteMdoc
PathToFramesInMdoc $(pwd)/$DATA_DIR/Frames
UseGPU 0
PairwiseFrames 7
GroupSize 1
AlignAndSumBinning 10 1
AntialiasFilter 4
RefineAlignment 0
StopIterationsAtShift 0.000000
ShiftLimit 20
MinForSplineSmoothing 0
FilterRadius2 0.050000
FilterSigma2 0.007145
VaryFilter 0.050000
ScalingOfSum 16.000000" > $MASTER_FRAMEWATCHER_COM


# ===================
# *** SUBMIT JOBS ***
# ===================

# Kill all processes to start from a clean slate
pkill -9 -f framewatcher
pkill -9 -f serieswatcher

FW_PIPELINE=fw_pipeline 
BRT_PIPELINE=brt_pipeline
DB_PIPELINE=db_pipeline

# SUBMIT FRAMEWATCHER JOB

if [[ $RUN_FRAMEWATCHER -eq 1 ]]; then
    typewriter "===== RUNNING FRAMEWATCHER =====" 0.02
    typewriter "See progress with 'tmux a -t $FW_PIPELINE'" 0.02

    tmux new-session -d -s $FW_PIPELINE

    if [[ $DO_MC_DOSEWEIGHT -eq 0 ]]; then 
        tmux send-keys "framewatcher -watch '$WATCH_DIR' -master '$MASTER_FRAMEWATCHER_COM' -o '$OUT_DIR' -pr '$PROC_DIR' -gpu 0 -bin 1 -po 1024 -thumb '$THUMB_DIR' -mdrop '$DROP_MEAN'" C-m 
    elif [[ $DO_MC_DOSEWEIGHT -eq 1 ]]; then 
        tmux send-keys "framewatcher -watch '$WATCH_DIR' -master '$MASTER_FRAMEWATCHER_COM' -o '$OUT_DIR' -pr '$PROC_DIR' -gpu 0 -bin 1 -po 1024 -thumb '$THUMB_DIR' -dmdoc -dnorm -volt 200 -ddrop 0.1,1.5" C-m
    else
        echo "$DO_MC_DOSEWEIGHT is not a valid choice for dose weighting during motion correction."
        exit 1
    fi 
fi

# SUBMIT SERIESWATCHER JOB 
typewriter "===== RUNNING SERIESWATCHER =====" 0.02
typewriter "===== See progress with 'tmux a -t $BRT_PIPELINE' =====" 0.02

tmux new-session -d -s $BRT_PIPELINE
tmux send-keys "python $PY_RECONSTRUCT $CPUS $GPUS $OUT_DIR $READ_MDOC $REMOVE_XRAYS $PREALIGN_BIN $TRACK_METHOD \
    $SIZE_GOLD $FINAL_BIN $DO_SIRT $DO_TRIMVOL $PIXEL_SIZE $TILTAXIS $DOSE_SYM $VOLTAGE $CS $REORIENT \
    $THICKNESS_BINNED $THICKNESS_UNBINNED $USE_SOBEL $NUM_BEADS $SOBEL_SIGMA $PATCH_SIZE_X $PATCH_SIZE_Y \
    $PATCH_OVERLAP_X $PATCH_OVERLAP_Y $DO_CTF $DEFOCUS_RANGE_LOW $DEFOCUS_RANGE_HIGH $AUTOFIT_RANGE $AUTOFIT_STEP \
    $TUNE_FITTING_SAMPLING $FAKE_SIRT_ITERS" C-m

# SUBMIT DATABASE TRANSFER JOB  
typewriter "===== WATCHING FOR COMPLETED JOBS TO TRANSFER TO DATABASE =====" 0.02
typewriter "===== See progress with 'tmux a -t $DB_PIPELINE' =====" 0.02

SUBFRAME_DIR=$DATA_DIR/Frames
tmux new-session -d -s $DB_PIPELINE
tmux send-keys "python $PY_TRANSFER $OUT_DIR $SUBFRAME_DIR $EXTENSION $USER_DB_ID" C-m

echo "See all active processes in tmux sessions: tmux ls"
echo "View a specific tmux session: tmux a -t session_name"
echo "Detach the session (to regain control of the terminal): Ctrl-b d"
echo "Kill a specific tmux session: tmux -t session_name kill-session"
echo "Kill all tmux sessions: tmux kill-server"


if [[ $DO_DENOISING -eq 1 ]]; then
    echo
    typewriter "Waiting to begin half-tomo generation for denoising..." 0.02
    echo

    # Watch for serieswatcher child processes. If no more after some time, move forward with denoising (if selected)
    num_sw=$(pgrep -c -f serieswatcher)

    # Don't want the denoising to start until some processes have began
    while [[ $num_sw -lt 2 ]]; do
        echo "No new reconstruction processes. Sleeping..."
        sleep 30s
        num_sw=$(pgrep -c -f serieswatcher)
    done

    # Now, want to wait until serieswatcher is done
    while [[ $num_sw -gt 1 ]]; do
        echo "Serieswatcher is running. Sleeping..."
        sleep 30s
        num_sw=$(pgrep -c -f serieswatcher)
    done

    DONE_DIR=$OUT_DIR/Done
    DN_PROC=dn_proc

    tmux new-session -d -s $DN_PROC
    tmux send-keys "python $PY_HALFTOMO 1 $DONE_DIR" C-m
fi


