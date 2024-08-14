# Matthew Martinez
# Sanofi - US Dept. of Large Molecules Research, Protein Engineering Group, Structural Biology

# Reconstruct tomograms during pipeline processing
# This script will be passed the data processing parameters by the main shell script 

import contextlib
import os
from pathlib import Path
import sys
import subprocess
import time
from typing import Generator


def setup_serieswatcher(**kwargs) -> tuple[Path, Path]:
    """
    Construct the COM and ADOC files to to run serieswatcher

    :param: **kwargs
        Pipeline arguments passed from shell script

    :return: tuple of files
        master com file, master adoc file
    """
    assert 'cpus' in kwargs
    assert 'gpus' in kwargs
    assert 'out_dir' in kwargs
    assert 'remove_xrays' in kwargs
    assert 'read_mdoc' in kwargs
    assert 'remove_xrays' in kwargs
    assert 'prealign_bin' in kwargs
    assert 'track_method' in kwargs
    assert 'size_gold' in kwargs
    assert 'final_bin' in kwargs
    assert 'do_sirt' in kwargs
    assert 'do_trimvol' in kwargs
    assert 'pixel_size' in kwargs
    assert 'tiltaxis' in kwargs
    assert 'dose_sym' in kwargs
    assert 'voltage' in kwargs
    assert 'cs' in kwargs
    assert 'reorient' in kwargs
    assert 'thickness_binned' in kwargs
    assert 'thickness_unbinned' in kwargs
    assert 'use_sobel' in kwargs
    assert 'num_beads' in kwargs
    assert 'sobel_sigma' in kwargs
    assert 'patch_size' in kwargs
    assert 'patch_overlap' in kwargs
    assert 'do_ctf' in kwargs
    assert 'defocus_range'  in kwargs
    assert 'autofit_range'  in kwargs
    assert 'autofit_step' in kwargs
    assert 'tune_fitting_sample' in kwargs
    assert 'fake_sirt_iters' in kwargs

    # Get those parameters which can be read from the mdoc - pixel size, exposure, tilt angles
    mdoc_info = read_mdoc(get_mdoc(kwargs['out_dir']))

    master_com = Path.cwd() / Path('coms/BRT_MASTER.com')
    master_com.parent.mkdir(parents=True, exist_ok=True)
    master_adoc = Path.cwd() / Path('coms/BRT_MASTER.adoc')
    master_adoc.parent.mkdir(parents=True, exist_ok=True)

    if kwargs['pixel_size'] == "EMPTY":
        kwargs['pixel_size'] == mdoc_info['Pixel Size']
    if kwargs['thickness_binned'] != "EMPTY":
        kwargs['thickness_unbinned'] = str(
            int(kwargs['thickness_binned']) * int(kwargs['final_bin'])
        )

    # Construct the com file
    master_com.write_text(
        '$batchruntomo -StandardInput\n'
        'NamingStyle     1\n'
        'MakeSubDirectory\n'
        f'CPUMachineList  localhost:{kwargs["cpus"]}\n'
        f'GPUMachineList  {kwargs["gpus"]}\n'
        'NiceValue       15\n'
        'EtomoDebug      0\n'
        f'DirectiveFile   {master_adoc}\n'
        f'CurrentLocation {kwargs["out_dir"]}\n'
        'BypassEtomo\n'
    )
    
    # Construct the adoc file
    master_adoc.write_text(
        'setupset.systemTemplate = /usr/local/IMOD/SystemTemplate/cryoSample.adoc\n'
        f'runtime.Preprocessing.any.removeXrays = {kwargs["remove_xrays"]}\n'
        f'comparam.prenewst.newstack.BinByFactor = {kwargs["prealign_bin"]}\n'
        f'runtime.Fiducials.any.trackingMethod = {kwargs["track_method"]}\n'
        f'setupset.copyarg.gold = {kwargs["size_gold"]}\n'
        f'runtime.AlignedStack.any.binByFactor = {kwargs["final_bin"]}\n'
        f'runtime.Reconstruction.any.useSirt = {kwargs["do_sirt"]}\n'
        'runtime.Trimvol.any.scaleFromZ = \n'
        f'runtime.Postprocess.any.doTrimvol = {kwargs["do_trimvol"]}\n'
        f'setupset.copyarg.pixel = {kwargs["pixel_size"]}\n'
        f'setupset.copyarg.rotation = {kwargs["tiltaxis"]}\n'
        f'setupset.copyarg.dosesym = {kwargs["dose_sym"]}\n'
        f'setupset.copyarg.voltage = {kwargs["voltage"]}\n'
        f'setupset.copyarg.Cs = {kwargs["cs"]}\n'
        'comparam.prenewst.newstack.AntialiasFilter = 4\n'
        'comparam.newst.newstack.AntialiasFilter = 4\n'
        f'runtime.Trimvol.any.reorient = {kwargs["reorient"]}\n'
        f'comparam.tilt.tilt.THICKNESS = {kwargs["thickness_unbinned"]}\n'
    )

    if int(kwargs['track_method']) == 0:
        # Fiducial tracking
        with master_adoc.open("a") as f:
            f.write(
                'runtime.Fiducials.any.seedingMethod = 1\n'
                f'comparam.track.beadtrack.SobelFilterCentering = {kwargs["use_sobel"]}\n'
                f'comparam.autofidseed.autofidseed.TargetNumberOfBeads = {kwargs["num_beads"]}\n'
            )
        
        if int(kwargs['use_sobel']) == 1:
            with master_adoc.open("a") as f:
                f.write(
                    f'comparam.track.beadtrack.KernelSigmaForSobel = {kwargs["sobel_sigma"]}\n'
                )
    elif int(kwargs['track_method']) == 1:
        # Patch tracking
        with master_adoc.open("a") as f:
            f.write(
                f'comparam.xcorr_pt.tiltxcorr.SizeOfPatchesXandY = {kwargs["patch_size"][0]},{kwargs["patch_size"][1]}\n'
                f'comparam.xcorr_pt.tiltxcorr.OverlapOfPatchesXandY = {kwargs["patch_overlap"][0]},{kwargs["patch_overlap"][1]}\n'
            )
    else:
        raise ValueError(f"Tracking method of {kwargs['track_method']} is not supported")
    
    if int(kwargs['do_ctf']) == 1:
        with master_adoc.open("a") as f:
            f.write(
                f'runtime.AlignedStack.any.correctCTF = {kwargs["do_ctf"]}\n'
                f'comparam.ctfplotter.ctfplotter.ScanDefocusRange = {kwargs["defocus_range"][0]},{kwargs["defocus_range"][1]}\n'
                f'runtime.CTFplotting.any.autoFitRangeAndStep = {kwargs["autofit_range"]},{kwargs["autofit_step"]}\n'
                'comparam.ctfplotter.ctfplotter.BaselineFittingOrder = 4\n'
                'comparam.ctfplotter.ctfplotter.SearchAstigmatism = 1\n'
            )

    if int(kwargs['do_sirt']) == 0:
        with master_adoc.open("a") as f:
            f.write(
                f'comparam.tilt.tilt.FakeSIRTiterations = {kwargs["fake_sirt_iters"]}'
            )
            
    return master_com, master_adoc

def get_mdoc(p: Path) -> Path:
    """ Get mdoc files for each processing directory """
    while True: 
        mdocs = [x for x in list(p.rglob('*.mdoc'))]
        if mdocs:
            break
        time.sleep(60)
    return mdocs[0]


def read_mdoc(mdoc: Path) -> dict[str, any]:
    """ 
    Get additional information from mdocs, including:
        mag, pixel size, defocus, tilt angles, tilt increment

    For every mrc file, return a dict like:
    mdoc_dict = {
        "mag": [int] mag,
        "pixel size": [float] pixSize,
        "defocus": [list[float]] defocus,
        "defocus avg": [float] avg_defocus,
        "tilt angles": [list[float]] [tilt angles],
        "tilt min": [float] tilt_min,
        "tilt max": [float] tilt_max
        "tilt increment": [int] tilt increment
    }
    """
    assert mdoc.exists()

    header_info = {}
    # Open each file and extract the relevant information
    with open(mdoc, 'r') as f:
        header_info['Tilt Angles'] = []
        header_info['Defocus'] = []
        for line in f:
            strip_line = line.rstrip()
            if 'TiltAngle' in strip_line:
                index = strip_line.find('=')
                angle = round(float(strip_line[index+2:]))
                header_info['Tilt Angles'].append(angle)
            if 'Defocus' in strip_line and 'Target' not in strip_line:
                index = strip_line.find('=')
                header_info['Defocus'].append(float(strip_line[index+2:]))
            if 'Magnification' in strip_line:
                index = strip_line.find('=')
                header_info['Magnification'] = strip_line[index+2:]
            if 'PixelSpacing' in strip_line:
                index = strip_line.find('=')
                header_info['Pixel Size'] = str(round(float(strip_line[index+2:]),2)/10)

    header_info['Tilt Min'] = min(header_info['Tilt Angles'])
    header_info['Tilt Max'] = max(header_info['Tilt Angles'])
    header_info['Tilt Step'] = round(abs(
        (header_info['Tilt Max'] - header_info['Tilt Min'])
        / len(header_info['Tilt Angles'])
    )) 
    header_info['Defocus Avg'] = round(
        sum(header_info['Defocus']) / float(
            len(header_info['Defocus'])
        ), 2)
    
    return header_info


@contextlib.contextmanager
def chdir(path: str | Path) -> Generator[None, None, None]:
    """ Changed working directory and returns to the previous on exit """
    prev_cwd = Path.cwd()
    os.chdir(path)
    try:
        yield
    finally:
        os.chdir(prev_cwd)
    

if len(sys.argv) < 34:
    raise ValueError("Not enough arguments passed to DB_RECONSTRUCT")

pipeline_args = {
    'cpus': sys.argv[1], 
    'gpus': sys.argv[2],
    'out_dir': sys.argv[3], 
    'read_mdoc': sys.argv[4], 
    'remove_xrays': sys.argv[5],
    'prealign_bin': sys.argv[6],
    'track_method': sys.argv[7], 
    'size_gold': sys.argv[8],
    'final_bin': sys.argv[9], 
    'do_sirt': sys.argv[10],
    'do_trimvol': sys.argv[11], 
    'pixel_size': sys.argv[12], 
    'tiltaxis': sys.argv[13],
    'dose_sym': sys.argv[14], 
    'voltage': sys.argv[15], 
    'cs': sys.argv[16], 
    'reorient': sys.argv[17], 
    'thickness_binned': sys.argv[18], 
    'thickness_unbinned': sys.argv[19],
    'use_sobel': sys.argv[20],
    'num_beads': sys.argv[21],
    'sobel_sigma': sys.argv[22],
    'patch_size': [sys.argv[23], sys.argv[24]],
    'patch_overlap': [sys.argv[25], sys.argv[26]],
    'do_ctf': sys.argv[27],
    'defocus_range': [sys.argv[28], sys.argv[29]],
    'autofit_range': sys.argv[30],
    'autofit_step': sys.argv[31],
    'tune_fitting_sample': sys.argv[32],
    'fake_sirt_iters': sys.argv[33]
}
assert 'out_dir' in pipeline_args
pipeline_args['out_dir'] = Path.cwd() / Path(pipeline_args['out_dir'])

print(pipeline_args)

com, adoc = setup_serieswatcher(**pipeline_args)
brt_pipeline = "brt_pipeline"
with chdir(pipeline_args['out_dir']):
    # cmd = f'tmux new-session -d -s {brt_pipeline}'
    # subprocess.run(cmd, shell=True)

    # cmd = f'tmux send-keys "serieswatcher -com {com} -adoc {adoc}" C-m'
    cmd = f'serieswatcher -com {com} -adoc {adoc}'
    subprocess.run(cmd, shell=True)
    
    print(f"CHECK STATUS OF PIPELINE RECONSTRUCTION WITH COMMAND: tmux a -t {brt_pipeline}")
    print("DETACH SESSION (i.e. still running, but now longer watching it) with: control-b d")
    print(f"KILL SESSION: tmux -t {brt_pipeline} kill-session")
    print("KILL ALL SESSIONS: tmux kill-server")