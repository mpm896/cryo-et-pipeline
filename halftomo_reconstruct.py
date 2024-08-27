# Matthew Martinez
# Sanofi US - Dept. of Large Molecules Research, Protein Engineering Group, Structural Biology

# Create tomograms composed of even and odd tilts for noise2noise denoising models
# such as Cryo-CARE and DeepDeWedge

# If running in standalone mode, must copy over appropriate files from database:
#   Original tilt series
#   Transform file .xf
#   Tilt file .tlt
#   Xtiltfile .xtilt

# Constructing tilt file may really need to first be copied from the OG tilt file - 
# there are some specifics that should be kept the same as the OG tilt file

import contextlib
import logging
import os
from pathlib import Path
import sys
import subprocess
import time
from typing import Generator

import mrcfile


DB_DIR = '/root/cloud-data/its-cmo-darwin-magellan-workspaces-folders/WS_Cryoem/CX_LMR/Project_directories/cryo-et-data'
EXT = '*_rec.mrc'  # Suffix of completed tomogram
BIN = 6  # Bin factor
GPU = 0  # Number GPU device - 0 for best

# Setup the logger
logging.basicConfig(
    level=logging.INFO,
    format="{asctime} - {levelname} - {message}",
    style="{",
    datefmt="%Y-%m-%d %H:%M",
    encoding="utf-8"
)
formatter = logging.Formatter(
    fmt="{asctime} - {levelname} - {message}",
    style="{",
    datefmt="%Y-%m-%d %H:%M"
)

@contextlib.contextmanager
def chdir(path: str | Path) -> Generator[None, None, None]:
    """ Changed working directory and returns to the previous on exit """
    prev_cwd = Path.cwd()
    os.chdir(path)
    try:
        yield
    finally:
        os.chdir(prev_cwd)


def get_datasets(p: str | Path) -> list[Path]:
    """ Get all the datasets in the given directory for half set generation """
    return [x for x in p.iterdir() if x.is_dir() and list(x.glob(EXT))]


def check_metadata(dataset: Path) -> int:
    """ 
    Check the metadata required for generation of half sets 
    
    :returns: int
        0 - no metadata present
        1 - all metadata present, can just run tilt
        2 - partial metadata, need to run both newstack and tilt
    """
    db_set = dataset.name.split('/')[0]
    name = [x for x in dataset.glob('*_rec.mrc') if 'full' not in x.name][0].name.split('_rec')[0]
    
    xf = list(dataset.glob(f'{name}.xf'))
    ali = list(dataset.glob(f'{name}_ali.mrc'))
    tlt = list(dataset.glob(f'{name}.tlt'))
    xtilt = list(dataset.glob(f'{name}.xtilt'))

    if (not xf and not ali) or not tlt or not xtilt:
        logger.info('Not enough metadata for dataset %s - %s' % (db_set, name))
        return 0
    elif ali and tlt and xtilt:
        logger.info(
            'Dataset %s - %s -- Procede with halfset generation via tilt'
            % (db_set, name)
        )
        return 1
    elif xf and tlt and xtilt and not ali:
        logger.info(
            '\nDataset %s - %s -- Procede with halfset generation via newstack'
            % (db_set, name)
        )
        return 2
    return 0


def construct_coms(dataset: Path, status: int) -> None:
    """ Construct the newstack.com and tilt.com files """
    name = [x for x in dataset.glob('*_rec.mrc') if 'full' not in x.name][0].name.split('_rec')[0]
    mdoc = list(dataset.glob(f'{name}*.mdoc'))[0]

    logger.info('Identified MDOC file %s' % mdoc.name)

    if status == 2:
        logger.info('Constructing com file for newstack -> newst.com')
        
        # Construct newst.com
        newst_com = dataset / Path('newst.com')
        newst_com.write_text(
            '$setenv IMOD_OUTPUT_FORMAT MRC\n'
            '$newstack -StandardInput\n'
            'AntialiasFilter	4\n'
            f'InputFile	{name}.mrc\n'
            f'OutputFile	{name}_ali.mrc\n'
            f'TransformFile	{name}.xf\n'
            'LinearInterpolation	0\n'
            f'BinByFactor	{BIN}\n'
            'TaperAtFill	1,1\n'
            'AdjustOrigin\n'
            'OffsetsInXandY	0,0\n'
            '#DistortionField	.idf\n'
            'ImagesAreBinned	1\n'
            f'#GradientFile	{name}.maggrad\n'
            '$if (-e ./savework) ./savework'
        )

    logger.info('Constructing com files for tilt -> tilt_evens.com and tilt_odds.com')
    
    # Construct tilt_evens.com and tilt_odds.com
    tilt_com = dataset / Path('tilt.com')
    tilt_evens = dataset / Path('tilt_evens.com')
    tilt_odds = dataset / Path('tilt_odds.com')
    
    if tilt_evens.exists():
       tilt_evens.unlink()
    if tilt_odds.exists():
        tilt_odds.unlink() 


    mrc = dataset / Path(f'{name}.mrc')
    rec = dataset / Path(f'{name}_rec.mrc')

    with mrcfile.open(mrc) as f:
        full_mage_size = [f.header.nx, f.header.ny, f.header.nz]
        pix = f.voxel_size
        if round(pix.x.item(), 2) == round(pix.y.item(), 2) == round(pix.z.item(), 2):
            mrc_pix = round(pix.x.item(), 2)
    with mrcfile.open(rec) as f:
        rec_image_size = [f.header.nx, f.header.ny, f.header.nz]
        pix = f.voxel_size
        if round(pix.x.item(), 2) == round(pix.y.item(), 2) == round(pix.z.item(), 2):
            rec_pix = round(pix.x.item(), 2)

    evens = f'INCLUDE {",".join(str(i+1) for i in range(full_mage_size[2]) if i % 2 == 0)}'
    odds = f'INCLUDE {",".join(str(i+1) for i in range(full_mage_size[2]) if i % 2 == 1)}'

    # Write the files based on the preexisting tilt.com
    if tilt_com.exists():
        logger.info('Creating tilt coms based off of existing com files.')
        with tilt_com.open('r') as f:
            for line in f:
                if 'InputProjections' in line:
                    line = f'InputProjections {name}_ali.mrc\n'
                if 'OutputFile' in line:
                    with tilt_evens.open('a') as t:
                        line = f'OutputFile	{name}_full_rec_evens.mrc\n'
                        t.write(line)
                    with tilt_odds.open('a') as t:
                        line = f'OutputFile	{name}_full_rec_odds.mrc\n'
                        t.write(line)
                    continue
                if 'IMAGEBINNED' in line:
                    line = f'IMAGEBINNED	{BIN}\n'
                if 'TILTFILE' in line:
                    line = f'TILTFILE {name}.tlt\n'
                if 'XTILTFILE' in line:
                    line = f'XTILTFILE {name}.xtilt\n'
                if 'useGPU' in line:
                    line = f'UseGPU	    {GPU}\n'

                with tilt_evens.open('a') as t:
                    t.write(line)
                with tilt_odds.open('a') as t:
                    t.write(line)
            
        with tilt_evens.open('a') as f:
            f.write(f'{evens}')
        with tilt_odds.open('a') as f:
            f.write(f'{odds}')
    else:
        logger.info('Creating brand new tilt com files.')
        rec_bin = full_mage_size[0] // rec_image_size[0]
        rec_thickness = int(rec_bin * rec_image_size[2])

        tilt_evens.write_text(
            '$setenv IMOD_OUTPUT_FORMAT MRC\n'
            '$tilt -StandardInput\n'
            'FakeSIRTiterations	10\n'
            f'InputProjections {name}_ali.mrc\n'
            f'OutputFile	{name}_full_rec_evens.mrc\n'
            f'IMAGEBINNED	{BIN}\n'
            f'TILTFILE {name}.tlt\n'
            f'XTILTFILE {name}.xtilt\n'
            f'UseGPU	{GPU}\n'
            f'THICKNESS	{rec_thickness}\n'
            'RADIAL .35 .035\n'
            'FalloffIsTrueSigma 1\n'
            'SCALE 0 0.00144\n'
            'PERPENDICULAR\n'
            'MODE 1\n'
            f'FULLIMAGE	 {full_mage_size[0]} {full_mage_size[1]}\n'
            'SUBSETSTART	0 0\n'
            'AdjustOrigin 1\n'
            f'{evens}\n'
            '$if (-e ./savework) ./savework'
        )

        tilt_odds.write_text(
            '$setenv IMOD_OUTPUT_FORMAT MRC\n'
            '$tilt -StandardInput\n'
            'FakeSIRTiterations	10\n'
            f'InputProjections {name}_ali.mrc\n'
            f'OutputFile	{name}_full_rec_odds.mrc\n'
            f'IMAGEBINNED	{BIN}\n'
            f'TILTFILE {name}.tlt\n'
            f'XTILTFILE {name}.xtilt\n'
            f'UseGPU	{GPU}\n'
            f'THICKNESS	{rec_thickness}\n'
            'RADIAL .35 .035\n'
            'FalloffIsTrueSigma 1\n'
            'SCALE 0 0.00144\n'
            'PERPENDICULAR\n'
            'MODE 1\n'
            f'FULLIMAGE	 {full_mage_size[0]} {full_mage_size[1]}\n'
            'SUBSETSTART	0 0\n'
            'AdjustOrigin 1\n'
            f'{odds}\n'
            '$if (-e ./savework) ./savework'
        )


def newstack(dataset: Path) -> bool:
    """ 
    Run newstack on a dataset with its newst.com file
     
    :returns: bool
        True if completed successfully
        False if error  
     """
    db_set = dataset.name.split('/')[0]
    name = [x for x in dataset.glob('*_rec.mrc') if 'full' not in x.name][0].name.split('_rec')[0]

    cmd = 'subm newst.com'
    with chdir(dataset):
        logger.info(
            'Beginning newstack to generate aligned stack for dataset %s -- %s'
            % (db_set, name)
        )

        out = subprocess.run(cmd, shell=True, capture_output=True, encoding='utf-8')
        if 'finished successfully' in out.stdout:
            logger.info(
                '%s -- Successfully generated aligned stack %s_ali.mrc'
                % (db_set, name)
            )
            return True
        elif 'ERROR' in out.stdout:
            logger.error(
                '%s -- Error in generating aligned stack'
                % db_set
            )
            return False
        return False


def tilt(dataset: Path) -> bool:
    """ 
    Run tilt on a dataset to generate half tomos with its tilt.com files

    :returns: bool
        True if completed successfully
        False is error 
     """
    db_set = dataset.name.split('/')[0]
    name = [x for x in dataset.glob('*_rec.mrc') if 'full' not in x.name][0].name.split('_rec')[0]
    
    with chdir(dataset):

        cmd = 'subm tilt_evens.com'
        logger.info(
            'tilt - Begining generation of tomogram from even tilts for dataset %s -- %s'
            % (db_set, name)
        )
        out = subprocess.run(cmd, shell=True, capture_output=True, encoding='utf-8')

        success = 0
        if 'finished successfully' in out.stdout:
            logger.info(
                '%s -- Successfully generated even tilt tomogram'
                % db_set
            )
            success += 1
        elif 'ERROR' in out.stdout:
            logger.error(
                '%s -- Error in generating even tilt tomogram'
                % db_set
            )

        cmd = 'subm tilt_odds.com'
        logger.info(
            'tilt - Begining generation of tomogram from odd tilts for dataset %s -- %s'
            % (db_set, name)
        )
        out = subprocess.run(cmd, shell=True, capture_output=True, encoding='utf-8')

        if 'finished successfully' in out.stdout:
            logger.info(
                '%s -- Successfully generated odd tilt tomogram'
                % db_set
            )
            success += 1
        elif 'ERROR' in out.stdout:
            logger.error(
                '%s -- Error in generating odd tilt tomogram'
                % db_set
            )

        return True if success == 2 else False


def trimvol(dataset: Path) -> bool:
    """ 
    Run trimvol on a dataset to rotate the final tomogram
    
    :returns: bool
        True if completed successfully
        False if error
    """
    db_set = dataset.name.split('/')[0]
    name = [x for x in dataset.glob('*_rec.mrc')if 'full' not in x.name][0].name.split('_rec')[0]
    with chdir(dataset):
        full_rec_evens = f'{name}_full_rec_evens.mrc'
        rec_evens = f'{name}_rec_evens.mrc'
        full_rec_odds = f'{name}_full_rec_odds.mrc'
        rec_odds = f'{name}_rec_odds.mrc'

        cmd = f'trimvol -rx {full_rec_evens} {rec_evens}'
        logger.info(
            'trimvol - rotating evens tilt dataset %s -- %s around the X axis'
            % (db_set, name)
        )
        out = subprocess.run(cmd, shell=True, capture_output=True, encoding='utf-8')
        
        success = 0
        if 'finished successfully' in out.stdout:
            logger.info(
                '%s -- Successfully rotated even tilt tomogram'
                % db_set
            )
            success += 1
        elif 'ERROR' in out.stdout:
            logger.error(
                '%s -- Error in rotating even tilt tomogram'
                % db_set
            )

        cmd = f'trimvol -rx {full_rec_odds} {rec_odds}'
        logger.info(
            'trimvol - rotating odds tilt dataset %s -- %s around the X axis'
            % (db_set, name)
        )
        out = subprocess.run(cmd, shell=True, capture_output=True, encoding='utf-8')

        if 'finished successfully' in out.stdout:
            logger.info(
                '%s -- Successfully rotated odd tilt tomogram'
                % db_set
            )
            success += 1
        elif 'ERROR' in out.stdout:
            logger.error(
                '%s -- Error in rotating odd tilt tomogram'
                % db_set
            )

        return True if success == 2 else False


def generate_halfsets(p: Path) -> None:
    """ Generate the half tomograms from even and odd tilts """
    # Begin the processing
    datasets: list[Path] = get_datasets(p)
    for d in datasets:
        # Check the metadata
        prog: int = check_metadata(d)

        if prog == 0:
            raise Exception('Must first transfer data from database and try again')
        else:
            # Construct the appropriate coms files
            construct_coms(d, prog)

            # Run newstack and tilt to generate the aligned tilt series and tomogram, respectively
            if prog == 2:
                if not newstack(d):
                    logging.error('%s -- Error in newstack. Terminating this dataset...' % d.name)
                    continue
            if not tilt(d):
                logging.error('%s -- Error in tilt. Terminating this dataset...' % d.name)
                continue
            if not trimvol(d):
                logging.error('%s -- error in trimvol. Terminating this dataset...' % d.name)
                continue


def subtransfer(p: Path) -> None:
    """ Move all halfset tomograms into a subdirectory """
    files_evens = [f for f in p.rglob('*_rec_evens.mrc') if 'full' not in f.name]
    files_odds = [f for f in p.rglob('*_rec_odds.mrc') if 'full' not in f.name]
    
    logger.info('\nTransferring halfset reconstructions into subdirectory "halfsets"')

    for i, f in enumerate(files_evens):
        subdir = Path(f'{f.parent}/halfsets')
        Path.mkdir(subdir, exist_ok=True)

        # Move files with Path.rename()
        f.rename(subdir / f.name)
        files_odds[i].rename(subdir / files_odds[i].name)

    logger.info('Completed transfer of halfset reconstructions into "halfsets" subdirectory')

def sync_db(p: Path) -> None:
    """ Sync halfset subdirs with the database S3 location """
    subdirs = p.rglob('halfsets')
    logger.info('\n%s -- Synching the halfset reconstructions to the database' % p.name)
    for d in subdirs:
        src = d.as_posix()
        subdir = src.split('/')[-2]
        dest = f'{DB_DIR}/{subdir}'

        if Path(dest).exists():
            cmd = f'rsync --progress -avhr {src} {dest}'
            subprocess.run(cmd, shell=True)
            logger.info('%s -- Finished synching halfset reconstructions to the database' % p.name)
        else:
            logger.warning('%s -- Does not yet exist in database. Skipping sync to database' % p.name)


# Determine if this is running during the pipeline or standalone
if len(sys.argv) < 2:
    raise ValueError("Not enough arguments. Add '1' if running during the pipeline, '0' if running standalone.")
RUNNING_PIPELINE = sys.argv[1]

if len(sys.argv) == 2:
    p = Path.cwd()
elif len(sys.argv) == 3:
    p = Path(sys.argv[2])
    print(p)
    assert p.exists(), f"Given path {p} does not exist"
    assert p.is_dir(), f"Given path {p} is not a valid directory"

# Setup logger
logger = logging.getLogger(__name__)
filename = f'{p.as_posix()}/{time.strftime("%Y%m%d_%H%M", time.localtime())}_HALF-TOMO.log'
handler = logging.FileHandler(filename)
handler.setFormatter(formatter)
logger.setLevel(logging.INFO)
logger.addHandler(handler)

# Generate halfsets and transfer to a subdirectory
generate_halfsets(p)
subtransfer(p)

# At this point, halfsets are generated and moved into the appropriate locations.
# Setup denoising in a new tmux session
pyDDW = "/cloud-data/its-cmo-darwin-magellan-workspaces-folders/WS_Cryoem/CX_LMR/Project_directories/cryo-et-pipeline/pyDDW.py"
cmd = f'tmux send-keys "python {pyDDW} 1 {p.as_posix()}" C-m'
subprocess.run(cmd, shell=True)

# Sync to the database S3 location
sync_db(p)

