#!/usr/bin/env python3
# A Python wrapper to run DeepDeWedge
# Integrated into the cryo-ET automated data processing pipeline at Sanofi

# Matthew Martinez
# Sanodi US - Dept. of Large Molecules Research, Protein Engineering Group, Structural Biology

# ARGUMENTS TO BE PASSED:
#   If running as pipeline (0, 1)
#   Basename of files for training (otherwise selects randomly)
#   Number of files for training
#   [Optional] Path to parent folder of files

import contextlib
import logging
import os
from pathlib import Path
import random
import subprocess
import sys
import time
from typing import Generator

import yaml


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

logger = logging.getLogger(__name__)
filename = f'{time.strftime("%Y%m%d_%H%M", time.localtime())}_DDW.log'
handler = logging.FileHandler(filename)
handler.setFormatter(formatter)
logger.setLevel(logging.INFO)
logger.addHandler(handler)


@contextlib.contextmanager
def chdir(path: str | Path) -> Generator[None, None, None]:
    """ Changed working directory and returns to the previous on exit """
    prev_cwd = Path.cwd()
    os.chdir(path)
    try:
        yield
    finally:
        os.chdir(prev_cwd)


def locate_halfsets(
        p: Path, 
        evens_ext: str='evens', 
        odds_ext: str='odds', 
        ext: str='mrc'
    ) -> tuple[tuple[Path, Path]]:
    """
    Find all the halfsets

    :param Path p: Path within which all the halfsets will be located
    :param str evens_ext: Added extension to the evens halfset
    :param str odds_ext: Added extension to the odds halfset
    :param str ext: File extension
    
    :return tuple(tuple(Path, Path)): Pairs of halfsets
    :raises FileNotFoundError: If no files found
    """
    evens = [f for f in p.rglob(f'*{evens_ext}.{ext}') if 'full' not in f.name]
    odds = [f for f in p.rglob(f'*{odds_ext}.{ext}') if 'full' not in f.name]

    logger.info('Found %d evens and %d odds' % (len(evens), len(odds)))

    if len(evens) == 0 and len(odds) == 0:
        logger.error('HALFSETS NOT FOUND WITHIN PATH %s' % p.name)
        raise FileNotFoundError

    halfsets = tuple(zip(evens, odds, strict=True))
    return halfsets


def get_random_halfsets(halfsets: tuple[tuple[Path, Path]], n: int=5) -> tuple[tuple[Path, Path]]:
    """
    Get random sample of halfsets for training

    :param tuple(tuple(Path, Path)) halfsets: Halfset pairs
    :param int n: Number of halfsets to randomly select
    
    return tuple(tuple(Path, Path)): pairs of randomly selected halfsets
    """
    samples = random.sample(range(0, len(halfsets)), n)
    evens = [halfsets[i][0] for i in samples]
    odds = [halfsets[i][1] for i in samples]

    logger.info('Added to the training set: {}'.format(', '.join([f.name for f in evens])))
    logger.info('Added to the training set: {}'.format(', '.join([f.name for f in odds])))
    return tuple(zip(evens, odds))


def construct_config(p: Path, proc: str, halfsets: tuple[tuple[Path, Path]], **kwargs) -> Path:
    """
    Construct the config files for fitting and refining

    :param Path p: Path to save the files
    :param str proc: 'fit' or 'refine', whether the config file is for model fitting or tomogram refinement
    :param tuple(tuple(Path, Path)): Halfset pairs
    :param **kwargs: Additional args for the config.yaml files

    :return Path: config file
    """
    if proc not in ['fit', 'refine']:
        logger.error('VALUE ERROR: When constructing config file, must be either fit or refine')
        raise ValueError
    
    if proc == 'fit':
        filename = p / Path('fit_config.yaml')
    elif proc == 'refine':
        assert 'model_checkpoint_file' in kwargs, "Must provide model checkpoint file"
        assert Path(kwargs['model_checkpoint_file']).exists()
        filename = p / Path('refine_config.yaml')
    
    # Shared
    subtomo_size = kwargs['subtomo_size'] if 'subtomo_size' in kwargs else 96
    mw_angle = kwargs['mw_angle'] if 'mw_angle' in kwargs else 60
    num_workers = kwargs['num_workers'] if 'num_workers' in kwargs else 4
    gpu = kwargs['gpu'] if 'gpu' in kwargs else 0
    seed = kwargs['seed'] if 'seed' in kwargs else 42
    overwrite = kwargs['overwrite'] if 'overwrite' in kwargs else True

    # Prepare data
    mask_files = kwargs['mask_files'] if 'mask_files' in kwargs else None
    min_nonzero_mask_fraction_in_subtomo = (
        kwargs['min_nonzero_mask_fraction_in_subtomo'] 
        if 'min_nonzero_mask_fraction_in_subtomo' in kwargs else 0.3
    )
    subtomo_extraction_strides = (
        kwargs['subtomo_extraction_strides'] 
        if 'subtomo_extraction_strides' in kwargs else [64, 80, 80]
    )
    val_fraction = kwargs['val_fraction'] if 'val_fraction' in kwargs else 0.2

    # Fit model
    unet_params_dict = {
        'chans': kwargs['chans'] if 'chans' in kwargs else 64,
        'num_downsample_layers': kwargs['num_downsample_layers'] if 'num_downsample_layers' in kwargs else 3,
        'drop_prob': kwargs['drop_prob'] if 'drop_prob' in kwargs else 0.0
    }
    adam_params_dict = {
        'lr': kwargs['lr'] if 'lr' in kwargs else 0.0004
    }
    num_epochs = kwargs['num_epochs'] if 'num_epochs' in kwargs else 1000
    batch_size_fit = kwargs['batch_size_fit'] if 'batch_size_fit' in kwargs else 5
    update_subtomo_missing_wedges_every_n_epochs = (
        kwargs['update_subtomo_missing_wedges_every_n_epochs']
        if 'update_subtomo_missing_wedges_every_n_epochs' in kwargs else 10
    )
    check_val_every_n_epochs = (
        kwargs['check_val_every_n_epochs'] if 'check_val_every_n_epochs' in kwargs
        else 10
    )
    save_n_models_with_lowest_val_loss = (
        kwargs['save_n_models_with_lowest_val_loss']
        if 'save_n_models_with_lowest_val_loss' in kwargs else 5
    )
    save_n_models_with_lowest_fitting_loss = (
        kwargs['save_n_models_with_lowest_fitting_loss'] 
        if 'save_n_models_with_lowest_fitting_loss' in kwargs else 5
    )
    save_model_every_n_epochs = (
        kwargs['save_model_every_n_epochs'] if 'save_model_every_n_epochs' in kwargs
        else 50
    )
    model_logger = kwargs['logger'] if 'logger' in kwargs else "csv"

    # Refine
    model_checkpoint_file = kwargs['model_checkpoint_file'] if 'model_checkpoint_file' in kwargs else None
    subtomo_overlap = kwargs['subtomo_overlap'] if 'subtomo_overlap' in kwargs else 32
    batch_size_refine = kwargs['batch_size_refine'] if 'batch_size_refine' in kwargs else 10

    # Split halfsets into evens and odds
    tomo_0_files, tomo_1_files = zip(*halfsets)

    # Construct the file data
    conf = {
        'shared': {
            'project_dir': (p / 'DDW').as_posix(),
            'tomo0_files': [f.as_posix() for f in tomo_0_files],
            'tomo1_files': [f.as_posix() for f in tomo_1_files],
            'subtomo_size': subtomo_size,
            'mw_angle': mw_angle,
            'num_workers': num_workers,
            'gpu': gpu,
            'seed': seed,
            'overwrite': overwrite
        },
        'prepare_data': {
            'mask_files': mask_files,
            'min_nonzero_mask_fraction_in_subtomo': min_nonzero_mask_fraction_in_subtomo,
            'subtomo_extraction_strides': subtomo_extraction_strides,
            'val_fraction': val_fraction
        },
        'fit_model': {
            'unet_params_dict': unet_params_dict,
            'adam_params_dict': adam_params_dict,
            'num_epochs': num_epochs,
            'batch_size': batch_size_fit,
            'update_subtomo_missing_wedges_every_n_epochs': update_subtomo_missing_wedges_every_n_epochs,
            'check_val_every_n_epochs': check_val_every_n_epochs,
            'save_n_models_with_lowest_val_loss': save_n_models_with_lowest_val_loss,
            'save_n_models_with_lowest_fitting_loss': save_n_models_with_lowest_fitting_loss,
            'save_model_every_n_epochs': save_model_every_n_epochs,
            'logger': model_logger
        },
        'refine_tomogram': {
            'model_checkpoint_file': model_checkpoint_file,
            'subtomo_overlap': subtomo_overlap,
            'batch_size': batch_size_refine
        }
    }

    # Save the YAML file
    with open(filename, 'w') as f:
        yaml.safe_dump(conf, f, sort_keys=False)

    logger.info('Saved config file %s' % filename)
    return Path(filename)


def get_best_model(p: Path, mode: str='val') -> Path:
    """
    Get the best model, according the the metric 'mode'
    
    :param Path p: project directory
    :param str mode: Get best val, fit model or latest model
    
    :return Path: Path to the model
    """
    assert mode in ['val', 'fit', 'latest']

    if mode == 'val':
        d = list(p.rglob('val_loss'))
        d = d[-1] if len(d) > 1 else d[0]
        ckpts = [x for x in d.iterdir() if x.suffix == '.ckpt']
        loss = [float(x.name.split('=')[-1].split('.ckpt')[0]) for x in ckpts]
        best = sorted(tuple(zip(loss, ckpts)))[0][1]

        logger.info('Identified best model by validation loss -- %s' % best.name)

        return best
    elif mode == 'fit':
        d = list(p.rglob('fitting_loss'))
        d = d[-1] if len(d) > 1 else d[0]
        ckpts = [x for x in d.iterdir() if x.suffix == '.ckpt']
        loss = [float(x.name.split('=')[-1].split('.ckpt')[0]) for x in ckpts]
        best = sorted(tuple(zip(loss, ckpts)))[0][1]

        logger.info('Identified best model by fitting loss -- %s' % best.name)

        return best
    else:
        # By lates epoch
        d = list(p.rglob('epoch'))
        d = d[-1] if len(d) > 1 else d[0]
        ckpts = [x for x in d.iterdir() if x.suffix == '.ckpt']
        loss = [float(x.name.split('=')[-1].split('.ckpt')[0]) for x in ckpts]
        best = sorted(tuple(zip(loss, ckpts)), reverse=True)[0][1]

        logger.info('Identified model by latest epoch -- %s' % best.name)

        return best



def prepare(config: Path) -> None:
    """
    run 'ddw prepare-data --config ./config.yaml
    
    :param Path config: Path to config.yaml file
    """
    assert config.exists()
    logger.info('\nBEGINNING TO PREPARE THE DATA')

    cmd = f'source activate DDW && ddw prepare-data --config {config}'
    if (exit_code := subprocess.run(cmd, shell=True).returncode) != 0:
        logger.error('Error in data preparation. Return code %d' % exit_code)
        sys.exit(exit_code)

    logger.info('COMPLETED PREPARING DATA')


def fit(config: Path) -> None:
    """
    Run 'ddw fit-model --config ./config.yaml'

    :param Path config: Path to config.yaml file
    """
    assert config.exists()
    logger.info('\nBEGINNING FIT-MODEL')

    cmd = f'source activate DDW && ddw fit-model --config {config}'
    if (exit_code := subprocess.run(cmd, shell=True).returncode) != 0:
        logger.error('Error in model fitting. Return code %d' % exit_code)
        sys.exit(exit_code)
    
    logger.info('TRAINING COMPLETE!')



def refine(config: Path) -> None:
    """
    Run 'ddw refine-tomogram --config ./config.yaml
    
    :param Path config: Path to config.yaml file
    """
    assert config.exists()

    logger.info('\nBEGINNING TO REGINE ALL TOMOGRAMS')

    cmd = f'source activate DDW && ddw refine-tomogram --config {config}'
    if (exit_code := subprocess.run(cmd, shell=True).returncode) != 0:
        logger.error('Error in refining the tomograms. Return code %d' % exit_code)
        sys.exit(exit_code)

    logger.info('REFINEMENT COMPLETE')
    

# Determine if this is running during the pipeline or standalone
if len(sys.argv) < 2:
    raise ValueError("Not enough arguments. Add '1' if running during the pipeline, '0' if running standalone.")
RUNNING_PIPELINE = sys.argv[1]

if len(sys.argv) == 2:
    p = Path.cwd()
elif len(sys.argv) == 3:
    p = Path(sys.argv[2])
    assert p.exists(), f"Given path {p} does not exist"
    assert p.is_dir(), f"Given path {p} is not a valid directory"

# For now, just choose randomly for training among the datasets for prototyping
halfsets = locate_halfsets(p)
training_sets = get_random_halfsets(halfsets)

# Prepare data and fit model
config_file = construct_config(p, 'fit', training_sets)
# prepare(config=config_file)
fit(config=config_file)

# Refine the tomograms
with open(config_file, 'r') as f:
    config_data = yaml.safe_load(f)
proj_dir = Path(config_data['shared']['project_dir'])
proj_dir = Path('/common/workdir/DDW_test')
model = get_best_model(proj_dir, mode='val').as_posix()
config_file = construct_config(p, 'refine', halfsets, model_checkpoint_file=model)
refine(config=config_file)

# Sync to proper workdir location and to database S3 location

