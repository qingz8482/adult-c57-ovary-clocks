#!/usr/bin/env python3
"""Run all prespecified stages, resume checkpoints, and summarize completed analyses."""
import argparse
import concurrent.futures
import json
import os
from pathlib import Path
import subprocess
import sys
import time

DEV = Path(__file__).resolve().parent.parent
parser = argparse.ArgumentParser()
parser.add_argument('--config', type=Path, default=DEV/'config/production.json')
parser.add_argument('--rscript', default=os.environ.get('RSCRIPT', 'Rscript'))
parser.add_argument('--phases', default='A;B,C;final;sensitivities;A_repeated')
parser.add_argument('--workers', type=int, help='Execution-only concurrency override; each R fit remains single threaded')
args = parser.parse_args()
args.config = args.config.resolve()
cfg = json.loads(args.config.read_text())
base = Path(cfg['output_dir'])
if not base.is_absolute():
    base = DEV/base
base.mkdir(parents=True, exist_ok=True)
logs = base/'logs'
logs.mkdir(exist_ok=True)
lock = base/'RUNNING.lock'
try:
    fd = os.open(lock, os.O_CREAT | os.O_EXCL | os.O_WRONLY)
except FileExistsError:
    raise SystemExit(f'Run lock exists: {lock}. Check the recorded PID before resuming.')
with os.fdopen(fd, 'w') as f:
    json.dump({'pid': os.getpid(), 'started': time.strftime('%Y-%m-%dT%H:%M:%S%z'),
               'config': str(args.config)}, f, indent=2)
env = dict(os.environ, OPENBLAS_NUM_THREADS='1', OMP_NUM_THREADS='1', MKL_NUM_THREADS='1', TZ='Asia/Singapore')
recipes = ['HVG-ElasticNet', 'VersionB-ElasticNet', 'VersionB-PCA']
previous_state = json.loads((base/'run_status.json').read_text()) if (base/'run_status.json').exists() else {}
actual_workers = args.workers if args.workers is not None else cfg['workers']
if actual_workers < 1:
    raise ValueError('workers must be positive')
state = {'status': 'running', 'completed_phases': previous_state.get('completed_phases', []),
         'config': str(args.config), 'pid': os.getpid(), 'actual_workers': actual_workers,
         'one_thread_per_R_fit': True}

def write_state():
    state['updated'] = time.strftime('%Y-%m-%dT%H:%M:%S%z')
    tmp = base/'run_status.json.tmp'
    tmp.write_text(json.dumps(state, indent=2)+'\n')
    tmp.replace(base/'run_status.json')

def run_recipe(recipe, phase, fold=None, repeat=None):
    suffix = f"_r{repeat:02d}_f{fold}" if fold is not None else ""
    log = logs/f'{recipe}_{phase.replace(",", "_")}{suffix}.log'
    command = [args.rscript, str(DEV/'scripts/run_models.R'), f'--recipe={recipe}',
               f'--config={args.config}', f'--only={phase}']
    if fold is not None:
        command += [f'--fold={fold}', f'--repeat={repeat}']
    with log.open('a', buffering=1) as f:
        f.write('\nSTART '+time.strftime('%Y-%m-%dT%H:%M:%S%z')+'\n')
        subprocess.run(command, env=env, cwd=DEV, stdout=f, stderr=subprocess.STDOUT, check=True)
    return recipe

def summarize():
    subprocess.run([args.rscript, str(DEV/'scripts/summarize_results.R'), f'--config={args.config}'],
                   env=env, cwd=DEV, check=True)

try:
    write_state()
    for phase in args.phases.split(';'):
        if not any(a in cfg['analyses'] for a in phase.split(',')):
            continue
        state['current_phase'] = phase
        write_state()
        print('START PHASE', phase, flush=True)
        # Prime shared run headers once per recipe before concurrent task writers.
        if not state['completed_phases']:
            for r in recipes:
                subprocess.run([args.rscript, str(DEV/'scripts/run_models.R'), f'--recipe={r}',
                                f'--config={args.config}', '--max-tasks=0'], env=env, cwd=DEV, check=True)
        jobs = []
        if phase == 'A':
            jobs = [(r, phase, k, 0) for k in range(1,5) for r in recipes]
        elif phase == 'A_repeated':
            jobs = [(r, phase, k, rep) for rep in range(1,cfg['repeats']+1) for k in range(1,5) for r in recipes]
        elif phase == 'sensitivities':
            jobs = [(r, 'A_'+v, k, 0) for v in cfg['sensitivities'] for k in range(1,5) for r in recipes]
        else:
            jobs = [(r, phase, None, None) for r in recipes]
        with concurrent.futures.ThreadPoolExecutor(max_workers=actual_workers) as pool:
            futures = {pool.submit(run_recipe, *job): job for job in jobs}
            for future in concurrent.futures.as_completed(futures):
                print('COMPLETED', future.result(), phase, flush=True)
        summarize()
        if phase not in state['completed_phases']:
            state['completed_phases'].append(phase)
        write_state()
    state['status'] = 'complete'
    state.pop('current_phase', None)
    print('ALL PRESPECIFIED PHASES COMPLETE', flush=True)
except BaseException as exc:
    state['status'] = 'failed'
    state['error'] = str(exc)
    raise
finally:
    write_state()
    lock.unlink(missing_ok=True)
