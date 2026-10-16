"""Relocatable wrapper for the frozen Roman 007 inputs and plotting functions."""
from pathlib import Path
import argparse, datetime, importlib.util, json, re, shlex, subprocess

REPO = Path(__file__).resolve().parents[1]
WORK = REPO / 'workflow'
SNAPSHOT = WORK / 'snapshot'
CHANNELS = ['omega', 'phi', 'rho', 'eta', 'etap', 'xphsp']
FLAGS = '-O -w -ffixed-line-length-132 -ff2c -fno-automatic -fdefault-real-8 -fallow-argument-mismatch'

def checked(args, cwd=REPO, **kwargs):
    return subprocess.run(args, cwd=cwd, check=True, **kwargs)

def root_flags(option):
    return shlex.split(subprocess.check_output(['root-config', option], text=True))

def compile_converter():
    cwd = REPO / 'util/root_tree'
    checked(['gfortran', '-O3', '-fno-automatic', '-c', 'make_root_tree.f'], cwd)
    checked(['g++', '-O3', *root_flags('--cflags'), '-DMYROOT', '-c', 'froot.c', '-o', 'froot.co'], cwd)
    checked(['gfortran', '-o', 'make_root_tree', 'make_root_tree.o', 'froot.co', *root_flags('--libs'), '-lstdc++'], cwd)

def analyze(root_paths, output):
    spec = importlib.util.spec_from_file_location('frozen_analysis', SNAPSHOT / 'analyze_root.py')
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    report = {}
    for channel, path in root_paths.items():
        report[channel] = module.analyze(path, 'h10', output/(channel+'_Roman_007')/'plots', channel+' SIMC', True)
    (output/'analysis_summary.json').write_text(json.dumps(report, indent=2)+'\n')

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('action', choices=['build', 'run', 'analyze'])
    parser.add_argument('--output', type=Path, help='New, non-existing output directory')
    args = parser.parse_args()
    if args.action == 'build':
        if not (REPO/'CTP/O.Linux/Linux/lib/libctp.a').is_file():
            parser.error('First run make -C CTP on the WSL host with RPC development headers.')
        checked(['make', 'FFLAGSA='+FLAGS])
        compile_converter()
        return
    stamp = datetime.datetime.now().strftime('%Y%m%d_%H%M%S_%f')
    output = (args.output or WORK/'output'/(args.action+'_'+stamp)).resolve()
    if output.exists(): parser.error('Refusing to overwrite an existing output directory.')
    output.mkdir(parents=True)
    paths = {}
    if args.action == 'run':
        import ROOT
        import shutil
        ROOT.gROOT.SetBatch(True)
        for directory in ['worksim','outfiles','runout']: (REPO/directory).mkdir(exist_ok=True)
        for channel in CHANNELS:
            name = channel+'_Roman_007'
            destination = output/name
            destination.mkdir()
            expected = (SNAPSHOT/'infiles'/(name+'.inp')).read_bytes()
            actual = (REPO/'infiles'/(name+'.inp')).read_bytes()
            if actual != expected: raise RuntimeError('Input differs from frozen version: '+name)
            shutil.copy2(REPO/'infiles'/(name+'.inp'), destination)
            with (destination/'simc.log').open('w') as log:
                checked([str(REPO/'simc')], input=name+'\n', text=True, stdout=log, stderr=subprocess.STDOUT, timeout=180)
            if not re.search(r'Last Event.*1000\s+successes', (destination/'simc.log').read_text()):
                raise RuntimeError('SIMC did not report 1000 successes: '+name)
            with (destination/'conversion.log').open('w') as log:
                checked(['./make_root_tree'], REPO/'util/root_tree', input=name+'\n', text=True, stdout=log, stderr=subprocess.STDOUT, timeout=60)
            root_path = REPO/'worksim'/(name+'.root')
            f = ROOT.TFile.Open(str(root_path))
            if not f or f.IsZombie() or not f.Get('h10') or f.Get('h10').GetEntries()!=1000:
                raise RuntimeError('Invalid ROOT output: '+name)
            f.Close()
            shutil.copy2(root_path, destination)
            for path in (REPO/'outfiles').glob(name+'*'):
                if path.is_file(): shutil.copy2(path,destination)
            paths[channel] = destination/(name+'.root')
    else:
        paths = {c:SNAPSHOT/'runs'/(c+'_Roman_007')/(c+'_Roman_007.root') for c in CHANNELS}
    analyze(paths, output)
    print('Results:', output)

if __name__ == '__main__': main()
