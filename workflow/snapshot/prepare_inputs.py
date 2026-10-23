"""Preserve each upstream reaction block; copy provisional Henry apparatus blocks."""
from pathlib import Path
import re, json, hashlib, subprocess
BASE = Path('/path/to/EIC/SIMC_RUNS/ACPOSTUMA_SIMC_GFORTRAN')
REPO = Path('/path/to/simc')
HENRY = Path('/path/to/EIC/SIMC_RUNS/HENRY_SIMC_FILE_INPUT/phi_baseline_10p60.inp')
TEMPLATES = dict(omega='omega_245_54_-3.inp', phi='phi_600_320_17_+0590.inp',
    rho='rho_600_320_17_+0590.inp', eta='eta_600_320_17_+0590.inp',
    etap='etap_600_320_17_+0590.inp', xphsp='xphsp_245_54_00.inp')
def block(text):
    return re.search(r'(?ims)^begin parm experiment\b.*?^end parm experiment[^\n]*', text)
def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()
def main():
    text = HENRY.read_text()
    rest = text[block(text).end():]
    # Correct only the misspelled closing block label, preserving all values.
    rest = rest.replace('end parm beamandtergetinfo', 'end parm beamandtargetinfo')
    dest = BASE/'infiles'; dest.mkdir(exist_ok=True)
    audit = {'status':'PROVISIONAL technical test, not production physics',
        'commit':subprocess.check_output(['git','rev-parse','HEAD'],cwd=REPO,text=True).strip(),
        'henry_input':str(HENRY), 'henry_sha256':digest(HENRY),
        'apparatus_change':'Only corrected beamandtergetinfo closing-label typo.', 'channels':{}}
    for channel, name in TEMPLATES.items():
        src = REPO/'infiles'/name
        reaction = block(src.read_text()).group()
        reaction, n = re.subn(r'(?m)^(\s*ngen\s*=\s*)\d+',r'\g<1>1000',reaction)
        assert n == 1, (name,n)
        output = '; PROVISIONAL: apparatus copied from local Henry baseline, NOT confirmed production input.\n'+reaction+rest
        target = channel+'_Roman_007.inp'
        (dest/target).write_text(output)
        (REPO/'infiles'/target).write_text(output)
        assert output[block(output).end():] == rest
        audit['channels'][channel] = {'template':name,'template_sha256':digest(src),'input':target,
            'sha256':digest(dest/target),'ngen_requested_successes':1000}
    (BASE/'input_audit.json').write_text(json.dumps(audit,indent=2)+'\n')
    print(json.dumps(audit,indent=2))
if __name__ == '__main__': main()
