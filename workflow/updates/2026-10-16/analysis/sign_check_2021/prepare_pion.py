"""Modern input field spelling; retain upstream physical settings."""
from pathlib import Path
import re, difflib
src=Path('/opt/simc/bill-master-8b43f9d/infiles')
dest=Path('/path/to/simc-workspace/work/sign_check_2021')
for old,new in [('test_pion_h','test_pion_h_sign'),('extra_nucpi','extra_nucpi_sign')]:
    s=(src/f'{old}.inp').read_text()
    lines=[]
    for line in s.splitlines(keepends=True):
        if '=' in line and not line.lstrip().startswith(';'):
            left,right=line.split('=',1)
            line=left.replace('.','%')+'='+right
        lines.append(line)
    result=''.join(lines).replace('end parm beamandtergetinfo','end parm beamandtargetinfo')
    if old=='test_pion_h':
        result=result.replace("extra_dbase_file='extra_nucpi'", "extra_dbase_file='extra_nucpi_sign'")
        result=re.sub(r'(?m)^(\s*ngen\s*=\s*)10000\b',r'\g<1>1000',result)
    (dest/f'{new}.inp').write_text(result)
    (src/f'{new}.inp').write_text(result)
    (dest/f'{new}.diff').write_text(''.join(difflib.unified_diff(s.splitlines(True),result.splitlines(True),fromfile=f'{old}.inp',tofile=f'{new}.inp')))
