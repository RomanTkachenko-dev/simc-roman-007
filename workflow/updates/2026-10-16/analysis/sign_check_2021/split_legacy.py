"""Compile legacy procedures separately, preserving their exact source text."""
from pathlib import Path
import re, hashlib
p = Path('/opt/simc/bill-master-8b43f9d/trg_track.f')
s = p.read_bytes()
pattern = rb'(?im)^      SUBROUTINE trgDeriv\([^\n]*\n.*?^      END[ \t]*\r?\n'
m = re.search(pattern, s, flags=re.S)
assert m is not None
body, rest = m.group(), s[:m.start()]+s[m.end():]
assert rest[:m.start()] + body + rest[m.start():] == s
build = p.parent / 'legacy_build'
build.mkdir(exist_ok=True)
(build/'trg_deriv.f').write_bytes(body)
(build/'trg_track_rest.f').write_bytes(rest)
print('Canonical trg_track.f unchanged:', hashlib.sha256(s).hexdigest())
