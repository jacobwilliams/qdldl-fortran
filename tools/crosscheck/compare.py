"""Compares the outputs of the C QDLDL and of the Fortran port (see run.sh)."""
import sys

def read(path):
    lines = open(path).read().split('\n')
    out, i = [], 0
    while i + 6 < len(lines) and lines[i].strip():
        n, s, r = map(int, lines[i].split())
        rec = dict(n=n, sumLnz=s, r=r,
                   etree=[int(v) for v in lines[i+1].split()],
                   Lp=[int(v) for v in lines[i+2].split()],
                   Li=[int(v) for v in lines[i+3].split()],
                   Lx=[float(v) for v in lines[i+4].split()],
                   D=[float(v) for v in lines[i+5].split()],
                   x=[float(v) for v in lines[i+6].split()])
        out.append(rec)
        i += 7
    return out

def reldiff(a, b):
    scale = max([abs(v) for v in a] + [1e-300])
    return max([abs(u - v) for u, v in zip(a, b)] + [0.0]) / scale

c, f = read('c_out.txt'), read('fortran_out.txt')
ok = len(c) == len(f) and len(c) > 0
for k, (a, b) in enumerate(zip(c, f), 1):
    same_int = all(a[key] == b[key] for key in ('n', 'sumLnz', 'r', 'etree', 'Lp', 'Li'))
    d = {key: reldiff(a[key], b[key]) for key in ('Lx', 'D', 'x')}
    good = same_int and all(v <= 1e-13 for v in d.values())
    ok = ok and good
    print(f"matrix {k}: n = {a['n']:5d}, nnz(L) = {a['sumLnz']:7d}, positive pivots = {a['r']:5d}, "
          f"structure {'identical' if same_int else 'DIFFERENT'}, "
          f"rel. diff L {d['Lx']:.1e}, D {d['D']:.1e}, x {d['x']:.1e}  {'ok' if good else 'FAIL'}")
print('cross-check passed' if ok else 'cross-check FAILED')
sys.exit(0 if ok else 1)
