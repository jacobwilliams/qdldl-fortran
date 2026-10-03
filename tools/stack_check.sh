#!/usr/bin/env bash
#
# List what the library could put on the stack: its local arrays with
# explicit bounds that aren't constants (automatic arrays), and the array
# temporaries the compiler creates for its expressions (gfortran's
# -Warray-temporaries). Some compilers and flags (gfortran -fstack-arrays or
# -Ofast, the Intel compilers without -heap-arrays) put both on the stack,
# which a large problem then overflows. The library must have neither: the
# script fails if it finds one of either.
#
# usage (from the repository root):
#
#    pixi run stack-check
#
# (Adapted from the same script of sqpopt.)

set -euo pipefail

FPM=${FPM:-fpm}

echo "automatic arrays in src/:"
python3 - <<'PYTHON'
import glob, re, sys
found = 0
for path in sorted(glob.glob('src/*.[fF]90')):
    source = open(path).read()
    # (the file's named constants: bounds made only of them and of numbers are constant)
    constants = set(m.lower() for m in re.findall(r'^\s*integer\s*,\s*parameter[^:]*::\s*(\w+)\s*=', source,
                                                  re.I | re.M))
    for i, line in enumerate(source.split('\n'), 1):
        code = line.split('!')[0]
        if '::' not in code:
            continue
        left, right = code.split('::', 1)
        if re.search(r'intent|allocatable|parameter|pointer', left, re.I):
            continue
        if not re.match(r'\s*(real|integer|logical|character|complex)', left, re.I):
            continue
        m = re.search(r'dimension\s*\((.*)\)', left, re.I)
        dims = [m.group(1)] if m else re.findall(r'\w+\s*\(([^()]*(?:\([^()]*\)[^()]*)*)\)', right.split('=')[0])
        for d in dims:
            if re.fullmatch(r'[\s:,]+', d):
                continue   # (deferred shape)
            if all(w.isdigit() or w.lower() in constants for w in re.findall(r'\w+', d)):
                continue   # (constant bounds)
            found += 1
            print(f'  {path}:{i}: {line.strip()[:110]}')
print(f'  {found} found')
sys.exit(1 if found else 0)
PYTHON

# a build with the warning turned on: fpm keeps each object's compiler output
# in a `.log` file next to it, and the commands in build/compile_commands.json
$FPM build --flag "-O1 -Warray-temporaries -fdiagnostics-plain-output" > /dev/null 2>&1

echo
echo "array temporaries in src/:"
python3 - <<'PYTHON'
import json, re, sys
seen = set()
for entry in json.load(open('build/compile_commands.json')):
    args = entry['arguments']
    if '-Warray-temporaries' not in args or not re.search(r'(^|/)src/', entry['file']):
        continue
    try:
        log = open(args[args.index('-o') + 1] + '.log').read()
    except OSError:
        sys.exit('no compiler log for ' + entry['file'])
    for m in re.finditer(r'^\S*?(src/\w+\.[fF]90):(\d+):\d+: Warning: Creating array temporary', log, re.M):
        seen.add((m.group(1), int(m.group(2))))
for path, line in sorted(seen):
    text = open(path).read().split('\n')[line - 1].strip()
    print(f'  {path}:{line}: {text[:110]}')
print(f'  {len(seen)} found')
sys.exit(1 if seen else 0)
PYTHON
