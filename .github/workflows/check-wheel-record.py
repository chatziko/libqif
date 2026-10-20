#!/usr/bin/env python
# Check that the files in the given wheels match their RECORD (required by PyPI)

import sys, glob, zipfile, csv, hashlib, base64, io

ok = True
for path in [p for arg in sys.argv[1:] for p in glob.glob(arg)]:
    with zipfile.ZipFile(path) as z:
        files = {n for n in z.namelist() if not n.endswith('/')}
        record = next(n for n in files if n.endswith('.dist-info/RECORD'))
        rows = {r[0]: r for r in csv.reader(io.TextIOWrapper(z.open(record), 'utf-8')) if r}

        problems = [f'in RECORD but missing: {n}' for n in rows.keys() - files]
        for n in sorted(files - {record}):
            if n not in rows:
                problems.append(f'not in RECORD: {n}')
                continue
            algo, _, digest = rows[n][1].partition('=')
            if base64.urlsafe_b64encode(hashlib.new(algo, z.read(n)).digest()).rstrip(b'=').decode() != digest:
                problems.append(f'hash mismatch: {n}')

    print(f'{path}: {"OK" if not problems else "RECORD mismatch"}')
    for p in problems:
        print(f'    {p}')
    ok = ok and not problems

sys.exit(0 if ok else 1)
