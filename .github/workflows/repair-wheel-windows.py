#!/usr/bin/env python

# Kostas
# Script from: https://github.com/Wandercraft/jiminy/blob/648c0ec918ca2c2f734a32bada1dbfaf6226de48/build_tools/wheel_repair.py
# With the following changes:
# - The dlls of the original wheel is used if -d is not provided
# - The original dlls from the source wheel are _not_ copied in the target wheel
# - <name>.libs is prefixes also in the dependencies of dlls
# - msvcp140.dll filename is not hashed, to use the system's version, if available



# This tool has been copied from https://github.com/vinayak-mehta/pdftopng/blob/main/scripts/wheel_repair.py
# and extended to supported hierarchical folder architecture with mulitple .pyd
# to update, and to move all the DLL in a common folder *package*.lib installed
# jointly with the package itself, similarly to auditwheel on Linux platform.
#(see also https://discuss.python.org/t/delocate-auditwheel-but-for-windows/2589/9).

import os
import shutil
import pathlib
import hashlib
import zipfile
import argparse
import tempfile

import pefile
from machomachomangler.pe import redll


def hash_filename(filepath, blocksize=65536):
    # don't hash msvcp140.dll, to use the system's version, if available
    # print("hash_filename", filepath.lower())
    if "msvcp140.dll" in filepath.lower():
        return os.path.basename(filepath)

    hasher = hashlib.sha256()

    with open(filepath, "rb") as afile:
        buf = afile.read(blocksize)
        while len(buf) > 0:
            hasher.update(buf)
            buf = afile.read(blocksize)

    root, ext = os.path.splitext(filepath)
    return f"{os.path.basename(root)}-{hasher.hexdigest()[:8]}{ext}"


def find_dll_dependencies(dll_filepath, lib_dir, tabs=""):
    print(tabs, "find_dll_dependencies", dll_filepath, lib_dir)
    dlls = [x.lower() for x in os.listdir(lib_dir)]
    dll_deps = {}
    for entry in pefile.PE(dll_filepath).DIRECTORY_ENTRY_IMPORT:
        entry_name = entry.dll.decode("utf-8")
        print(tabs, "    entry", entry_name,  entry_name.lower() in dlls )
        if entry_name.lower() in dlls:
            dll_deps.setdefault(
                os.path.basename(dll_filepath), set()).add(entry_name)
            nested_dll_deps = find_dll_dependencies(
                os.path.join(lib_dir, entry_name), lib_dir, tabs+"    ")
            dll_deps.update(nested_dll_deps)
    print(tabs, "   find res", dll_deps)
    return dll_deps


def mangle_filename(old_filename, new_filename, mapping):
    print("mangle", old_filename, new_filename, mapping)
    with open(old_filename, "rb") as f:
        buf = f.read()
    new_buf = redll(buf, mapping)
    with open(new_filename, "wb") as f:
        f.write(new_buf)


parser = argparse.ArgumentParser(
    description="Vendor in external shared library dependencies of a wheel."
)

parser.add_argument("WHEEL_FILE", type=str, help="Path to wheel file")
parser.add_argument(
    "-d", "--dll-dir", dest="DLL_DIR", type=str, help="Directory to find the DLLs"
)
parser.add_argument(
    "-w",
    "--wheel-dir",
    dest="WHEEL_DIR",
    type=str,
    help=('Directory to store delocated wheels (default: "wheelhouse/")'),
    default="wheelhouse/",
)

args = parser.parse_args()

wheel_name = os.path.basename(args.WHEEL_FILE)
repaired_wheel = os.path.join(os.path.abspath(args.WHEEL_DIR), wheel_name)

old_wheel_dir = tempfile.mkdtemp()
new_wheel_dir = tempfile.mkdtemp()
package_name = wheel_name.split("-")[0]
bundle_name = package_name + "/.libs"
bundle_path = os.path.join(new_wheel_dir, bundle_name)
os.makedirs(bundle_path)

if(args.DLL_DIR is None):
    args.DLL_DIR = old_wheel_dir + "/" + package_name

with zipfile.ZipFile(args.WHEEL_FILE, "r") as wheel:
    print("==", wheel.namelist())
    wheel.extractall(old_wheel_dir)
    wheel.extractall(new_wheel_dir, [name for name in wheel.namelist() if not name.lower().endswith(".dll")])
    pyd_rel_paths = [os.path.normpath(path)
                     for path in wheel.namelist() if path.endswith(".pyd")]

dll_dependencies = {}
for rel_path in pyd_rel_paths:
    print("rel_path", rel_path)
    abs_path = os.path.join(old_wheel_dir, rel_path)
    dll_dependencies.update(find_dll_dependencies(abs_path, args.DLL_DIR))

for dll, dependencies in dll_dependencies.items():
    print("dll", dll)
    mapping = {}

    if dll.endswith(".pyd"):
        rel_path = next(path for path in pyd_rel_paths if path.endswith(dll))

    for dep in dependencies:
        print("    dll dep", dep)
        hashed_name = hash_filename(os.path.join(args.DLL_DIR, dep))  # already basename
        
        # Kostas: the original version, only the dependencies of .pyd files were set to <module>.libs/<dllname>
        # But the dependecies of .dll files also need the <module>.libs/ prefix to be loaded
        #
        #if dll.endswith(".pyd"):
        bundle_rel_path = os.path.join(
            "..\\" * rel_path.count(os.path.sep), bundle_name)
        mapping[dep.encode("ascii")] = os.path.join(
            bundle_rel_path, hashed_name).encode("ascii")
        #else:
        #    mapping[dep.encode("ascii")] = hashed_name.encode("ascii")

        shutil.copy(
            os.path.join(args.DLL_DIR, dep),
            os.path.join(bundle_path, hashed_name))

    if dll.endswith(".pyd"):
        old_name = os.path.join(old_wheel_dir, rel_path)
        new_name = os.path.join(new_wheel_dir, rel_path)
    else:
        old_name = os.path.join(args.DLL_DIR, dll)
        hashed_name = hash_filename(os.path.join(args.DLL_DIR, dll))  # already basename
        new_name = os.path.join(bundle_path, hashed_name)

    mangle_filename(old_name, new_name, mapping)

pathlib.Path(os.path.dirname(repaired_wheel)).mkdir(parents=True, exist_ok=True)
with zipfile.ZipFile(repaired_wheel, "w", zipfile.ZIP_DEFLATED) as new_wheel:
    for root, dirs, files in os.walk(new_wheel_dir):
        new_root = os.path.relpath(root, new_wheel_dir)
        for file in files:
            new_wheel.write(
                os.path.join(root, file), os.path.join(new_root, file))
