#!/bin/bash
#
# Builds the docs using the qif package from PyPI, without compiling the library. Used by Cloudflare Pages:
#   build command:     bash misc/docs/build.sh
#   output directory:  misc/docs/_build/html
#   env:               PYTHON_VERSION=3.13
#
# Signatures are taken from the .pyi stubs in the repository, so they follow this branch even if the
# installed qif is older. Functions not yet released, however, cannot be documented until the next release.

set -e
cd "$(dirname "$0")"

pip install -r requirements.txt
rm -rf _build
sphinx-build -b html . _build/html
