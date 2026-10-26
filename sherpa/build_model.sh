#!/usr/bin/env bash
# Build and install the Sherpa model library for ChiPT_UFO, with the
# form-factor plugin FormFactors.C compiled in.
#
# Usage: sherpa/build_model.sh [path/to/Sherpa-generate-model]
#   defaults to Sherpa-generate-model on $PATH; NCORE sets the build jobs (4).
set -euo pipefail

here=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
ufo=$(dirname "$here")/ChiPT_UFO
out=$ufo/.sherpa
gen=${1:-Sherpa-generate-model}
ncore=${NCORE:-4}

# Sherpa-generate-model writes param_ChiPT_UFO.dat to the working directory;
# keep it inside the (ignored) .sherpa directory.
mkdir -p "$out"
(cd "$out" && "$gen" --overwrite --ncore "$ncore" "$ufo")

cp "$here/FormFactors.C" "$out/"
cmake_file=$out/CMakeLists.txt
if ! grep -q 'FormFactors.C' "$cmake_file"; then
    sed -i 's/^add_library(SherpaChiPT_UFO SHARED /&FormFactors.C /' "$cmake_file"
    grep -q 'SHARED FormFactors.C ' "$cmake_file" \
        || { echo "build_model.sh: add_library(SherpaChiPT_UFO ...) not found in $cmake_file" >&2; exit 1; }
    # No -ffast-math: FormFactors.C caches with a NaN sentinel.
    echo 'set_source_files_properties(FormFactors.C PROPERTIES COMPILE_OPTIONS "-g;-O2")' >> "$cmake_file"
fi

# Sherpa-generate-model ignores build failures, so rebuild and install here.
cmake -S "$out" -B "$out/build"
cmake --build "$out/build" -j"$ncore"
cmake --install "$out/build"
