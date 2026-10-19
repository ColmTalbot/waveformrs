#!/bin/bash
# Build liblalsimulation 6.2.0 with the fix of lalsuite merge request !2593
# ("Fixing inconsistency in antisymmetric phase for XO4a/XPNR that affects
# parity symmetry", merged into lalsuite master as 64d2004def) for the
# imrphenomxpnr_parityfix* reference suites of generate_reference_data.py.
#
# The library is built without autotools or SWIG from `git archive
# lalsimulation-v6.2.0` plus `git diff 64d2004def^1 64d2004def --
# lalsimulation/lib`, against the lal of the conda environment, with the
# flags of the conda-forge osx-arm64 build of lalsimulation 6.2.0 (-O2
# -ftree-vectorize, OpenMP). With conda-forge's clang 19.1.7 for
# arm64-apple-darwin20.0.0 (CC below) the unpatched 6.2.0 sources reproduce
# the conda-forge liblalsimulation 6.2.0 bit for bit on the IMRPhenomXPNR
# reference suites; Apple clang 15 differs from it by ~1e-9.
#
# usage: CC=/path/to/clang-19 scripts/build_patched_lalsimulation.sh LALSUITE_GIT_DIR OUT_DIR [--unpatched]
#
#   CONDA_PREFIX  the environment with lal 7.7.0, gsl and llvm-openmp
#   CC            the C compiler (default clang)
set -euo pipefail
LALSUITE=$1
OUT=$2
PATCHED=1
if [ "${3:-}" = "--unpatched" ]; then PATCHED=0; fi
E=${CONDA_PREFIX:?set CONDA_PREFIX to the environment with lal}
CC=${CC:-clang}
TAG=lalsimulation-v6.2.0
FIX=64d2004def

mkdir -p "$OUT/src" "$OUT/obj" "$OUT/include/lal"
OUT=$(cd "$OUT" && pwd)
git -C "$LALSUITE" archive --format=tar "$TAG" lalsimulation | tar -x -C "$OUT/src"
SRC=$OUT/src/lalsimulation
TAG_COMMIT=$(git -C "$LALSUITE" rev-parse "$TAG^{commit}")
if [ "$PATCHED" = 1 ]; then
  git -C "$LALSUITE" diff "$FIX^1" "$FIX" -- lalsimulation/lib > "$OUT/lalsuite-2593.patch"
  patch -d "$OUT/src" -p1 < "$OUT/lalsuite-2593.patch"
  CLEAN=UNCLEAN
  STATUS="UNCLEAN: $TAG + lalsuite!2593 (git diff $FIX^1 $FIX -- lalsimulation/lib)"
else
  CLEAN=CLEAN
  STATUS="CLEAN: $TAG rebuilt"
fi

# the installed headers from the source tree, the generated ones of 6.2.0
for h in $(awk '/^pkginclude_HEADERS = /{f=1;next} f&&/^$/{f=0} f{print $1}' "$SRC/lib/Makefile.am"); do
  if [ -f "$SRC/lib/$h" ]; then ln -sf "$SRC/lib/$h" "$OUT/include/lal/$h"; fi
done
cp "$E/include/lal/LALSimulationConfig.h" "$E/include/lal/LALSimulationVCSInfo.h" "$OUT/include/lal/"
echo '#define HAVE_UNISTD_H 1' > "$OUT/config.h"

COMPILER=$("$CC" --version | head -1)
CFLAGS="-target arm64-apple-darwin20.0.0 -ftree-vectorize -fPIC -fstack-protector-strong -O2 -pipe"
if [ -n "${SDKROOT:-}" ] || SDKROOT=$(xcrun --show-sdk-path 2>/dev/null); then
  CFLAGS="$CFLAGS -isysroot $SDKROOT"
fi
CPPFLAGS="-DLAL_FALLBACK_DATA_PATH=\"$E/share/lalsimulation\" -D_FORTIFY_SOURCE=2 -mmacosx-version-min=11.0"
CPPFLAGS="$CPPFLAGS -Xclang -fopenmp -I$OUT -I$OUT/include -I$SRC/lib -isystem $E/include"

cat > "$OUT/LALSimulationVCSInfo.c" <<EOF
#include <stdlib.h>
#include <lal/LALVCSInfoType.h>
#include <lal/LALVCSInfoHeader.h>
#include <lal/LALSimulationVCSInfo.h>
#include <lal/LALSimulationConfig.h>
const LALVCSInfo lalSimulationVCSInfo = {
  .name = "LALSimulation",
  .version = LALSIMULATION_VERSION,
  .vcsId = "$TAG_COMMIT",
  .vcsDate = "",
  .vcsBranch = "None",
  .vcsTag = "$TAG",
  .vcsAuthor = "",
  .vcsCommitter = "",
  .vcsClean = "$CLEAN",
  .vcsStatus = "$STATUS",
  .configureArgs = "scripts/build_patched_lalsimulation.sh: $COMPILER -O2 -ftree-vectorize -Xclang -fopenmp",
  .configureDate = "",
  .buildDate = "",
};
const LALVCSInfo lalSimulationVCSIdentInfo = {
  .name = "\$LALSimulationName: LALSimulation \$",
  .vcsStatus = "\$LALSimulationVCSStatus: $STATUS \$",
};
const LALVCSInfoList lalSimulationVCSInfoList = { &lalVCSInfo, &lalSimulationVCSInfo, NULL };
EOF

SOURCES=$(awk '/^liblalsimulation_la_SOURCES = /{f=1;next} f&&/END_OF_LIST/{f=0} f{print $1}' "$SRC/lib/Makefile.am" | grep '\.c$')
cd "$OUT/obj"
export CC CFLAGS CPPFLAGS SRC
printf '%s\n' $SOURCES | xargs -P 8 -n 1 sh -c '$CC $CFLAGS $CPPFLAGS -w -c "$SRC/lib/$0" -o "$(basename "$0" .c).o"'
"$CC" $CFLAGS $CPPFLAGS -c "$OUT/LALSimulationVCSInfo.c" -o LALSimulationVCSInfo.o
"$CC" $CFLAGS -dynamiclib -mmacosx-version-min=11.0 -o "$OUT/liblalsimulation.37.dylib" ./*.o \
  -install_name "$OUT/liblalsimulation.37.dylib" \
  -L"$E/lib" -Wl,-rpath,"$E/lib" -llal -llalsupport -lgsl -lomp -lm
echo "built $OUT/liblalsimulation.37.dylib ($STATUS)"
