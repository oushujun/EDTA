#!/bin/bash
# Install parallelized ProcessRepeats into the current conda env's RepeatMasker.
#
# What it does:
#   1. Backs up ProcessRepeats and RepeatMasker as .ori
#   2. Installs the parallelized ProcessRepeats (adds -pa N option)
#   3. Patches RepeatMasker to pass its -pa value through to ProcessRepeats
#
# Usage (portable: run from wherever this folder lives — the script installs
#        the ProcessRepeats sitting next to it):
#   conda activate <env>        # e.g. EDTA2.3
#   bash /path/to/install_ProcessRepeats.sh
#
# To restore originals:
#   cp $CONDA_PREFIX/share/RepeatMasker/ProcessRepeats.ori $CONDA_PREFIX/share/RepeatMasker/ProcessRepeats
#   cp $CONDA_PREFIX/share/RepeatMasker/RepeatMasker.ori   $CONDA_PREFIX/share/RepeatMasker/RepeatMasker

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
SRC="$SCRIPT_DIR/ProcessRepeats"

# --- Detect current conda env ---
if [ -z "${CONDA_PREFIX:-}" ]; then
    echo "Error: no conda env active. Run 'conda activate <env>' first." >&2
    exit 1
fi

RM_DIR="$CONDA_PREFIX/share/RepeatMasker"
PR="$RM_DIR/ProcessRepeats"
RM="$RM_DIR/RepeatMasker"

if [ ! -f "$PR" ]; then
    echo "Error: ProcessRepeats not found at $PR" >&2
    exit 1
fi
if [ ! -f "$SRC" ]; then
    echo "Error: source ProcessRepeats not found at $SRC" >&2
    exit 1
fi

ENV_NAME="$(basename "$CONDA_PREFIX")"
echo "Target env: $ENV_NAME ($CONDA_PREFIX)"

# --- Check RM version ---
RM_VER="unknown"
if [ -f "$RM_DIR/RepeatMaskerConfig.pm" ]; then
    RM_VER=$(grep -m1 -oP 'VERSION\s*=\s*"\K[^"]+' "$RM_DIR/RepeatMaskerConfig.pm" || echo "unknown")
fi
echo "RepeatMasker version: $RM_VER"

# --- Check ProcessRepeats compatibility ---
# Our modified ProcessRepeats is based on 4.2.2. Check if the target is close enough.
# Compare against original backup if one exists (from a prior install).
ORIG_PR="$PR"
[ -f "$PR.ori" ] && ORIG_PR="$PR.ori"
[ -f "$PR.orig" ] && ORIG_PR="$PR.orig"

DIFF_LINES=$(diff "$SRC" "$ORIG_PR" 2>/dev/null | grep -c '^[<>]' || true)
# Our patch adds ~130 changed lines vs 4.2.2. If diff is much larger, the base version differs.
if [ "$DIFF_LINES" -gt 250 ]; then
    echo "Warning: ProcessRepeats in $ENV_NAME (RM $RM_VER) differs significantly from the"
    echo "  parallelized version ($DIFF_LINES changed lines). This version was built for RM 4.2.x."
    if [ "${FORCE:-0}" = "1" ]; then
        echo "  FORCE=1 set — installing anyway."
    elif { exec 3</dev/tty; } 2>/dev/null; then
        read -r -p "  Install anyway? [y/N] " answer <&3 || answer=""
        exec 3<&-
        [ "$answer" = "y" ] || [ "$answer" = "Y" ] || { echo "Aborted."; exit 1; }
    else
        echo "  No TTY for confirmation; aborting. Re-run with FORCE=1 to install anyway." >&2
        exit 1
    fi
fi

# --- Step 1: Back up ProcessRepeats ---
if [ ! -f "$PR.ori" ]; then
    cp -p "$PR" "$PR.ori"
    echo "Backed up ProcessRepeats -> ProcessRepeats.ori"
else
    echo "ProcessRepeats.ori already exists, skipping backup"
fi

# --- Step 2: Install parallelized ProcessRepeats ---
# --remove-destination breaks any conda hardlink first, so we never write through
# a shared inode into the pkgs cache / sibling envs (their originals stay intact).
cp --remove-destination "$SRC" "$PR"
chmod +x "$PR"
echo "Installed parallelized ProcessRepeats"

# --- Step 3: Patch RepeatMasker to pass -pa to ProcessRepeats ---
if [ ! -f "$RM" ]; then
    echo "Warning: RepeatMasker script not found at $RM, skipping patch"
    exit 0
fi

# Check if already patched
if grep -q 'prOptions.*-pa.*parallel' "$RM"; then
    echo "RepeatMasker already patched to pass -pa to ProcessRepeats"
else
    # Back up RepeatMasker
    if [ ! -f "$RM.ori" ]; then
        cp -p "$RM" "$RM.ori"
        echo "Backed up RepeatMasker -> RepeatMasker.ori"
    fi

    # Insert -pa passthrough right before the compressCatFile block.
    # Anchor: the last $prOptions assignment before ".cat" is appended.
    # Works across RM 4.1.8 - 4.2.3 (all have this exact pattern).
    perl -i -pe '
        if (/\$prOptions\s*\.=\s*"-xsmall/ && !$done) {
            $_ .= "    \$prOptions .= \"-pa \" . \$options{'\''parallel'\''} . \" \" if (\$options{'\''parallel'\''});\n";
            $done = 1;
        }
    ' "$RM"

    # Verify the patch took effect
    if grep -q 'prOptions.*-pa.*parallel' "$RM"; then
        echo "Patched RepeatMasker to pass -pa to ProcessRepeats"
    else
        echo "Error: failed to patch RepeatMasker. Restoring backup." >&2
        cp -p "$RM.ori" "$RM"
        exit 1
    fi
fi

echo ""
echo "Done. ProcessRepeats will now run in parallel when RepeatMasker is called with -pa N."
echo "  Example: RepeatMasker -pa 8 -lib library.fa genome.fa"
