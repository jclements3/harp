#!/usr/bin/env python3
"""Batch render Mutopia SATB hymns with Aiken 7-shape noteheads.

Finds all SATB LilyPond files in the Mutopia clone, injects \aikenHeads
into each Staff context, updates syntax with convert-ly, and renders PDFs.
"""

import os
import re
import subprocess
import sys
import shutil

LILYPOND = os.path.expanduser("~/local/lilypond-2.24.4/bin/lilypond")
CONVERT_LY = os.path.expanduser("~/local/lilypond-2.24.4/bin/convert-ly")
MUTOPIA = os.path.expanduser("~/projects/harp/mutopia/ftp")
OUTDIR = os.path.expanduser("~/projects/harp/mutopia_shapes")

# All 141 SATB pieces found by the search
SATB_FILES = []

def find_satb_files():
    """Find all .ly files in Mutopia that contain SATB voices."""
    satb_patterns = [
        r'\\partCombine|\\partcombine',
        r'Soprano\s*=|Alto\s*=|Tenor\s*=|Bass\s*=',
        r'soprano\s*=|alto\s*=|tenor\s*=|bass\s*=',
        r'\\new Staff.*\\new Staff.*\\new Staff.*\\new Staff',
    ]
    found = []
    for root, dirs, files in os.walk(MUTOPIA):
        for f in files:
            if not f.endswith('.ly'):
                continue
            path = os.path.join(root, f)
            try:
                with open(path, 'r', encoding='latin-1') as fh:
                    content = fh.read()
                # Check for SATB indicators
                has_satb = False
                for pat in satb_patterns:
                    if re.search(pat, content, re.IGNORECASE):
                        has_satb = True
                        break
                # Also check for at least 2 Staff contexts
                staff_count = len(re.findall(r'\\context\s+Staff|\\new\s+Staff', content))
                if has_satb and staff_count >= 2:
                    found.append(path)
            except Exception:
                continue
    return sorted(found)


def inject_aiken_heads(content):
    """Inject \\aikenHeads into Staff contexts in a LilyPond file."""
    # Strategy: Add \aikenHeads after each \context Staff or \new Staff line
    # that doesn't already have it
    if '\\aikenHeads' in content:
        return content

    # Add to Staff contexts: insert \aikenHeads after << or { following Staff
    # Pattern: \context Staff = "name" << or \new Staff <<
    lines = content.split('\n')
    result = []
    for i, line in enumerate(lines):
        result.append(line)
        stripped = line.strip()
        # If this line defines a Staff context, add \aikenHeads
        if re.search(r'\\(context|new)\s+Staff', stripped):
            # Check if next non-empty line has << or {
            # Insert \aikenHeads after the staff definition
            # Find the right place - after << or after the Staff line
            if '<<' in stripped or '{' in stripped:
                result.append('        \\aikenHeads')
            else:
                # The << might be on the next line, we'll add after it
                pass

    content = '\n'.join(result)

    # Fallback: if we didn't manage to insert it contextually,
    # try a more targeted approach for common patterns
    if '\\aikenHeads' not in content:
        # Insert after \context Staff = XXX <<
        content = re.sub(
            r'(\\(?:context|new)\s+Staff\s*(?:=\s*\w+\s*)?<<)',
            r'\1\n        \\aikenHeads',
            content
        )

    # Another fallback for patterns like \new Staff { or \context Staff {
    if '\\aikenHeads' not in content:
        content = re.sub(
            r'(\\(?:context|new)\s+Staff\s*(?:=\s*\w+\s*)?\{)',
            r'\1\n        \\aikenHeads',
            content
        )

    return content


def fix_obsolete_spacing(content):
    """Remove obsolete spacing directives that crash modern LilyPond."""
    # Remove obsolete-between-system-space lines
    content = re.sub(
        r'obsolete-between-system-space\s*=\s*[^\n]*',
        '',
        content
    )
    # Fix references to staff-space variable
    content = re.sub(
        r'system-system-spacing\.basic-distance\s*=\s*#\(/ obsolete-between-system-space staff-space\)',
        'system-system-spacing.basic-distance = #8',
        content
    )
    content = re.sub(
        r'score-system-spacing\.basic-distance\s*=\s*#\(/ obsolete-between-system-space staff-space\)',
        'score-system-spacing.basic-distance = #8',
        content
    )
    return content


def process_file(src_path, outdir):
    """Process a single LilyPond file: convert syntax, inject shapes, render."""
    # Create a meaningful output name
    rel = os.path.relpath(src_path, MUTOPIA)
    name = rel.replace('/', '_').replace('.ly', '')
    work_dir = os.path.join(outdir, 'work')
    os.makedirs(work_dir, exist_ok=True)

    work_file = os.path.join(work_dir, name + '.ly')

    # Check if this file uses \include for other files
    with open(src_path, 'r', encoding='latin-1') as f:
        content = f.read()

    # If file includes other files, copy the whole directory
    src_dir = os.path.dirname(src_path)
    includes = re.findall(r'\\include\s+"([^"]+)"', content)
    if includes:
        # Copy all files from source directory to work directory
        for inc_file in includes:
            inc_src = os.path.join(src_dir, inc_file)
            if os.path.exists(inc_src):
                inc_dst = os.path.join(work_dir, inc_file)
                os.makedirs(os.path.dirname(inc_dst), exist_ok=True)
                shutil.copy2(inc_src, inc_dst)

    # Step 1: Copy and convert syntax
    shutil.copy2(src_path, work_file)
    subprocess.run([CONVERT_LY, '-e', work_file],
                   capture_output=True, text=True, timeout=30)

    # Step 2: Read converted file, inject \aikenHeads, fix spacing
    with open(work_file, 'r', encoding='latin-1') as f:
        content = f.read()

    content = inject_aiken_heads(content)
    content = fix_obsolete_spacing(content)

    with open(work_file, 'w', encoding='latin-1') as f:
        f.write(content)

    # Step 3: Render
    pdf_path = os.path.join(outdir, 'pdf', name + '.pdf')
    os.makedirs(os.path.join(outdir, 'pdf'), exist_ok=True)

    result = subprocess.run(
        [LILYPOND, '-o', os.path.join(outdir, 'pdf', name), work_file],
        capture_output=True, text=True, timeout=120,
        cwd=work_dir
    )

    if result.returncode == 0 and os.path.exists(pdf_path):
        return True, name
    else:
        return False, f"{name}: {result.stderr[:200] if result.stderr else 'unknown error'}"


def main():
    os.makedirs(OUTDIR, exist_ok=True)

    print("Finding SATB files in Mutopia...")
    files = find_satb_files()
    print(f"Found {len(files)} SATB LilyPond files")

    success = 0
    failed = 0
    errors = []

    for i, path in enumerate(files):
        rel = os.path.relpath(path, MUTOPIA)
        print(f"[{i+1}/{len(files)}] {rel}...", end=' ', flush=True)
        try:
            ok, msg = process_file(path, OUTDIR)
            if ok:
                print("OK")
                success += 1
            else:
                print(f"FAIL")
                errors.append(msg)
                failed += 1
        except Exception as e:
            print(f"ERROR: {e}")
            errors.append(f"{rel}: {e}")
            failed += 1

    print(f"\n{'='*60}")
    print(f"Results: {success} succeeded, {failed} failed out of {len(files)}")
    print(f"PDFs in: {OUTDIR}/pdf/")

    if errors:
        print(f"\nFailed files:")
        for e in errors[:20]:
            print(f"  {e}")


if __name__ == '__main__':
    main()
