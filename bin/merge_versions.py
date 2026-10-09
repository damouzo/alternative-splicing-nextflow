#!/usr/bin/env python3
"""
Merge every per-process versions.yml emitted during a run into a single
run_info/software_versions.yml.

Nextflow versions.yml files look like:

    "RMATS_PREP":
        rmats: 4.3.0

Input paths are passed as a single ':'-separated list (work-dir paths contain
no ':' on the supported filesystems). Missing or empty files are skipped.

Stdlib only.
"""

import argparse
import sys


def parse_versions(path):
    """Return {process: {tool: version}} for one versions.yml."""
    entries = {}
    process = None
    with open(path, 'r', encoding='utf-8', errors='replace') as handle:
        for line in handle:
            if not line.strip() or line.lstrip().startswith('#'):
                continue
            stripped = line.strip()
            indent = len(line) - len(line.lstrip())
            if stripped.endswith(':'):
                name = stripped[:-1].strip().strip('"')
                # Top-level process keys have no leading indent
                process = name if indent == 0 else process
                if process is not None:
                    entries.setdefault(process, {})
            elif ':' in stripped and process is not None:
                key, value = stripped.split(':', 1)
                entries[process][key.strip()] = value.strip()
    return entries


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--inputs', required=True,
                        help="':' -separated list of versions.yml paths")
    parser.add_argument('--output', default='software_versions.yml')
    args = parser.parse_args()

    merged = {}
    for path in [p for p in args.inputs.split(':') if p]:
        try:
            entries = parse_versions(path)
        except OSError:
            continue
        for process, values in entries.items():
            merged.setdefault(process, {}).update(values)

    lines = ['# Consolidated software versions (one entry per process)',
             'software_versions:']
    for process in sorted(merged):
        lines.append('  %s:' % process)
        for tool in sorted(merged[process]):
            lines.append('    %s: %s' % (tool, merged[process][tool]))
    with open(args.output, 'w', encoding='utf-8') as handle:
        handle.write('\n'.join(lines) + '\n')

    print('[merge_versions] %d processes -> %s' % (len(merged), args.output))
    return 0


if __name__ == '__main__':
    sys.exit(main())
