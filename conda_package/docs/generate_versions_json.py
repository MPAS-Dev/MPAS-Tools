#!/usr/bin/env python
import argparse
import json
import os

# Mode: local or production
parser = argparse.ArgumentParser(
    description='Generate versions.json for MPAS-Tools documentation.')
parser.add_argument(
    '--local',
    action='store_true',
    help='Generate versions.json for local build.'
)
args = parser.parse_args()
local = args.local
base_dir = '_build/html' if local else 'gh-pages'
shared_dir = os.path.join(base_dir, 'shared')

entries = []

if not os.path.exists(base_dir) or not os.listdir(base_dir):
    raise FileNotFoundError(
        f"Base directory '{base_dir}' does not exist or is empty.")

versions = sorted(os.listdir(base_dir), reverse=True)
if 'master' in versions:
    versions.insert(0, versions.pop(versions.index('master')))

for name in versions:
    path = os.path.join(base_dir, name)
    if os.path.isdir(path) and name not in ('shared', '.git'):
        entries.append({
            'version': name,
            'url': f'../{name}/' if local else f'/MPAS-Tools/{name}/'
        })

os.makedirs(shared_dir, exist_ok=True)
with open(os.path.join(shared_dir, 'versions.json'), 'w') as f:
    json.dump(entries, f, indent=2)
