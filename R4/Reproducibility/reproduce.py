"""Replay R4 from the compact input snapshot and verify the resulting tables."""
import argparse
import json
from pathlib import Path
import subprocess
import sys

import pandas as pd

SKIP = {'icgc_sample_crosswalk.csv', 'original_module_definitions.csv'}
HERE = Path(__file__).resolve().parent

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--package', type=Path, default=HERE.parent)
    parser.add_argument('--out', type=Path, required=True, help='A new output directory')
    args = parser.parse_args()
    package = args.package.resolve()
    output = args.out.resolve()
    if output.exists():
        parser.error('Choose a new output directory')
    subprocess.run([sys.executable,str(HERE/'analyze.py'),'--snapshot',
                    str(package/'Analysis_Inputs'),'--out',str(output)],check=True)
    checked = []
    for original in sorted((package/'Supplementary_Data').rglob('*.csv')):
        if original.name in SKIP:
            continue
        generated = output/original.name
        pd.testing.assert_frame_equal(
            pd.read_csv(original,float_precision='round_trip'),
            pd.read_csv(generated,float_precision='round_trip'),
            check_dtype=False,check_exact=False,rtol=1e-8,atol=1e-10)
        checked.append(original.name)
    expected = json.loads((package/'Reproducibility'/'summary.json').read_text())
    actual = json.loads((output/'summary.json').read_text())
    if actual != expected:
        raise ValueError('Recomputed summary differs from the supplied summary')
    (output/'replay_verification.json').write_text(json.dumps(
        {'verified_numeric_tables':len(checked),'summary_matches':True,'files':checked},indent=2))
    print(f'Verified {len(checked)} numerical tables and the analysis summary.')

if __name__=='__main__':
    main()
