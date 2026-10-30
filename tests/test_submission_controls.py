"""Checks for the new duration and paired-analysis evidence."""
import json
import hashlib
from pathlib import Path
import sys
import unittest
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'code'))
from bd_core import _simulate_from_field

ROOT = Path(__file__).resolve().parents[1]


class SubmissionControls(unittest.TestCase):
    def test_figure_source_manifests_cover_all_five_figures(self):
        outputs=set()
        for name,inputs in [
                ('figure_sources.json',[('summary','summary_sha256'),('noise_validation','noise_sha256')]),
                ('fig5_sources.json',[('summary','summary_sha256'),('paired_sources','paired_sources_sha256')])]:
            record=json.loads((ROOT/'figs/revision_20261026'/name).read_text())
            for source,key in inputs:
                self.assertFalse(Path(record[source]).is_absolute())
                self.assertEqual(hashlib.sha256((ROOT/record[source]).read_text(encoding='utf-8').encode('utf-8')).hexdigest(),record[key])
            for output in record['outputs']:
                self.assertFalse(Path(output).is_absolute())
                self.assertTrue((ROOT/output).is_file())
                outputs.add(output)
        expected={p.relative_to(ROOT).as_posix() for p in (ROOT/'figs/revision_20261026').iterdir()
                  if p.suffix in ['.pdf','.svg','.png']}
        self.assertEqual(outputs,expected)

    def test_nested_prefix_matches_direct_parallel_update(self):
        eta = np.random.default_rng(43).integers(0,2,size=(2,24,8),dtype=np.uint8)
        for T in [6,12,24]:
            profiles=np.zeros((2,1,8),dtype=np.int32)
            widths=np.zeros((2,T),dtype=np.float32)
            _simulate_from_field(eta[:,:T,:],0,8,T,1,np.array([T]),profiles,widths)
            h=np.zeros((2,8),dtype=np.int64)
            expected=[]
            for t in range(T):
                h=np.maximum.reduce([h+eta[:,t,:],np.roll(h,1,axis=1),np.roll(h,-1,axis=1)])
                expected.append(h.var(axis=1))
            np.testing.assert_array_equal(profiles[:,0,:],h)
            np.testing.assert_allclose(widths,np.asarray(expected).T,rtol=1e-6)

    def test_paired_central_values_match_production_fits(self):
        original=json.loads((ROOT/'data/analysis/revision_20261026/summary.json').read_text())
        new=json.loads((ROOT/'data/analysis/submission_20261030/summary.json').read_text())
        for row in new['paired']:
            fits={r['generator']:r['alpha'] for r in original['exponents'] if abs(r['theta']-row['theta'])<1e-4}
            self.assertAlmostEqual(row['delta_alpha'],fits['fgn']-fits['powerlaw'],places=12)
        for row in new['duration']:
            with np.load(ROOT/row['raw_file']) as z:
                self.assertEqual(z['squared_width'].shape,(3,32))
                self.assertEqual(int(z['record_length']),4*int(z['T']))


if __name__=='__main__':
    unittest.main()
