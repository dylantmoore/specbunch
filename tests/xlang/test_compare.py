"""Focused tests of the selected-answer contract and exhaustive reporting."""
import copy
import gzip
import json
import struct
import tempfile
import unittest
from pathlib import Path
from compare import compare, CrossReport


def number(x):
    x=float(x)
    return dict(hex=x.hex(), bits=struct.pack('>d', x).hex())


def result():
    return dict(status=0, strings={'reason':'', 'flags_str':''}, scalars={
        'B':number(.2), 'sigma':number(.04), 'mu':number(0), 'tau':number(.3),
        'tier':number(1), 'flags':number(0), 'window_lo':number(10),
        'window_hi':number(20), 'band_lo':number(10), 't_lo':number(10),
        'fit_p':number(.1), 'rank':number(1)}, vectors={'windows.B':[number(.2)]})


class ContractTests(unittest.TestCase):
    def test_protected_fields(self):
        expected=result()
        for key in ('B','sigma','mu','tau','tier','flags','window_lo','window_hi','band_lo','t_lo'):
            actual=copy.deepcopy(expected);actual['scalars'][key]=number(.99)
            with self.subTest(key=key), self.assertRaises(AssertionError): compare(actual,expected)
        for group,key,value in [('strings','reason','no_coherent_window'),('strings','flags_str','poor_fit')]:
            actual=copy.deepcopy(expected);actual[group][key]=value
            with self.assertRaises(AssertionError):compare(actual,expected)
        actual=copy.deepcopy(expected);actual['status']=1
        with self.assertRaises(AssertionError):compare(actual,expected)

    def test_estimate_tolerance_and_strict_local(self):
        expected=result();actual=copy.deepcopy(expected)
        actual['scalars']['B']=number(.2+1e-9);actual['scalars']['mu']=number(9e-13)
        compare(actual,expected)
        with self.assertRaises(AssertionError):compare(actual,expected,exact=True)
        actual['scalars']['mu']=number(2e-12)
        with self.assertRaises(AssertionError):compare(actual,expected)

    def test_diagnostics_and_lengths_reported(self):
        expected=result();actual=copy.deepcopy(expected)
        actual['scalars']['fit_p']=number(.99);actual['scalars']['rank']=number(4)
        actual['vectors']['windows.B']=[number(.3),number(.4)]
        compare(actual,expected)
        with tempfile.TemporaryDirectory() as tmp:
            report=CrossReport(tmp)
            compare(actual,expected,report=report,fixture='diagnostic')
            # A later selected flip must not suppress earlier diagnostic evidence.
            actual['scalars']['flags']=number(1)
            compare(actual,expected,report=report,fixture='flip')
            compare(actual,expected,report=report,fixture='flip')
            with self.assertRaises(AssertionError):report.assert_pass()
            counts=json.loads((Path(tmp)/'counts.json').read_text())
            self.assertEqual(counts['selected_answer_flips'],1)
            self.assertEqual(counts['fixtures'],2)
            with gzip.open(Path(tmp)/'differences.jsonl.gz','rt') as f:rows=[json.loads(l) for l in f]
            self.assertEqual(counts['differences'],len(rows))
            self.assertTrue(any(r['field']=='vectors.windows.B[1]' for r in rows))
            self.assertIn('Selected-answer flips: 1',(Path(tmp)/'summary.md').read_text())


if __name__=='__main__':unittest.main()
