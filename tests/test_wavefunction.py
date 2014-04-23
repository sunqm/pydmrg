#!/usr/bin/env python

import unittest
import numpy
import pydmrg

class KnowValues(unittest.TestCase):
    def test_load(self):
        wfn = pydmrg.Wavefunction()
        wfn.load(0, 3, 1, prefix='/dev/shm/')
        self.assertEqual(wfn.deltaQuantum.particleNumber, 6)
        self.assertEqual(wfn.deltaQuantum.totalSpin, 0)
        self.assertEqual(wfn.stateInfo.totalStates, 55)
        mat = wfn.stateInfo.allowedQuanta.copy()
        mat[1,1] = mat[2,5] = mat[3,0] = mat[3,6] = False
        self.assertTrue(not mat.any())

        diff = wfn.stateInfo.quantaStates - numpy.array([6, 36, 7, 6],dtype=int)
        self.assertEqual(abs(diff).sum(), 0)
        diff = wfn.stateInfo.leftUnMapQuanta - numpy.array([2, 6, 7, 9],dtype=int)
        self.assertEqual(abs(diff).sum(), 0)
        diff = wfn.stateInfo.rightUnMapQuanta - numpy.array([2, 6, 7, 9],dtype=int)
        self.assertEqual(abs(diff).sum(), 0)

        spinquanta = wfn.stateInfo.get_quanta(0)
        self.assertEqual(spinquanta.particleNumber, 6)
        self.assertEqual(spinquanta.totalSpin, 0)
        self.assertEqual(spinquanta.irrep, 0)
        self.assertEqual(len(wfn.stateInfo.get_quantaMap(1,1)), 0)


if __name__ == "__main__":
    print "Full Tests for Lattice DMET"
    unittest.main()

