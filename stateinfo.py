#
# File: stateinfo.py
# Author: Qiming Sun <qimings@princeton.edu>
#

import numpy
import _dmrg
import quanta

class StateInfo(object):
    def __init__(self):
        self.totalStates = 0
        self.quantaStates = 0
        self.allowedQuanta = numpy.ones((0,0),dtype=bool)
        self.leftUnMapQuanta = []
        self.rightUnMapQuanta = []

    def refresh_by(self, rawstateinfo):
        assert(isinstance(rawstateinfo, _dmrg.RawStateInfo))
        self._raw = rawstateinfo
        self._sync_raw2self()

    def _sync_raw2self(self):
        self.totalStates = self._raw.totalStates
        self.quantaStates = self._raw.quantaStates
        self.allowedQuanta = self._raw.get_whole_allowedQuanta()
        self.leftUnMapQuanta = self._raw.leftUnMapQuanta
        self.rightUnMapQuanta = self._raw.leftUnMapQuanta

    def _sync_self2raw(self):
        pass

    def get_quanta(self, quanta_id):
        spinquanta = quanta.SpinQuantum()
        spinquanta.refresh_by(self._raw.get_quanta(quanta_id))
        return spinquanta

    def get_quantaMap(self, lquanta_id, rquanta_id):
        '''return a list'''
        return self._raw.get_quantaMap(lquanta_id, rquanta_id)

    def get_allowedQuanta(self, lquanta_id, rquanta_id):
        '''return True/False'''
        return self._raw.get_allowedQuanta(lquanta_id, rquanta_id)
