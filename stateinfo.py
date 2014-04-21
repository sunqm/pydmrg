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
        self.quantaStates = []
        self.allowedQuanta = numpy.ones((0,0),dtype=bool)
        self.unBlockedIndex = []
        self.leftUnMapQuanta = []
        self.rightUnMapQuanta = []

        # point to another StateInfo, (not RawStateInfo)
        self.leftStateInfo = None
        self.rightStateInfo = None
        self.unCollectedStateInfo = None

    def refresh_by(self, rawstateinfo):
        assert(isinstance(rawstateinfo, _dmrg.RawStateInfo))
        self._raw = rawstateinfo
        self._sync_raw2self()
        #how to initialize:
        # self.unCollectedStateInfo self.leftStateInfo self.rightStateInfo

    def init_by_spinquantum(self, n, sq, qs_lst):
        if n > 1 or len(qs_lst) > 1:
            raise ValueError('TODO: initialize StateInfo')
        self._raw = _dmrg.NewRawStateInfo()
        self._raw.init_by_a_spinquantum(sq._raw)

    def _sync_raw2self(self):
        self.totalStates = sum(self._raw.quantaStates)
        self.quantaStates = self._raw.quantaStates
        self.unBlockedIndex = [sum(self.quantaStates[:i]) \
                               for i,_ in enumerate(self.quantaStates)]
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

def CollectQuanta(old_stateinfo):
    # make self as an unCollectedStateInfo for a new StateInfo
    newsi = StateInfo()
    newsi._raw = _dmrg.NewRawStateInfo()
    _dmrg.Pyunion_StateInfo_quanta(newsi._raw, old_stateinfo._raw)
    newsi.unCollectedStateInfo = old_stateinfo
    newsi.set_unCollectedStateInfo(old_stateinfo._raw)
    newsi._sync_raw2self()

    newsi.leftStateInfo = old_stateinfo.leftStateInfo
    newsi.rightStateInfo = old_stateinfo.rightStateInfo
    # NOTE avoid the following two being overwritten by _sync_raw2self
    newsi.leftUnMapQuanta = old_stateinfo.leftUnMapQuanta;
    newsi.rightUnMapQuanta = old_stateinfo.rightUnMapQuanta;
    return newsi


def TensorProduct(a, b, constraint=0):
    '''return StateInfo c, c = (StateInfo a) * (StateInfo b)

    constraint = 0 for NO_PARTICLE_SPIN_NUMBER_CONSTRAINT
    constraint = 1 for PARTICLE_SPIN_NUMBER_CONSTRAINT
    '''
    c = StateInfo()
    c.refresh_by(_dmrg.PyTensorProduct(a._raw, b._raw, constraint))
    return c
