import os
import glob
import FinalStateAnalysis.TagAndProbe.HetauCorrection as HetauCorrection
import FinalStateAnalysis.TagAndProbe.MuonPOGCorrections as MuonPOGCorrections
import FinalStateAnalysis.TagAndProbe.PileupWeight as PileupWeight
from FinalStateAnalysis.PlotTools.decorators import memo, memo_last

@memo
def getVar(name, var):
    return name+var

is7TeV = bool('7TeV' in os.environ['jobid'])
pu_distributions  = {
    'singlee'  : glob.glob(os.path.join( 'inputs', os.environ['jobid'], 'data_SingleElectron*pu.root')),
    'singlem'  : glob.glob(os.path.join( 'inputs', os.environ['jobid'], 'data_SingleMuon*pu.root'))}
pu_distributionsUp  = {
    'singlee'  : glob.glob(os.path.join( 'inputs', os.environ['jobid'], 'data_SingleElectron*pu_up.root')),
    'singlem'  : glob.glob(os.path.join( 'inputs', os.environ['jobid'], 'data_SingleMuon*pu_up.root'))}
pu_distributionsDown  = {
    'singlee'  : glob.glob(os.path.join( 'inputs', os.environ['jobid'], 'data_SingleElectron*pu_down.root')),
    'singlem'  : glob.glob(os.path.join( 'inputs', os.environ['jobid'], 'data_SingleMuon*pu_down.root'))}
mc_pu_tag                  = 'Asympt25ns'#'S6' if is7TeV else 'S10'


def make_puCorrector(dataset, kind=None):
    'makes PU reweighting according to the pu distribution of the reference data and the MC, MC distribution can be forced'
    if not kind:
        kind = mc_pu_tag
    weights = []
    if dataset in pu_distributions:# and dataset in pu_distributionsUp and dataset in pu_distributionsDown:
        #print 'pile up corrector', dataset, mc_pu_tag
        #print PileupWeight.PileupWeight( mc_pu_tag, *(pu_distributions[dataset]))
        return PileupWeight.PileupWeight( mc_pu_tag , *(pu_distributions[dataset]))
#        weights = (PileupWeight.PileupWeight( 'S6' if is7TeV else 'S10', *(pu_distributions[dataset])), PileupWeight.PileupWeight( 'S6' if is7TeV else 'S10', *(pu_distributionsUp[dataset])), PileupWeight.PileupWeight( 'S6' if is7TeV else 'S10', *(pu_distributionsDown[dataset])))
#        return weights
    else:
        raise KeyError('dataset not present. Please check the spelling or add it to mcCorrectors.py')

def make_puCorrectorUp(dataset, kind=None):
    'makes PU reweighting according to the pu distribution of the reference data and the MC, MC distribution can be forced'
    if not kind:
        kind = mc_pu_tag
    pileupid = 'Asympt25ns'
    if dataset in pu_distributions:
        #print 'pile up corrector up', dataset
        return PileupWeight.PileupWeight(  mc_pu_tag , *(pu_distributionsUp[dataset]))
        #return PileupWeight.PileupWeight( mc_pu_tag , *(pu_distributions[dataset]))
    else:
        raise KeyError('dataset not present. Please check the spelling or add it to mcCorrectors.py')

def make_puCorrectorDown(dataset, kind=None):
    'makes PU reweighting according to the pu distribution of the reference data and the MC, MC distribution can be forced'
    if not kind:
        kind = mc_pu_tag
    if dataset in pu_distributions:
        #print 'pile up corrector down', dataset
        return PileupWeight.PileupWeight(  mc_pu_tag , *(pu_distributionsDown[dataset]))
        #return PileupWeight.PileupWeight( mc_pu_tag, *(pu_distributions[dataset]))
    else:
        raise KeyError('dataset not present. Please check the spelling or add it to mcCorrectors.py')

##put here the trigger correction as in https://github.com/mverzett/UWHiggs/blob/WH_At_Paper/wh/mcCorrectors.py
"""
correct_e = HetauCorrection.correct_hamburg_e
correct_eid13_mva = HetauCorrection.correct_eid13_mva
correct_eiso13_mva = HetauCorrection.correct_eiso13_mva
correct_eid13_p1s_mva = HetauCorrection.correct_eid13_p1s_mva
correct_eiso13_p1s_mva = HetauCorrection.correct_eiso13_p1s_mva
correct_eid13_m1s_mva = HetauCorrection.correct_eid13_m1s_mva
correct_eiso13_m1s_mva = HetauCorrection.correct_eiso13_m1s_mva
correct_eid_mva = HetauCorrection.scale_eleId_hww
correct_eReco_mva = HetauCorrection.scale_elereco_hww
correct_eIso_mva = HetauCorrection.scale_eleIso_hww
correct_trigger_mva = HetauCorrection.single_ele_mva
"""

### 

def get_electron_corrections(row,*args):
    'makes corrections to iso and id of electrons'
    ret = 1.
    for arg in args:
        abseta = abs(getattr(row, '%sEta' % arg))
        pt     = getattr(row, '%sPt'  % arg)
        ret   *= correct_e(pt,abseta)
    return ret

def get_electronId_corrections_MVA(row, *args):
    ret = 1.
    for arg in args:
        eta  = getattr(row, '%sEta' % arg)
        pt   = getattr(row, '%sPt'  % arg)
        ret *= correct_eid_mva(pt,eta)[0]
    return ret
def get_electronReco_corrections_MVA(row, *args):
    ret = 1.
    for arg in args:
        eta  = getattr(row, '%sEta' % arg)
        pt   = getattr(row, '%sPt'  % arg)
        ret *= correct_eReco_mva(pt,eta)[0]
    return ret
def get_electronIso_corrections_MVA(row, *args):
    ret = 1.
    for arg in args:
        eta  = getattr(row, '%sEta' % arg)
        pt   = getattr(row, '%sPt'  % arg)
        ret *= correct_eIso_mva(pt,eta)[0]
    return ret

def get_trigger_corrections_MVA(row, *args):
    ret = 1.
    for arg in args:
        eta  = getattr(row, '%sEta' % arg)
        pt   = getattr(row, '%sPt'  % arg)
        ret *= correct_trigger_mva(pt,eta)[0]
    return ret
def get_trigger_corrections_p1s_MVA(row, *args):
    ret = 1.
    for arg in args:
        eta  = getattr(row, '%sEta' % arg)
        pt   = getattr(row, '%sPt'  % arg)
        ret *= correct_trigger_mva(pt,eta)[0]+correct_trigger_mva(pt,eta)[1]
    return ret
def get_trigger_corrections_m1s_MVA(row, *args):
    ret = 1.
    for arg in args:
        eta  = getattr(row, '%sEta' % arg)
        pt   = getattr(row, '%sPt'  % arg)
        ret *= correct_trigger_mva(pt,eta)[0]-correct_trigger_mva(pt,eta)[1]

    return ret
        
def get_electronId_corrections13_MVA(row, *args):
    ret = 1.
    for arg in args:
        eta  = getattr(row, '%sEta' % arg)
        pt   = getattr(row, '%sPt'  % arg)
        ret *= correct_eid13_mva(pt,eta)
    return ret
def get_electronIso_corrections13_MVA(row, *args):
    ret = 1.
    for arg in args:
        eta  = getattr(row, '%sEta' % arg)
        pt   = getattr(row, '%sPt'  % arg)
        ret *= correct_eiso13_mva(pt,eta)
    return ret
def get_electronId_corrections13_p1s_MVA(row, *args):
    ret = 1.
    for arg in args:
        eta  = getattr(row, '%sEta' % arg)
        pt   = getattr(row, '%sPt'  % arg)
        ret *= correct_eid13_p1s_mva(pt,eta)
    return ret
def get_electronIso_corrections13_p1s_MVA(row, *args):
    ret = 1.
    for arg in args:
        eta  = getattr(row, '%sEta' % arg)
        pt   = getattr(row, '%sPt'  % arg)
        ret *= correct_eiso13_p1s_mva(pt,eta)
    return ret
def get_electronId_corrections13_m1s_MVA(row, *args):
    ret = 1.
    for arg in args:
        eta  = getattr(row, '%sEta' % arg)
        pt   = getattr(row, '%sPt'  % arg)
        ret *= correct_eid13_m1s_mva(pt,eta)
    return ret
def get_electronIso_corrections13_m1s_MVA(row, *args):
    ret = 1.
    for arg in args:
        eta  = getattr(row, '%sEta' % arg)
        pt   = getattr(row, '%sPt'  % arg)
        ret *= correct_eiso13_m1s_mva(pt,eta)
    return ret

def make_shifted_weights(default, shifts, functors):
    '''make_shifted_weights(default, shifts, functors) --> functor
    takes as imput the central value functor and two lists
    the name of the shifts and the shifted functors
    the returned functor takes one additional string to
    select the shifted functor. If the shift kwarg is missing
    or does not match any shift tag the central (default)
    fuctor output is returned'''
    #make default fetching faster
    default = default 
    def functor(*args, **kwargs):
        shift = ''
        if 'shift' in kwargs:
            print "Shift ",shift
            shift = kwargs['shift']
            del kwargs['shift']

            #check if to apply shifts
            for tag, fcn in zip(shifts, functors):
                if tag == shift:
                    #print 'check if to apply shifts',shift
                    return fcn(*args, **kwargs)
        #print 'def functor',kwargs, default
        return default(*args, **kwargs)
    #print 'functor', shifts, functor
    return functor

def make_multiple(fcn,mu_id="not_muon_or_not_isocorr", indexed=False, shift=0):
    '''make_multiple(fcn, indexed=True, shift=0) --> functor
    takes as imput a weight correction function of pt and eta
    and returns a functor multiple(row,*args) --> weight
    where *args are the base name of the objects upon which 
    compute the correction.

    If indexed is true means that fcn returns a tuple 
    (weight, error) shift +/-1 makes the functor return the 
    weight+/-error in this case'''
    def multiple(row,*args):
        ret = 1.
        for arg in args:
            abseta = getattr(
                row,
                getVar(arg,'Eta')
            ) 
            pt     = getattr(row, getVar(arg,'Pt'))
            if pt<30: pt =30 #only fakerate checks allow pt < 30. This is an approximation to not re-run the Tag and Probe
            if mu_id=="not_muon_or_not_isocorr":
                fcn_ret = fcn(pt,abseta)
            else:
                fcn_ret=fcn(mu_id,pt,abseta)
            if indexed:
                value, err = fcn_ret
                if shift == 1:
                    ret *= (value + err)
                elif shift == -1:
                    ret *= (value - err)
                else:
                    ret *= value
            else:
                ret   *= fcn_ret
        return ret
    return multiple


##put here the trigger correction as in https://github.com/mverzett/UWHiggs/blob/WH_At_Paper/wh/mcCorrectors.py
"""
correct_e             = make_multiple(HetauCorrection.correct_hamburg_e    )
correct_eid13_mva     = make_multiple(HetauCorrection.correct_eid13_mva    )
correct_eid13_p1s_mva = make_multiple(HetauCorrection.correct_eid13_p1s_mva)
correct_eid13_m1s_mva = make_multiple(HetauCorrection.correct_eid13_m1s_mva)


correct_muid15_tight_mva     = make_multiple(HetauCorrection.correct_muid15_tight_mva)
correct_muid15_tight_p1s_mva = make_multiple(HetauCorrection.correct_muid15_tight_p1s_mva)
correct_muid15_tight_m1s_mva = make_multiple(HetauCorrection.correct_muid15_tight_m1s_mva)


#iso tight  and id tight
correct_muiso15_mva     = make_multiple(HetauCorrection.correct_muiso15_mva)
correct_muiso15_p1s_mva = make_multiple(HetauCorrection.correct_muiso15_p1s_mva)
correct_muiso15_m1s_mva = make_multiple(HetauCorrection.correct_muiso15_m1s_mva)

correct_eid15_mva     = make_multiple(HetauCorrection.correct_eid15_mva    )
correct_eid15_p1s_mva = make_multiple(HetauCorrection.correct_eid15_p1s_mva)
correct_eid15_m1s_mva = make_multiple(HetauCorrection.correct_eid15_m1s_mva)


correct_mutrig15_mva_p3     = make_multiple(HetauCorrection.correct_mutrig15_mva_p3    )
correct_mutrig15_p1s_mva_p3 = make_multiple(HetauCorrection.correct_mutrig15_p1s_mva_p3)
correct_mutrig15_m1s_mva_p3 = make_multiple(HetauCorrection.correct_mutrig15_m1s_mva_p3)



correct_mutrig15_mva_p2     = make_multiple(HetauCorrection.correct_mutrig15_mva_p2    )
correct_mutrig15_p1s_mva_p2 = make_multiple(HetauCorrection.correct_mutrig15_p1s_mva_p2)
correct_mutrig15_m1s_mva_p2 = make_multiple(HetauCorrection.correct_mutrig15_m1s_mva_p2)


correct_mutrig15_mva_combined     = make_multiple(HetauCorrection.correct_mutrig15_mva_combined)
correct_mutrig15_p1s_mva_combined = make_multiple(HetauCorrection.correct_mutrig15_p1s_mva_combined)
correct_mutrig15_m1s_mva_combined = make_multiple(HetauCorrection.correct_mutrig15_m1s_mva_combined)


correct_eiso13_mva     = make_multiple(HetauCorrection.correct_eiso13_mva    )
correct_eiso13_p1s_mva = make_multiple(HetauCorrection.correct_eiso13_p1s_mva)
correct_eiso13_m1s_mva = make_multiple(HetauCorrection.correct_eiso13_m1s_mva)


#correct_eid_mva = make_multiple(HetauCorrection.scale_eleId_hww)
#correct_eReco_mva = make_multiple(HetauCorrection.scale_elereco_hww)
#correct_eIso_mva = make_multiple(HetauCorrection.scale_eleIso_hww)
correct_trigger_mva    = make_multiple(HetauCorrection.single_ele_mva, indexed=True)
correct_trigger_mva_up = make_multiple(HetauCorrection.single_ele_mva, indexed=True, shift=1)
correct_trigger_mva_dw = make_multiple(HetauCorrection.single_ele_mva, indexed=True, shift=-1)

efficiency_trigger_mva    = make_multiple(HetauCorrection.single_ele_eff_mva, indexed=True)
efficiency_trigger_mva_up = make_multiple(HetauCorrection.single_ele_eff_mva, indexed=True, shift=1)
efficiency_trigger_mva_dw = make_multiple(HetauCorrection.single_ele_eff_mva, indexed=True, shift=-1)

correct_eEmb     = make_multiple(HetauCorrection.correct_eEmb,indexed=True)
correct_eEmb_p1s = make_multiple(HetauCorrection.correct_eEmb,indexed=True, shift=1)
correct_eEmb_m1s = make_multiple(HetauCorrection.correct_eEmb,indexed=True, shift=-1)
"""


correct_muid15_tight_mva= make_multiple(MuonPOGCorrections.make_muon_pog_PFTight_2015CD())
correct_muiso15_mva = make_multiple(MuonPOGCorrections.make_muon_pog_TightIso_2015CD(),"Tight")
correct_mutrig15_mva_combined= make_multiple(MuonPOGCorrections.make_muon_pog_IsoMu20oIsoTkMu20_2015())

#correction for muon POG https://twiki.cern.ch/twiki/bin/viewauth/CMS/MuonReferenceEffsRun2
muid15_tight_correction = make_shifted_weights(
    correct_muid15_tight_mva,
    [],
    []
)
#iso tight id tight
muiso15_correction = make_shifted_weights(
    correct_muiso15_mva,
    [],
    []
)
#single muon trigger combined c and d with the lumi weight 
mutrig15_mva_combined = make_shifted_weights(
    correct_mutrig15_mva_combined,
    [],
    []
)

"""
#for dataC
mutrig15_mva_p2_correction=make_shifted_weights(
    correct_mutrig15_mva_p2,
    ['mutrigp1s_2','mutrigm1s_2'],
    [correct_mutrig15_p1s_mva_p2,correct_mutrig15_m1s_mva_p2]
)
#for dataD
mutrig15_mva_p3_correction = make_shifted_weights(
    correct_mutrig15_mva_p3,
    ['mutrigp1s_3','mutrigm1s_3'],
    [correct_mutrig15_p1s_mva_p3,correct_mutrig15_m1s_mva_p3]
)

eiso_correction = make_shifted_weights(
    correct_eiso13_mva, 
    ['eisop1s','eisom1s'], 
    [correct_eiso13_p1s_mva, correct_eiso13_m1s_mva],
)

eid_correction = make_shifted_weights(
    correct_eid13_mva,
    ['eidp1s','eidm1s'],
    [correct_eid13_p1s_mva, correct_eid13_m1s_mva]
)
#it is id and iso from https://indico.cern.ch/event/370511/contribution/3/attachments/1168717/1687113/tnP_EGM_Oct_12.pdf
eid15_correction = make_shifted_weights(
    correct_eid15_mva,
    ['eidp1s','eidm1s'],
    [correct_eid15_p1s_mva, correct_eid15_m1s_mva]
)

eEmb_correction = make_shifted_weights(
    correct_eEmb,
    ['eEmbp1s','eEmbm1s'],
    [correct_eEmb_p1s, correct_eEmb_m1s]
)


trig_correction = make_shifted_weights(
    correct_trigger_mva,
    ['trp1s', 'trm1s'],
    [correct_trigger_mva_up, correct_trigger_mva_dw]
)

trig_efficiency = make_shifted_weights(
    efficiency_trigger_mva,
    ['trp1s', 'trm1s'],
    [efficiency_trigger_mva_up, efficiency_trigger_mva_dw]
)
"""
