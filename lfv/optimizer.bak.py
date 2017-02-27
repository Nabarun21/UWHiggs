2 #author Mauro Verzetti
'small interface module to deal with optimizization'

import os
import itertools

os.environ['RUN_OPTIMIZATION']='1'

RUN_OPTIMIZATION = ('RUN_OPTIMIZATION' in os.environ) and eval(os.environ['RUN_OPTIMIZATION'])
print RUN_OPTIMIZATION
#RUN_OPTIMIZATION = True
grid_search = {}

_0jets = {
#   'mPt'  : range(25,60,5)+[43],
#   'ePt'  : range(10,35,10),
#   'dphi' : [1.9,2.1,2.3, 2.5,2.7],
#   'mMtToPfMet' : range(10,100,20),
#   'eMtToPfMet' : range(45,95,20),

   'dphiMetToE':[0.3,0.5,0.7,0.9,1.2,1.5,1.8,2.1,2.4]
   
}

_0jets_default = {
   'mPt' : 30,
   'ePt' : 50,
   'dphi': 2.5,
   'mMtToPfMet' :50,
   
}
_0jet_region_templates = ['mPt%i', 'ePt%i', 'dphi%.2f', 'mMtToPfMet%i','eMtToPfMet%i','dphiMetToE%.2f'] #'mPt%i_ePt%i_dphi%.2f_mMtToPfMet%i'

def _get_0jet_regions(mPt, ePt, dphi, mMtToPfMet,eMtToPfMet,dphiMetToE):
#   pass_mPt        = [i for i in _0jets['mPt'       ] if mPt        > i] 
#   pass_ePt        = [i for i in _0jets['ePt'       ] if ePt        > i] 
#   pass_dphi       = [i for i in _0jets['dphi'      ] if dphi       > i] 
#   pass_mMtToPfMet = [i for i in _0jets['mMtToPfMet'] if mMtToPfMet > i]
#   pass_eMtToPfMet = [i for i in _0jets['eMtToPfMet'] if eMtToPfMet < i]
   pass_dphiMetToE = [i for i in _0jets['dphiMetToE'] if dphiMetToE < i]

#   cuts = [pass_mPt, pass_ePt,pass_dphi, pass_mMtToPfMet,pass_eMtToPfMet,pass_dphiMetToE]
#   pass_default_mPt        = mPt        > _0jets_default['mPt'       ] 
#   pass_default_ePt        = ePt        > _0jets_default['ePt'       ] 
#   pass_default_dphi       = dphi       > _0jets_default['dphi'      ] 
#   pass_default_mMtToPfMet = mMtToPfMet < _0jets_default['mMtToPfMet']
#
#   defaults = [pass_default_mPt, pass_default_ePt, pass_default_dphi, pass_default_mMtToPfMet]
   ret = []
   for cut_idx, opts in enumerate(cuts):
#       if all(j for i,j in enumerate(defaults) if i != cut_idx):
      ret.extend([_0jet_region_templates[cut_idx] % i for i in opts])

   return ret
_1jets = {
#   'mPt'  : range(25,60,5)+[43],
#   'ePt'  : range(10,35,10),
#   'dphi' : [1.7,1.9,2.1,2.3,2.5,2.7,2.8],
#   'mMtToPfMet' : range(10,100,20),
#   'eMtToPfMet' : range(45,95,20),
   'dphiMetToE':[0.3,0.4,0.5,0.7,0.9,1.2,1.5,1.8,2.1,2.4]
}


_1jets_default = {
    'mPt' : 30,
    'ePt' : 40,
    'mMtToPfMet' :35,
}

_1jet_region_templates =  ['mPt%i', 'ePt%i', 'dphi%.2f', 'mMtToPfMet%i','eMtToPfMet%i','dphiMetToE%.2f'] 
def _get_1jet_regions(mPt, ePt, dphi, mMtToPfMet,eMtToPfMet,dphiMetToE):
   pass_mPt        = [i for i in _1jets['mPt'       ] if mPt        > i] 
   pass_ePt        = [i for i in _1jets['ePt'       ] if ePt        > i] 
   pass_dphi       = [i for i in _1jets['dphi'      ] if dphi       > i] 
   pass_mMtToPfMet = [i for i in _1jets['mMtToPfMet'] if mMtToPfMet > i]
   pass_eMtToPfMet = [i for i in _1jets['eMtToPfMet'] if eMtToPfMet < i]
   pass_dphiMetToE = [i for i in _1jets['dphiMetToE'] if dphiMetToE < i]

   cuts = [pass_mPt, pass_ePt, pass_dphi, pass_mMtToPfMet,pass_eMtToPfMet,pass_dphiMetToE]

#   pass_default_mPt        = mPt        > _1jets_default['mPt'       ] 
#   pass_default_ePt        = ePt        > _1jets_default['ePt'       ] 
#   pass_default_mMtToPfMet = mMtToPfMet < _1jets_default['mMtToPfMet']
#
 #  defaults = [pass_default_mPt, pass_default_ePt,  pass_default_mMtToPfMet]
   ret = []
   for cut_idx, opts in enumerate(cuts):
#       if all(j for i,j in enumerate(defaults) if i != cut_idx):
      ret.extend([_1jet_region_templates[cut_idx] % i for i in opts])
           
   return ret

    
    
_2jets = {
#    'mPt'  : [25,30,40],
#    'ePt'  : [10,20],
#    'mMtToPfMet' : [15,30,45],
#    'eMtToPfMet' : [25,40,60],
    'dphiMetToE':[0.3,0.4,0.5,0.7,0.9,1.2,1.5,1.8,2.1,2.4]
#    'vbf_mass' : range(400, 550, 50),
#    'vbf_deta' : [1, 1.5, 2],
}
_2jets_default = {
    'mPt' : 30,
    'ePt' : 40,
    'mMtToPfMet' : 35,
    'vbf_mass' : 400,
    'vbf_deta' : 2.5,
}

_2jet_region_templates =  ['mPt%i', 'ePt%i','mMtToPfMet%i','eMtToPfMet%i','dphiMetToE%.2f', 'vbf_mass%i', 'vbf_deta%.1f' ] 

def _get_2jet_regions(mPt, ePt, mMtToPfMet,eMtToPfMet,dphiMetToE,vbf_mass, vbf_deta):
   pass_mPt        = [i for i in _0jets['mPt'       ] if mPt        > i] 
   pass_ePt        = [i for i in _0jets['ePt'       ] if ePt        > i] 
   pass_mMtToPfMet = [i for i in _0jets['mMtToPfMet'] if mMtToPfMet > i]
   pass_eMtToPfMet = [i for i in _0jets['eMtToPfMet'] if eMtToPfMet < i]
   pass_dphiMetToE = [i for i in _0jets['dphiMetToE'] if dphiMetToE < i]
   pass_vbf_mass = [i for i in _2jets['vbf_mass'] if vbf_mass > i]
   pass_vbf_deta = [i for i in _2jets['vbf_deta'] if vbf_deta > i]

   cuts = [pass_mPt, pass_ePt, pass_mMtToPfMet,pass_eMtToPfMet,pass_dphiMetToE,pass_vbf_mass, pass_vbf_deta]
#    pass_default_mPt        = mPt        > _2jets_default['mPt'       ] 
#    pass_default_ePt        = ePt        > _2jets_default['ePt'       ] 
#    pass_default_mMtToPfMet = mMtToPfMet < _2jets_default['mMtToPfMet']
#    pass_default_vbf_mass = vbf_mass >  _2jets_default['vbf_mass'] 
#    pass_default_vbf_deta = vbf_deta > _2jets_default['vbf_deta']
#
  #  defaults = [pass_default_mPt, pass_default_ePt,  pass_default_mMtToPfMet, pass_default_vbf_mass, pass_default_vbf_deta]
   ret = []
   for cut_idx, opts in enumerate(cuts):
#        if all(j for i,j in enumerate(defaults) if i != cut_idx):
      ret.extend([_2jet_region_templates[cut_idx] % i for i in opts])
            
   return ret
  


def empty(*args):
    return []
##
compute_regions_0jet = _get_0jet_regions if RUN_OPTIMIZATION else empty
compute_regions_1jet = _get_1jet_regions if RUN_OPTIMIZATION else empty
compute_regions_2jet = _get_2jet_regions if RUN_OPTIMIZATION else empty


#print compute_regions_0jet(100000,100000,100000,1000,-1000000,-100000)
#print
#print compute_regions_1jet(100000,100000,100000,1000,-1000000,-100000)
#print
#rint compute_regions_2jet(100000,100000,100000,-1000,-1000000,100000,1000000)
#



if __name__ == "__main__":
    from pdb import set_trace
    set_trace()
    #print '\n'.join(grid_search.keys())
else:
    print "Running optimization: %s" % RUN_OPTIMIZATION
