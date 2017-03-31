

def GetSF(WP,x,flavour,syst):
    if (WP==1):
        if (abs(flavour)==4 or abs(flavour)==5):
            if (syst==0):
                return 0.498094*((1.+(0.422991*x))/(1.+(0.210944*x)))
      
        else: 
            if (syst==0):
                return 1.0589+0.000382569*x+-2.4252e-07*x*x+2.20966e-10*x*x*x;
   
    else: return 0;


def  bTagEventWeight( nBtaggedJets, bjetpt_1,  bjetflavour_1,  bjetpt_2,  bjetflavour_2,  WP,  syst, nBTags):
    if (nBtaggedJets > 2): return -10000
    if ( nBTags > 2 ): return -10000

 
    ##################################################################
 #   Event weight matrix:
 #   ------------------------------------------------------------------
 #   nBTags\b-tagged jets  |    0        1             2
 #   ------------------------------------------------------------------
 #     0                   |    1      1-SF      (1-SF1)(1-SF2)
 #                         |
 #     1                   |    0       SF    SF1(1-SF2)+(1-SF1)SF2
 #                         |
 #     2                   |    0        0           SF1SF2
    ##################################################################
 
  
    if ( nBTags > nBtaggedJets): return 0
    if ( nBTags==0 and nBtaggedJets==0): return 1

    weight = 0
    if (nBtaggedJets==1):
        SF = GetSF(WP,bjetpt_1,bjetflavour_1,syst)
        for  i in range(2):
            if ( i != nBTags ): continue
            weight += pow(SF,i)*pow(1-SF,1-i)
    
  
    elif (nBtaggedJets==2):
        SF1 = GetSF(WP,bjetpt_1,bjetflavour_1,syst)
        SF2 = GetSF(WP,bjetpt_2,bjetflavour_2,syst)
    
        for i in range(2):
            for j in range(2):
                if( (i+j) != nBTags ): continue
                weight += pow(SF1,i)*pow(1-SF1,1-i)*pow(SF2,j)*pow(1-SF2,1-j)

    return weight


