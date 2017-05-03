listi=["",'1','2','3','4']
for i in listi:
#    target="DY"+i+"JetsToLNu_TuneCUETP8M1_13TeV-madgraphMLM-pythia8.root.txt"
    target="DY"+i+"JetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8.root.txt"
    fyal=open(target,"r")
    jambu=0
    jambi=0
    for line in fyal.readlines():
        count=line.split("  ")[0]
        lost=line.split("  ")[1]
        jambu+=int(count)
        jambi+=int(lost)

    print target,float(jambi)/float(jambu)


