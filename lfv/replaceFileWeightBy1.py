import os 
import sys
for filename in os.listdir("inputs/"+sys.argv[1]+"/"):
    if "lumicalc.sum" in filename and "LFV" not in filename and "data" not in filename:
        with open("inputs/"+sys.argv[1]+"/"+filename,"wb") as lumifile:
#            lumistr=float(lumifile.readline().strip())
            lumifile.write(sys.argv[2])
            lumifile.close()
        with open("inputs/"+sys.argv[1]+"/"+filename,"r") as lumifile:
            lumistr2=float(lumifile.readline().strip())
        print filename,":  ",1/lumistr2,"   ",lumistr2
