## redone on 12/13/13 to run mask50 for the prk non-penal methods ONLY.

import re, os, sys, subprocess, fnmatch
import shutil

def GetFiles(ext, dirpath):
        files=[]
        dir=os.listdir(dirpath)
        for file in dir:
                if fnmatch.fnmatch(file, str(ext)+'.+'):
                        files.append(file)
        return files
########

n=int(sys.argv[1]) - 1

alndir='cp /home/sjs3495/alntree/nucguided_linsi_'+gene+'_seqs_real/'
treedir='cp /home/sjs3495/alntree/aatrees_linsi_'+gene+'_seqs_real/aatree'

copy=treedir+str(n)+'.txt tree.tre'
subprocess.call(copy, shell=True)

name=alg+mask+str(n)+'.fasta'
copy = alndir+name+' temp.fasta'
subprocess.call(copy, shell=True)

cline='/home/sjs3495/bin/codeml codeml.ctl'
runit=subprocess.call(cline, shell=True)
final_file='rst'
shutil.move(final_file, outdir+name+'.rst')
