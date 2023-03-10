#!/usr/bin/env ccp4-python

import os, sys, shutil
import argparse

def make_best_fasta(scratch_directory=None, input_fasta=None, output_fasta="output.fasta", seqsdb=None, nseqs=None, maxhits=10, dbtype=None, nproc=1, debug=False):

    phmmerEXE=os.path.join(os.environ["CCP4"], "libexec", "phmmer")
    workingDIR=os.getcwd()

    scratch=os.path.abspath(scratch_directory)
    if not os.path.isdir(scratch):
        os.mkdir(scratch)
    os.chdir(scratch)
    
    splits_list=[]
    if nseqs is None:
        splits=0
        for x in os.listdir(seqsdb):
            if x.endswith(".fasta"):
                splits+=1
                splits_list.append(int(x.split("_")[-1].replace(".fasta", "")))
    else:
        splits=int(nseqs)
        for x in range(1, nseqs+1):
            splits_list.append(x)
    
    ssub="#!/bin/bash\n"
    ssub+="#SBATCH --job-name=phmmer_%a     # Job name\n"
    ssub+="#SBATCH --array=1-%d%%%d             # set array\n" % (max(splits_list), nproc)
    ssub+="#SBATCH --ntasks=1                                   # Run on a single CPU\n"
    ssub+="#SBATCH --mem=1gb                                    # Job memory request\n"
    ssub+="#SBATCH --time=00:60:00                              # Time limit hrs:min:sec\n"
    ssub+="#SBATCH --output=job_%s_%%a.log   # Standard output and error log\n" % dbtype 
    
    ssub+="\n"
    ssub+="source %s\n" % os.path.join(os.environ["CCP4"], "bin", "ccp4.setup-sh")
    ssub+="\n"

    #ssub+="\n"
    #ssub+="sleep $((RANDOM%30+1))\n"
    #ssub+="\n"
    
    #ssub+="time cp -v %s/seq_${SLURM_ARRAY_TASK_ID}.fasta /dev/shm/seq_${SLURM_ARRAY_TASK_ID}.fasta\n" % (seqsdb)  
    #ssub+="%s -o phmmer_${SLURM_ARRAY_TASK_ID}.log %s /dev/shm/seq_${SLURM_ARRAY_TASK_ID}.fasta\n" % (phmmerEXE, input_fasta)  
    #ssub+="%s -o phmmer_${SLURM_ARRAY_TASK_ID}.log %s %s/seq_${SLURM_ARRAY_TASK_ID}.fasta\n" % (phmmerEXE, input_fasta, seqsdb)  
    ssub+="if [[ -f %s/seq_${SLURM_ARRAY_TASK_ID}.fasta ]]\n" % seqsdb
    ssub+="then\n"
    ssub+="%s --cpu 1 --F1 1e-15 --F2 1e-15 -o phmmer_%s_${SLURM_ARRAY_TASK_ID}.log %s %s/seq_${SLURM_ARRAY_TASK_ID}.fasta\n" % (phmmerEXE, dbtype, input_fasta, seqsdb)  
    ssub+="time ccp4-python -m mrparse.phmmer_cl -g -l phmmer_%s_${SLURM_ARRAY_TASK_ID}.log -m %d -b temp_%s_${SLURM_ARRAY_TASK_ID}.fasta -s %s\n" % (dbtype, maxhits, dbtype, seqsdb)
    #ssub+="time ccp4-python -m mrparse.phmmer_cl -g -l phmmer_${SLURM_ARRAY_TASK_ID}.log -m %d -b temp_${SLURM_ARRAY_TASK_ID}.fasta -s /dev/shm\n" % (maxhits)
    #ssub+="time rm -v /dev/shm/seq_${SLURM_ARRAY_TASK_ID}.fasta\n"
    ssub+="touch FINISHED_%s_${SLURM_ARRAY_TASK_ID}.txt\n" % dbtype
    ssub+="fi\n"
    
    subscript="afdb.sub" 
    
    subfile=open(subscript, "w")
    subfile.write(ssub)
    subfile.close()
    
    os.system("sbatch %s" % subscript)
    
    bestfasta=""
    pcount=0
    foundlist=[]
    while pcount < splits:
        for x in splits_list: #range(1, (splits+1)):
            if os.path.isfile("FINISHED_%s_%d.txt" % (dbtype, x)) and x not in foundlist:
                if os.path.isfile("temp_%s_%d.fasta" % (dbtype, x)):
                    myseqs=open("temp_%s_%d.fasta" % (dbtype,  x))
                    bestfasta+="".join(myseqs.readlines())
                    myseqs.close()
                    pcount+=1
                    foundlist.append(x)
                    if not debug:
                        os.remove("temp_%s_%d.fasta" % (dbtype, x))
                else:
                    print("Warning - sequence list file not found for %d" % x)
                if not debug:
                    os.remove("FINISHED_%s_%d.txt" % (dbtype, x))
                    if os.path.isfile("job_%d.log" % x):
                        os.remove("job_%d.log" % x)
                    if os.path.isfile("phmmer_%d.log" % x):
                        os.remove("phmmer_%d.log" % x)

    if not debug:
        if os.path.isfile("afdb.sub"):
            os.remove("afdb.sub")

    out=open(output_fasta,"w")
    out.write(bestfasta)
    out.close()

    os.chdir(workingDIR)
    import time
    time.sleep(1)
    if not debug:
        shutil.rmtree(scratch_directory, ignore_errors=True)

def getmytophits(logfile=None, maxhits=10, seqsdb=None, output=None, debug=False):

    myfastasequences=""
    hitscount=0
    
    if os.path.isfile(logfile):
        plog=open(logfile)
        plines=plog.readlines()
        plog.close()
    
        myid=int(logfile.split("_")[2].replace(".log",""))
    
        myhitcount=0
        myhitlist=[]
    
        for line in plines:
            #if ">> AFDB" in line:
            if ">> " in line:
                if ">> AFDB" in line: 
                    name=line.split()[1].replace("AFDB:","")
                elif ">> MGY" in line:
                    name=line.split()[1]
                myhitcount+=1
                myhitlist.append(name)
                if myhitcount == maxhits:
                    myfasta=open(os.path.join(seqsdb, "seq_%d.fasta" % myid), "r")
                    myfline=myfasta.readline()
                    mysequences=""
                    while myfline:
                        if ">" in myfline:
                            if myfline.split()[0].replace(">AFDB:","") in myhitlist or myfline.split()[0].replace(">","") in myhitlist:
                                myfastasequences+=myfline + myfasta.readline()
                                hitscount+=1
                                if hitscount >= maxhits:
                                    break
                        myfline=myfasta.readline()
                    myfasta.close()
                    break

        if myhitcount < maxhits:
            myfasta=open(os.path.join(seqsdb, "seq_%d.fasta" % myid), "r")
            myfline=myfasta.readline()
            mysequences=""
            while myfline:
                if ">" in myfline:
                    if myfline.split()[0].replace(">AFDB:","") in myhitlist or myfline.split()[0].replace(">","") in myhitlist:
                        myfastasequences+=myfline + myfasta.readline()
                        hitscount+=1
                        if hitscount >= myhitcount:
                            break
                myfline=myfasta.readline()
            myfasta.close()
    
    else:
       print("logfile %s not found" % logfile)
       sys.exit()
    
    outfasta=open(output,"w")
    outfasta.write(myfastasequences)
    outfasta.close()

if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument('-f', '--seqin', help="Input sequence file")
    parser.add_argument('-t', '--scratch', help="scratch directory (default=scratch)", default="scratch")
    parser.add_argument('-n', '--nseqs', help="number of sequence files to use (defaults to number in seqsdb folder)", default=None)
    parser.add_argument('-o', '--output', help="output best fasta file (default output.fasta)", default="output.fasta")
    parser.add_argument('-s', '--seqsdb', help="location of fasta files", default=None)
    parser.add_argument('-y', '--dbtype', help="database type (afdb or esmatlas)", default=None)
    
    parser.add_argument('-l', '--logfile', help="phmmer log file", default=None)
    parser.add_argument('-b', '--bestfasta', help="best hits fasta file for this job", default=None)
    parser.add_argument('-g', '--gethits', action='store_true', help="get best hits for this job")
    
    parser.add_argument('-m', '--maxhits', help="maximum number of hits (default=10)", type=int, default=10)
    parser.add_argument('-d', '--debug', action='store_true', help="debug output")
    args = parser.parse_args()
    
    seqsdb=args.seqsdb
    workingDIR=os.getcwd()

    
    if not args.gethits:
        seqin=os.path.abspath(args.seqin)
        scratch=os.path.abspath(args.scratch)
        outfasta=os.path.abspath(args.output)
    
        if os.path.isdir(scratch) == False:
            os.mkdir(scratch)
    
        if args.nseqs is None:
            nseqs=0
            for x in os.listdir(seqsdb):
                if x.endswith(".fasta"):
                    nseqs+=1
        else:
            nseqs=int(args.nseqs)
    
        make_best_fasta(scratch_directory=scratch, input_fasta=seqin, output_fasta=outfasta, seqsdb=seqsdb, nseqs=nseqs, maxhits=args.maxhits, dbtype=args.dbtype, debug=args.debug)
    else:
        getmytophits(logfile=args.logfile, maxhits=args.maxhits, seqsdb=seqsdb, output=args.bestfasta, debug=args.debug)
