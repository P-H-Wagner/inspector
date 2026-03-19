import os
from glob import glob
from datetime import datetime
import sys 
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('channel') # sig or hb or data
parser.add_argument('nFiles')  # nr of miniAOD files to process
parser.add_argument('inspector')  # which inspector to run (reco, gen) 
args = parser.parse_args()

######################################

#500 jobs per user on std queue
nMaxJobs = 800

#default
filesPerJob = 1
queue = 'short'
time =45
print("========> processing ", filesPerJob, " files per job")

######################################

#queue = 'standard'
#time = 60
nevents = 500
nFiles = int(args.nFiles)
now = datetime.now()
dt_string = now.strftime("%d_%m_%Y_%H_%M_%S")

######################################

def filesFromFolder(direc):
  filenames = os.listdir(direc)
  return ['file:' + direc + filename for filename in filenames ]

def filesFromTxt(txtFile):
  with open(txtFile) as dataFiles: 
    filenames = [line[0:-1] for line in dataFiles] #-2 to remove the newline character \n
  return ['file:root://cms-xrd-global.cern.ch//' + filename for filename in filenames ]

#######################################

# Input source

if args.channel == 'sig':

    #directory = '/pnfs/psi.ch/cms/trivcat/store/user/manzoni/all_signals_HbToDsPhiKKPiMuNu_MT_MINI_21jan23_v1/' #old signal from MA thesis
    #inputfiles = filesFromFolder(directory)
    directory = '/pnfs/psi.ch/cms/trivcat/store/user/pahwagne/mc/signals/all_signals_request_21_11_23.txt' # new signals!!
    inputfiles = filesFromTxt(directory)
    naming = 'all_signals'

if args.channel == 'hb':

  if args.inspector == "gen":
    #directory = '/pnfs/psi.ch/cms/trivcat/store/user/pahwagne/miniAOD/hb_fragment_11_06_2024_18_06_45/'
    #directory = '/pnfs/psi.ch/cms/trivcat/store/user/pahwagne/miniAOD/hb_fragment_11_06_2024_18_56_33/' #higher stats
    directory = '/pnfs/psi.ch/cms/trivcat/store/user/pahwagne/mc/hb/inclusive/private_prod/hb_inclusive_crab_20250605_133227.txt' #private production
    #inputfiles = filesFromFolder(directory)
    inputfiles = filesFromTxt(directory)


  if args.inspector == "reco":
    directory = '/pnfs/psi.ch/cms/trivcat/store/user/manzoni/inclusive_HbToDsPhiKKPiMuNu_MINI_25mar21_v1/' #hb 
    inputfiles = filesFromFolder(directory)
  naming = 'hb_inclusive'

if args.channel == 'b0':

  if args.inspector == "gen":
    directory = '/pnfs/psi.ch/cms/trivcat/store/user/pahwagne/miniAOD/b0_fragment_03_06_2024_22_52_56/' 
    inputfiles = filesFromFolder(directory)

  if args.inspector == "reco":
    directory = '/pnfs/psi.ch/cms/trivcat/store/user/pahwagne/mc/hb/b0/b0.txt' #b0 
    inputfiles = filesFromTxt(directory)

  naming = 'b0'

if args.channel == 'bplus':

  if args.inspector == "gen":
    directory = '/pnfs/psi.ch/cms/trivcat/store/user/pahwagne/miniAOD/bplus_fragment_03_06_2024_19_15_19/' 
    inputfiles = filesFromFolder(directory)

  if args.inspector == "reco":
    directory = '/pnfs/psi.ch/cms/trivcat/store/user/pahwagne/mc/hb/bplus/bplus.txt' #bplus 
    inputfiles = filesFromTxt(directory)
  naming = 'bplus'

if args.channel == 'bs':

  if args.inspector == "gen":
    directory = '/pnfs/psi.ch/cms/trivcat/store/user/pahwagne/miniAOD/bs_fragment_03_06_2024_19_15_47/'  # old filter
    directory = '/pnfs/psi.ch/cms/trivcat/store/user/pahwagne/miniAOD/bs_fragment_26_08_2024_16_13_34/'  # old filte but different channel
    directory = '/pnfs/psi.ch/cms/trivcat/store/user/pahwagne/miniAOD/bs_fragment_26_08_2024_19_27_38/'  # old filte but different channel and more stats!
    #directory = '/pnfs/psi.ch/cms/trivcat/store/user/pahwagne/miniAOD/bs_fragment_06_06_2024_09_52_40'  # new filter
    #directory = '/pnfs/psi.ch/cms/trivcat/store/user/pahwagne/miniAOD/bs_fragment_06_06_2024_20_38_57'  # new filter
    inputfiles = filesFromFolder(directory)

  if args.inspector == "reco":
    directory = '/pnfs/psi.ch/cms/trivcat/store/user/pahwagne/mc/hb/bs/bs.txt' #bs 
    inputfiles = filesFromTxt(directory)
  naming = 'bs'

if args.channel == 'lambdab':

  if args.inspector == "gen":
    directory = '/pnfs/psi.ch/cms/trivcat/store/user/pahwagne/miniAOD/lambdab_fragment_04_06_2024_09_40_54' 
    inputfiles = filesFromFolder(directory)

  if args.inspector == "reco":
    directory = '/pnfs/psi.ch/cms/trivcat/store/user/pahwagne/mc/hb/lambdab/lambdab.txt' #lambdab 
    inputfiles = filesFromTxt(directory)
  naming = 'lambdab'

if args.channel == 'data':
  directory = '/pnfs/psi.ch/cms/trivcat/store/user/pahwagne/data/BPark_2018_D/BPark_2018D_part1.txt' #data bParking 2018 part D
  #txtFile = '/pnfs/psi.ch/cms/trivcat/store/user/pahwagne/data/dataTest/single.txt' # test
  inputfiles = filesFromTxt(directory)
  naming = 'data'

if args.channel == "dstau":
  directory = '/pnfs/psi.ch/cms/trivcat/store/user/pahwagne/miniAOD/dstau_fragment_12_11_2024_14_56_51/'  
  directory = '/pnfs/psi.ch/cms/trivcat/store/user/pahwagne/miniAOD/dstau_fragment_14_11_2024_13_17_18/' #more stats!  
  inputfiles = filesFromFolder(directory)
  naming = 'dstau'

if args.channel == "dsstartau":
  directory = '/pnfs/psi.ch/cms/trivcat/store/user/pahwagne/miniAOD/dsstartau_fragment_12_11_2024_14_57_15/' 
  directory = '/pnfs/psi.ch/cms/trivcat/store/user/pahwagne/miniAOD/dsstartau_fragment_14_11_2024_15_26_04/' #more stats! 
  inputfiles = filesFromFolder(directory)
  naming = 'dsstartau'

if args.channel == "dsmu":
 
  directory = '/pnfs/psi.ch/cms/trivcat/store/user/pahwagne/miniAOD/dsmu_fragment_12_11_2024_15_21_32/' 
  directory = '/pnfs/psi.ch/cms/trivcat/store/user/pahwagne/miniAOD/dsmu_fragment_13_11_2024_18_58_31/' #more stats! 
  inputfiles = filesFromFolder(directory)
  naming = 'dsmu'

if args.channel == "dsmu_isgw2":
 
  directory = '/pnfs/psi.ch/cms/trivcat/store/user/pahwagne/miniAOD/dsmu_isgw2_fragment_12_11_2024_22_14_08/' 
  directory = '/pnfs/psi.ch/cms/trivcat/store/user/pahwagne/miniAOD/dsmu_isgw2_fragment_16_11_2024_12_18_09/' #more stats! 
  inputfiles = filesFromFolder(directory)
  naming = 'dsmu_isgw2'

if args.channel == "dsstarmu":
  directory = '/pnfs/psi.ch/cms/trivcat/store/user/pahwagne/miniAOD/dsstarmu_fragment_12_11_2024_15_21_41/'  
  directory = '/pnfs/psi.ch/cms/trivcat/store/user/pahwagne/miniAOD/dsstarmu_fragment_14_11_2024_09_38_51/' #more stats!  
  inputfiles = filesFromFolder(directory)
  naming = 'dsstarmu'

if args.channel == "dsstarmu_isgw2":
  directory = '/pnfs/psi.ch/cms/trivcat/store/user/pahwagne/miniAOD/dsstarmu_isgw2_fragment_18_11_2024_15_42_17/' #more stats!  
  inputfiles = filesFromFolder(directory)
  naming = 'dsstarmu_isgw2'



if nFiles != -1:
  #process not the whole dataset but only nFiles
  inputfiles = inputfiles[0:nFiles] #50 files give ca 200k events

os.makedirs("/pnfs/psi.ch/cms/trivcat/store/user/pahwagne/nanoAOD/"+dt_string)
os.makedirs(dt_string+"/logs")
os.makedirs(dt_string+"/errs")

for i,j in enumerate(range(0, len(inputfiles), filesPerJob)):

  fin = inputfiles[j:j+filesPerJob]
  #print(fin)
  #template
  temp = open("temp_cfg.py", "rt")
  #file to write to
  cfg = open(dt_string+"/cfg_chunk_{0}.py".format(i),"wt")
  #file to save things
  fout = "/scratch/pahwagne/{0}/{1}_chunk_{2}.root".format(dt_string,naming,i)  

  for line in temp:
    if   "HOOK_CHANNEL" in line: cfg.write(line.replace("HOOK_CHANNEL", args.channel))
    elif "HOOK_INPUT" in line: cfg.write(line.replace("HOOK_INPUT", directory))
    elif "HOOK_N_EVENTS" in line: cfg.write(line.replace("HOOK_N_EVENTS", str(nevents)))
    elif "HOOK_FILE_IN" in line: cfg.write(line.replace("HOOK_FILE_IN", str(fin)))
    elif "HOOK_FILE_OUT" in line: cfg.write(line.replace("HOOK_FILE_OUT", fout))
    elif "HOOK_INSP" in line: cfg.write(line.replace("HOOK_INSP", args.inspector))
    else: cfg.write(line)

  temp.close()
  cfg.close()

  to_write = '\n'.join([
         '#!/bin/bash',
         # --- create scratch dir and create temp .sh file in scratch ---
         'mkdir -p /scratch/pahwagne/'+dt_string,
         'payload=/scratch/pahwagne/'+dt_string+'/apptainer-payload-{0}.sh'.format(i),
         # --- write into temp ---
         'cat > "$payload" << EOF',
         'cd /work/pahwagne/inspector/CMSSW_10_6_37/src/rds/inspector/',
         'source $VO_CMS_SW_DIR/cmsset_default.sh',
         'export SCRAM_ARCH=slc7_amd64_gcc700', #export new arch
         'cmsenv',
         'echo ">>>> cmsenv activated"',
                 
         'cd /work/pahwagne/CMSSW_10_6_37/src/' ,
         'scramv1 runtime -sh',
         'cd /work/pahwagne/inspector/CMSSW_10_6_37/src/rds/inspector/test/'+dt_string ,

         'ls /scratch/pahwagne/',
         'cmsRun cfg_chunk_{1}.py'.format(dt_string,i),
         'xrdcp /scratch/pahwagne/{0}/{1}_chunk_{2}.root root://t3dcachedb03.psi.ch:1094///pnfs/psi.ch/cms/trivcat/store/user/pahwagne/inspector/{0}/{1}_chunk_{2}.root'.format(dt_string,naming,i),
         'rm /scratch/pahwagne/{0}/{1}_chunk_{2}.root'.format(dt_string,naming,i),
         '',
         'EOF',
         # --- close tmp file ---
         'echo "printing payload content:"',
         'cat /scratch/pahwagne/'+dt_string+'/apptainer-payload-{0}.sh'.format(i),
         'rm /scratch/pahwagne/{0}/{1}_chunk_{2}.root'.format(dt_string,naming,i),
         'echo ">>>> Done, launching cfg "',
         # --- make payload executable and run it in el7 singularity ---
         'chmod u+x "$payload"',
         '/cvmfs/cms.cern.ch/common/cmssw-el7  --bind /scratch,/work --command-to-run $payload'
         

     ])

  with open("{0}/submitter_chunk_{1}.sh".format(dt_string,i), "wt") as flauncher:
    flauncher.write(to_write)


  command_sh_batch = ' '.join([

        'sbatch',
        '-p '+queue,
        '--account=t3',
        '-o {0}/logs/chunk_{1}.log'.format(dt_string,i),
        '-e {0}/errs/chunk_{1}.err'.format(dt_string,i),
        #'--mem=1200M',
        '--job-name=MINI_{0}_{1}'.format(i,args.channel),
        '--time={0}'.format(time),
        '{0}/submitter_chunk_{1}.sh'.format(dt_string,i),
     ])

  print(command_sh_batch)
  os.system(command_sh_batch)







