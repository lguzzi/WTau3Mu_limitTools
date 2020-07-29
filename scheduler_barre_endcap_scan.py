import os
from itertools import product
import numpy

import argparse
parser = argparse.ArgumentParser()
parser.add_argument('--maxjobs'     , default = 5       , type = int)
parser.add_argument('--outdir'      , default = 'result', type = str)
parser.add_argument('--workdir'     , default = '.'     , type = str)
parser.add_argument('--run_local'   , action = 'store_true')
args = parser.parse_args()

CMD = 'python run_limit_cfg_barrel_endcap.py --batch --bdt_cut_barrel {BWP} --bdt_cut_endcap {EWP} --outdir {OUT} --workdir {WKD}'
JOB = '''\
Universe = vanilla          \n\
Executable = condor_exe.sh  \n\
use_x509userproxy = true    \n\
Log        = condor.log     \n\
Output     = condor.out     \n\
Error      = condor.error   \n\
getenv      = True          \n\
request_memory = 1024       \n\
queue 1
'''

barrel_wps = numpy.arange(0.99, 1, 0.001)
endcap_wps = barrel_wps
#barrel_wps = [0.99, 0.991, 0.992, 0.997, 0.999]
#endcap_wps = [0.99, 0.992, 0.993, 0.994, 0.999]

task_queue  = [CMD.format(BWP = bwp, EWP = ewp, OUT = args.outdir, WKD = args.workdir) for bwp, ewp in product(barrel_wps, endcap_wps)]
schedule    = [' & '.join(task_queue[jj:jj+args.maxjobs]) for jj in range(0, len(task_queue), args.maxjobs)]

if args.run_local:
    for sch in schedule: 
        exit_code = os.system(sch)

else:
    print '[ERR] not implemented'
    import sys
    sys.exit()
    import htcondor
    schedd = htcondor.Schedd()

    print '[INFO] Submitting jobs'

    for tsk in task_queue:
        with open('condor_exe.sh', 'w') as executable:
            executable.write('#!/bin/bash\n')
            executable.write(tsk)
        os.system('chmod +xwr condor_exe.sh')
        
        sub = htcondor.Submit(JOB)
        with schedd.transaction() as txn:
            sub.queue(txn)
    
    print '[INFO]', len(task_queue), 'jobs submitted'
    