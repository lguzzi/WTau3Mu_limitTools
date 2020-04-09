import os
from itertools import product
import numpy

import argparse
parser = argparse.ArgumentParser()
parser.add_argument('--maxjobs', default = 5       , type = int)
parser.add_argument('--outdir' , default = 'result', type = str)
args = parser.parse_args()

CMD = 'python run_limit_cfg_barrel_endcap.py --batch --bdt_cut_barrel {BWP} --bdt_cut_endcap {EWP} --outdir {OUT}'

barrel_wps = numpy.arange(0.990, 0.999, 0.001)
endcap_wps = barrel_wps

task_queue  = [CMD.format(BWP = bwp, EWP = ewp, OUT = args.outdir) for bwp, ewp in product(barrel_wps, endcap_wps)]
schedule    = [' & '.join(task_queue[jj:jj+args.maxjobs]) for jj in range(0, len(task_queue), args.maxjobs)]

for sch in schedule: 
    exit_code = os.system(sch)
