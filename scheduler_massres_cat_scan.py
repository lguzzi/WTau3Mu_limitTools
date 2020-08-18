import os
from itertools import product
import numpy

import argparse
parser = argparse.ArgumentParser()
parser.add_argument('--maxjobs', default = 5       , type = int)
parser.add_argument('--outdir' , default = 'result', type = str)
args = parser.parse_args()

CMD = 'python run_limit_cfg_massres_cat_sing.py --batch --bdt_lo {CUT} --outdir {OUT}'

eta_wps = numpy.arange(0.0, 0.2, 0.01)

task_queue = [CMD.format(CUT = cut, OUT = args.outdir) for cut in eta_wps]
schedule   = [' & '.join(task_queue[jj:jj+args.maxjobs]) for jj in range(0, len(task_queue), args.maxjobs)]

for sch in schedule: 
    exit_code = os.system(sch)
