import ROOT
import os
from itertools import product
from random import uniform
from time import sleep

## TODO
## - should I apply the correct signal normalization at ntuple level or right here by default? [DONE]
## - check difference in ROOT and RooFit Integral() value (difference of about 1e-4) [DONE]

COMBINE_CMD = 'combineCards.py {IN} > {OUT}'
CREATE_MODEL_CMD = 'text2workspace.py --X-assign-flatParam-prior {IN} -o {OUT}'
RUN_COMBINE_HYBRID_CMD  = "combine -M {METHOD} --X-rtd MINIMIZER_freezeDisassociatedParams --testStat={STAT} --fitNuisances={BLIND} --frequentist {MODEL} -T {NTOYS} --expectedFromGrid {GRID} -C {CL}  --plot='{PDF}/limit_combined_hybridnew_{CL}_WP_{LABEL}.pdf' --rMin -1 --rMax 10 --setParameterRanges {PARAMETERS} {DISCRETE} | grep Limit >  '{RES}/limit_combined_hybridnew_CL_{CL}_central_WP_{LABEL}.txt'"
#RUN_COMBINE_HYBRID_CMD = "combine -M HybridNew {WSP} {BLI} --LHCmode LHC-limits --X-rtd MINIMIZER_freezeDisassociatedParams -T {TOY} -C {CL}  --plot='limit_combined_hybridnew_{CL}.pdf' --rMin {r} --rMax {R} {PAR} | grep Limit >  '{RES}/limit_combined_hybridnew_CL_{CL}_central_WP_{LABEL}.txt'""
RUN_COMBINE_ASYMPTOTIC_CMD = 'combineTool.py -M AsymptoticLimits  --run blind  -d  {CARD} --cl {CL} | grep Expected > "{RES}/limit_combined_asymptotic_CL_{CL}_WP_{LABEL}.txt"'
class Configuration:
    def __init__(self,  baseline, sig_file_path, bkg_file_path, tree_name,
                        result_dir = 'central', work_dir = '.',
                        sig_norm = 1):
        self.baseline      = baseline
        self.sig_file_path = sig_file_path
        self.bkg_file_path = bkg_file_path
        self.tree_name     = tree_name
        self.sig_norm      = sig_norm
        self.categories    = []

        self.result_dir    = '/'.join(['.', work_dir, result_dir])
        self.datacard_dir  = '/'.join([self.result_dir, 'datacards'])
        self.workspace_dir = '/'.join([self.result_dir, 'workspaces'])
        self.pdf_dir       = '/'.join([self.result_dir, 'pdf'])

        ## mutlithread makedirs safety
        for _ in range(100):
            sleep(uniform(1, 3./1000))
        if not os.path.exists(self.pdf_dir)         : os.makedirs(self.pdf_dir)
        if not os.path.exists(self.datacard_dir)    : os.makedirs(self.datacard_dir)
        if not os.path.exists(self.result_dir)      : os.makedirs(self.result_dir)
        if not os.path.exists(self.workspace_dir)   : os.makedirs(self.workspace_dir)

    def add_category(self, category):
        print '[INFO] adding category', category.name, '(%s)' %"BLINDED" if category.blind else "NOT BLINDED"
        print '[INFO] signal normalization coefficient is', category.sig_norm

        category.set_sig_file_path(self.sig_file_path)
        category.set_bkg_file_path(self.bkg_file_path)
        category.set_tree_name    (self.tree_name    )
        category.set_pdf_dir      (self.pdf_dir      )
        category.update_selection (self.baseline     )

        category.datacard_dir   = self.datacard_dir
        category.result_dir     = self.result_dir
        category.workspace_dir  = self.workspace_dir

        self.categories.append(category)
    
    def fit_model(self):
        for cc in self.categories: 
            cc.fit_model()
    
    def write_datacards(self):
        for cc in self.categories:
            cc.write_datacard()

    def combine_datacards(self):
        global COMBINE_CMD
        global CREATE_MODEL_CMD

        print '[INFO] merging datacards'
        
        inputs = ['{NAM}={DTC}/CMS_T3MSignal_13TeV_W_{CAT}.txt'.format(DTC = self.datacard_dir, NAM = cc.name, CAT = cc.name) for cc in self.categories]
        input_str = ' '.join(inputs)
        
        self.output_datacard_path   = '/'.join([self.datacard_dir , 'CMS_T3MSignal_13TeV_W_Combined.txt'])
        self.output_model_path = '/'.join([self.workspace_dir, 'CMS_T3M_13TeV_W_Combined.root'])
        #output_datacard_path   = self.datacard_dir + '/datacard_comb_' + '_'.join([cc.label for cc in self.categories]) + '.txt'
        #self.output_model_path = self.datacard_dir + '/model_comb_'    + '_'.join([cc.label for cc in self.categories]) + '.root'
        os.system(COMBINE_CMD.format(IN = input_str, OUT = self.output_datacard_path))
        print '[INFO] creating combined model'
        os.system(CREATE_MODEL_CMD.format(IN = self.output_datacard_path, OUT = self.output_model_path))

    def run_combine(self,   ntoys, cl = 0.90, method = 'HybridNew', stat = 'LHC', grid = 0.5, asymptotic = False):
        global RUN_COMBINE_HYBRID_CMD, RUN_COMBINE_ASYMPTOTIC_CMD

        print '[INFO] Running combine'

        label = '_'.join([cc.name for cc in self.categories])

        parameters  = [ROOT.RooArgList(cat.out_wspace.pdf('bkg').getParameters(ROOT.RooArgSet(cat.wspace.var("cand_refit_tau_mass"))))[i] for cat in self.categories for i in range(ROOT.RooArgList(cat.out_wspace.pdf('bkg').getParameters(ROOT.RooArgSet(cat.wspace.var("cand_refit_tau_mass")))).getSize())] # What do a zoo owner and a Python data analyst have in common? They both import panda
        parameters += [ROOT.RooArgList(cat.out_wspace.pdf('sig').getParameters(ROOT.RooArgSet(cat.wspace.var("cand_refit_tau_mass"))))[i] for cat in self.categories for i in range(ROOT.RooArgList(cat.out_wspace.pdf('sig').getParameters(ROOT.RooArgSet(cat.wspace.var("cand_refit_tau_mass")))).getSize())] # What was a python's first words? print("s"*10)
        parameters  = [p for p in parameters if not p.isConstant() and not 'roomultipdf_cat_W' in p.GetName()]
        
        parameters_selection  = ['bkgNorm_{}=0,1000000'.format(cat.name) for cat in self.categories]
        parameters_selection += ['{P}={m},{M}'.format(P=p.GetName(), m=p.getMin(), M=p.getMax()) for p in parameters]
        parameters_selection  =  ':'.join(parameters_selection)

        discrete_profiling =  ','.join(['{}={}'.format("roomultipdf_cat_W_{}".format(cat.name), cat.out_wspace.pdf("bkg").getCurrentIndex()) for cat in self.categories if cat.discp])
        isblind=all(cc.blind for cc in self.categories)
        RUN_CMD = RUN_COMBINE_HYBRID_CMD.format(    LABEL       = label,
                                                    MODEL       = self.output_model_path,
                                                    NTOYS       = ntoys,
                                                    BLIND       = int(all(not cc.blind for cc in self.categories)),
                                                    CL          = cl,
                                                    METHOD      = method,
                                                    STAT        = stat,
                                                    GRID        = grid,
                                                    PARAMETERS  = parameters_selection,
                                                    RES         = self.result_dir,
                                                    PDF         = self.pdf_dir,
                                                    DISCRETE    = '--setParameters '+discrete_profiling if discrete_profiling else '',
        ) if not asymptotic else RUN_COMBINE_ASYMPTOTIC_CMD.format( CL    = cl, 
                                                                    CARD  = self.output_datacard_path,
                                                                    RES   = self.result_dir, 
                                                                    LABEL = label,
        )
        print ">"*50
        print RUN_CMD
        print ">"*50
        os.system(RUN_CMD)
