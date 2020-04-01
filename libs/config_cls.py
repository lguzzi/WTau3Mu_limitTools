import ROOT
import os

## TODO
## - should I apply the correct signal normalization at ntuple level or right here by default?

COMBINE_CMD = 'combineCards.py {IN} > {OUT}'
CREATE_MODEL_CMD = 'text2workspace.py {IN} -o {OUT}'
RUN_COMBINE_HYBRID_CMD  = "combine -n {LABEL} -M {METHOD} --testStat={STAT} --frequentist {MODEL} -T {NTOYS} --expectedFromGrid {GRID} -C {CL}  --plot='{PDF}/limit_combined_hybridnew_scan_central_{CL}_WP_{LABEL}.pdf' --rMin 0 --rMax 10 --setParameterRanges {PARAMETERS} | grep Limit >  '{RES}/limit_combined_hybridnew_CL_{CL}_central_WP_{LABEL}.txt'"

class Configuration:
    def __init__(self,  baseline, sig_file_path, bkg_file_path, tree_name,
                        pdf_dir = 'pdf', datacard_dir = 'datacards', result_dir = 'result',
                        sig_norm = 90480./(1.E6 + 902.E3)*(8580+11370)*0.1138/0.1063*1E-7):
        self.baseline      = baseline
        self.sig_file_path = sig_file_path
        self.bkg_file_path = bkg_file_path
        self.tree_name     = tree_name
        self.sig_norm      = sig_norm
        self.pdf_dir       = pdf_dir             
        self.categories    = []

        self.pdf_dir       = pdf_dir 
        self.datacard_dir  = datacard_dir 
        self.result_dir    = result_dir 

        if not os.path.exists(self.pdf_dir)     : os.mkdir(self.pdf_dir)
        if not os.path.exists(self.datacard_dir): os.mkdir(self.datacard_dir)
        if not os.path.exists(self.result_dir)  : os.mkdir(self.result_dir)

    def add_category(self, category):
        print '[INFO] adding category', category.name, '(%s)' %"BLINDED" if category.blind else "NOT BLINDED"
        print '[INFO] signal normalization coefficient is', self.sig_norm

        category.set_sig_file_path(self.sig_file_path)
        category.set_bkg_file_path(self.bkg_file_path)
        category.set_tree_name    (self.tree_name    )
        category.set_pdf_dir      (self.pdf_dir      )
        category.update_selection (self.baseline     )
        category.set_sig_norm     (self.sig_norm     )

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
        
        inputs = ['{NAM}={DTC}/datacard_{NAM}.txt'.format(DTC = self.datacard_dir, NAM = cc.name) for cc in self.categories]
        input_str                   = ' '.join(inputs)
        
        output_datacard_path   = self.datacard_dir + '/datacard_comb_' + '_'.join([cc.label for cc in self.categories]) + '.txt'
        
        self.output_model_path = self.datacard_dir + '/model_comb_'    + '_'.join([cc.label for cc in self.categories]) + '.root'

        os.system(COMBINE_CMD.format(IN = input_str, OUT = output_datacard_path))

        print '[INFO] creating combined model'
        os.system(CREATE_MODEL_CMD.format(IN = output_datacard_path, OUT = self.output_model_path))

    def run_combine(self, ntoys, cl = 0.90, method = 'HybridNew', stat = 'LHC', grid = 0.5):
        global RUN_COMBINE_HYBRID_CMD

        print '[INFO] Running combine'

        label = '_'.join([cc.label for cc in self.categories])

        parameters_selection =  ['bkgNorm_barrel$wpb=0,1000000'] + \
                                ['bkgNorm_{CAT}=0,1000000:a0_{CAT}=-100,100'.format(CAT = cc.label) for cc in self.categories]
        parameters_selection =  ':'.join(parameters_selection)

        os.system(RUN_COMBINE_HYBRID_CMD.format(    LABEL       = label,
                                                    MODEL       = self.output_model_path,
                                                    NTOYS       = ntoys,
                                                    CL          = cl,
                                                    METHOD      = method,
                                                    STAT        = stat,
                                                    GRID        = grid,
                                                    PARAMETERS  = parameters_selection,
                                                    RES         = self.result_dir,
                                                    PDF         = self.pdf_dir,
        ))
