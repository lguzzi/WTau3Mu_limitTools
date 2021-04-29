import ROOT
import math

class Systematics:
    def __init__(self, name, distribution, value):
        self.name = name
        self.distribution = distribution
        self.value = value

class Category:
    def __init__(self,  name, working_point, selection, wspace,
                        tree_name = None, sig_file_path = None, bkg_file_path = None, blind = True, pdf_dir = None, sig_norm = 1, useLabel = True):
        self.selection     = selection
        self.name          = name
        self.wspace        = wspace
        self.working_point = working_point
        self.label         = ''.join([self.name, self.working_point]) if useLabel else self.name

        self.blind = blind
        self.mass_range = (self.wspace.var('cand_refit_tau_mass').getMin(), self.wspace.var('cand_refit_tau_mass').getMax())
        self.nbins      = self.wspace.var('cand_refit_tau_mass').getBins()
        self.tree_name = tree_name

        self.sig_file_path = sig_file_path
        self.bkg_file_path = bkg_file_path
        self.pfd_dir       = pdf_dir
        self.sig_norm      = sig_norm

        self.systematics = []

        self.datacard_dir  = '/'.join(['./', 'datacards'])
        self.result_dir    = '/'.join(['./', 'results'])
        self.workspace_dir = '/'.join(['./', 'workspaces'])

    def add_systematics(self, syst):
        for se in syst:
            self.systematics.append(se)
    def set_sig_file_path(self, file_path): self.sig_file_path = file_path
    def set_bkg_file_path(self, file_path): self.bkg_file_path = file_path
    def set_tree_name    (self, tree_name): self.tree_name     = tree_name
    def set_pdf_dir      (self, pdf_dir  ): self.pdf_dir       = pdf_dir
    def update_selection (self, baseline ): self.selection     = ' & '.join([self.selection, baseline])

    def load_files(self):
        '''
        load the signal and bakground datasets into roofit
        '''
        print '[INFO] loading files'
        self.sig_file = ROOT.TFile.Open(self.sig_file_path, 'READ')
        self.bkg_file = ROOT.TFile.Open(self.bkg_file_path, 'READ')

        self.sig_tree = self.sig_file.Get(self.tree_name)
        self.bkg_tree = self.bkg_file.Get(self.tree_name)

        self.sig_dataset = ROOT.RooDataSet('sig_dataset_%s' %self.label, '', self.sig_tree, self.wspace.allVars(), self.selection, 'mcweight').reduce(
            ROOT.RooArgSet(self.wspace.var("cand_refit_tau_mass")), 
            "cand_refit_tau_mass > %s & cand_refit_tau_mass < %s" %(self.mass_range[0], self.mass_range[1])
        )
        #import pdb; pdb.set_trace()

        blinder = 'abs(cand_refit_tau_mass-1.78)>0.04'
        bkg_selection = ' & '.join([self.selection, blinder]) if self.blind else self.selection
        self.bkg_dataset = ROOT.RooDataSet('bkg_dataset_%s' %self.label, '', self.bkg_tree, self.wspace.allVars(), bkg_selection).reduce(
            ROOT.RooArgSet(self.wspace.var("cand_refit_tau_mass")))
        self.norm_sig_integ = self.sig_dataset.sumEntries() * self.sig_norm
        self.norm_bkg_integ = self.bkg_dataset.sumEntries()

    def fit_model(self):
        '''
        fit the signal and background parameters to their samples
        '''
        print '[INFO] fitting model for category', self.name
        self.frame = self.wspace.var('cand_refit_tau_mass').frame()

        self.load_files()
        
        self.fit_bkg()
        self.fit_sig()

        self.frame.GetYaxis().SetRangeUser(0., 12.)

        self.frame.Draw()

        ROOT.gPad.SaveAs('%s/mass_%s_%dbins.pdf'%(self.pdf_dir, self.label, self.nbins))

        self.save_workspace()

    def save_workspace(self):
        '''
        create a new workspace to use as combine datacard
        '''
        self.out_wspace = ROOT.RooWorkspace('t3m_shapes')

        self.output = ROOT.TFile.Open('%s/CMS_T3M_13TeV_W_%s.root' %(self.workspace_dir, self.label), 'RECREATE')
        self.output.cd()

        self.out_wspace.factory('cand_refit_tau_mass[{RLO},{RHI}]'.format(  RLO = self.mass_range[0]    , 
                                                                            RHI = self.mass_range[1]   ))
        self.out_wspace.factory('mean[{VAL}]' .format(  VAL = self.wspace.var('mean' ).getValV()))
        self.out_wspace.factory('sigma[{VAL}]'.format(  VAL = self.wspace.var('sigma').getValV()))

        self.out_wspace.factory("Exponential::bkg(cand_refit_tau_mass, a0_{CAT}[{VAL},{ERR},{ERR}])".format(CAT = self.label,
                                                                                                            VAL = self.wspace.var('slope').getValV()    ,
                                                                                                            ERR = self.wspace.var('slope').getError()  ))
        self.out_wspace.factory('RooGaussian::sig(cand_refit_tau_mass, mean, sigma)')
        
        self.out_wspace.var('mean' ).setConstant()
        self.out_wspace.var('sigma').setConstant()

        data = ROOT.RooDataSet('data_obs', 'data_obs', self.bkg_dataset, ROOT.RooArgSet(self.wspace.var('cand_refit_tau_mass')))
        getattr(self.out_wspace, 'import')(data)

        ## praise https://root-forum.cern.ch/t/ownership-in-python/9573/4
        ROOT.SetOwnership(self.wspace    , False)
        ROOT.SetOwnership(self.out_wspace, False)

        self.out_wspace.Write()
        self.output.Close()

    def fit_sig(self):
        self.sig_fit_results =  self.wspace.pdf('sig').fitTo(self.sig_dataset,  ROOT.RooFit.Range('sig_region') ,
                                                                                ROOT.RooFit.Save()              , 
                                                                                ROOT.RooFit.SumW2Error(True)    )
        self.sig_dataset.plotOn(self.frame              , 
            #ROOT.RooFit.Binning(self.nbins)             ,
            ROOT.RooFit.Rescale(self.sig_norm)          ,
            ROOT.RooFit.DrawOption('B')                 , 
            ROOT.RooFit.DataError(ROOT.RooAbsData.None) , 
            ROOT.RooFit.XErrorSize(0)                   , 
            ROOT.RooFit.LineWidth(2)                    ,
            ROOT.RooFit.FillColor(ROOT.kRed)            ,
            ROOT.RooFit.FillStyle(3003)                 ) 
        self.wspace.pdf('sig').plotOn(self.frame, ROOT.RooFit.LineColor(ROOT.kRed))

    def fit_bkg(self):

        self.bkg_fit_results = self.wspace.pdf('bkg').fitTo(self.bkg_dataset,   ROOT.RooFit.Range('left,right') , 
                                                                                ROOT.RooFit.Save()              , 
                                                                                ROOT.RooFit.SumW2Error(True)    )
        self.bkg_dataset.plotOn(self.frame, #ROOT.RooFit.Binning(self.nbins)  , 
                                            ROOT.RooFit.MarkerSize(1.)  )
        self.wspace.pdf('bkg').plotOn(self.frame, ROOT.RooFit.LineColor(ROOT.kBlue))

    def write_datacard(self):
        '''
        write the combine datacard from the pre-fitted roofit model
        '''
        ## FIXED
        ## bypass the incorrect mcweights given at ntuple production. Should be fixed
        #BYPASS = (90480./(1.E6 + 902.E3)*(8580+11370)*0.1138/0.1063*1E-7)
        #import pdb; pdb.set_trace()

        category_systematics = '\n'.join([
            '{NAM}          {DIS}                       {VAL}               -'.format(
                NAM = sys.name, DIS = sys.distribution, VAL = sys.value
            ) for sys in self.systematics
        ])

        with open('%s/CMS_T3MSignal_13TeV_W_%s.txt' %(self.datacard_dir, self.label), 'w') as card:
            card.write(
'''
imax 1 number of bins
jmax * number of processes minus 1
kmax * number of nuisance parameters
--------------------------------------------------------------------------------
shapes background    Wtau3mu_{cat}       {wdr}/CMS_T3M_13TeV_W_{cat}.root t3m_shapes:bkg
shapes signal        Wtau3mu_{cat}       {wdr}/CMS_T3M_13TeV_W_{cat}.root t3m_shapes:sig
shapes data_obs      Wtau3mu_{cat}       {wdr}/CMS_T3M_13TeV_W_{cat}.root t3m_shapes:data_obs
--------------------------------------------------------------------------------
bin               Wtau3mu_{cat}
observation       {obs:d}
--------------------------------------------------------------------------------
bin                                     Wtau3mu_{cat}       Wtau3mu_{cat}
process                                 signal              background
process                                 0                   1
rate                                    {signal:.4f}        {bkg:.4f}
--------------------------------------------------------------------------------
{SYS}
mc_stat_{cat} lnN                       {mcstat:.4f}        -   
--------------------------------------------------------------------------------
bkgNorm_{cat} rateParam                 Wtau3mu_{cat}        background      1.
a0_{cat}      param   {slopeval:.4f} {slopeerr:.4f}
'''.format(
         cat      = self.label,
         wdr      = self.workspace_dir,
         obs      = self.norm_bkg_integ if self.blind==False else -1,
         signal   = self.norm_sig_integ,
         bkg      = self.wspace.var("nbkg").getVal(),
         mcstat   = 1. + math.sqrt(self.norm_sig_integ / self.sig_norm) / (self.norm_sig_integ / self.sig_norm),
         slopeval = self.wspace.var("slope").getVal(),  
         slopeerr = self.wspace.var("slope").getError(),
         SYS  = category_systematics,
         )
)
'''
mu_hlt{cat}   lnN                       {mu_hlt:.4f}        -   
trk_hlt{cat}  lnN                       {trk_hlt:.4f}       -   
hlt_extrap    lnN                       1.05                -   
'''
        