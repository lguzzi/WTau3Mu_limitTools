import ROOT

class Category:
    def __init__(self, name, selection, wspace, tree_name = None, sig_file_path = None, bkg_file_path = None, blind = True):
        self.selection  = selection
        self.name       = name
        self.wspace     = wspace

        self.blind = blind

        self.tree_name = tree_name

        self.wspace.var('a0{CAT}').SetName( wspace.var('a0{CAT}').GetName().format(CAT = self.name) )

        self.sig_file_path = sig_file_path
        self.bkg_file_path = bkg_file_path

    def set_sig_file_path(self, file_path): self.sig_file_path = file_path
    def set_bkg_file_path(self, file_path): self.bkg_file_path = file_path
    def set_tree_name    (self, tree_name): self.tree_name     = tree_name
    def set_pdf_dir      (self, pdf_dir  ): self.pdf_dir       = pdf_dir
    def update_selection (self, baseline ): self.selection     = ' & '.join([self.selection, baseline])

    def load_files(self):
        print '[INFO] loading files'
        self.sig_file = ROOT.TFile.Open(self.sig_file_path, 'READ')
        self.bkg_file = ROOT.TFile.Open(self.bkg_file_path, 'READ')

        self.sig_tree = self.sig_file.Get(self.tree_name)
        self.bkg_tree = self.bkg_file.Get(self.tree_name)

        self.sig_dataset = ROOT.RooDataSet('sig_dataset_%s' %self.name, '', self.sig_tree, self.wspace.allVars(), self.selection, 'mcweight').reduce(
            ROOT.RooArgSet(self.wspace.var("cand_refit_tau_mass")))

        blinder = 'abs(cand_refit_tau_mass-1.78)>0.06'
        bkg_selection = ' & '.join([self.selection, blinder]) if self.blind else self.selection
        self.bkg_dataset = ROOT.RooDataSet('bkg_dataset_%s' %self.name, '', self.bkg_tree, self.wspace.allVars(), bkg_selection).reduce(
            ROOT.RooArgSet(self.wspace.var("cand_refit_tau_mass")))

    def fit_model(self):
        print '[INFO] fitting model for category', self.name
        self.frame = self.wspace.var('cand_refit_tau_mass').frame()

        self.load_files()
        
        self.fit_sig()
        self.fit_bkg()

        self.frame.Draw()

        ROOT.gPad.SaveAs('%s/mass_%s_%dbins.pdf'%(self.pdf_dir, self.name, self.wspace.var("cand_refit_tau_mass").getBins()))

        self.save_workspace()

    def save_workspace(self):
        self.out_wspace = ROOT.RooWorkspace('t3m_shapes')

        self.output = ROOT.TFile.Open('datacards/datacard%s.root' %self.name, 'RECREATE')
        self.output.cd()

        self.out_wspace.factory('cand_refit_tau_mass[{RLO},{RHI}]'.format(  RLO = self.wspace.var('cand_refit_tau_mass').getMin(), 
                                                                            RHI = self.wspace.var('cand_refit_tau_mass').getMax()))
        self.out_wspace.factory('mean[{VAL}]' .format(  VAL = self.wspace.var('mean' ).getVal()))
        self.out_wspace.factory('sigma[{VAL}]'.format(  VAL = self.wspace.var('sigma').getVal()))

        self.out_wspace.factory("Exponential::bkg(cand_refit_tau_mass, a0{CAT}[{VAL},{ERR},{ERR}])".format( CAT = self.name,
                                                                                                            VAL = self.wspace.var('slope').getVal() ,
                                                                                                            ERR = self.wspace.var('slope').getValE()))
        self.out_wspace.factory('RooGaussian::sig(cand_refit_tau_mass, mean, sigma)')
        
        self.out_wspace.var('mean' ).setConstant()
        self.out_wspace.var('sigma').setConstant()

        data = ROOT.RooDataSet('data_obs', 'data_obs', self.bkg_dataset, ROOT.RooArgSet(self.wspace.var('cand_refit_tau_mass')))
        getattr(self.out_wspace, 'import')(data)
        import pdb; pdb.set_trace()
        self.out_wspace.Write()
        self.output.Close()


    def fit_sig(self):
        self.sig_fit_results =  self.wspace.pdf('sig').fitTo(self.sig_dataset,  ROOT.RooFit.Range('sig_region') ,
                                                                                ROOT.RooFit.Save()              , 
                                                                                ROOT.RooFit.SumW2Error(True)    )
        self.sig_frame = self.wspace.var('cand_refit_tau_mass').frame()
        self.sig_dataset.plotOn(self.frame              , 
            #ROOT.RooFit.Binning(self.nbins)             , 
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
        