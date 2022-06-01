import ROOT
import math

class Systematics:
    def __init__(self, name, distribution, value):
        self.name = name
        self.distribution = distribution
        self.value = value

class Category:
    def __init__(self,  name, working_point, selection, wspace,
                        tree_name = None, sig_file_path = None, bkg_file_path = None, blind = True, pdf_dir = None, sig_norm = 1, useLabel = True,
                        Z_frac = 0., sigma = 0.01, dp = False):
        self.selection     = selection
        self.name          = name
        self.wspace        = wspace
        self.working_point = working_point
        self.label         = ''.join([self.name, self.working_point]) if useLabel else self.name
        self.Z_frac = Z_frac
        self.sigma = sigma
        self.discp = dp

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
        self.sig_integ      = self.sig_dataset.sumEntries()
        self.norm_sig_integ = self.sig_integ * self.sig_norm
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
        self.out_wspace.factory('mean[{VAL}]'  .format(  VAL = self.wspace.var('mean').getValV()))
        #self.out_wspace.factory('sigma[{VAL}]' .format(  VAL = self.wspace.var('sigma' ).getValV()))
        self.out_wspace.factory('sigma_{CAT}[{VAL},0.005,0.1]'.format(CAT = self.label, VAL = self.wspace.var('sigma').getValV(), ERR = self.wspace.var('sigma').getError()))

        if not self.discp:
            self.out_wspace.factory("Exponential::bkg(cand_refit_tau_mass, a0_{CAT}[{VAL},-100, 100])".format(CAT = self.label,
                                                                                                                VAL = self.wspace.var('slope').getValV()    ,
                                                                                                                ERR = self.wspace.var('slope').getError()  ))
        else:
            getattr(self.out_wspace, 'import')(self.wspace.pdf('bkg'))
        #self.out_wspace.factory('RooGaussian::sig(cand_refit_tau_mass, mean, sigma)')
        self.out_wspace.factory('RooGaussian::sig(cand_refit_tau_mass, mean, sigma_{})'.format(self.label))
        
        self.out_wspace.var('mean' ).setConstant()
        #self.out_wspace.var('sigma').setConstant()

        data = ROOT.RooDataSet('data_obs', 'data_obs', self.bkg_dataset, ROOT.RooArgSet(self.wspace.var('cand_refit_tau_mass')))
        getattr(self.out_wspace, 'import')(data)

        ## praise https://root-forum.cern.ch/t/ownership-in-python/9573/4
        ROOT.SetOwnership(self.wspace    , False)
        ROOT.SetOwnership(self.out_wspace, False)

        self.out_wspace.Write()
        self.output.Close()

    def fit_sig(self):
        self.sig_fit_results =  self.wspace.pdf('sig').fitTo(self.sig_dataset,  ROOT.RooFit.Range('prefit_range') ,
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
        
        can = ROOT.TCanvas()
        frame = self.wspace.var("cand_refit_tau_mass").frame(40)
        self.sig_dataset.plotOn(frame)
        self.wspace.pdf('sig').plotOn(frame, ROOT.RooFit.LineColor(ROOT.kRed))
        frame.Draw()
        can.SaveAs("{}/prefit_signal_{}.pdf".format(self.pdf_dir, self.name), "pdf")

    def fit_bkg(self):
        frame = self.wspace.var("cand_refit_tau_mass").frame(40)
        leg = ROOT.TLegend()
        self.bkg_dataset.plotOn(frame,  #ROOT.RooFit.Binning(self.nbins)  , 
                                        ROOT.RooFit.MarkerSize(1.)  )
        self.bkg_dataset.plotOn(self.frame, #ROOT.RooFit.Binning(self.nbins)  , 
                                            ROOT.RooFit.MarkerSize(1.)  )

        if not self.discp:
            self.bkg_fit_results = self.wspace.pdf('bkg').fitTo(self.bkg_dataset,   ROOT.RooFit.Range('left,right') , 
                                                                                    ROOT.RooFit.Save()              , 
                                                                                    ROOT.RooFit.SumW2Error(True)    ,
                                                                                    ROOT.RooFit.Extended(True)      )
        else:
            bkg_pdf = self.wspace.pdf("bkg").getCurrentPdf()
            getattr(self.wspace, 'import')(ROOT.RooRealVar("nbkg", "",  2000, 0, 550000))
            ext_bkg = ROOT.RooAddPdf('bkg', '', ROOT.RooArgList(bkg_pdf), ROOT.RooArgList(self.wspace.var("nbkg")))
            self.bkg_fit_results = ext_bkg.fitTo(self.bkg_dataset,  ROOT.RooFit.Range('left,right') , 
                                                                    ROOT.RooFit.Save()              , 
                                                                    ROOT.RooFit.SumW2Error(True)    ,
                                                                    ROOT.RooFit.Extended(True)      )
        #for i in range(getattr(self.wspace.pdf('bkg'), 'getNumPdfs', lambda: 0)()):
        #    pdf = self.wspace.pdf('bkg').getPdf(i)
        #    num = ROOT.RooRealVar("num"+str(i) if i > 0 else 'nbkg', '', 2000, 0, 550000)
        #    epdf = ROOT.RooAddPdf("ext."+pdf.GetName(), '',  ROOT.RooArgList(pdf),  ROOT.RooArgList(num))
        #    getattr(self.wspace, 'import')(epdf)
        #    res = self.wspace.pdf(epdf.GetName()).fitTo(self.bkg_dataset,   ROOT.RooFit.Range('left,right') , 
        #                                        ROOT.RooFit.Save()              , 
        #                                        ROOT.RooFit.SumW2Error(True)    ,
        #                                        ROOT.RooFit.Extended(True)      )
        #    self.wspace.pdf(epdf.GetName()).plotOn(frame, ROOT.RooFit.LineColor(i+1), ROOT.RooFit.Name(self.wspace.pdf(epdf.GetName()).GetName()))
        #    leg.AddEntry(frame.findObject(self.wspace.pdf(epdf.GetName()).GetName()), self.wspace.pdf(epdf.GetName()).GetName(), 'l')
        
        self.wspace.pdf('bkg').plotOn(self.frame, ROOT.RooFit.LineColor(ROOT.kBlue))

        can = ROOT.TCanvas()
        self.wspace.pdf('bkg').plotOn(frame, ROOT.RooFit.LineColor(2), ROOT.RooFit.Name('prefit'), ROOT.RooFit.LineStyle(ROOT.kDashed))
        leg.AddEntry(frame.findObject('prefit'), 'bkg pre-fit', 'l')
        frame.Draw()
        leg.Draw("SAME")
        can.SaveAs("{}/prefit_background_{}.pdf".format(self.pdf_dir, self.name), "pdf")

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
        aset = ROOT.RooArgSet(self.wspace.var("cand_refit_tau_mass"))
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
sigma_{cat} param {sigmaval:.4f} {sigmaerr:.4f}
{bkg_pdf}
'''.format(
        cat      = self.label,
        wdr      = self.workspace_dir,
        obs      = int(self.norm_bkg_integ) if self.blind==False else -1, # NOTE: norm_bkg_integ is actually the obsrved entries in the full range when running unblinded
        signal   = self.norm_sig_integ*(1.+self.Z_frac),
        bkg      = self.wspace.var("nbkg").getVal(),
        # if testing non-extended prefit, extrapolate the background from the pdf integral ratio
        #bkg      = self.bkg_dataset.sumEntries()*self.wspace.pdf('bkg').createIntegral(aset, ROOT.RooFit.NormSet(aset), ROOT.RooFit.Range("mass_range")).getValV()/self.wspace.pdf('bkg').createIntegral(aset, ROOT.RooFit.NormSet(aset), ROOT.RooFit.Range("left,right")).getValV(),
        # mcstat is the relative uncertainty associated to 1 / N_GEN. Do the math and find sigma_{1/N_GEN} / (1/N_GEN) = 1 / sqrt(N_GEN)
        mcstat   = 1. + (1.-self.Z_frac) /  math.sqrt(self.sig_integ),
        sigmaval = self.wspace.var('sigma').getVal(),
        sigmaerr = self.wspace.var('sigma').getVal()*self.sigma,
        SYS  = category_systematics,
        bkg_pdf = 'a0_{cat}      flatParam   {slopeval:.4f}'.format(
            slopeval = self.wspace.var("slope").getVal(),
            slopeerr = self.wspace.var("slope").getError(),
            cat      = self.name,
           ) if not self.discp else 'roomultipdf_cat_W_{} discrete'.format(self.name)
        )
)
        