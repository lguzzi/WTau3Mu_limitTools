import ROOT
import os
import math
ROOT.gInterpreter.Declare('#include "/gwpool/users/lguzzi/Tau3Mu/2017_2018/Fancy/CMS_lumi.h"')
ROOT.gInterpreter.Declare('#include "/gwpool/users/lguzzi/Tau3Mu/2017_2018/Fancy/TDR_style.h"')
ROOT.setTDRStyle()
class Systematics:
    def __init__(self, name, distribution, value):
        self.name = name
        self.distribution = distribution
        self.value = value

class Category:
    def __init__(self,  name, working_point, selection, wspace,
                        tree_name = None, sig_file_path = None, bkg_file_path = None, blind = True, pdf_dir = None, sig_norm = 1, useLabel = True,
                        Z_frac = 0., sigma = 0.01, dp = False, norm=1.0, flat=False):
        self.selection     = selection
        self.name          = name
        self.wspace        = wspace
        self.working_point = working_point
        self.label         = ''.join([self.name, self.working_point]) if useLabel else self.name
        self.Z_frac = Z_frac
        self.sigma = sigma
        self.discp = dp
        self.norm = norm
        self.blind = blind
        self.mass_range = (self.wspace.var('cand_refit_tau_mass').getMin(), self.wspace.var('cand_refit_tau_mass').getMax())
        self.nbins      = self.wspace.var('cand_refit_tau_mass').getBins()
        self.tree_name = tree_name

        self.sig_file_path = sig_file_path
        self.bkg_file_path = bkg_file_path
        self.pfd_dir       = pdf_dir
        self.sig_norm      = sig_norm
        self.flat          = flat

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

        self.sig_dataset = ROOT.RooDataSet('sig_dataset_%s' %self.name, '', self.sig_tree, self.wspace.allVars(), ' && '.join([self.selection, 'tau_sv_ls>2.0']), 'mcweight').reduce(
            ROOT.RooArgSet(self.wspace.var("cand_refit_tau_mass")), 
            "cand_refit_tau_mass > %s & cand_refit_tau_mass < %s" %(self.mass_range[0], self.mass_range[1])
        )
        #import pdb; pdb.set_trace()

        blinder = 'abs(cand_refit_tau_mass-1.78)>0.04'
        bkg_selection = ' & '.join([self.selection, blinder]) if self.blind else self.selection
        self.bkg_dataset = ROOT.RooDataSet('bkg_dataset_%s' %self.name, '', self.bkg_tree, self.wspace.allVars(), bkg_selection).reduce(
            ROOT.RooArgSet(self.wspace.var("cand_refit_tau_mass")))
        self.sig_integ      = self.sig_dataset.sumEntries()
        self.norm_sig_integ = self.sig_integ * self.sig_norm
        self.norm_bkg_integ = self.bkg_dataset.sumEntries()

    def fit_model(self):
        '''
        fit the signal and background parameters to their samples
        '''
        print '[INFO] fitting model for category', self.name
        self.frame = self.wspace.var('cand_refit_tau_mass').frame(40)

        self.load_files()
        
        self.fit_bkg()
        self.fit_sig()

        self.frame.GetYaxis().SetRangeUser(0., 12.)
        self.frame.GetXaxis().SetTitle('m(3mu) [GeV]')
        self.frame.GetYaxis().SetTitle('entries / 10 MeV')

        self.frame.SetTitle("")
        self.frame.Draw()

        ROOT.CMS_lumi(ROOT.gPad.GetPad(0), 5 if '17' in self.name else 6, 10)
        ROOT.gPad.SaveAs('%s/mass_%s_%dbins.pdf'%(self.pdf_dir, self.name, self.nbins))

        self.save_workspace()

    def save_workspace(self):
        '''
        create a new workspace to use as combine datacard
        '''
        self.out_wspace = ROOT.RooWorkspace('t3m_shapes')

        self.output = ROOT.TFile.Open('%s/CMS_T3M_13TeV_W_%s.root' %(self.workspace_dir, self.name), 'RECREATE')
        self.output.cd()

#        getattr(self.out_wspace, 'import')(ROOT.RooRealVar("bkgNorm_{}".format(self.name),"bkgNorm_{}".format(self.name),0,1e+6))
        self.out_wspace.factory('cand_refit_tau_mass[{RLO},{RHI}]'.format(  RLO = self.mass_range[0]    , 
                                                                            RHI = self.mass_range[1]   ))
        self.out_wspace.var('cand_refit_tau_mass').setBins(self.wspace.var("cand_refit_tau_mass").getBins())
        self.out_wspace.factory('mean[{VAL}]'  .format(  VAL = self.wspace.var('mean').getValV()))
        #self.out_wspace.factory('sigma[{VAL}]' .format(  VAL = self.wspace.var('sigma' ).getValV()))
        self.out_wspace.factory('sigma_{CAT}[{VAL},0.005,0.1]'.format(CAT = self.name, VAL = self.wspace.var('sigma').getValV(), ERR = self.wspace.var('sigma').getError()))

        if not self.discp:
           # self.out_wspace.factory("Exponential::bkg(cand_refit_tau_mass, a0_{CAT}[{VAL},-100, 100])".format(  CAT = self.name,
           #                                                                                                     VAL = self.wspace.var('slope').getValV()    ,
           #                                                                                                     ERR = self.wspace.var('slope').getError()  ))
            #self.out_wspace.factory("Bernstein::brn(cand_refit_tau_mass, b0_{CAT}[{VAL},-100, 100])".format(    CAT = self.name,
            #                                                                                                    VAL = self.wspace.var('bern1').getValV()    ,
            #                                                                                                    ERR = self.wspace.var('bern1').getError()  ))
            #self.out_wspace.factory('SUM::bkg(f0_{CAT}[{VAL},0,1]*exp, brn)'.format(    CAT=self.name,
            #                                                                            VAL=1.*self.wspace.var('cexpo').getValV()))
#            pol0 = ROOT.RooPolynomial("pol0", "", self.out_wspace.var("cand_refit_tau_mass"), ROOT.RooArgSet())
#            norm = ROOT.RooRealVar("bkgNorm_{}".format(self.name),"bkgNorm_{}".format(self.name),0,1e+3)
#            bkgf = ROOT.RooAddPdf('bkg', '', ROOT.RooArgList(pol0), ROOT.RooArgList(norm))
#            getattr(self.out_wspace, 'import')(bkgf)
            self.out_wspace.factory("Polynomial::bkg(cand_refit_tau_mass, {})")
#            self.out_wspace.factory("Bernstein::bkg(cand_refit_tau_mass, {{c0_{C}[{N},0,100]}})".format(C=self.name, N=self.wspace.var("nbkg").getVal()))
        else:
            getattr(self.out_wspace, 'import')(self.wspace.pdf('bkg'))
        #self.out_wspace.factory('RooGaussian::sig(cand_refit_tau_mass, mean, sigma)')
        self.out_wspace.factory('RooGaussian::sig(cand_refit_tau_mass, mean, sigma_{})'.format(self.name))
        
        self.out_wspace.var('mean' ).setConstant()
        #self.out_wspace.var('sigma').setConstant()

        if not self.blind:
            data = ROOT.RooDataSet('data_obs', '', self.bkg_dataset, ROOT.RooArgSet(self.wspace.var('cand_refit_tau_mass')))
        else:
            thepdf = self.wspace.pdf("bkg") if not self.discp else self.extend_from_multipdf()
            self.wspace.var("nbkg").setVal(0.01 if self.wspace.var("nbkg").getVal() < 1 else self.wspace.var("nbkg").getVal())
            self.wspace.var("nbkg").setVal(self.wspace.var("nbkg").getVal()*self.norm)
            data = ROOT.RooStats.AsymptoticCalculator.GenerateAsimovData(thepdf, ROOT.RooArgSet(self.wspace.var("cand_refit_tau_mass")))
            data.SetName("data_obs")

        getattr(self.out_wspace, 'import')(data)

        ## praise https://root-forum.cern.ch/t/ownership-in-python/9573/4
        ROOT.SetOwnership(self.wspace    , False)
        ROOT.SetOwnership(self.out_wspace, False)

        self.out_wspace.Write()
        self.output.Close()

    def extend_from_multipdf(self):
        bkg_pdf = self.wspace.pdf("bkg").getCurrentPdf()
        if not self.wspace.var("nbkg"):
            getattr(self.wspace, 'import')(ROOT.RooRealVar("nbkg", "",  2000, 0, 550000))
        return ROOT.RooAddPdf('bkg', '', ROOT.RooArgList(bkg_pdf), ROOT.RooArgList(self.wspace.var("nbkg")))

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
        #import pdb; pdb.set_trace()

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
                                                                                    ROOT.RooFit.SumW2Error(False)    ,
                                                                                    ROOT.RooFit.Extended(True)      )
        else:
            ext_bkg = self.extend_from_multipdf()
            self.bkg_fit_results = ext_bkg.fitTo(self.bkg_dataset,  ROOT.RooFit.Range('left,right') , 
                                                                    ROOT.RooFit.Save()              , 
                                                                    ROOT.RooFit.SumW2Error(False)    ,
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

        ## derive the correct expected number of events and its error
        #ctest   = ROOT.TCanvas()
        #slope   = ROOT.RooRealVar("slope", "", -10, 10)
        #expo    = ROOT.RooExponential("expo", "", self.wspace.var("cand_refit_tau_mass"), slope)
        #fres    = expo.fitTo(self.bkg_dataset, ROOT.RooFit.Range('left,right'), ROOT.RooFit.Save(), ROOT.RooFit.SumW2Error(True), ROOT.RooFit.Extended(False))
        #frame   = self.wspace.var("cand_refit_tau_mass").frame(40)
        #self.bkg_dataset.plotOn(frame)
        #expo.plotOn(frame)
        #expo.createIntegral(ROOT.RooArgSet(self.wspace.var("cand_refit_tau_mass")), ROOT.RooFit.NormSet(ROOT.RooArgSet(self.wspace.var("cand_refit_tau_mass"))), ROOT.RooFit.Range("mass_range")).getVal()
        can = ROOT.TCanvas()
        self.wspace.pdf('bkg').plotOn(frame, ROOT.RooFit.LineColor(2), ROOT.RooFit.Name('prefit'), ROOT.RooFit.LineStyle(ROOT.kDashed))
        leg.AddEntry(frame.findObject('prefit'), 'bkg pre-fit', 'l')
        frame.Draw()
        leg.Draw("SAME")
        if self.blind:
            can.SaveAs("{}/prefit_background_{}.pdf".format(self.pdf_dir, self.name), "pdf")

    def write_datacard(self):
        '''
        write the combine datacard from the pre-fitted roofit model
        '''
        ## FIXED
        ## bypass the incorrect mcweights given at ntuple production. Should be fixed
        #BYPASS = (90480./(1.E6 + 902.E3)*(8580+11370)*0.1138/0.1063*1E-7)
        #import pdb; pdb.set_trace()
        CLOPPER={
            'A17': (0.05,4.64),
            'B17': (0.29,2.45),
            'C17': (0.34,2.23),
            'A18': (0.05,4.64),
            'B18': (0.00,500.),
            'C18': (0.67,1.46),
        }
        category_systematics = '\n'.join([
            '{NAM}          {DIS}                       {VAL}               -'.format(
                NAM = sys.name, DIS = sys.distribution, VAL = sys.value
            ) for sys in self.systematics
        ])
        aset = ROOT.RooArgSet(self.wspace.var("cand_refit_tau_mass"))
        with open('%s/CMS_T3MSignal_13TeV_W_%s.txt' %(self.datacard_dir, self.name), 'w') as card:
            card.write(
'''
imax 1 number of bins
jmax * number of processes minus 1
kmax * number of nuisance parameters
--------------------------------------------------------------------------------
shapes bkg          Wtau3mu_{cat}       {wdr}/CMS_T3M_13TeV_W_{cat}.root t3m_shapes:bkg
shapes sig          Wtau3mu_{cat}       {wdr}/CMS_T3M_13TeV_W_{cat}.root t3m_shapes:sig
shapes data_obs     Wtau3mu_{cat}       {wdr}/CMS_T3M_13TeV_W_{cat}.root t3m_shapes:data_obs
--------------------------------------------------------------------------------
bin               Wtau3mu_{cat}
observation       {obs:d}
--------------------------------------------------------------------------------
bin                                     Wtau3mu_{cat}       Wtau3mu_{cat}
process                                 sig                 bkg
process                                 0                   1
rate                                    {signal:.4f}        {bkg:.4f}
--------------------------------------------------------------------------------
{SYS}
mc_stat_{cat} lnN                       {mcstat:.4f}        -   
bkgNorm_{cat} rateParam                 Wtau3mu_{cat}        bkg      1.    {FLAT}
--------------------------------------------------------------------------------
sigma_{cat} param {sigmaval:.4f} {sigmaerr:.4f}
'''.format(
        cat      = self.name,
        wdr      = os.path.abspath(self.workspace_dir),
        obs      = int(self.norm_bkg_integ) if self.blind==False else -1, # NOTE: norm_bkg_integ is actually the obsrved entries in the full range when running unblinded
        signal   = self.norm_sig_integ*(1.+self.Z_frac),
        bkg      = self.wspace.var("nbkg").getVal() if self.wspace.var("nbkg").getVal() > 0.01 else 0.01,
        # if testing non-extended prefit, extrapolate the background from the pdf integral ratio
        #bkg      = self.bkg_dataset.sumEntries()*self.wspace.pdf('bkg').createIntegral(aset, ROOT.RooFit.NormSet(aset), ROOT.RooFit.Range("mass_range")).getValV()/self.wspace.pdf('bkg').createIntegral(aset, ROOT.RooFit.NormSet(aset), ROOT.RooFit.Range("left,right")).getValV(),
        # mcstat is the relative uncertainty associated to 1 / N_GEN. Do the math and find sigma_{1/N_GEN} / (1/N_GEN) = 1 / sqrt(N_GEN)
        mcstat   = 1. + (1.-self.Z_frac) /  math.sqrt(self.sig_integ),
        sigmaval = self.wspace.var('sigma').getVal(),
        sigmaerr = self.wspace.var('sigma').getVal()*self.sigma,
        SYS  = category_systematics,
        FLAT = "[{LOW:.4f},{HIGH:.4f}]\nbkgNorm_{cat} flatParam".format(
#            HIGH=1+1./self.wspace.var("nbkg").getError() if self.wspace.var("nbkg").getVal() > 0.01    else 500 ,
#            LOW =1-1./self.wspace.var("nbkg").getError() if self.wspace.var("nbkg").getVal() > 1       else 0   ,
            HIGH=CLOPPER[self.name][1],
            LOW =CLOPPER[self.name][0],
            cat = self.name,
        )# if self.flat else ""
#fitbias_{cat} lnN                       {FB:.4f}            -
 #       bkg_pdf = 'a0_{cat} flatParam'.format(
 #           slopeval = self.wspace.var("slope").getVal(),
 #           slopeerr = self.wspace.var("slope").getError(),
 #           cat      = self.name,
 #          ) if not self.discp else 'roomultipdf_cat_W_{} discrete'.format(self.name)
        )
)
#{bkg_pdf}
#a0_{cat} flatParam
        
#b0_{cat} flatParam
#bkgNorm_{cat} flatParam                 Wtau3mu_{cat}        bkg      1.
#bkgNorm_{cat} rateParam                 Wtau3mu_{cat}        bkg      1.    [{LOW:.4f},{HIGH:.4f}]
#bkgNorm_{cat} flatParam
