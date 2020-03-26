import ROOT
import sys

sys.path.append('./libs/')
from config_cls     import Configuration
from category_cls   import Category

## GLOBAL PARAMETERS
##
mass_range  = (1.6, 2.0)
sig_range   = (1.72, 1.84)
right_range = (1.84, 2.0)
left_range  = (1.6, 1.72)
nbins       = 40

baseline =  "   (   ((abs(cand_refit_mass12-1.020)<0.02)*(cand_charge12==0))    + \
                    ((abs(cand_refit_mass13-1.020)<0.02)*(cand_charge13==0))    + \
                    ((abs(cand_refit_mass23-1.020)<0.02)*(cand_charge23==0))    ) == 0"

baseline += " & (   ((abs(cand_refit_mass12-0.782)<0.02)*(cand_charge12==0))    +\
                    ((abs(cand_refit_mass13-0.782)<0.02)*(cand_charge13==0))    +\
                    ((abs(cand_refit_mass23-0.782)<0.02)*(cand_charge23==0))    ) == 0"

baseline += " & abs(cand_charge) == 1 & abs(cand_refit_tau_mass - 1.8) < 0.2"

path_data = "samples/background_10jan2020.root"
path_mc   = "samples/signal_10jan2020.root"

## CREATE THE WORKSPACE
##
wspace = ROOT.RooWorkspace('wspace')

getattr(wspace, 'import')(ROOT.RooRealVar('cand_refit_tau_mass', '3-#mu mass'          , mass_range[0], mass_range[1], 'GeV'))
getattr(wspace, 'import')(ROOT.RooRealVar('bdt'                , 'bdt'                 , -1 , 1))
getattr(wspace, 'import')(ROOT.RooRealVar('cand_charge'        , 'charge'              , -4 , 4))
getattr(wspace, 'import')(ROOT.RooRealVar('mcweight'           , 'mcweight'            ,  0., 5))
getattr(wspace, 'import')(ROOT.RooRealVar('cand_refit_tau_eta' , 'cand_refit_tau_eta'  , -5., 5))
getattr(wspace, 'import')(ROOT.RooRealVar('cand_refit_mass12'  , 'cand_refit_mass12'   ,  0., 2))
getattr(wspace, 'import')(ROOT.RooRealVar('cand_refit_mass13'  , 'cand_refit_mass13'   ,  0., 2))
getattr(wspace, 'import')(ROOT.RooRealVar('cand_refit_mass23'  , 'cand_refit_mass23'   ,  0., 2))
getattr(wspace, 'import')(ROOT.RooRealVar('cand_mass12'        , 'cand_mass12'         ,  0., 2))
getattr(wspace, 'import')(ROOT.RooRealVar('cand_mass13'        , 'cand_mass13'         ,  0., 2))
getattr(wspace, 'import')(ROOT.RooRealVar('cand_mass23'        , 'cand_mass23'         ,  0., 2))
getattr(wspace, 'import')(ROOT.RooRealVar('cand_charge12'      , 'cand_charge12'       , -3 , 3))
getattr(wspace, 'import')(ROOT.RooRealVar('cand_charge13'      , 'cand_charge13'       , -3 , 3))
getattr(wspace, 'import')(ROOT.RooRealVar('cand_charge23'      , 'cand_charge23'       , -3 , 3))
getattr(wspace, 'import')(ROOT.RooRealVar('mu1_hlt_doublemu3_trk_tau3mu_type', 'mu1_hlt_doublemu3_trk_tau3mu_type', 0 , 100))
getattr(wspace, 'import')(ROOT.RooRealVar('mu2_hlt_doublemu3_trk_tau3mu_type', 'mu2_hlt_doublemu3_trk_tau3mu_type', 0 , 100))
getattr(wspace, 'import')(ROOT.RooRealVar('mu3_hlt_doublemu3_trk_tau3mu_type', 'mu3_hlt_doublemu3_trk_tau3mu_type', 0 , 100))
getattr(wspace, 'import')(ROOT.RooRealVar('mu1_refit_muonid_tight', 'mu1_refit_muonid_tight', 0 , 1))
getattr(wspace, 'import')(ROOT.RooRealVar('mu2_refit_muonid_tight', 'mu2_refit_muonid_tight', 0 , 1))
getattr(wspace, 'import')(ROOT.RooRealVar('mu3_refit_muonid_tight', 'mu3_refit_muonid_tight', 0 , 1))

wspace.var("cand_refit_tau_mass").setBins(nbins)
wspace.var("cand_refit_tau_mass").setRange("left"  , *left_range )
wspace.var("cand_refit_tau_mass").setRange("right" , *right_range)
wspace.var("cand_refit_tau_mass").setRange("sig_region", *sig_range)
wspace.factory("RooGaussian::sig(cand_refit_tau_mass, mean[1.77, 1.7, 1.8], sigma[0.02, 0.001, 0.5])")

# see this https://root-forum.cern.ch/t/fit-only-the-sidebands-yield-on-full-range-using-rooextendpdf/31868
slope = ROOT.RooRealVar('a0{CAT}', '', -0.001, -1e3, 1e3)
nbkg  = ROOT.RooRealVar('nbkg', 'nbkg', 2000, 0, 550000)
expo  = ROOT.RooPolynomial('bkg_expo', 'bkg_expo', wspace.var('cand_refit_tau_mass'), ROOT.RooArgList(slope))

expomodel = ROOT.RooAddPdf('bkg', '', ROOT.RooArgList(expo), ROOT.RooArgList(nbkg))
getattr(wspace, 'import')(expomodel)

## CREATE THE CFG OBJECT
##
cfg = Configuration(baseline = baseline, bkg_file_path = path_data, sig_file_path = path_mc, tree_name = 'tree')

barrel = Category(name = 'barrel', selection = "bdt > 0.996 & abs(cand_refit_tau_eta) <  1.6", wspace = wspace.Clone())
endcap = Category(name = 'endcap', selection = "bdt > 0.996 & abs(cand_refit_tau_eta) >= 1.6", wspace = wspace.Clone())

cfg.add_category(barrel)
cfg.add_category(endcap)

## RUN COMBINE
##
cfg.fit_model()
#cfg.write_datacards()
#cfg.combine_datacards()
#cfg.run_combine()