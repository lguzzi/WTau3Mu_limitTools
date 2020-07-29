import ROOT
import sys

import argparse
parser = argparse.ArgumentParser()
parser.add_argument('--inner'   , default = 0.007   , type = float)
parser.add_argument('--outer'   , default = 0.012   , type = float)

parser.add_argument('--bdt_low' , default = 0.996  , type = float)
parser.add_argument('--bdt_mid' , default = 0.996  , type = float)
parser.add_argument('--bdt_hig' , default = 0.996  , type = float)

parser.add_argument('--outdir'  , default = 'result', type = str  )
parser.add_argument('--workdir' , default = '.'     , type = str  )
parser.add_argument('--batch'   , action = 'store_true')
args = parser.parse_args()

ROOT.gROOT.SetBatch(args.batch)

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
MC_NORM = 90480./(390032 + 381627)*(8580+11370)*0.1138/0.1063*1E-7
#MC_NORM = 1

baseline =  "   (   ((abs(cand_refit_mass12-1.020)<0.02)*(cand_charge12==0))    + \
                    ((abs(cand_refit_mass13-1.020)<0.02)*(cand_charge13==0))    + \
                    ((abs(cand_refit_mass23-1.020)<0.02)*(cand_charge23==0))    ) == 0"

#baseline += " & (   ((abs(cand_refit_mass12-0.782)<0.02)*(cand_charge12==0))    +\
#                    ((abs(cand_refit_mass13-0.782)<0.02)*(cand_charge13==0))    +\
#                    ((abs(cand_refit_mass23-0.782)<0.02)*(cand_charge23==0))    ) == 0"

baseline += " & abs(cand_charge) == 1 & abs(cand_refit_tau_mass - 1.8) < 0.2"
baseline = ' '.join(baseline.split())


path_data = "samples/background_26may2020.root"
path_mc   = "samples/signal_26may2020.root"

## CREATE THE WORKSPACE
##
wspace = ROOT.RooWorkspace('wspace')

## NOTE careful with the ranges
getattr(wspace, 'import')(ROOT.RooRealVar('cand_refit_tau_mass', '3-#mu mass'          , mass_range[0], mass_range[1], 'GeV'))
getattr(wspace, 'import')(ROOT.RooRealVar('bdt'                , 'bdt'                 , -1 , 1))
getattr(wspace, 'import')(ROOT.RooRealVar('cand_charge'        , 'charge'              , -4 , 4))
getattr(wspace, 'import')(ROOT.RooRealVar('mcweight'           , 'mcweight'            ,  0., 1000))    ## MC PU and SFs weights
#getattr(wspace, 'import')(ROOT.RooRealVar('weight'             , 'weight'              ,  0., 1000))   ## NOTE NOT these! These are the mass-flattening weights!
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
wspace.factory("RooGaussian::sig(cand_refit_tau_mass, mean[1.78, -1.7, 1.9], sigma[0.02, 0, 0.1])")

# see this https://root-forum.cern.ch/t/fit-only-the-sidebands-yield-on-full-range-using-rooextendpdf/31868
slope = ROOT.RooRealVar('slope', '', -0.001, -1e3, 1e3)
nbkg  = ROOT.RooRealVar('nbkg', 'nbkg', 2000, 0, 550000)
#expo  = ROOT.RooPolynomial('bkg_expo', 'bkg_expo', wspace.var('cand_refit_tau_mass'), ROOT.RooArgList(slope))
expo  = ROOT.RooExponential('bkg_expo', 'bkg_expo', wspace.var('cand_refit_tau_mass'), slope)

expomodel = ROOT.RooAddPdf('bkg', '', ROOT.RooArgList(expo), ROOT.RooArgList(nbkg))
getattr(wspace, 'import')(expomodel)

## CREATE THE CFG OBJECT
##
cfg = Configuration(baseline = baseline, bkg_file_path = path_data, sig_file_path = path_mc, tree_name = 'tree',
    result_dir = args.outdir, work_dir = args.workdir, sig_norm = MC_NORM,
)

low_sel = ' & '.join([
    '( bdt_ > {BDT} )'.format(BDT = args.bdt_low),
    '( (sqrt(cand_refit_tau_massE) / cand_refit_tau_mass) <= {IWP} )'.format(IWP = args.inner),
])
mid_sel = ' & '.join([
    '( bdt_ > {BDT} )'.format(BDT = args.bdt_mid),
    '( (sqrt(cand_refit_tau_massE) / cand_refit_tau_mass) >  {IWP} )'.format(IWP = args.inner),
    '( (sqrt(cand_refit_tau_massE) / cand_refit_tau_mass) <= {OWP} )'.format(OWP = args.outer),
])
hig_sel = ' & '.join([
    '( bdt_ > {BDT} )'.format(BDT = args.bdt_hig),
    '( (sqrt(cand_refit_tau_massE) / cand_refit_tau_mass) >  {OWP} )'.format(OWP = args.outer),
])

low = Category(name = 'low', working_point = 'bdt{BDT}_i{IWP}_o{OWP}'.format(BDT = args.bdt_low, IWP = args.inner, OWP = args.outer), selection = low_sel, wspace = wspace.Clone())
mid = Category(name = 'mid', working_point = 'bdt{BDT}_i{IWP}_o{OWP}'.format(BDT = args.bdt_mid, IWP = args.inner, OWP = args.outer), selection = mid_sel, wspace = wspace.Clone())
hig = Category(name = 'hig', working_point = 'bdt{BDT}_i{IWP}_o{OWP}'.format(BDT = args.bdt_hig, IWP = args.inner, OWP = args.outer), selection = hig_sel, wspace = wspace.Clone())

cfg.add_category(low)
cfg.add_category(mid)
cfg.add_category(hig)

## RUN COMBINE
##
cfg.fit_model()
cfg.write_datacards()
cfg.combine_datacards()
cfg.run_combine(ntoys = 5000)