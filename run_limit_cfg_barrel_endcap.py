import ROOT
import sys

import argparse
parser = argparse.ArgumentParser()
parser.add_argument('--bdt_cut_barrel', default = 0.996    , type = float)
parser.add_argument('--bdt_cut_endcap', default = 0.996    , type = float)
parser.add_argument('--outdir'        , default = 'result' , type = str  )
parser.add_argument('--workdir'       , default = '.'      , type = str  )
parser.add_argument('--batch', action = 'store_true')
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

baseline =  "   ((  ((abs(cand_refit_mass12-1.020)<0.02)*(cand_charge12==0))    + \
                    ((abs(cand_refit_mass13-1.020)<0.02)*(cand_charge13==0))    + \
                    ((abs(cand_refit_mass23-1.020)<0.02)*(cand_charge23==0))    ) == 0)"

#baseline += " & ((  ((abs(cand_refit_mass12-0.782)<0.02)*(cand_charge12==0))    +\
#                    ((abs(cand_refit_mass13-0.782)<0.02)*(cand_charge13==0))    +\
#                    ((abs(cand_refit_mass23-0.782)<0.02)*(cand_charge23==0))    ) == 0)"

baseline += " & ((  ((abs(cand_refit_mass12-0.547)<0.02)*(cand_charge12==0))    +\
                    ((abs(cand_refit_mass13-0.547)<0.02)*(cand_charge13==0))    +\
                    ((abs(cand_refit_mass23-0.547)<0.02)*(cand_charge23==0))    ) == 0)"

baseline += " & (abs(cand_charge) == 1 & abs(cand_refit_tau_mass - 1.8) < 0.2)"
#baseline = "abs(cand_charge) == 1 & abs(cand_refit_tau_mass - 1.8) < 0.2"
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

## NOTE match exactly Riccardo's script (round norm factor to float)
#FACTOR = float('%f' %(90480./(1.E6 + 902.E3)*(8580+11370)*0.1138/0.1063*1E-7))  / (90480./2.e6*(8580+11370)*0.1138/0.1063*1E-7)

cfg = Configuration(baseline = baseline, bkg_file_path = path_data, sig_file_path = path_mc, tree_name = 'tree',
    result_dir = args.outdir, work_dir = args.workdir, sig_norm = MC_NORM,
)

barrel = Category(name = 'barrel', working_point = str(args.bdt_cut_barrel), selection = "bdt > {CUT} & abs(cand_refit_tau_eta) <  1.6".format(CUT = args.bdt_cut_barrel), wspace = wspace.Clone())
endcap = Category(name = 'endcap', working_point = str(args.bdt_cut_endcap), selection = "bdt > {CUT} & abs(cand_refit_tau_eta) >= 1.6".format(CUT = args.bdt_cut_endcap), wspace = wspace.Clone())

cfg.add_category(barrel)
cfg.add_category(endcap)

## RUN COMBINE
##
cfg.fit_model()
cfg.write_datacards()
cfg.combine_datacards()
cfg.run_combine(ntoys = 5000, grid = 0.5)