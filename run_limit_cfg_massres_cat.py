import ROOT
import sys

import argparse
parser = argparse.ArgumentParser()
parser.add_argument('--bdt_low' , default = 0.992  , type = float)
parser.add_argument('--bdt_mid' , default = 0.988  , type = float)
parser.add_argument('--bdt_hig' , default = 0.995  , type = float)

parser.add_argument('--outdir'  , default = 'central', type = str  )
parser.add_argument('--workdir' , default = '.'      , type = str  )
parser.add_argument('--batch'   , action = 'store_true')

parser.add_argument('--asymptotic', action = 'store_true')

args = parser.parse_args()

ROOT.gROOT.SetBatch(args.batch)

sys.path.append('/gwpool/users/lguzzi/Tau3Mu/2017_2018/combine/multicategory/libs/')
from config_cls     import Configuration
from category_cls   import Category, Systematics

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

baseline += " & (abs(cand_charge) == 1 & abs(cand_refit_tau_mass - 1.8) < 0.2)"
baseline += " & (HLT_Tau3Mu_Mu5_Mu1_TkMu1_IsoTau10_Charge1_matched || HLT_Tau3Mu_Mu7_Mu1_TkMu1_IsoTau15_matched)"
baseline = ' '.join(baseline.split())

## UltraLegacy
MC_NORM   = 90369./(492e+3 + 500e+3)*(8580+11370)*0.1138/0.1063*1E-7
MC_NORM17 = 30541./(500e+3)*(8580+11370)*0.1138/0.1063*1E-7
MC_NORM18 = 59828./(492e+3)*(8580+11370)*0.1138/0.1063*1E-7
path_mc   = '/gwpool/users/lguzzi/Tau3Mu/2017_2018/BDT/singleclass/ntuples/signal_29mar2021.root'
path_data = '/gwpool/users/lguzzi/Tau3Mu/2017_2018/BDT/singleclass/ntuples/background_29mar2021.root'

## CREATE THE WORKSPACE
##
wspace = ROOT.RooWorkspace('wspace')

## NOTE careful with the ranges
getattr(wspace, 'import')(ROOT.RooRealVar('cand_refit_tau_mass' , '3-#mu mass'          , mass_range[0], mass_range[1], 'GeV'))
getattr(wspace, 'import')(ROOT.RooRealVar('cand_refit_tau_massE', '3-#mu mass error ^2' ,  0., 1, 'GeV^{2}'))
getattr(wspace, 'import')(ROOT.RooRealVar('bdt'                 , 'bdt'                 , -1 , 1))
getattr(wspace, 'import')(ROOT.RooRealVar('cand_charge'         , 'charge'              , -4 , 4))
getattr(wspace, 'import')(ROOT.RooRealVar('mcweight'            , 'mcweight'            ,  0., 1000))    ## MC PU and SFs weights
#getattr(wspace, 'import')(ROOT.RooRealVar('weight'              , 'weight'              ,  0., 1000))   ## NOTE NOT these! These are the mass-flattening weights!
getattr(wspace, 'import')(ROOT.RooRealVar('cand_refit_tau_eta'  , 'cand_refit_tau_eta'  , -5., 5))
getattr(wspace, 'import')(ROOT.RooRealVar('cand_refit_mass12'   , 'cand_refit_mass12'   ,  0., 2))
getattr(wspace, 'import')(ROOT.RooRealVar('cand_refit_mass13'   , 'cand_refit_mass13'   ,  0., 2))
getattr(wspace, 'import')(ROOT.RooRealVar('cand_refit_mass23'   , 'cand_refit_mass23'   ,  0., 2))
getattr(wspace, 'import')(ROOT.RooRealVar('cand_mass12'         , 'cand_mass12'         ,  0., 2))
getattr(wspace, 'import')(ROOT.RooRealVar('cand_mass13'         , 'cand_mass13'         ,  0., 2))
getattr(wspace, 'import')(ROOT.RooRealVar('cand_mass23'         , 'cand_mass23'         ,  0., 2))
getattr(wspace, 'import')(ROOT.RooRealVar('cand_charge12'       , 'cand_charge12'       , -3 , 3))
getattr(wspace, 'import')(ROOT.RooRealVar('cand_charge13'       , 'cand_charge13'       , -3 , 3))
getattr(wspace, 'import')(ROOT.RooRealVar('cand_charge23'       , 'cand_charge23'       , -3 , 3))
getattr(wspace, 'import')(ROOT.RooRealVar('mu1_hlt_doublemu3_trk_tau3mu_type', 'mu1_hlt_doublemu3_trk_tau3mu_type', 0 , 100))
getattr(wspace, 'import')(ROOT.RooRealVar('mu2_hlt_doublemu3_trk_tau3mu_type', 'mu2_hlt_doublemu3_trk_tau3mu_type', 0 , 100))
getattr(wspace, 'import')(ROOT.RooRealVar('mu3_hlt_doublemu3_trk_tau3mu_type', 'mu3_hlt_doublemu3_trk_tau3mu_type', 0 , 100))
getattr(wspace, 'import')(ROOT.RooRealVar('mu1_refit_muonid_tight', 'mu1_refit_muonid_tight', 0 , 1))
getattr(wspace, 'import')(ROOT.RooRealVar('mu2_refit_muonid_tight', 'mu2_refit_muonid_tight', 0 , 1))
getattr(wspace, 'import')(ROOT.RooRealVar('mu3_refit_muonid_tight', 'mu3_refit_muonid_tight', 0 , 1))
getattr(wspace, 'import')(ROOT.RooRealVar('HLT_Tau3Mu_Mu5_Mu1_TkMu1_IsoTau10_Charge1_matched', 'HLT_Tau3Mu_Mu5_Mu1_TkMu1_IsoTau10_Charge1_matched', 0 , 1))
getattr(wspace, 'import')(ROOT.RooRealVar('HLT_Tau3Mu_Mu7_Mu1_TkMu1_IsoTau15_matched', 'HLT_Tau3Mu_Mu7_Mu1_TkMu1_IsoTau15_matched', 0 , 1))
getattr(wspace, 'import')(ROOT.RooRealVar('year', 'year', 0 , 20))

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
    result_dir = args.outdir, work_dir = args.workdir,
)

category_selection = ' & '.join([
    '( bdt > {BDT} )',
    '( (sqrt(cand_refit_tau_massE) / cand_refit_tau_mass) >  {IWP} )',
    '( (sqrt(cand_refit_tau_massE) / cand_refit_tau_mass) <= {OWP} )',
    'year == {YAR}',
])

low17_sel = category_selection.format(BDT = 0.98, IWP = 0.000, OWP = 0.007, YAR = 17) ## 0.998
mid17_sel = category_selection.format(BDT = 0.98, IWP = 0.007, OWP = 0.012, YAR = 17) ## 0.994    
hig17_sel = category_selection.format(BDT = 0.98, IWP = 0.012, OWP = 10000, YAR = 17) ## 0.992

low18_sel = category_selection.format(BDT = 0.98, IWP = 0.000, OWP = 0.007, YAR = 18) ## 0.988    
mid18_sel = category_selection.format(BDT = 0.98, IWP = 0.007, OWP = 0.012, YAR = 18) ## 0.982        
hig18_sel = category_selection.format(BDT = 0.98, IWP = 0.012, OWP = 10000, YAR = 18) ## 0.992    

low17 = Category(name = 'A17', sig_norm = MC_NORM17, working_point = 'bdt{BDT}'.format(BDT = args.bdt_low), selection = low17_sel, wspace = wspace.Clone())
mid17 = Category(name = 'B17', sig_norm = MC_NORM17, working_point = 'bdt{BDT}'.format(BDT = args.bdt_mid), selection = mid17_sel, wspace = wspace.Clone())
hig17 = Category(name = 'C17', sig_norm = MC_NORM17, working_point = 'bdt{BDT}'.format(BDT = args.bdt_hig), selection = hig17_sel, wspace = wspace.Clone())

low18 = Category(name = 'A18', sig_norm = MC_NORM18, working_point = 'bdt{BDT}'.format(BDT = args.bdt_low), selection = low18_sel, wspace = wspace.Clone())
mid18 = Category(name = 'B18', sig_norm = MC_NORM18, working_point = 'bdt{BDT}'.format(BDT = args.bdt_mid), selection = mid18_sel, wspace = wspace.Clone())
hig18 = Category(name = 'C18', sig_norm = MC_NORM18, working_point = 'bdt{BDT}'.format(BDT = args.bdt_hig), selection = hig18_sel, wspace = wspace.Clone())

## SYSTEMATICS
##
common_systematics = [ 
    Systematics(name = 'lumi'     , distribution = 'lnN', value = '1.017' ),
    Systematics(name = 'xs_W'     , distribution = 'lnN', value = '1.037' ),
    Systematics(name = 'br_Wtaunu', distribution = 'lnN', value = '1.0021'),
    Systematics(name = 'br_Wmunu' , distribution = 'lnN', value = '1.0015'),
]
low17.add_systematics(common_systematics+[
    Systematics(name = 'muonID_A17', distribution = 'lnN', value = '1.011'),
])
mid17.add_systematics(common_systematics+[
    Systematics(name = 'muonID_B17', distribution = 'lnN', value = '1.014'),
])
hig17.add_systematics(common_systematics+[
    Systematics(name = 'muonID_C17', distribution = 'lnN', value = '1.017'),
])
low18.add_systematics(common_systematics+[
    Systematics(name = 'muonID_A18', distribution = 'lnN', value = '1.032'),
])
mid18.add_systematics(common_systematics+[
    Systematics(name = 'muonID_B18', distribution = 'lnN', value = '1.045'),
])
hig18.add_systematics(common_systematics+[
    Systematics(name = 'muonID_C18', distribution = 'lnN', value = '1.055'),
])

cfg.add_category(low17)
cfg.add_category(mid17)
cfg.add_category(hig17)
cfg.add_category(low18)
cfg.add_category(mid18)
cfg.add_category(hig18)

## RUN COMBINE
##
cfg.fit_model()
cfg.write_datacards()
cfg.combine_datacards()
cfg.run_combine(ntoys = 2000, grid = 0.5, asymptotic = args.asymptotic)
