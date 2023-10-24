import ROOT
import sys
from itertools import product
import argparse
parser = argparse.ArgumentParser()

parser.add_argument('--outdir'  , default = 'central', type = str, help = 'output sub-dir')
parser.add_argument('--cl'      , default = 0.9, type = float, help = 'Confidence level')
parser.add_argument('--quantile', default = 0.5, type = float, help = "\
    -3 sigma: 0.002 ;\n\
    -2 sigma: 0.023 ;\n\
    -1 sigma: 0.159 ;\n\
    MEDIAN: 0.5 ;\n\
    +1 sigma: 0.841 ;\n\
    +2 sigma: 0.977 ;\n\
    +3 sigma: 0.998 ;\n"
)
parser.add_argument('--workdir' , default = '.'      , type = str, help = 'main output dir')
parser.add_argument('--unblind' , action = 'store_true', help = 'run UNBLINDED limits')
parser.add_argument('--batch'   , action = 'store_true')
parser.add_argument('--discrete-profiling', action = 'store_true')

parser.add_argument('--asymptotic', action = 'store_true')

args = parser.parse_args()

ROOT.gROOT.SetBatch(args.batch)

sys.path.append('/gwpool/users/lguzzi/Tau3Mu/2017_2018/combine/multicategory/libs/')
from config_cls     import Configuration
from category_cls   import Category, Systematics

## GLOBAL PARAMETERS
##
mass_range  = (1.6, 2.0)
sig_range   = (1.74, 1.82)
right_range = (1.82, 2.0)
left_range  = (1.6, 1.74)
prefit_range= (1.7, 1.85)
nbins       = 100

baseline =  "   (   ((abs(cand_refit_mass12-1.020)<0.02)*(cand_charge12==0))    + \
                    ((abs(cand_refit_mass13-1.020)<0.02)*(cand_charge13==0))    + \
                    ((abs(cand_refit_mass23-1.020)<0.02)*(cand_charge23==0))    ) == 0"
baseline += " & (   ((abs(cand_refit_mass12-0.782)<0.02)*(cand_charge12==0))    + \
                    ((abs(cand_refit_mass13-0.782)<0.02)*(cand_charge13==0))    + \
                    ((abs(cand_refit_mass23-0.782)<0.02)*(cand_charge23==0))    ) == 0"

# muon POG selection 
baseline += " & (abs(cand_charge) == 1 & abs(cand_refit_tau_mass - 1.8) < 0.2)"
baseline += " & (HLT_Tau3Mu_Mu5_Mu1_TkMu1_IsoTau10_Charge1_matched==1 || HLT_Tau3Mu_Mu7_Mu1_TkMu1_IsoTau15_Charge1_matched==1)"
baseline += " & ((mu1_refit_pt > 3.5 & abs(mu1_refit_eta) < 1.2) || (mu1_refit_pt > 2.0 & abs(mu1_refit_eta) >= 1.2 & abs(mu1_refit_eta) < 2.4))"
baseline += " & ((mu2_refit_pt > 3.5 & abs(mu2_refit_eta) < 1.2) || (mu2_refit_pt > 2.0 & abs(mu2_refit_eta) >= 1.2 & abs(mu2_refit_eta) < 2.4))"
baseline += " & ((mu3_refit_pt > 3.5 & abs(mu3_refit_eta) < 1.2) || (mu3_refit_pt > 2.0 & abs(mu3_refit_eta) >= 1.2 & abs(mu3_refit_eta) < 2.4))"

baseline += " & (mu1_refit_muonid_medium==1 & mu2_refit_muonid_medium==1 & mu3_refit_muonid_medium==1)"
# HLT reinforcement
baseline += " & (mu1_refit_pt > 7 & mu2_refit_pt > 1 & mu3_refit_pt > 1)"
baseline += " & (cand_refit_dR12 < 0.5 || cand_refit_dR13 < 0.5 || cand_refit_dR23 < 0.5)"
baseline += " & (cand_refit_mass12 < 1.9 || cand_refit_mass13 < 1.9 || cand_refit_mass23 < 1.9)"
baseline += " & (cand_refit_tau_pt > 15)"
baseline += " & (abs(cand_refit_tau_eta) < 2.5)"
baseline += " & (tau_sv_ls>2)"

baseline = ' '.join(baseline.split())

## UltraLegacy
MC_NORM   = 90369./(492e+3 + 500e+3)*(8580+11370)*0.1138/0.1063*1E-7
#MC_NORM17 = (0.918/0.89)*30541./(500e+3)*(8580+11370)*0.1138/0.1063*1E-7
MC_NORM17 = 30541./(500e+3)*(8580+11370)*0.1138/0.1063*1E-7
#MC_NORM17 = 27120./(500e+3)*(8580+11370)*0.1138/0.1063*1E-7
MC_NORM18 = 59828./(492e+3)*(8580+11370)*0.1138/0.1063*1E-7

#path_mc   = '/gwpool/users/lguzzi/Tau3Mu/2017_2018/BDT/singleclass/ntuples/signal_10jun2021-deepMET-threeMedium-v2.root'
#path_mc   = '/gwpool/users/lguzzi/Tau3Mu/2017_2018/BDT/singleclass/ntuples/signal_threeMedium_WEIGHTS_11OCT2021.root'
#path_data = '/gwpool/users/lguzzi/Tau3Mu/2017_2018/BDT/singleclass/ntuples/background_10jun2021-deepMET-threeMedium-v2.root'

path_mc   = '/gwpool/users/lguzzi/Tau3Mu/2017_2018/BDT/singleclass/ntuples/signal_threeMedium_weighted_16Mar2022.root'
path_data = '/gwpool/users/lguzzi/Tau3Mu/2017_2018/BDT/singleclass/ntuples/background_threeMedium-UNBLINDED.root'
import pdb; pdb.set_trace()
#path_mc   = '/gwpool/users/lguzzi/Tau3Mu/2017_2018/BDT/singleclass/ntuples/signal_threeMedium_weighted_6Apr2022-catA.root'
#path_data = '/gwpool/users/lguzzi/Tau3Mu/2017_2018/BDT/singleclass/ntuples/background_threeMedium_weighted_6Apr2022-catA.root'

## CREATE THE WORKSPACE
##
wspace = ROOT.RooWorkspace('wspace')

## NOTE careful with the ranges
getattr(wspace, 'import')(ROOT.RooRealVar('cand_refit_tau_mass' , '3-#mu mass'          , mass_range[0], mass_range[1], 'GeV'))
getattr(wspace, 'import')(ROOT.RooRealVar('cand_refit_tau_massE', '3-#mu mass error ^2' ,  0., 100, 'GeV^{2}'))
getattr(wspace, 'import')(ROOT.RooRealVar('bdt'                 , 'bdt'                 , -1 , 2))
getattr(wspace, 'import')(ROOT.RooRealVar('cand_charge'         , 'charge'              , -4 , 4))
getattr(wspace, 'import')(ROOT.RooRealVar('mcweight'            , 'mcweight'            ,  0., 1e+9))    ## MC PU and SFs weights
#getattr(wspace, 'import')(ROOT.RooRealVar('weight'              , 'weight'              ,  0., 1000))   ## NOTE NOT these! These are the mass-flattening weights!

getattr(wspace, 'import')(ROOT.RooRealVar('tau_sv_ls'           , 'tau_sv_ls'           , -1000, 1000))
getattr(wspace, 'import')(ROOT.RooRealVar('tau_sv_prob'         , 'tau_sv_prob'         , -1000, 1000))
getattr(wspace, 'import')(ROOT.RooRealVar('cand_refit_tau_eta'  , 'cand_refit_tau_eta'  , -5, 5 ))
getattr(wspace, 'import')(ROOT.RooRealVar('cand_refit_tau_pt'   , 'cand_refit_tau_pt'   , 0, 1000))
getattr(wspace, 'import')(ROOT.RooRealVar('cand_refit_mass12'   , 'cand_refit_mass12'   ,  0., 100))
getattr(wspace, 'import')(ROOT.RooRealVar('cand_refit_mass13'   , 'cand_refit_mass13'   ,  0., 100))
getattr(wspace, 'import')(ROOT.RooRealVar('cand_refit_mass23'   , 'cand_refit_mass23'   ,  0., 100))
#getattr(wspace, 'import')(ROOT.RooRealVar('cand_mass12'         , 'cand_mass12'         ,  0., 100))
#getattr(wspace, 'import')(ROOT.RooRealVar('cand_mass13'         , 'cand_mass13'         ,  0., 100))
#getattr(wspace, 'import')(ROOT.RooRealVar('cand_mass23'         , 'cand_mass23'         ,  0., 100))
getattr(wspace, 'import')(ROOT.RooRealVar('cand_charge12'       , 'cand_charge12'       , -3 , 3))
getattr(wspace, 'import')(ROOT.RooRealVar('cand_charge13'       , 'cand_charge13'       , -3 , 3))
getattr(wspace, 'import')(ROOT.RooRealVar('cand_charge23'       , 'cand_charge23'       , -3 , 3))
getattr(wspace, 'import')(ROOT.RooRealVar('mu1_hlt_doublemu3_trk_tau3mu_type', 'mu1_hlt_doublemu3_trk_tau3mu_type', 0 , 100))
getattr(wspace, 'import')(ROOT.RooRealVar('mu2_hlt_doublemu3_trk_tau3mu_type', 'mu2_hlt_doublemu3_trk_tau3mu_type', 0 , 100))
getattr(wspace, 'import')(ROOT.RooRealVar('mu3_hlt_doublemu3_trk_tau3mu_type', 'mu3_hlt_doublemu3_trk_tau3mu_type', 0 , 100))

getattr(wspace, 'import')(ROOT.RooRealVar('mu1_refit_eta', 'mu1_refit_eta', -4 , 4))
getattr(wspace, 'import')(ROOT.RooRealVar('mu2_refit_eta', 'mu2_refit_eta', -4 , 4))
getattr(wspace, 'import')(ROOT.RooRealVar('mu3_refit_eta', 'mu3_refit_eta', -4 , 4))
getattr(wspace, 'import')(ROOT.RooRealVar('mu1_refit_pt', 'mu1_refit_pt', 0 , 1000))
getattr(wspace, 'import')(ROOT.RooRealVar('mu2_refit_pt', 'mu2_refit_pt', 0 , 1000))
getattr(wspace, 'import')(ROOT.RooRealVar('mu3_refit_pt', 'mu3_refit_pt', 0 , 1000))

getattr(wspace, 'import')(ROOT.RooRealVar('cand_refit_dR12', 'cand_refit_dR12', 0 , 100))
getattr(wspace, 'import')(ROOT.RooRealVar('cand_refit_dR13', 'cand_refit_dR13', 0 , 100))
getattr(wspace, 'import')(ROOT.RooRealVar('cand_refit_dR23', 'cand_refit_dR23', 0 , 100))

getattr(wspace, 'import')(ROOT.RooRealVar('mu1_refit_muonid_tight', 'mu1_refit_muonid_tight', -1 , 2))
getattr(wspace, 'import')(ROOT.RooRealVar('mu2_refit_muonid_tight', 'mu2_refit_muonid_tight', -1 , 2))
getattr(wspace, 'import')(ROOT.RooRealVar('mu3_refit_muonid_tight', 'mu3_refit_muonid_tight', -1 , 2))
getattr(wspace, 'import')(ROOT.RooRealVar('mu1_refit_muonid_soft', 'mu1_refit_muonid_soft', -1 , 2))
getattr(wspace, 'import')(ROOT.RooRealVar('mu2_refit_muonid_soft', 'mu2_refit_muonid_soft', -1 , 2))
getattr(wspace, 'import')(ROOT.RooRealVar('mu3_refit_muonid_soft', 'mu3_refit_muonid_soft', -1 , 2))
getattr(wspace, 'import')(ROOT.RooRealVar('mu1_refit_muonid_medium', 'mu1_refit_muonid_medium', -1 , 2))
getattr(wspace, 'import')(ROOT.RooRealVar('mu2_refit_muonid_medium', 'mu2_refit_muonid_medium', -1 , 2))
getattr(wspace, 'import')(ROOT.RooRealVar('mu3_refit_muonid_medium', 'mu3_refit_muonid_medium', -1 , 2))

getattr(wspace, 'import')(ROOT.RooRealVar('HLT_Tau3Mu_Mu5_Mu1_TkMu1_IsoTau10_Charge1_matched', 'HLT_Tau3Mu_Mu5_Mu1_TkMu1_IsoTau10_Charge1_matched', -1 , 2))
getattr(wspace, 'import')(ROOT.RooRealVar('HLT_Tau3Mu_Mu7_Mu1_TkMu1_IsoTau15_Charge1_matched', 'HLT_Tau3Mu_Mu7_Mu1_TkMu1_IsoTau15_Charge1_matched', -1 , 2))
getattr(wspace, 'import')(ROOT.RooRealVar('year', 'year', 0 , 20))

wspace.var("cand_refit_tau_mass").setBins(nbins)
wspace.var("cand_refit_tau_mass").setRange("left"  , *left_range )
wspace.var("cand_refit_tau_mass").setRange("right" , *right_range)
wspace.var("cand_refit_tau_mass").setRange("sig_region", *sig_range)
wspace.var("cand_refit_tau_mass").setRange("prefit_range", *prefit_range)
wspace.var("cand_refit_tau_mass").setRange("mass_range", *mass_range)
wspace.factory("RooGaussian::sig(cand_refit_tau_mass, mean[1.78, -1.7, 1.9], sigma[0.02, 0, 0.1])")

# see this https://root-forum.cern.ch/t/fit-only-the-sidebands-yield-on-full-range-using-rooextendpdf/31868
slope = ROOT.RooRealVar('slope', '', -0.001, -1e3, 1e3)
poly1 = ROOT.RooRealVar("poly1", '', 0, -3, 3)
poly2 = ROOT.RooRealVar("poly2", '', 0, -3, 3)
bern1 = ROOT.RooRealVar("bern1", '', 0.1, 0, 1)
pow1  = ROOT.RooRealVar("pow1", '',-3 ,-6 ,-0.0001)
nbkg  = ROOT.RooRealVar('nbkg', 'nbkg', 2000, 0, 550000)
# multidim pdf 
# http://cms-analysis.github.io/HiggsAnalysis-CombinedLimit/part3/nonstandard/#discrete-profiling
#expo  = ROOT.RooPolynomial('bkg_expo', 'bkg_expo', wspace.var('cand_refit_tau_mass'), ROOT.RooArgList(slope))
expo = ROOT.RooExponential('bkg_expo', 'bkg_expo', wspace.var('cand_refit_tau_mass'), slope)
poly = ROOT.RooChebychev  ('bkg_poly', 'bkg_poly', wspace.var('cand_refit_tau_mass'), ROOT.RooArgList(poly1))
pol0 = ROOT.RooPolynomial ('bkg_pol0', 'bkg_pol0', wspace.var('cand_refit_tau_mass'), ROOT.RooArgList())
powl = ROOT.RooGenericPdf ('bkg_powl', "TMath::Power(@0,@1)", ROOT.RooArgList(wspace.var('cand_refit_tau_mass'), pow1))
bern = ROOT.RooBernstein('bkg_bern', 'bkg_bern', wspace.var('cand_refit_tau_mass'), ROOT.RooArgList(bern1))

nexpo = ROOT.RooRealVar('nexpo', 'nexpo', 2000, 0, 550000)
npoly = ROOT.RooRealVar('npoly', 'npoly', 2000, 0, 550000)
npol0 = ROOT.RooRealVar('npol0', 'npol0', 2000, 0, 550000)
nbern = ROOT.RooRealVar('nbern', 'nbern', 2000, 0, 550000)
npowl = ROOT.RooRealVar('npowl', 'npowl', 2000, 0, 550000)

eexpo = ROOT.RooAddPdf('eexpo', '', ROOT.RooArgList(expo), ROOT.RooArgList(nexpo))
epoly = ROOT.RooAddPdf('epoly', '', ROOT.RooArgList(poly), ROOT.RooArgList(npoly))
epowl = ROOT.RooAddPdf('epowl', '', ROOT.RooArgList(powl), ROOT.RooArgList(npowl))
ebern = ROOT.RooAddPdf('ebern', '', ROOT.RooArgList(bern), ROOT.RooArgList(nbern))
epol0 = ROOT.RooAddPdf('epol0', '', ROOT.RooArgList(pol0), ROOT.RooArgList(npol0))

cexpo = ROOT.RooRealVar('cexpo', '', 0.5, 0, 1)
bexpo = ROOT.RooAddPdf('bexp', '', ROOT.RooArgList(expo, bern), ROOT.RooArgList(cexpo))

wspace_A17 = wspace.Clone()
wspace_B17 = wspace.Clone()
wspace_C17 = wspace.Clone()
wspace_A18 = wspace.Clone()
wspace_B18 = wspace.Clone()
wspace_C18 = wspace.Clone()

def load_dp(path):
    fil = ROOT.TFile.Open(path, "READ")
    wsp = fil.Get("ospace")
    rmp = wsp.pdf("multipdf").Clone("bkg")
    return rmp

if args.discrete_profiling:
    getattr(wspace_A17, 'import')(load_dp('/gwpool/users/lguzzi/Tau3Mu/2017_2018/combine_test/T3MuCombine/python/MultiPdfWorkspaces/W_A17.root'))
    getattr(wspace_B17, 'import')(load_dp('/gwpool/users/lguzzi/Tau3Mu/2017_2018/combine_test/T3MuCombine/python/MultiPdfWorkspaces/W_B17.root'))
    getattr(wspace_C17, 'import')(load_dp('/gwpool/users/lguzzi/Tau3Mu/2017_2018/combine_test/T3MuCombine/python/MultiPdfWorkspaces/W_C17.root'))
    getattr(wspace_A18, 'import')(load_dp('/gwpool/users/lguzzi/Tau3Mu/2017_2018/combine_test/T3MuCombine/python/MultiPdfWorkspaces/W_A18.root'))
    getattr(wspace_B18, 'import')(load_dp('/gwpool/users/lguzzi/Tau3Mu/2017_2018/combine_test/T3MuCombine/python/MultiPdfWorkspaces/W_B18.root'))
    getattr(wspace_C18, 'import')(load_dp('/gwpool/users/lguzzi/Tau3Mu/2017_2018/combine_test/T3MuCombine/python/MultiPdfWorkspaces/W_C18.root'))
else:
#    background = ROOT.RooAddPdf('bkg', '', ROOT.RooArgList(bexpo), ROOT.RooArgList(nbkg))
    background = ROOT.RooAddPdf('bkg', '', ROOT.RooArgList(pol0), ROOT.RooArgList(nbkg))
    #background = ROOT.RooExponential('bkg', 'bkg_expo', wspace.var('cand_refit_tau_mass'), slope)
    getattr(wspace_A17, 'import')(background.Clone('bkg'))
    getattr(wspace_B17, 'import')(background.Clone('bkg'))
    getattr(wspace_C17, 'import')(background.Clone('bkg'))
    getattr(wspace_A18, 'import')(background.Clone('bkg'))
    getattr(wspace_B18, 'import')(background.Clone('bkg'))
    getattr(wspace_C18, 'import')(background.Clone('bkg'))

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

low17_sel = category_selection.format(BDT = 0.991, IWP = 0.000, OWP = 0.007, YAR = 17) ## 0.991
mid17_sel = category_selection.format(BDT = 0.994, IWP = 0.007, OWP = 0.012, YAR = 17) ## 0.994
hig17_sel = category_selection.format(BDT = 0.992, IWP = 0.012, OWP = 10000, YAR = 17) ## 0.992

low18_sel = category_selection.format(BDT = 0.995, IWP = 0.000, OWP = 0.007, YAR = 18) ## 0.995
mid18_sel = category_selection.format(BDT = 0.998, IWP = 0.007, OWP = 0.012, YAR = 18) ## 0.998
hig18_sel = category_selection.format(BDT = 0.994, IWP = 0.012, OWP = 10000, YAR = 18) ## 0.994

low17 = Category(norm=1.0, flat=True , blind = not args.unblind, name = 'A17', sig_norm = MC_NORM17, working_point = 'bdt{BDT}'.format(BDT = 0.991), selection = low17_sel, wspace = wspace_A17, sigma = 0.02, dp = args.discrete_profiling)#, Z_frac = 0.117)
mid17 = Category(norm=1.0, flat=False, blind = not args.unblind, name = 'B17', sig_norm = MC_NORM17, working_point = 'bdt{BDT}'.format(BDT = 0.994), selection = mid17_sel, wspace = wspace_B17, sigma = 0.06, dp = args.discrete_profiling)#, Z_frac = 0.085)
hig17 = Category(norm=1.0, flat=False, blind = not args.unblind, name = 'C17', sig_norm = MC_NORM17, working_point = 'bdt{BDT}'.format(BDT = 0.992), selection = hig17_sel, wspace = wspace_C17, sigma = 0.02, dp = args.discrete_profiling)#, Z_frac = 0.090)
low18 = Category(norm=1.0, flat=True , blind = not args.unblind, name = 'A18', sig_norm = MC_NORM18, working_point = 'bdt{BDT}'.format(BDT = 0.995), selection = low18_sel, wspace = wspace_A18, sigma = 0.02, dp = args.discrete_profiling)#, Z_frac = 0.052)
mid18 = Category(norm=1.0, flat=False, blind = not args.unblind, name = 'B18', sig_norm = MC_NORM18, working_point = 'bdt{BDT}'.format(BDT = 0.998), selection = mid18_sel, wspace = wspace_B18, sigma = 0.06, dp = args.discrete_profiling)#, Z_frac = 0.043)
hig18 = Category(norm=1.0, flat=False, blind = not args.unblind, name = 'C18', sig_norm = MC_NORM18, working_point = 'bdt{BDT}'.format(BDT = 0.994), selection = hig18_sel, wspace = wspace_C18, sigma = 0.02, dp = args.discrete_profiling)#, Z_frac = 0.048)
## SYSTEMATICS
##
Z_frac = 0#0.84
W_frac = 1. - Z_frac

common_systematics = [ 
    Systematics(name = 'xs_W'     , distribution = 'lnN', value = str(1.037*W_frac)),
    Systematics(name = 'br_Wtaunu', distribution = 'lnN', value = str(1.018*W_frac)),
    Systematics(name = 'br_Wmunu' , distribution = 'lnN', value = str(1.014*W_frac)),
    Systematics(name = 'WNLO'     , distribution = 'lnN', value = str(1.040*W_frac)),
]
low17.add_systematics(common_systematics+[
    Systematics(name = 'Lumi17'    , distribution = 'lnN', value = '1.023'),
    Systematics(name = 'muonID_A17', distribution = 'lnN', value = '1.013'),
    Systematics(name = 'HLT_Mu_A17', distribution = 'lnN', value = '1.019'),
    Systematics(name = 'HLT_iso17' , distribution = 'lnN', value = '1.12'),
    #Systematics(name = 'frac_Z_A17', distribution = 'lnN', value = str(1.007*Z_frac)),
    Systematics(name = 'HLT_TkMu_A17', distribution = 'lnN', value = '1.11'),
])
mid17.add_systematics(common_systematics+[
    Systematics(name = 'Lumi17'    , distribution = 'lnN', value = '1.023'),
    Systematics(name = 'muonID_B17', distribution = 'lnN', value = '1.014'),
    Systematics(name = 'HLT_Mu_B17', distribution = 'lnN', value = '1.021'),
    Systematics(name = 'HLT_iso17' , distribution = 'lnN', value = '1.12'),
    #Systematics(name = 'frac_Z_B17', distribution = 'lnN', value = str(1.008*Z_frac)),
    Systematics(name = 'HLT_TkMu_B17', distribution = 'lnN', value = '1.1'),
])
hig17.add_systematics(common_systematics+[
    Systematics(name = 'Lumi17'    , distribution = 'lnN', value = '1.023'),
    Systematics(name = 'muonID_C17', distribution = 'lnN', value = '1.015'),
    Systematics(name = 'HLT_Mu_C17', distribution = 'lnN', value = '1.022'),
    Systematics(name = 'HLT_iso17' , distribution = 'lnN', value = '1.12'),
    #Systematics(name = 'frac_Z_C17', distribution = 'lnN', value = str(1.014*Z_frac)),
    Systematics(name = 'HLT_TkMu_C17', distribution = 'lnN', value = '1.15'),
])
low18.add_systematics(common_systematics+[
    Systematics(name = 'Lumi18'    , distribution = 'lnN', value = '1.025'),
    Systematics(name = 'muonID_A18', distribution = 'lnN', value = '1.039'),
    Systematics(name = 'HLT_Mu_A18', distribution = 'lnN', value = '1.01'),
    Systematics(name = 'HLT_iso18' , distribution = 'lnN', value = '1.07'),
    #Systematics(name = 'frac_Z_A18', distribution = 'lnN', value = str(1.01*Z_frac)),
    Systematics(name = 'HLT_TkMu_A18', distribution = 'lnN', value = '1.08'),
])
mid18.add_systematics(common_systematics+[
    Systematics(name = 'Lumi18'    , distribution = 'lnN', value = '1.025'),
    Systematics(name = 'muonID_B18', distribution = 'lnN', value = '1.047'),
    Systematics(name = 'HLT_Mu_B18', distribution = 'lnN', value = '1.01'),
    Systematics(name = 'HLT_iso18' , distribution = 'lnN', value = '1.07'),
    #Systematics(name = 'frac_Z_B18', distribution = 'lnN', value = str(1.012*Z_frac)),
    Systematics(name = 'HLT_TkMu_B18', distribution = 'lnN', value = '1.08'),
])
hig18.add_systematics(common_systematics+[
    Systematics(name = 'Lumi18'    , distribution = 'lnN', value = '1.025'),
    Systematics(name = 'muonID_C18', distribution = 'lnN', value = '1.052'),
    Systematics(name = 'HLT_Mu_C18', distribution = 'lnN', value = '1.01'),
    Systematics(name = 'HLT_iso18' , distribution = 'lnN', value = '1.07'),
    #Systematics(name = 'frac_Z_C18', distribution = 'lnN', value = str(1.014*Z_frac)),
    Systematics(name = 'HLT_TkMu_C18', distribution = 'lnN', value = '1.09'),
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
#cfg.combine_datacards()
#cfg.run_combine(ntoys = 5000, grid = args.quantile, asymptotic = args.asymptotic, cl = args.cl)
