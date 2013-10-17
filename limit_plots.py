#!/usr/bin/env python

################################################
#
# This script imports to plots.py
# and plots limits for the di-lepton
# and di-photon analyses, and a few supplementary plots
#
#   - input data from ASCII files
#
# Usage:
#        ./plots.py -p zplrime_limits -f pdf
#
#
# Gena Kukartsev, 2011
#
################################################

import ROOT
from plotter import *

def fit(plot_format = 'png'):
    #
    # Fit k-factor
    #    

    data={}
    
    ROOT.gStyle.SetHatchesLineWidth(2)
    ROOT.gStyle.SetLineStyleString(11,"60 30");

    data[1] = Data(9.9)

    _scale = 1.0
    _scale_theory = 1.0
    
    data[1].add('obs', 'theory/belyaev_tab_lhc_d_cteq_20dec2010v2.txt', 'CTEQ6L k-factor',
                _scale,
                1, 1, 3,
                1, 1.3, 29,
                dofit=True, fit_min = 700, fit_max = 2500)

    makeMultiGraph(data[1], 'kfactor_fit.'+plot_format,
                   100, 3000, 1.0, 1.5,
                   xlabel = "M [GeV]", ylabel = 'k-factor',
                   xLegend = .55, yLegend = .40,
                   legendWidth = 0.45, legendHeight = 0.50,
                   fillStyle = 1002,
                   drawOption = 'APC3', letter = '',
                   toplegend = '#font[42]{CTEQ6L k-factor}')






def theory_plots(plot_format = 'png'):
    #
    # Graviton theory curves
    #    

    data={}
    
    ROOT.gStyle.SetHatchesLineWidth(2)
    ROOT.gStyle.SetLineStyleString(11,"60 30");

    data[1] = Data(36.1)
#
#    _scale = 1.0
#    _scale_theory = 1.0
#    
#    data[1].add_expected('exp', 'expected_tmp.ascii', 'Expected 95% C.L. limit',
#                         _scale,
#                         ROOT.kRed-3, 0, 0,
#                         ROOT.kRed-3, 1.3, 0,
#                         value_type = 'median',
#                         error_type = 'quantile',
#                         fill_style = 1002,
#                         extra_scale_name = 'unit',
#                         fill_color = ROOT.kYellow-7,
#                         fill2_color = ROOT.kRed-3)
#    data[1].add('RS0.1', 'rsg01_mumu.txt', 'G_{#font[132]{KK}} k/#bar{M_{#font[132]{Pl}}}=0.1',
#                1.0e+9*_scale_theory,
#                ROOT.kBlue-1, 1, 1,
#                ROOT.kBlue-1, 0.1, 8,
#                var_scale = 'graviton',
#                #fill_style = 3013,
#                fill_style = 1002,
#                fill_color = ROOT.kBlue+2)
#    data[1].add('RS0.05', 'rsg005_mumu.txt', 'G_{#font[132]{KK}} k/#bar{M_{#font[132]{Pl}}}=0.05',
#                1.0e+9*_scale_theory,
#                ROOT.kGreen-6, 1, 1,
#                ROOT.kGreen -6, 0.1, 8,
#                var_scale = 'graviton',
#                fill_style = 1002)
#    data[1].add('obs', 'observed_tmp.ascii', '95% C.L. limit',
#                _scale,
#                1, 1, 3,
#                1, 1.3, 29)
#
#    makeMultiGraph(data[1], 'gkk_obs.'+plot_format,
#                   300, 1450, 0, 0.5,
#                   xlabel = "M [GeV]", ylabel = '#sigma [pb]',
#                   xLegend = .55, yLegend = .40,
#                   legendWidth = 0.45, legendHeight = 0.50,
#                   fillStyle = 1002,
#                   drawOption = 'APC3', letter = '',
#                   toplegend = '#font[42]{#gamma#gamma+#mu^{+}#mu^{-}+ee}')
#
#




def limit_plots(plot_format = 'png'):
    #
    # Dielectron limits
    #    

    data={}
    
    ROOT.gStyle.SetHatchesLineWidth(2)
    ROOT.gStyle.SetLineStyleString(11,"60 30");

    data[2] = Data(20.6)

    _scale = 0.001
    #_scale = 10.0*1.02
    _scale_theory = 1.0/1120.0
    
    #data[2].add_expected('exp', 'dielectron_expected_714pb.ascii', 'Expected 95% C.L. limit',
    #data[2].add_expected('exp', 'dielectron_mcmc_1100ipb_exp_10jul2011.ascii', 'Expected 95% C.L. limit',
    #data[2].add_expected('exp', 'dielectron_mcmc_expected_1100ipb_16jul2011.ascii', 'Expected 95% C.L. limit',
    #data[2].add_expected('exp', 'update_30sep2011/dielectron_expected_v1.ascii', 'Expected 95% C.L. limit',
    #data[2].add_expected('exp', 'update_22oct2011/dielectron_expected.ascii', 'Expected 95% C.L. limit',
    #
    #data[2].add_expected('exp', 'dielectron_expected_4600ipb_21dec2011v1.ascii', 'Expected 95% C.L. limit',
    #data[2].add_expected('exp', 'petyt_4700/expected_4700_scaled.ascii', 'Expected 95% C.L. limit',
    #data[2].add_expected('exp', 'dimuon_ratio_expected_20invfb.ascii', 'Expected 95% C.L. limit',
    data[2].add_expected('exp', 'dimuon_ratio_expected_20fb_GaussianSyst.ascii', 'Expected 95% C.L. limit',
                         _scale,
                         ROOT.kGreen, 0, 0,
                         ROOT.kYellow, 1.3, 0,
                         value_type = 'median',
                         error_type = 'quantile',
                         fill_style = 1002,
                         extra_scale_name = 'unit',
                         fill_color = ROOT.kGreen,
                         fill2_color = ROOT.kYellow,
                         smooth = None)
    data[2].add('Stu06', 'theory_8TeV/zprime_stu_ep060_cteq6l1_masswindow.txt', 'Z\'_{Stu} #epsilon = 0.06',
                1.0e+9*_scale_theory,
                ROOT.kMagenta+3, 1, 3,
                ROOT.kMagenta+3, 0.1, 8,
                var_scale = 'cteq6l1_fit')
    data[2].add('Stu05', 'theory_8TeV/zprime_stu_ep050_cteq6l1_masswindow.txt', 'Z\'_{Stu} #epsilon = 0.05',
                1.0e+9*_scale_theory,
                ROOT.kMagenta, 1, 3,
                ROOT.kMagenta, 0.1, 8,
                var_scale = 'cteq6l1_fit')
    data[2].add('Stu04', 'theory_8TeV/zprime_stu_ep040_cteq6l1_masswindow.txt', 'Z\'_{Stu} #epsilon = 0.04',
                1.0e+9*_scale_theory,
                ROOT.kMagenta-3, 1, 3,
                ROOT.kMagenta-3, 0.1, 8,
                var_scale = 'cteq6l1_fit')
    data[2].add('Stu03', 'theory_8TeV/zprime_stu_ep030_cteq6l1_masswindow.txt', 'Z\'_{Stu} #epsilon = 0.03',
                1.0e+9*_scale_theory,
                ROOT.kMagenta-7, 1, 3,
                ROOT.kMagenta-7, 0.1, 8,
                var_scale = 'cteq6l1_fit')
    data[2].add('Stu02', 'theory_8TeV/zprime_stu_ep020_cteq6l1_masswindow.txt', 'Z\'_{Stu} #epsilon = 0.02',
                1.0e+9*_scale_theory,
                ROOT.kMagenta-10, 1, 3,
                ROOT.kMagenta-10, 0.1, 8,
                var_scale = 'cteq6l1_fit')
    #data[2].add_median('obs', 'update_22oct2011/dielectron_mass_graviton.ascii', '95% C.L. limit',
    #data[2].add_median('obs', 'update_22oct2011/dielectron_mass_cteq.ascii', '95% C.L. limit',
    #data[2].add_median('obs', 'update_22oct2011/dielectron_observed.ascii', '95% C.L. limit',
    #data[2].add_median('obs', 'update_30sep2011/dielectron_observed_v1.ascii', '95% C.L. limit',
    #data[2].add_median('obs', 'dielectron_mcmc_mass_graviton_1100ipb_16jul2011.ascii', '95% C.L. limit',
    #data[2].add_median('obs', 'dielectron_mcmc_mass_cteq_1100ipb_16jul2011.ascii', '95% C.L. limit',
    #data[2].add_median('obs', 'dielectron_mcmc_observed_1100ipb_16jul2011.ascii', '95% C.L. limit',
    #data[2].add('obs', 'dielectron_mass_cteq_1100ipb.ascii', '95% C.L. limit',
    #data[2].add('obs', 'dielectron_mass_graviton_1100ipb.ascii', '95% C.L. limit',
    #data[2].add('obs', 'dielectron_mcmc_1100ipb_10jul2011.ascii', '95% C.L. limit',
    #data[2].add('obs', 'dielectron_mcmc_1100ipb_obs_10jul2011.ascii', '95% C.L. limit',
    #data[2].add('obs', 'dielectron_observed_714pb.ascii', '95% C.L. limit',
    #data[2].add('obs', 'dielectron_cteq_714ipb.ascii', '95% C.L. limit',
    #data[2].add('obs', 'dielectron_graviton_714ipb.ascii', '95% C.L. limit',
    #
    #data[2].add_median('obs', 'dielectron_observed_4600ipb_21dec2011v1.ascii', '95% C.L. limit',
    #data[2].add_median('obs', 'petyt_4700/observed_4700_scaled.ascii', '95% C.L. limit',
    #data[2].add_median('obs', 'dimuon_observed_20637ipb.ascii', '95% C.L. limit',
    #            _scale,
    #            1, 1, 3,
    #            1, 0.1, 29)
    makeMultiGraph(data[2], 'zprime_dielectron_20ifb.'+plot_format,
                   300, 3500, 4e-8, 0.0001,
                   #xlabel = "M [GeV]", ylabel = 'R_{#sigma} #upoint 10^{-4}',
                   xlabel = "M [GeV]", ylabel = 'R_{#sigma}',
                   xLegend = .59, yLegend = .40,
                   legendWidth = 0.45, legendHeight = 0.50,
                   fillStyle = 1002,
                   drawOption = 'APC3', letter = '',
                   toplegend = '#font[42]{ee}',
                   #find_intersection = True,
                   #draw_limit_line = True
                   )


    #########################################################
    #
    # Dimuon limits
    #    

    ROOT.gStyle.SetHatchesLineWidth(2)
    ROOT.gStyle.SetLineStyleString(11,"60 30");

    data[3] = Data(20.637)

    #data[3].add_expected('exp', 'dimuon_expected_746ipb.ascii', 'Expected 95% C.L. limit',
    #data[3].add_expected('exp', 'dimuon_expected_1125ipb.ascii', 'Expected 95% C.L. limit',
    #data[3].add_expected('exp', 'dimuon_mcmc_1100ipb_exp_10jul2011.ascii', 'Expected 95% C.L. limit',
    #data[3].add_expected('exp', 'dimuon_mcmc_expected_1100ipb_17jul2011.ascii', 'Expected 95% C.L. limit',
    #data[3].add_expected('exp', 'dimuon_ratio_expected_14TeV_20invfb.ascii', 'Expected 95% C.L. limit',
    #data[3].add_expected('exp', 'dimuon_ratio_expected_20invfb.ascii', 'Expected 95% C.L. limit',
    data[3].add_expected('exp', 'dimuon_ratio_expected_20fb_GaussianSyst.ascii', 'Expected 95% C.L. limit',
                         _scale,
                         ROOT.kGreen, 0, 0,
                         ROOT.kGreen, 1.3, 0,
                         value_type = 'median',
                         error_type = 'quantile',
                         fill_style = 1002,
                         extra_scale_name = 'unit',
                         fill_color = ROOT.kGreen,
                         fill2_color = ROOT.kYellow,
                         smooth = None)
    data[3].add('Stu06', 'theory_8TeV/zprime_stu_ep060_cteq6l1_masswindow.txt', 'Z\'_{Stu} #epsilon = 0.06',
                1.0e+9*_scale_theory,
                ROOT.kMagenta+3, 1, 3,
                ROOT.kMagenta+3, 0.1, 8,
                var_scale = 'cteq6l1_fit')
    data[3].add('Stu05', 'theory_8TeV/zprime_stu_ep050_cteq6l1_masswindow.txt', 'Z\'_{Stu} #epsilon = 0.05',
                1.0e+9*_scale_theory,
                ROOT.kMagenta, 1, 3,
                ROOT.kMagenta, 0.1, 8,
                var_scale = 'cteq6l1_fit')
    data[3].add('Stu04', 'theory_8TeV/zprime_stu_ep040_cteq6l1_masswindow.txt', 'Z\'_{Stu} #epsilon = 0.04',
                1.0e+9*_scale_theory,
                ROOT.kMagenta-3, 1, 3,
                ROOT.kMagenta-3, 0.1, 8,
                var_scale = 'cteq6l1_fit')
    data[3].add('Stu03', 'theory_8TeV/zprime_stu_ep030_cteq6l1_masswindow.txt', 'Z\'_{Stu} #epsilon = 0.03',
                1.0e+9*_scale_theory,
                ROOT.kMagenta-7, 1, 3,
                ROOT.kMagenta-7, 0.1, 8,
                var_scale = 'cteq6l1_fit')
    data[3].add('Stu02', 'theory_8TeV/zprime_stu_ep020_cteq6l1_masswindow.txt', 'Z\'_{Stu} #epsilon = 0.02',
                1.0e+9*_scale_theory,
                ROOT.kMagenta-10, 1, 3,
                ROOT.kMagenta-10, 0.1, 8,
                var_scale = 'cteq6l1_fit')
    #data[3].add_median('obs', 'update_22oct2011/dimuon_mass_graviton.ascii', '95% C.L. limit',
    #data[3].add_median('obs', 'update_22oct2011/dimuon_mass_cteq.ascii', '95% C.L. limit',
    #data[3].add_median('obs', 'update_22oct2011/dimuon_observed.ascii', '95% C.L. limit',
    #data[3].add_median('obs', 'update_30sep2011/dimuon_observed_v1.ascii', '95% C.L. limit',
    #data[3].add_median('obs', 'dimuon_mcmc_mass_graviton_1100ipb_16jul2011.ascii', '95% C.L. limit',
    #data[3].add_median('obs', 'dimuon_mcmc_mass_cteq_1100ipb_16jul2011.ascii', '95% C.L. limit',
    #data[3].add_median('obs', 'dimuon_mcmc_observed_1100ipb_16jul2011.ascii', '95% C.L. limit',
    #data[3].add('obs', 'dimuon_mass_cteq_1100ipb.ascii', '95% C.L. limit',
    #data[3].add('obs', 'dimuon_mass_graviton_1100ipb.ascii', '95% C.L. limit',
    #data[3].add('obs', 'dimuon_mcmc_1100ipb_10jul2011.ascii', '95% C.L. limit',
    #data[3].add('obs', 'dimuon_observed_1125ipb.ascii', '95% C.L. limit',
    #data[3].add('obs', 'dimuon_observed_746ipb.ascii', '95% C.L. limit',
    #data[3].add('obs', 'dimuon_cteq_746ipb.ascii', '95% C.L. limit',
    #data[3].add('obs', 'dimuon_graviton_746ipb.ascii', '95% C.L. limit',
#    data[3].add_median('exp_med', 'dimuon_ratio_expected_14TeV_20invfb.ascii', 'Expected 95% C.L. limit',
#   data[3].add_median('obs', 'dimuon_observed_20637ipb.ascii', '95% C.L. observed limit',
#                _scale,
#                1, 1, 3,
#                1, 0.1, 29)
    makeMultiGraph(data[3], 'zprime_dimuon_20ifb.'+plot_format,
                   #300, 2500, 0.005, 0.6,
                   200, 3500, 4e-8, 0.001,
                   #xlabel = "M [GeV]", ylabel = 'R_{#sigma} #upoint 10^{-4}',
                   xlabel = "M [GeV]", ylabel = 'R_{#sigma}',
                   xLegend = .59, yLegend = .40,
                   legendWidth = 0.45, legendHeight = 0.50,
                   fillStyle = 1002,
                   drawOption = 'APC3', letter = '',
                   toplegend = '#font[42]{#mu^{+}#mu^{-}}',
                   find_intersection = True,
                   #draw_limit_line = True
                   )


    ##########################################################################
    #
    # Dilepton limits
    #    

    #data={}
    """
    ROOT.gStyle.SetHatchesLineWidth(2)
    ROOT.gStyle.SetLineStyleString(11,"60 30");

    data[4] = Data() # no lumi

    #data[4].add_expected('exp', 'dilepton_expected_mu746ipb_e714ipb.ascii', 'Expected 95% C.L. limit',
    #data[4].add_expected('exp', 'dilepton_mcmc_1100ipb_exp_10jul2011.ascii', 'Expected 95% C.L. limit',
    #data[4].add_expected('exp', 'dilepton_mcmc_expected_1100ipb_16jul2011.ascii', 'Expected 95% C.L. limit',
    #data[4].add_expected('exp', 'update_30sep2011/dilepton_expected_v1.ascii', 'Expected 95% C.L. limit',
    #data[4].add_expected('exp', 'update_22oct2011/dilepton_expected.ascii', 'Expected 95% C.L. limit',
    data[4].add_expected('exp', 'dilepton_exp_16jan2012v3.txt', 'Expected 95% C.L. limit',
                         _scale,
                         ROOT.kGreen, 0, 0,
                         ROOT.kYellow, 1.3, 0,
                         value_type = 'median',
                         error_type = 'quantile',
                         fill_style = 1002,
                         extra_scale_name = 'unit',
                         fill_color = ROOT.kGreen,
                         fill2_color = ROOT.kYellow,
                         smooth = None)
    #data[4].add('SSM', 'theory/cs_ssm_belyaev_15jan2011v1.dat', 'Z\'_{SSM}',
    data[4].add('SSM', 'theory/zprime_ssm_cteq6l1_masswindow.txt', 'Z\'_{SSM}',
                1.0e+9*_scale_theory,
                ROOT.kGreen-6, 1, 3,
                ROOT.kGreen -6, 0.1, 8,
                var_scale = 'cteq6l1_fit')
    #data[4].add('Psi', 'theory/cs_psi_belyaev_15jan2011v1.dat', 'Z\'_{#Psi}',
    data[4].add('Psi', 'theory/zprime_psi_cteq6l1_masswindow.txt', 'Z\'_{#Psi}',
                1.0e+9*_scale_theory,
                ROOT.kBlue, 1, 3,
                ROOT.kBlue, 0.1, 8,
                var_scale = 'cteq6l1_fit')
    #data[4].add('RS0.1', 'rsg01_mumu.txt', 'G_{#font[132]{KK}} k/#bar{M_{#font[132]{Pl}}}=0.1',
    data[4].add('RS0.1', 'theory/rsg01_cteq6l1_masswindow.txt', 'G_{#font[132]{KK}} k/#bar{M_{#font[132]{Pl}}}=0.1',
                1.0e+9*_scale_theory,
                ROOT.kBlue-1, 1, 1,
                ROOT.kBlue-1, 0.1, 8,
                var_scale = 'graviton',
                #fill_style = 3013,
                fill_style = 1002,
                fill_color = ROOT.kBlue+2)
    #data[4].add('RS0.05', 'rsg005_mumu.txt', 'G_{#font[132]{KK}} k/#bar{M_{#font[132]{Pl}}}=0.05',
    data[4].add('RS0.05', 'theory/rsg005_cteq6l1_masswindow.txt', 'G_{#font[132]{KK}} k/#bar{M_{#font[132]{Pl}}}=0.05',
                1.0e+9*_scale_theory,
                ROOT.kGreen-6, 1, 1,
                ROOT.kGreen -6, 0.1, 8,
                var_scale = 'graviton',
                fill_style = 1002)
    #data[4].add_median('obs', 'update_22oct2011/dilepton_mass_graviton.ascii', '95% C.L. limit',
    #data[4].add_median('obs', 'update_22oct2011/dilepton_mass_cteq.ascii', '95% C.L. limit',
    #data[4].add_median('obs', 'update_22oct2011/dilepton_observed.ascii', '95% C.L. limit',
    #data[4].add_median('obs', 'update_30sep2011/dilepton_observed_v1.ascii', '95% C.L. limit',
    #data[4].add_median('obs', 'dilepton_mcmc_mass_graviton_1100ipb_16jul2011.ascii', '95% C.L. limit',
    #data[4].add_median('obs', 'dilepton_mcmc_mass_cteq_1100ipb_16jul2011.ascii', '95% C.L. limit',
    #data[4].add_median('obs', 'dilepton_mcmc_observed_1100ipb_16jul2011.ascii', '95% C.L. limit',
    #data[4].add('obs', 'dilepton_mass_cteq_1100ipb.ascii', '95% C.L. limit',
    #data[4].add('obs', 'dilepton_mass_graviton_1100ipb.ascii', '95% C.L. limit',
    #data[4].add('obs', 'dilepton_mcmc_1100ipb_10jul2011.ascii', '95% C.L. limit',
    #data[4].add('obs', 'dilepton_mcmc_1100ipb_obs_10jul2011.ascii', '95% C.L. limit',
    #data[4].add('obs', 'dilepton_observed_mu746ipb_e714ipb.ascii', '95% C.L. limit',
    #data[4].add('obs', 'dilepton_cteq_mu746ipb_e714ipb.ascii', '95% C.L. limit',
    #data[4].add('obs', 'dilepton_graviton_mu746ipb_e714ipb.ascii', '95% C.L. limit',
    data[4].add_median('obs', 'dilepton_obs_15jan2012v1.txt', '95% C.L. limit',
                       _scale,
                       1, 1, 3,
                       1, 0.1, 29)

    makeMultiGraph(data[4], 'zprime_dilepton_5ifb.'+plot_format,
                   300, 2500, 0.0000003, 0.00006,
                   #xlabel = "M [GeV]", ylabel = 'R_{#sigma} #upoint 10^{-4}',
                   xlabel = "M [GeV]", ylabel = 'R_{#sigma}',
                   xLegend = .59, yLegend = .40,
                   legendWidth = 0.45, legendHeight = 0.50,
                   fillStyle = 1002,
                   drawOption = 'APC3', letter = '',
                   toplegend = '#font[42]{ee(4.7fb^{-1})+#mu^{+}#mu^{-}(5.0fb^{-1})}',
                   #find_intersection = True,
                   #draw_limit_line = True
                   )

    """
