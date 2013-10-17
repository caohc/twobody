#!/usr/bin/env python

################################################
#
# This script plots limits for the di-lepton
# and di-photon analyses
#
#   - input data from ASCII files
#   
#
# Gena Kukartsev, 2011
#
################################################

import ROOT
from plotter import *

def mass_limits(plot_format = 'png'):
    #
    # Graviton theory curves
    #    

    data={}
    
    ROOT.gStyle.SetHatchesLineWidth(2)
    ROOT.gStyle.SetLineStyleString(11,"60 30");

    data[1] = Data(36.1)

    _scale = 1.0
    _scale_theory = 1.0
    
    data[1].add('RS0.1', 'rsg01_mumu.txt', 'G_{#font[132]{KK}} k/#bar{M_{#font[132]{Pl}}}=0.1',
                1.0e+9*_scale_theory,
                ROOT.kBlue-1, 1, 1,
                ROOT.kBlue-1, 0.1, 8,
                var_scale = 'graviton',
                #fill_style = 3013,
                fill_style = 1002,
                fill_color = ROOT.kBlue+2)
    data[1].add('RS0.05', 'rsg005_mumu.txt', 'G_{#font[132]{KK}} k/#bar{M_{#font[132]{Pl}}}=0.05',
                1.0e+9*_scale_theory,
                ROOT.kGreen-6, 1, 1,
                ROOT.kGreen -6, 0.1, 8,
                var_scale = 'graviton',
                fill_style = 1002)
    data[1].add('obs', 'gkk_mass_comb.ascii', '95% C.L. limit',
                _scale,
                1, 1, 3,
                1, 1.3, 29,
                )

    makeMultiGraph(data[1], 'gkk_obs.'+plot_format,
                   300, 1450, 0, 0.5,
                   xlabel = "M [GeV]", ylabel = '#sigma [pb]',
                   xLegend = .55, yLegend = .40,
                   legendWidth = 0.45, legendHeight = 0.50,
                   fillStyle = 1002,
                   drawOption = 'APC3', letter = '',
                   toplegend = '#font[42]{#gamma#gamma+#mu^{+}#mu^{-}+ee}',
                   find_intersection = True)


