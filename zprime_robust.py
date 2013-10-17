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

def robust_plots(plot_format = 'png'):

    data={}
    
    ROOT.gStyle.SetHatchesLineWidth(2)
    ROOT.gStyle.SetLineStyleString(11,"60 30");

    data[1] = Data(200)

    _scale = 1.0
    _scale_theory = 1.0
    
    data[1].add('low', 'dilepton_ratio_low_boundary_scan_27may2011v1.ascii', 'low boundary',
                1.0,
                2, 1, 3,
                2, 1.3, 29)

    data[1].add('high', 'dilepton_ratio_high_boundary_scan_27may2011v1.ascii', 'high boundary',
                1.0,
                4, 1, 3,
                4, 1.3, 29)

    data[1].add('flat_low', 'dilepton_ratio_low_boundary_scan_flatbg_27may2011v2.ascii', 'flat bg low boundary',
                1.0,
                1, 1, 3,
                1, 1.3, 29)

    data[1].add('flat_high', 'dilepton_ratio_high_boundary_scan_flatbg_27may2011v1.ascii', 'flat bg high boundary',
                1.0,
                6, 1, 3,
                6, 1.3, 29)

    makeMultiGraph(data[1], 'robust1.'+plot_format,
                   200, 2000, 0.0, 5.0,
                   xlabel = "dilepton mass, GeV", ylabel = 'N sigma',
                   xLegend = .42, yLegend = .10,
                   legendWidth = 0.45, legendHeight = 0.40,
                   fillStyle = 1002,
                   drawOption = 'APC3', letter = '',
                   toplegend = '#font[42]{#mu^{+}#mu^{-}+ee}')


