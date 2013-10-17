#!/usr/bin/env python

################################################
#
# This script calculates significance with LEE
# and plots the p-value plot
#
#   - input data from ASCII files
#   
# Usage:
#       ./plots.py -p lee -t '2.692381' -d zprime_significance_714ee_746mumu_22-Jun-2011_20-49-11/ -f pdf
#         
# Gena Kukartsev, 2011
#
################################################

import ROOT
import fileinput
import math
import os
import sys
from plotter import *
from ROOT import TMath
from ROOT import TCanvas
from ROOT import TH1F
from ROOT import TLine

import cms_prel as cms

def compare_columns(a, b):
    # sort on ascending index 1, then 0
    return cmp(a[1], b[1]) or cmp(a[0], b[0])


    
def plotSignificance(observed_test_stat, dir = None,
                     plot_format = 'png',
                     infile = None):
    #
    # Plots sampling distribution and test statistic from data
    # Used for the trials factor-corrected significance
    #
    # Input either from a file or from a dir with *.ascii files
    #

    legend = '[significance]:'

    intLumiEe = 4.7;
    intLumiMumu = 5.0;

    if observed_test_stat == None:
        print 'Observed test statistic value not specified (use -t option)'

    _limits =[]

    if dir:
        if infile:
            print legend, 'cannot process both dir and file, specify just one'
            sys.exit(-1)

        for _subdir, _dir, files in os.walk(dir):

            print legend, files

            for filename in files:

                if filename.rfind('ascii') == -1:
                    continue

                with open(dir+'/'+filename) as file:

                    lines=file.readlines()
                
                    for line in lines:
        
                        _words = line.strip().split()
        
                        if _words[0][0] == '#':
                            continue

                        _number = []
                        for word in _words:
                            _number.append(float(word))
            
                        _limits.append(_number)

    if infile:
        with open(infile) as file:

            lines=file.readlines()
                
            for line in lines:
        
                _words = line.strip().split()
                
                if _words[0][0] == '#':
                    continue

                _number = []
                for word in _words:
                    _number.append(float(word))
                    
                _limits.append(_number)

    _limits.sort(compare_columns)

    _h = TH1F("-2 ln L/L_{0}", "", 50, 0, 10)

    _p = observed_test_stat

    _count = 0
    for l in _limits:
        if float(l[1])<_p:
            _count += 1

        _h.Fill(l[1])

    _nup = len(_limits) - _count
    _pval = _nup/float(len(_limits))
    print 'test_stat:', _p
    print 'p-value:', _pval
    _z = math.sqrt(2.0)*TMath.ErfInverse(1.0-2.0*_pval)
    print 'significance:', _z

    c = TCanvas("c", "c", 600, 400)
    c.SetFrameLineWidth(2)
    c.SetRightMargin(0.05)
    c.cd()

    _h.SetLineWidth(2)
    _h.SetLineColor(ROOT.kBlue+2)
    _h.SetFillStyle(1001)
    _h.SetFillColor(ROOT.kBlue-9)
    _h.Draw("")

    _l = TLine(_p,0,_p,_h.GetMaximum() * 0.9)
    _l.SetLineWidth(2)
    _l.SetLineColor(ROOT.kRed)
    _l.Draw("same")
    latex = TLatex()
    latex.SetNDC()
    latex.SetTextSize(0.045)
    latex.SetTextAlign(5) # align left
    latex.DrawLatex(0.45,0.45,"#font[62]{p-value = "+'%.2f'%_pval+"}")
    latex.DrawLatex(0.45,0.4,"#font[62]{Z = "+'%.1f'%_z+"}")

    # x axis
    latex.SetTextAlign(31) # align left
    latex.DrawLatex(0.95,0.01, 'test statistic')

    # y axis
    latex.SetTextAlign(31) # align left
    latex.SetTextAngle(90)
    latex.DrawLatex(0.04,0.95, 'number of pseud-experiments')

    # legend
    xLegend = 0.65
    yLegend = 0.75
    legendWidth = 0.3
    legendHeight = 0.15
    tLegend = TLegend(xLegend, yLegend,
                      xLegend + legendWidth, yLegend + legendHeight)
    tLegend.AddEntry(_h, "pseudo-experiments", 'l')
    tLegend.AddEntry(_l, "observed data", 'l')
    tLegend.Draw()

    # CMS preliminary
    cms.CmsPrel()
    cms.Title('Z-score trials [600-2500] GeV, ee (4.7fb^{-1}) + #mu^{+}#mu^{-} (5.0fb^{-1`})')
    #cms.Title('Z-score trials [600-2500] GeV, ee (4.7fb^{-1})')
    #cms.Title('Z-score trials [600-2500] GeV, #mu^{+}#mu^{-} (5.0fb^{-1`})')
    #cms.Title('Z-score trials 965 GeV, ee (4.7fb^{-1}) + #mu^{+}#mu^{-} (5.0fb^{-1`})')
    #cms.Title('Z-score trials [600-2500] GeV, ee (4.7fb^{-1})')
    #cms.Title('Z-score trials [600-2500] GeV, #mu^{+}#mu^{-} (5.0fb^{-1`})')
    
    c.SaveAs("lee."+plot_format)
    
