#!/usr/bin/env python

################################################
#
# This script for making limit plots and such
#
#   - plot limits from ASCII file
#   
# Usage:
#        LEE plot: ./plots.py -p lee -t '2.692381' -i zprime_significance_dilepton.ascii
#        LEE plot: ./plots.py -p lee -t '2.692381' -d zprime_significance_714ee_746mumu_22-Jun-2011_20-49-11
#
#        Limit plots: ./plots.py --plot-number 1 -f pdf
#        Limit plots: ./plots.py --plot-number 2 -f pdf
#        Limit plots: ./plots.py --plot-number 3 -f pdf
#        Limit plots (old but working): ./plots.py -p zprime_limits -f pdf
#
#        k-factor fit: ./plots.py -p fit
#
# Gena Kukartsev, 2010
#
#           April 2011: modified for di-photon analysis
#                       and combination with di-lepton
#
################################################

from ROOT import TGraph
from ROOT import TMultiGraph
from ROOT import TCanvas
from ROOT import TLatex
from array import array
from ROOT import RooWorkspace
from ROOT import RooArgSet
from ROOT import RooFit
from ROOT import gROOT
from ROOT import TFile
from ROOT import TGraphAsymmErrors
from ROOT import TF1
from ROOT import RooNumIntConfig

from plotter import *

import ROOT
import math,random
from ROOT import TMath
from ROOT import TH1F
from ROOT import RooDataHist
from ROOT import RooArgList
from ROOT import RooRealVar
from ROOT import RooVoigtian

from operator import mod
import os

_legend = '[plots:]'

#
# command line parser
#
from optparse import OptionParser
add_help_option = "./metZplots.py -p PLOT_SET [other options]"

parser = OptionParser(add_help_option)

parser.add_option("-p", "--plot-set", dest="plot_set", default=None,
                  help="Which plot(s) to make", metavar="PLOTSET")

parser.add_option("-f", "--plot-format", dest="plot_format", default='png',
                  help="Plot file format", metavar="PLOTFORMAT")

parser.add_option("-l", "--lumi", dest="lumi", default=0.8,
                  help="Integrated luminosity", metavar="LUMI")

parser.add_option("--legend", dest="legend", default='channel',
                  help="Legend for the plot", metavar="LUMI")

parser.add_option("-d", "--dir", dest="dir", default=None,
                  help="Input directory", metavar="DIR")

parser.add_option("-i", "--input-file", dest="input_file", default=None,
                  help="Input file")

parser.add_option("-t", "--test-statistic", dest="teststat", default=None,
                  help="Value of the test statistic for look-elsewhere effect")

parser.add_option("--log-x", dest="logx",
                  default=False,
                  action="store_true",
                  help="Set log X axis scale")

parser.add_option("--log-y", dest="logy",
                  default=False,
                  action="store_true",
                  help="Set log Y axis scale")

# plots by number
parser.add_option("--plot-number", dest="plot_number", default=None,
                  action="append",
                  help="Make plot with given number as specified in tprime_plots.py")

print _legend, 'parsing command line options...',
(options, args) = parser.parse_args()
print 'done'

format = options.plot_format
lumi   = options.lumi

do_limits      = False
do_mass_limits = False
do_lee         = False
do_robust      = False
do_fit         = False
do_cls         = False
do_plot_ascii  = False

#tprime plots
do_tprime_limits = False

# define which plots to make based on params
if options.plot_set == 'zprime_limits':
    do_limits = True

if options.plot_set == 'mass':
    do_mass_limits = True

if options.plot_set == 'lee':
    do_lee = True

if options.plot_set == 'robust':
    do_robust = True

if options.plot_set == 'tprime_limits':
    do_tprime_limits = True

if options.plot_set == 'fit':
    do_fit = True

if options.plot_set == 'cls':
    do_cls = True

if options.plot_set == 'plot_ascii':
    do_plot_ascii = True



# end of command line parser

from limit_plots import * 
from mass_limits import * 
from look_elsewhere import * 
from zprime_robust import * 


############################################################
#
# Main
#
#

gROOT.SetStyle('Plain')

if do_limits:
    limit_plots(format)

if do_mass_limits:
    mass_limits(format)

if do_lee:
    ROOT.gStyle.SetOptStat(0)
    plotSignificance(float(options.teststat), options.dir, format,
                     options.input_file)

if do_robust:
    robust_plots(format)

if do_tprime_limits:
    from tprime_limits import *
    #tprime_limits(format)
    #tprime_limits_mcmc(format)
    tprime_cls_limits(format)

if do_cls:
    # - takes condor output dir as input
    # - read input root files from condor dir
    # - merge the results for each mass point and get CLs numbers
    # - make the xsec limit plot as function of mass
    #
    # Usage:
    #        ./plots.py -p cls -d <condor_dir_name> -l <lumi> --legend <#mu+jets> -f format
    #
    from cls_limits import *
    cls_limits(options.dir,format)
    tprime_cls_limits(options.dir,options.lumi, options.legend,format)

if do_fit:
    fit(format)

if do_plot_ascii:
    from plot_ascii import *
    plotAscii(options.input_file, 1, format)

if options.plot_number:
    #
    # Plots by number (from zprime_plots.py)
    #
    # Usage:
    #        ./plots.py --plot-number 2 --plot-number 5
    #
    import zprime_plots as zp
    zp.zprime_plots(format, options.plot_number,
                    options=options)

