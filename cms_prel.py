#!/usr/bin/env python
#########################################################################
#
# cms_prel.py
#
# Write CMS preliminary and other standard disclaimers on plots
#
# Usage:
#        import cms_prel
#        cms_prel.cmsPrel()
#
#
# Author: Gena Kukartsev, December 2011
#
#########################################################################
from __future__ import division

import ROOT


def Title(title = '', x=0.10, y=0.92):
    latex = ROOT.TLatex()
    latex.SetNDC()
    latex.SetTextSize(0.045)

    #latex.SetTextAlign(31) # align right
    latex.SetTextAlign(5) # align left
    latex.DrawLatex(x, y, '#font[62]{'+title+'}')

    return



def XLabel(label = '', x=0.95, y=0.04):
    latex = ROOT.TLatex()
    latex.SetNDC()
    latex.SetTextSize(0.045)

    latex.SetTextAlign(31) # align right
    #latex.SetTextAlign(5) # align left
    latex.DrawLatex(x, y, label)

    return



def YLabel(label = '', x=0.04, y=0.95):
    latex = ROOT.TLatex()
    latex.SetNDC()
    latex.SetTextSize(0.045)

    latex.SetTextAlign(31) # align right
    #latex.SetTextAlign(5) # align left
    latex.SetTextAngle(90)
    latex.DrawLatex(x, y, label)

    return



def CmsPrel(x=0.13, y=0.83):
    latex = ROOT.TLatex()
    latex.SetNDC()
    latex.SetTextSize(0.045)

    #latex.SetTextAlign(31) # align right
    latex.SetTextAlign(5) # align left
    latex.DrawLatex(x, y, "#font[62]{CMS preliminary}")

    return



def Cms(x=0.13, y=0.83):
    latex = ROOT.TLatex()
    latex.SetNDC()
    latex.SetTextSize(0.045)

    #latex.SetTextAlign(31) # align right
    latex.SetTextAlign(5) # align left
    latex.DrawLatex(x, y, "#font[62]{CMS}")

    return


