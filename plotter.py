#!/usr/bin/env python
##########################################################################################
######
###### plotter module
######
###### adopted from a similar script for MET performance study
######
######  - containers for data and plots
######  - plotting routines
###### 
###### Use: from plotter import *
######
###### Author: Gena Kukartsev, October 2010
######
##########################################################################################

from __future__ import division

import sys, ROOT
from ROOT import TFile
from ROOT import TTree
from ROOT import TH1F
from ROOT import TH2F
from ROOT import TCanvas
from ROOT import TLegend
from ROOT import SetOwnership

from ROOT import THStack
from ROOT import TLatex
from ROOT import TStyle
from ROOT import TProfile
from ROOT import TBox
from ROOT import TCut

import operator
from array import array

from tdrStyle import *
from ROOT import gROOT
from ROOT import TGraph
from ROOT import TMultiGraph
from ROOT import TMath
from ROOT import TGraphErrors
from ROOT import TGraphAsymmErrors
from ROOT import TF1

import math



setlogscale = 1



def get_scale(scale_name, x):

    if scale_name == 'cteq6l1_fit':
        # 2010
        a = 1.23524 - 0.132727*(x-1000.0)/1000.0
        # 2011 1.1/fb
        b = 1.23312 - 0.154901*(x-1000.0)/1000.0 +0.0516781*(x-1000.0)/1000.0*(x-1000.0)/1000.0
        #print a,b,a-b
        return b

    elif scale_name == 'graviton':
        return 1.6

    elif scale_name == 'rs_acceptance_ee':
        _a =   0.757
        _b = -82.021
        _c =  81.125
        _a2 = 0.739
        _b2 = -14.698
        _c2 = -43.930
        return (_a+_b/(x+_c)) / (_a2+_b2/(x+_c2))

    elif scale_name == 'rs_acceptance_mumu':
        _a =   0.878
        _b = -50.1
        _c =  22.8
        _p0 = 0.901
        _p1 = -830000.0
        _p2 = -0.0000188
        return (_a+_b/(x+_c)) / (_p0+_p1/x/x/x+_p2*x)

    elif scale_name == 'unit':
        return 1.0

    else:
        return None

def get_scale_err(scale_name, x):
    if scale_name == 'cteq6l1_fit':
        #return 0.0960395 + 0.0660697*math.log(x/1000.0)
        return 0.0 + 0.00996035*math.log(x/250.0)
    if scale_name == 'graviton':
        return 0.10
    else:
        return None

# container for files and trees
class Data:
    def __init__(self, lumi = None):
        self.legend = '[Data]:'
        self.lumi = lumi   # /pb  
        self.cat=[]        # master list of keys for all dictionaries
        self.file = {}
        self.x={}          # dict of arrays
        self.y={}
        self.exh={}
        self.exl={}
        self.eyh={}
        self.eyl={}
        # errors for the 2sigma quantile
        self.exh2={}
        self.exl2={}
        self.eyh2={}
        self.eyl2={}
        self.line_color={}
        self.line_style={}
        self.line_width={}
        self.marker_color={}
        self.marker_size={}
        self.marker_style={}
        self.fill_color={}
        self.fill2_color={}
        self.fill_style={}
        self.fill2_style={}
        self.tlegend={}        # None if do not want to be added
        self.legend_index={}   # order in the legend
        self.type={} 
        self.dofit={} 
        self.fit_min={} 
        self.fit_max={}
        self.frame_fill_color = 19

    def add(self, cat, filename, legend, scale = 1.0,
            line_color = 1, line_style = 1, line_width = 1,
            marker_color = 1, marker_size = 1.2, marker_style=8,
            dofit=False, fit_min = 100, fit_max = 3000,
            var_scale = None,
            fill_style = 1002,
            fill_color = None,
            fill2_style = 1002,
            fill2_color = ROOT.kRed,
            extra_scale_name = 'unit',
            legend_index = None):

        if fill_color == None:
            fill_color = line_color

        self.type[cat] = 'observed'

        self.dofit[cat] = dofit
        self.fit_min[cat] = fit_min
        self.fit_max[cat] = fit_max

        self.cat.append(cat)
        self.file = filename
        self.tlegend[cat] = legend
        #self.legend_index[cat] = legend_index
        if legend_index:
            self.legend_index[legend_index] = cat
        else:
            #self.legend_index[len(self.cat)+1] = cat
            _keylist = self.legend_index.keys() # keys are indices
            _keylist.sort()
            print 'keylist:', _keylist
            print 'keylist length:', len(_keylist)
            if len(_keylist)>0:
                _largest_index = _keylist[len(_keylist)-1]
            else:
                _largest_index = -1
            print 'largest_index:', _largest_index
            self.legend_index[_largest_index+1] = cat
        self.line_color[cat] = line_color
        self.line_style[cat] = line_style
        self.line_width[cat] = line_width
        self.marker_color[cat] = marker_color
        self.marker_size[cat] = marker_size
        self.marker_style[cat] = marker_style
        self.fill_style[cat] = fill_style
        self.fill2_style[cat] = fill2_style
        self.fill_color[cat] = fill_color
        self.fill2_color[cat] = fill2_color

        _x = array('d')
        _y = array('d')
        _exh = array('d')
        _exl = array('d')
        _eyh = array('d')
        _eyl = array('d')
        with open(filename, 'r') as _file:
            for line in _file:
                #print self.legend, line

                _buf = line.split()

                if len(_buf) == 0:
                    continue

                if _buf[0][0] == '#':
                    continue

                _mass = float(_buf[0])
                _scale = get_scale(var_scale, _mass)
                # get SSM k-factor
                #if cat=='SSM':
                #    print _mass, _scale,
                _scale_err = get_scale_err(var_scale, _mass)

                if _scale != None:
                    #print self.legend, 'mass =', _mass, 'k-factor =', _scale
                    _scale = _scale*scale
                else:
                    _scale = scale

                _xval = float(_buf[0])
                _x.append(_xval)
                # print SSM theory xsec
                #if cat=='SSM':
                #    print '            ', float(_buf[1])
                _yval = float(_buf[1])*_scale*get_scale(extra_scale_name, _xval)
                _y.append(_yval)

                # now fill errors (if any)
                if _scale_err != None:
                    #print self.legend, 'mass =', _mass, 'k-factor =', _scale
                    _ey = _yval * _scale_err * get_scale(extra_scale_name, _xval)
                else:
                    _ey = 0.0
                _ex = 0
                _eyh.append(_ey)
                _eyl.append(_ey)
                _exh.append(_ex)
                _exl.append(_ex)

        self.x[cat] = _x
        self.y[cat] = _y
        self.exh[cat] = _exh
        self.exl[cat] = _exl
        self.eyh[cat] = _eyh
        self.eyl[cat] = _eyl

        
    def add_function(self, cat, filename, legend, scale = 1.0,
                     line_color = 1, line_style = 1, line_width = 1,
                     marker_color = 1, marker_size = 1.2, marker_style=8,
                     fill_style = 1002,
                     fill_color = None,
                     fill2_style = 1002,
                     fill2_color = ROOT.kRed,
                     index_x = 0,
                     index_y = 1,
                     functor = None,
                     entry = None):

        if fill_color == None:
            fill_color = line_color

        self.type[cat] = 'observed'

        self.dofit[cat] = False
        self.fit_min[cat] = None
        self.fit_max[cat] = None

        self.cat.append(cat)
        self.file = filename
        self.tlegend[cat] = legend
        self.line_color[cat] = line_color
        self.line_style[cat] = line_style
        self.line_width[cat] = line_width
        self.marker_color[cat] = marker_color
        self.marker_size[cat] = marker_size
        self.marker_style[cat] = marker_style
        self.fill_style[cat] = fill_style
        self.fill2_style[cat] = fill2_style
        self.fill_color[cat] = fill_color
        self.fill2_color[cat] = fill2_color

        x_array = array('d')
        y_array = array('d')
        exh_array = array('d')
        exl_array = array('d')
        eyh_array = array('d')
        eyl_array = array('d')

        _nentries = {} # number of rows with a given mass (key:mass)
        _ydict = {} # key: mass, value: list yval

        _ydict_p = {}
        _ydict_m = {}

        with open(filename, 'r') as _file:
            for line in _file:
                #print self.legend, line

                _buf = line.split()

                if len(_buf) == 0:
                    continue

                if _buf[0][0] == '#':
                    continue

                _x_input = float( _buf[index_x].strip() )
                _y_input = float( _buf[index_y].strip() )

                _y_m = float( _buf[3].strip() )
                _y_p = float( _buf[5].strip() )

                _scale = scale

                _xval = _x_input
                
                _yval = _y_input * _scale

                #print 'debug:', _mass
                if _xval in _nentries:
                    _nentries[_xval] += 1
                else:
                    _nentries[_xval] = 0
                    _ydict[_xval] = []
                    _ydict_m[_xval] = []
                    _ydict_p[_xval] = []

                _ydict[_xval].append(_yval)

                _ydict_m[_xval].append(_y_m)
                _ydict_p[_xval].append(_y_p)


        _keys = _nentries.keys()
        _keys.sort()
        #print _keys
        #print _ydict
        for _x in _keys:

            x_array.append(_x)
            
            if functor:
                _y = functor(_x, _ydict)
            elif entry:
                # normalized residual
                #_y = math.fabs(1.0 - _ydict[_x][entry]/_ydict[_x][0])

                #residual as a fraction of expected band
                _y = 2.0*( _ydict[_x][0] - _ydict[_x][entry]) / (_ydict_p[_x][entry] - _ydict_m[_x][entry])
                
                #print _y
            else:
                _y = _ydict[_x][0]
            
            y_array.append(_y)
            eyh_array.append(0)
            eyl_array.append(0)
            exh_array.append(0)
            exl_array.append(0)

        self.x[cat] = x_array
        self.y[cat] = y_array
        self.exh[cat] = exh_array
        self.exl[cat] = exl_array
        self.eyh[cat] = eyh_array
        self.eyl[cat] = eyl_array



    def add_median(self, cat, filename, legend, scale = 1.0,
                   line_color = 1, line_style = 1, line_width = 1,
                   marker_color = 1, marker_size = 1.2, marker_style=8,
                   dofit=False, fit_min = 100, fit_max = 3000,
                   var_scale = None,
                   fill_style = 1002,
                   fill_color = None,
                   fill2_style = 1002,
                   fill2_color = ROOT.kRed,
                   extra_scale_name = 'unit',
                   value_type = 'median',
                   error_type = 'quantile',
                   legend_index = None):
        #
        # same as add() but for multiple values
        #
        

        if fill_color == None:
            fill_color = line_color

        self.type[cat] = 'observed'

        self.dofit[cat] = dofit
        self.fit_min[cat] = fit_min
        self.fit_max[cat] = fit_max

        self.cat.append(cat)
        self.file = filename
        self.tlegend[cat] = legend
        #self.legend_index[cat] = legend_index
        if legend_index:
            self.legend_index[legend_index] = cat
        else:
            self.legend_index[len(self.cat)+1] = cat
        self.line_color[cat] = line_color
        self.line_style[cat] = line_style
        self.line_width[cat] = line_width
        self.marker_color[cat] = marker_color
        self.marker_size[cat] = marker_size
        self.marker_style[cat] = marker_style
        self.fill_style[cat] = fill_style
        self.fill2_style[cat] = fill2_style
        self.fill_color[cat] = fill_color
        self.fill2_color[cat] = fill2_color

        _x = array('d')
        _y = array('d')
        _exh = array('d')
        _exl = array('d')
        _eyh = array('d')
        _eyl = array('d')
        _d = [] # list of tuples, one tuple per mass point
        with open(filename, 'r') as _file:
            _mass = None
            _arr = array('d')
            _prev_limit = 0.0
            for line in _file:
                #print self.legend, line

                _buf = line.split()

                if _buf[0][0] == '#':
                    continue

                #_mass = float(_buf[0])
                aMass = float(_buf[0])
                _scale = get_scale(var_scale, aMass)
                _scale_err = get_scale_err(var_scale, aMass)

                if _scale != None:
                    _scale = _scale*scale
                else:
                    _scale = scale

                #_xval = float(_buf[0])
                #_x.append(_xval)

                #_yval = float(_buf[1])*_scale*get_scale(extra_scale_name, _xval)
                aLim = float(_buf[1])*_scale*get_scale(extra_scale_name, aMass)
                #_y.append(_yval)

                if _mass == None:
                    _mass = aMass
                
                if abs(_mass-aMass) > 0.1:
                    _d.append((_mass, _arr))
                    #print self.legend, cat, 'mass =', _mass, 'PEs per point:', len(_arr)
                    _mass = aMass
                    _arr=array('d')

                    
                #if abs(aLim - _prev_limit) < 0.000001:
                #    #print self.legend, 'WARNING: duplicate expected limit, ignoring...'
                #    continue
                
                #_arr.append(float(_buf[1])*scale)
                _arr.append(aLim)
                _prev_limit = aLim
                
            _d.append((_mass, _arr))
                
        for t in _d:

            _x.append(t[0])
            _mean = TMath.Mean(len(t[1]), t[1])
            _median = TMath.Median(len(t[1]), t[1])
            _rms  = TMath.RMS(len(t[1]), t[1])
            _nrms = _rms/len(t[1])


            # cross check
            #CheckBimodal(t)
 
            #print self.legend, t[0], _mean, _rms, _nrms
            if value_type == 'mean':
                _value = _mean
                _y.append(_mean)
            elif value_type == 'median':
                _value = _median
                _y.append(_median)
            _exh.append(0)
            _exl.append(0)
            if error_type == 'rms':
                _eyh.append(_rms/2)
                _eyl.append(_rms/2)
            elif error_type == 'quantile':
                _prob = array('d')
                _prob.append(0.021)
                _prob.append(0.159)
                _prob.append(0.841)
                _prob.append(0.979)
                _nprob = 4
                _quantiles = array('d')
                _quantiles.append(0)
                _quantiles.append(0)
                _quantiles.append(0)
                _quantiles.append(0)
                TMath.Quantiles(len(t[1]), _nprob, t[1], _quantiles, _prob)
                #print self.legend, 'Quantiles:', _quantiles

                #_eyh.append(abs(_value-_quantiles[2]))
                #_eyl.append(abs(_value-_quantiles[1]))
                _eyh.append(0)
                _eyl.append(0)

            # dumping to text file
            #print t[0], _median

        self.x[cat] = _x
        self.y[cat] = _y
        self.exh[cat] = _exh
        self.exl[cat] = _exl
        self.eyh[cat] = _eyh
        self.eyl[cat] = _eyl

        
    def add_with_errors(self, cat, filename, legend, scale = 1.0,
                        line_color = 1, line_style = 1, line_width = 1,
                        marker_color = 1, marker_size = 1.2, marker_style=8,
                        dofit=False, fit_min = 100, fit_max = 3000,
                        var_scale = None,
                        fill_style = 1002,
                        fill_color = None,
                        fill2_style = 1002,
                        fill2_color = ROOT.kRed,
                        extra_scale_name = 'unit',
                        index_x = 0,
                        index_y = 1,
                        index_ex_up = None,
                        index_ey_up = None,
                        index_ex_down = None,
                        index_ey_down = None,
                        n_entry = 1, # number of entry if multiple entries with the same value of x
                        legend_index = None
                        ):

        if fill_color == None:
            fill_color = line_color

        self.type[cat] = 'observed'

        self.dofit[cat] = dofit
        self.fit_min[cat] = fit_min
        self.fit_max[cat] = fit_max

        self.cat.append(cat)
        self.file = filename
        self.tlegend[cat] = legend
        #self.legend_index[cat] = legend_index
        if legend_index:
            self.legend_index[legend_index] = cat
        else:
            self.legend_index[len(self.cat)+1] = cat
        self.line_color[cat] = line_color
        self.line_style[cat] = line_style
        self.line_width[cat] = line_width
        self.marker_color[cat] = marker_color
        self.marker_size[cat] = marker_size
        self.marker_style[cat] = marker_style
        self.fill_style[cat] = fill_style
        self.fill2_style[cat] = fill2_style
        self.fill_color[cat] = fill_color
        self.fill2_color[cat] = fill2_color

        _x = array('d')
        _y = array('d')
        _exh = array('d')
        _exl = array('d')
        _eyh = array('d')
        _eyl = array('d')

        _nentries = {} # number of rows with a given mass
        
        with open(filename, 'r') as _file:
            for line in _file:
                #print self.legend, line

                _buf = line.split()

                if _buf[0][0] == '#':
                    continue

                _mass = float(_buf[index_x])

                # only add if entry has the right number, .i.e. "third 600"
                #print 'debug:', _mass
                if _mass in _nentries:
                    _nentries[_mass] += 1
                else:
                    _nentries[_mass] = 1

                #print 'debug:', _nentries
                if _nentries[_mass] != n_entry:
                    continue

                _scale = None
                _scale_err = None
                if var_scale:
                    _scale = get_scale(var_scale, _mass)
                    _scale_err = get_scale_err(var_scale, _mass)

                if _scale != None:
                    #print self.legend, 'mass =', _mass, 'k-factor =', _scale
                    _scale = _scale*scale
                else:
                    _scale = scale

                _xval = float(_buf[index_x])
                _x.append(_xval)
                _yval = float(_buf[index_y])*_scale*get_scale(extra_scale_name, _xval)
                _y.append(_yval)

                # now fill errors (if any)
                if index_ex_up:
                    _ex_up = float(_buf[index_ex_up])
                else:
                    _ex_up = 0.0

                if index_ey_up:
                    _ey_up = float(_buf[index_ey_up])-_yval
                else:
                    _ey_up = 0.0
                    
                if index_ex_down:
                    _ex_down = float(_buf[index_ex_down])
                else:
                    _ex_down = 0.0
                    
                if index_ey_down:
                    _ey_down = _yval-float(_buf[index_ey_down])
                else:
                    _ey_down = 0.0

                # additional errors if any
                if _scale_err != None:
                    #print self.legend, 'mass =', _mass, 'k-factor =', _scale
                    _ey_extra = _yval * _scale_err * get_scale(extra_scale_name, _xval)
                    _ey_up = math.sqrt(_ey_extra*_ey_extra + _ey_up*_ey_up)
                    _ey_down = math.sqrt(_ey_extra*_ey_extra + _ey_down*_ey_down)

                _eyh.append(_ey_up)
                _eyl.append(_ey_down)
                _exh.append(_ex_up)
                _exl.append(_ex_down)

        self.x[cat] = _x
        self.y[cat] = _y
        self.exh[cat] = _exh
        self.exl[cat] = _exl
        self.eyh[cat] = _eyh
        self.eyl[cat] = _eyl


    def smooth_band(self, x, y, m, excl,
                    xmin = 300.0, xmax = 1000,
                    mode = 'low'):
        #
        # Smooth out a line defined by a set of points
        # Returns an array of "y" for the same x
        #

        legend = '[smooth]:'
        print legend, 'smoothing...'

        _y = array('d')
        print legend, y

        _shifted = array('d') # actual line for smoothing
        x_ = array('d') # x with excluded points
        _c = 0
        for y_ in y:
            #print 300+100*_c,m[_c],y[_c]
            if _c not in excl:
                x_.append(x[_c])
                if mode=='low':
                    _shifted.append(m[_c] - y[_c])
                else:
                    _shifted.append(m[_c] + y[_c])
            _c += 1
        
        _g = TGraph(len(x_), x_, _shifted)
        f1 = TF1("f1", "[0]+[1]*exp([2]*x+[3])", xmin, xmax)
        _fr = _g.Fit("f1", "MEWS", "", xmin, xmax)
        _fr.Print()

        c = TCanvas()
        _g.SetMarkerStyle(20)
        _g.SetMarkerSize(1.5)
        _g.Draw('AP')
        raw_input('press <enter> to continue...')

        _count = 0
        for _x in x:
            if _x > (xmin-0.1) and _x < (xmax+0.1):
                print legend, y[_count], f1.Eval(_x)
                if mode=='low':
                    _y.append( m[_count] - f1.Eval(_x) )
                else:
                    _y.append( f1.Eval(_x) - m[_count] )
            else:
                _y.append( y[_count] )
            _count += 1

        return _y

        
    def add_expected(self, cat, filename, legend, scale = 1.0,
                     line_color = 1, line_style = 1, line_width = 1,
                     marker_color = 1, marker_size = 1.2, marker_style=8,
                     value_type = 'median',
                     error_type = 'rms',
                     fill_style = 1002,
                     fill_color = ROOT.kYellow,
                     fill2_style = 1002,
                     fill2_color = ROOT.kRed,
                     extra_scale_name = 'unit',
                     smooth = None, # 'single' or 'comb' or None: special tricks for smoothing expected bands
                     legend_index = None):

        self.type[cat] = 'expected'

        #legend_index = None # safety - not implemented right yet

        #self.cat.append(cat)
        self.cat.append(cat)
        #self.cat.append(cat+'1sig')
        #self.cat.append(cat+'2sig')

        self.file = filename

        #self.tlegend[cat] = legend
        self.tlegend[cat] = 'median expected'
        self.tlegend[cat+'1sig'] = '68% expected'
        self.tlegend[cat+'2sig'] = '95% expected'

        #self.legend_index[cat] = legend_index
        if legend_index:
            self.legend_index[legend_index] = cat
            self.legend_index[legend_index+1] = cat+'1sig'
            self.legend_index[legend_index+2] = cat+'2sig'
            #print len(self.legend_index), self.legend_index
        else:
            self.legend_index[len(self.cat)+1] = cat
            self.legend_index[len(self.cat)+2] = cat+'1sig'
            self.legend_index[len(self.cat)+3] = cat+'2sig'
        self.line_color[cat] = line_color
        self.line_style[cat] = line_style
        self.line_width[cat] = line_width
        self.marker_color[cat] = marker_color
        self.marker_size[cat] = marker_size
        self.marker_style[cat] = marker_style
        self.fill_style[cat] = fill_style
        self.fill2_style[cat] = fill2_style
        self.fill_color[cat] = fill_color
        self.fill2_color[cat] = fill2_color

        _x  = array('d')
        _y  = array('d')
        _exh = array('d')
        _exl = array('d')
        _eyh = array('d')
        _eyl = array('d')
        _exh2 = array('d')
        _exl2 = array('d')
        _eyh2 = array('d')
        _eyl2 = array('d')
        _d = [] # list of tuples, one tuple per mass point
        with open(filename, 'r') as _file:
            _mass = None
            _arr = array('d')
            _prev_limit = 0.0
            for line in _file:
                #print self.legend, line

                _buf = line.split()

                aMass = float(_buf[0])
                aLim = float(_buf[1])*scale*get_scale(extra_scale_name, aMass)

                if _mass == None:
                    _mass = aMass
                
                if abs(_mass-aMass) > 0.1:
                    _d.append((_mass, _arr))
                    #print self.legend, cat, 'mass =', _mass, 'PEs per point:', len(_arr)
                    _mass = aMass
                    _arr=array('d')

                    
                #if abs(aLim - _prev_limit) < 0.000001:
                #    #print self.legend, 'WARNING: duplicate expected limit, ignoring...'
                #    continue
                
                #_arr.append(float(_buf[1])*scale)
                _arr.append(aLim)
                _prev_limit = aLim
                
            _d.append((_mass, _arr))

        #print self.legend, _d

        for t in _d:

            _x.append(t[0])
            _mean = TMath.Mean(len(t[1]), t[1])
            _median = TMath.Median(len(t[1]), t[1])
            _rms  = TMath.RMS(len(t[1]), t[1])
            _nrms = _rms/len(t[1])

            # cross check
            #CheckBimodal(t, 'exp')

            #print self.legend, t[0], _mean, _rms, _nrms
            if value_type == 'mean':
                _value = _mean
                _y.append(_mean)
            elif value_type == 'median':
                _value = _median
                _y.append(_median)
            _exh.append(0)
            _exl.append(0)
            _exh2.append(0)
            _exl2.append(0)
            if error_type == 'rms':
                _eyh.append(_rms/2)
                _eyl.append(_rms/2)
                _eyh2.append(_rms)
                _eyl2.append(_rms)
            elif error_type == 'quantile':
                _prob = array('d')
                _prob.append(0.021)
                _prob.append(0.159)
                _prob.append(0.841)
                _prob.append(0.979)
                _nprob = 4
                _quantiles = array('d')
                _quantiles.append(0)
                _quantiles.append(0)
                _quantiles.append(0)
                _quantiles.append(0)
                TMath.Quantiles(len(t[1]), _nprob, t[1], _quantiles, _prob)
                #print self.legend, 'Quantiles:', _quantiles

                _eyh.append(abs(_value-_quantiles[2]))
                _eyl.append(abs(_value-_quantiles[1]))
                _eyh2.append(abs(_value-_quantiles[3]))
                _eyl2.append(abs(_value-_quantiles[0]))

            # dumping to text file
            #print t[0], _median

        self.x[cat] = _x
        self.y[cat] = _y
        self.exh[cat] = _exh
        self.exl[cat] = _exl
        self.exh2[cat] = _exh2
        self.exl2[cat] = _exl2

        # smoothed bands
        if (smooth == 'single'):
            # ee, mumu
            self.eyh[cat] = self.smooth_band(_x,_eyh,_y,[1],300,1500,'high')
            self.eyh2[cat] = self.smooth_band(_x,_eyh2,_y,[1],300,1500,'high')
            self.eyl[cat] = self.smooth_band(_x,_eyl,_y,[],300,1500)
            self.eyl2[cat] = self.smooth_band(_x,_eyl2,_y,[1],300,1500)
        elif (smooth == 'comb'):
            # comb
            self.eyh[cat] = self.smooth_band(_x,_eyh,_y,[],300,600,'high')
            self.eyh2[cat] = self.smooth_band(_x,_eyh2,_y,[1],300,600,'high')
            self.eyl[cat] = self.smooth_band(_x,_eyl,_y,[1],300,600)
            self.eyl2[cat] = self.smooth_band(_x,_eyl2,_y,[1],300,600)
        else:
            # non-smoothed bands
            self.eyh[cat] = _eyh
            self.eyl[cat] = _eyl
            self.eyh2[cat] = _eyh2
            self.eyl2[cat] = _eyl2



    


class PlotGraphs:
    def __init__(self, data,
                 xlow, xhigh, ylow, yhigh,
                 xlabel = "", ylabel = "",
                 xLegend = .45, yLegend = .60,
                 legendWidth = 0.20, legendHeight = 0.45,
                 fillStyle = 3395,
                 drawOption = 'APL3',
                 make_tfile = False):

        self.graph_index = {}
        self.ngraphs = 0

        self.data = data
        self.drawOption = drawOption
        self.xlow=xlow
        self.xhigh=xhigh
        self.ylow=ylow
        self.yhigh=yhigh
        self.xlabel = xlabel
        self.ylabel = ylabel
        self.multigraph=TMultiGraph("MultiGraph", "")
        self.legend_type = {}

        if make_tfile:
            self.tfile = TFile('tgraphs.root', 'recreate')

        self.legend = TLegend(xLegend, yLegend,
                              xLegend + legendWidth, yLegend + legendHeight)

        self.g_ = {}
        self.g2_ = {}
        self.g3_ = {}
        for cat in data.cat:
            if data.type[cat] == 'observed':
                #self.g_[cat] = TGraph(len(data.x[cat]), data.x[cat], data.y[cat])
                self.g_[cat] = TGraphAsymmErrors(len(data.x[cat]), data.x[cat], data.y[cat], data.exl[cat], data.exh[cat], data.eyl[cat], data.eyh[cat])
                self.g_[cat].SetName(cat)
                if make_tfile:
                    _gcopy = self.g_[cat].Clone()
                    self.tfile.Append(_gcopy)
                if data.dofit[cat]:
                    #self.g_[cat].Fit("pol3", "M", "", data.fit_min[cat], data.fit_max[cat])
                    f1 = TF1("f1", "[0]+[1]*(x-1000.0)/1000.0+[2]*(x-1000.0)*(x-1000.0)/1000000.0+[3]*(x-1000.0)*(x-1000.0)*(x-1000.0)/1000000000.0", data.fit_min[cat], data.fit_max[cat])
                    f2 = TF1("f2", "[0]+[1]*(x-1000.0)/1000.0+[2]*(x-1000.0)*(x-1000.0)/1000000.0", data.fit_min[cat], data.fit_max[cat])
                    f3 = TF1("f3", "[0]+[1]*(x-1000.0)/1000.0", data.fit_min[cat], data.fit_max[cat])
                    _fr = self.g_[cat].Fit("f2", "MEWS", "", data.fit_min[cat], data.fit_max[cat])
                    _fr.Print()
                    
            elif data.type[cat] == 'expected':
                self.g_[cat] = TGraphAsymmErrors(len(data.x[cat]), data.x[cat], data.y[cat], data.exl[cat], data.exh[cat], data.eyl[cat], data.eyh[cat])
                self.g2_[cat] = TGraphAsymmErrors(len(data.x[cat]), data.x[cat], data.y[cat], data.exl2[cat], data.exh2[cat], data.eyl2[cat], data.eyh2[cat])

                self.g3_[cat] = TGraphAsymmErrors(len(data.x[cat]), data.x[cat], data.y[cat], data.exl2[cat], data.exh2[cat], data.y[cat], data.exl2[cat])
                self.g3_[cat].SetFillStyle(1002)
                #self.g3_[cat].SetFillStyle(3008)

                # 95% quantile 
                self.g2_[cat].SetMarkerColor(data.fill2_color[cat])
                self.g2_[cat].SetMarkerStyle(data.marker_style[cat])
                self.g2_[cat].SetMarkerSize(data.marker_size[cat])
                self.g2_[cat].SetLineColor(data.fill2_color[cat])
                self.g2_[cat].SetLineStyle(data.line_style[cat])
                self.g2_[cat].SetLineWidth(data.line_width[cat])
                #self.g2_[cat].SetFillColor(data.marker_color[cat]+2)
                self.g2_[cat].SetFillColor(data.fill2_color[cat])
                #self.g2_[cat].SetFillStyle(3008)
                #self.g2_[cat].SetFillStyle(3003)
                self.g2_[cat].SetFillStyle(data.fill2_style[cat])

                
            self.g_[cat].SetMarkerColor(data.marker_color[cat])
            self.g_[cat].SetMarkerStyle(data.marker_style[cat])
            self.g_[cat].SetMarkerSize(data.marker_size[cat])
            self.g_[cat].SetLineColor(data.line_color[cat])
            self.g_[cat].SetLineStyle(data.line_style[cat])
            self.g_[cat].SetLineWidth(data.line_width[cat])
            self.g_[cat].SetFillColor(data.fill_color[cat])
            if data.fill_style[cat] == None:
                self.g_[cat].SetFillStyle(fillStyle)
            else:
                self.g_[cat].SetFillStyle(data.fill_style[cat])

            if data.type[cat] == 'observed':
                # only the main observed limit is a line,
                # everything else is a filled area
                if cat[0:3] == 'obs':
                    _legend_type = 'lp'
                    if 'PC' in drawOption:
                        _draw_option = 'PC'
                    else:
                        _draw_option = 'PL'
                    self.multigraph.Add(self.g_[cat], _draw_option)
                elif cat[0:3] == 'SSM':
                    _legend_type = 'lp'
                    _draw_option = '3'
                    self.multigraph.Add(self.g_[cat], _draw_option)
                    # print theory curve
                    #self.g_[cat].Print()
                elif cat[0:3] == 'Psi':
                    _legend_type = 'lp'
                    _draw_option = '3'
                    self.multigraph.Add(self.g_[cat], _draw_option)
                elif cat[0:3] == 'Stu':
                    _legend_type = 'lp'
                    _draw_option = '3'
                    self.multigraph.Add(self.g_[cat], _draw_option)
                elif cat[0:2] == 'RS':
                    _legend_type = 'f'
                    _draw_option = '3'
                    self.multigraph.Add(self.g_[cat], _draw_option)
                else:
                    _legend_type = 'f'
                    self.multigraph.Add(self.g_[cat])

                #self.multigraph.Add(self.g_[cat])
                self.graph_index[cat] = self.ngraphs
                self.ngraphs += 1

                #self.legend . AddEntry( self.g_[cat], data.tlegend[cat], _legend_type);
                self.legend_type[cat] = _legend_type
                
            elif data.type[cat] == 'expected':
                self.g3_[cat].SetFillColor(0)
                self.g3_[cat].SetMinimum(self.ylow)
                self.g3_[cat].SetMaximum(self.yhigh)
                self.g3_[cat].GetXaxis().SetLimits(self.xlow, self.xhigh)
                self.multigraph.Add(self.g3_[cat], 'C3')
                self.ngraphs += 1

                # 95 expected band
                self.multigraph.Add(self.g2_[cat], 'C4L')
                self.ngraphs += 1

                self.multigraph.Add(self.g_[cat], 'C4L')
                self.graph_index[cat] = self.ngraphs
                self.ngraphs += 1

                # median line
                gline = self.g_[cat].Clone()
                gline.SetLineWidth(3)
                gline.SetLineColor(ROOT.kBlue)
                gline.SetLineStyle(2)
                # median expected marker
                gline.SetMarkerStyle(8)
                gline.SetMarkerSize(1)
                # print graph contents
                #gline.Print()
                #self.multigraph.Add(gline, 'LXC')
                self.multigraph.Add(gline, 'C3X')

                self.ngraphs += 1

                # median expected legend
                #self.legend . AddEntry( gline, 'median expected', "l");
                self.legend_type[cat] = "l"
                _g = self.g_[cat]
                self.g_[cat]=gline

                # 68% expected band legend
                #self.legend . AddEntry( self.g_[cat], '68% expected', "f");
                self.legend_type[cat+'1sig'] = "f"
                self.g_[cat+'1sig']=_g

                # 95% expected band legend
                #self.legend . AddEntry( self.g2_[cat], '95% expected', "f");
                self.legend_type[cat+'2sig'] = "f"
                self.g_[cat+'2sig']=self.g2_[cat]


        keylist = data.legend_index.keys() # keys are indices
        keylist.sort()
        for key in keylist:
            print key, data.legend_index[key]
            self.legend . AddEntry( self.g_[data.legend_index[key]],
                                    data.tlegend[data.legend_index[key]],
                                    self.legend_type[data.legend_index[key]]);

        self.legend . SetShadowColor(0)
        self.legend . SetFillColor(0)
        self.legend . SetLineColor(0)

        if make_tfile:
            self.tfile.Write()
            self.tfile.Close()



    def draw(self, yLabelSize = 0.055):

        self.multigraph.SetMinimum(self.ylow)
        self.multigraph.SetMaximum(self.yhigh)
        self.multigraph . Draw(self.drawOption)

        self.multigraph.GetXaxis().SetNdivisions(405)
        self.multigraph.GetYaxis().SetNdivisions(405)

        self.multigraph.GetXaxis().SetLimits(self.xlow, self.xhigh)
        #self.multigraph.GetYaxis().SetTitle("")

        self.multigraph.GetYaxis().SetLabelSize(yLabelSize)
        self.multigraph.GetXaxis().SetLabelSize(yLabelSize)

        latex = TLatex()
        latex.SetNDC()
        #latex.SetTextSize(0.04)
        latex.SetTextSize(yLabelSize)
        latex.SetTextAlign(31) # align right

        latex.DrawLatex(0.95,0.01, self.xlabel)

        latex.SetTextAngle(90)

        latex.DrawLatex(0.03,0.9, self.ylabel)

        self.legend.SetFillStyle(0)
        self.legend.SetBorderSize(0)
        #self.legend.SetTextSize(0.04)
        #self.legend.SetTextSize(yLabelSize)
        self.legend.SetTextSize(yLabelSize*0.7)
        self.legend.SetTextFont(42)
        self.legend.SetTextAlign(11)
        self.legend . Draw()

        return self.multigraph


    def draw_line(self, xline, ymin = -0.02, ymax = 0.20):

        if xline == None:
            return
        
        #print 'XXXX', xline
        _x = array('d')
        _y = array('d')
        _x.append(xline)
        _y.append(ymin)
        _x.append(float(xline))
        _y.append(ymax)
        g_ = TGraph(len(_x), _x, _y)
        self.multigraph.Add(g_, 'L')



    def print_values(self, cat1, cat2):

        legend = 'PlotGraphs::print_values():'

        if cat1 in self.data.cat:
            graph1 = self.multigraph.GetListOfGraphs().At(self.graph_index[cat1])
        else:
            return None

        if cat2 in self.data.cat:
            graph2 = self.multigraph.GetListOfGraphs().At(self.graph_index[cat2])
        else:
            return None

        for x in range (750, 1150, 50):
            dmax = 0.0
            drelmax = 0.0
            v1 = graph1.Eval(x)
            v2 = graph2.Eval(x)
            d = math.fabs(v2-v1)
            drel = math.fabs((v2-v1)/v1)
            if d>dmax:
                dmax = d
            if drel>drelmax:
                drelmax = drel
            print legend, 'x =', x
            print legend, cat1, 'value =', v1
            print legend, cat2, 'value =', v2
            print legend, 'abs diff =', d
            print legend, 'abs relative diff =', drel
            
        print legend, 'max abs diff =', dmax
        print legend, 'max abs relative diff =', drelmax


    def find_intersection(self, cat1, cat2,
                          xmin = 350, xmax = 3500,
                          precision = 0.000000001):
        legend = 'PlotGraphs::find_intersection(%s,%s):'%(cat1,cat2)

	#always put "exp" or "obs" in cat1 
	if cat1 is "exp" or cat1 is "obs":
	    if cat1 in self.data.cat:
            	graph1 = self.multigraph.GetListOfGraphs().At(self.graph_index[cat1])
            else:
	    	print "Failed to get ", cat1
	 	return ["None"]
	else:
	    print cat1, " is not obs/exp limit."
            return ["None"]
        
        if cat2 in self.data.cat:
            graph2 = self.multigraph.GetListOfGraphs().At(self.graph_index[cat2])
            #graph2.Print()
        else:
            print "Failed to get ", cat2
            return ["None"]

        _x = xmin

        print legend, graph1.Eval(_x) - graph2.Eval(_x)
        
        if graph1.Eval(_x) > graph2.Eval(_x):
            result = ["lowlimit"]
        else:
            result = [] 
            
        _interval = 1
        
        _sign1 = graph1.Eval(_x) > graph2.Eval(_x)
        while _x < xmax:
            _x += _interval
            _sign2 = graph1.Eval(_x) > graph2.Eval(_x)
            if _sign1 != _sign2:
                _d = 100000
                _step = -0.5*_interval
                _crosspt=_x
                
                while abs(_d) > precision:

                    _crosspt += _step
                    _d = graph1.Eval(_crosspt) - graph2.Eval(_crosspt)
                    _sign2 = _d > 0
                    if _sign1 != _sign2:
                        _step = -0.5*abs(_step)
                    else:
                        _step = 0.5*abs(_step)
                _sign1 = not(_sign1)
                result.append(_crosspt)
        if result:
            print_limits(result)
        else:
            print "[",xmin,"~",xmax,"] is excluded."
        return result

def cmsPrel(intLumi = None,
            x=0.12, y=0.93,
            letter = None,
            prel = True):
    latex = TLatex()
    latex.SetNDC()
    latex.SetTextSize(0.045)

    latex.SetTextAlign(31) # align right
    if (intLumi > 0.):
        latex.SetTextAlign(31) # align right
    
    latex.SetTextAlign(5) # align left

    if intLumi:
        if prel:
            #latex.DrawLatex(x,y,"#font[62]{CMS preliminary, #int #font[42]{L}dt = "+'%.1f'%intLumi+"fb^{-1}}")
            latex.DrawLatex(x,y,"#font[62]{CMS preliminary, #font[42]{L}dt = "+'%.1f'%intLumi+"fb^{-1}}")
        else:
            #latex.DrawLatex(x,y,"#font[62]{CMS, #int #font[42]{L}dt = "+'%.1f'%intLumi+"fb^{-1}}")
            latex.DrawLatex(x,y,"#font[62]{CMS, #font[42]{L}dt = "+'%.1f'%intLumi+"fb^{-1}}")
    else:
        if prel:
            latex.DrawLatex(x,y,"#font[62]{CMS preliminary}")
        else:
            latex.DrawLatex(x,y,"#font[62]{CMS}")

    latex.SetTextSize(0.055)

    if letter:
        latex.DrawLatex(0.13,0.15,"#font[62]{"+letter+")}")



def top_legend(text, x=0.892771, y=0.929, letter = None):
    latex = TLatex()
    latex.SetNDC()
    latex.SetTextSize(0.055)

    latex.SetTextAlign(31) # align right
    #latex.SetTextAlign(5) # align left
    latex.DrawLatex(x,y,text)




def legend(legend):
    latex = TLatex()
    latex.SetNDC()
    latex.SetTextSize(0.04)
    latex.SetTextColor(14)

    #latex.SetTextAlign(31) # align right
    latex.SetTextAlign(5) # align left
    latex.DrawLatex(0.5, 0.05, legend)


def print_limits(result):
    print "Allowed region:",
    
    _x=0
    while _x < len(result):
        if _x+1 < len(result):
            print "[",result[_x],",",result[_x+1],"];",
        else:
            if _x+1 == len(result):
                print "[",result[_x],", inf ];"
            else:
                print
        _x+=2
            
def makeMultiGraph(ds, filename,
                   xlow, xhigh, ylow, yhigh,
                   xlabel = "", ylabel = "",
                   xLegend = .65, yLegend = .30,
                   legendWidth = 0.30, legendHeight = 0.35,
                   fillStyle = 3395,
                   drawOption = 'APL3',
                   make_tfile = False,
                   find_intersection = False,
                   draw_limit_line = False,
                   letter = None,
                   print_values = None,
                   toplegend = None,
                   prel = True,
                   plotlegend = None):

#
# make_tfile - if True, will generate a file tgraphs.root with every plot
#              as a separate TGraph object named by its category
#
    _legend = '[makeMultiGraph]:'

    gROOT.SetStyle('Plain')

    c = TCanvas("c", "c", 600, 400)
    c.SetFrameLineWidth(1)
    c.SetRightMargin(0.05)
    c.SetLeftMargin(0.12)
    c.SetBottomMargin(0.12)
    c.cd()

    c.SetFrameFillColor(ds.frame_fill_color)
    #c.SetFrameBorderMode(0)
    c.SetTickx(1)
    c.SetTicky(1)
    
    ROOT.gPad.SetLogy(setlogscale)

    pMass = PlotGraphs(ds,
                       xlow, xhigh, ylow, yhigh,
                       xlabel, ylabel,
                       xLegend, yLegend,
                       legendWidth, legendHeight,
                       fillStyle,
                       drawOption,
                       make_tfile)


    if find_intersection:
        _name='exp'
#        _xssm   = pMass.find_intersection('obs', 'SSM')
#        _xpsi   = pMass.find_intersection('obs', 'Psi')
#        _xrs01  = pMass.find_intersection('obs', 'RS0.1')
#        _xrs005 = pMass.find_intersection('obs', 'RS0.05')
        _stu06 = pMass.find_intersection(_name, 'Stu06')
        _stu05 = pMass.find_intersection(_name, 'Stu05')
        _stu04 = pMass.find_intersection(_name, 'Stu04')
        _stu03 = pMass.find_intersection(_name, 'Stu03')
        _stu02 = pMass.find_intersection(_name, 'Stu02')
#        print _legend, filename, 'SSM:', _xssm
#        print _legend, filename, 'Psi:', _xpsi
#        print _legend, filename, 'RS 0.1:', _xrs01
#        print _legend, filename, 'RS 0.05:', _xrs005

    if draw_limit_line:
        pMass.draw_line(_xssm)
        pMass.draw_line(_xpsi)
        pMass.draw_line(_xrs01)
        pMass.draw_line(_xrs005)

    if print_values:
        pMass.print_values(_name, 'obs2')

    mg = pMass.draw()

    c.Modified()

    cmsPrel(ds.lumi, 0.12, 0.93, letter, prel = prel)

    if toplegend:
        top_legend(toplegend, 0.89, 0.92)

    #frame = c.GetFrame()
    #frame.SetFillStyle(0)
    #frame.Draw()
    _h = mg.GetHistogram()
    _h.SetTitle("")
    _h.GetXaxis().SetNdivisions(507)
    _h.GetYaxis().SetNdivisions(405)
    # for log, diff for diff plots
    _h.SetMaximum(yhigh)
    _h.SetMinimum(ylow)
    _h.Draw("axis same")

    c.SaveAs(filename)

 
    #return c




def makePosteriorPlot(infile,
                      outfile = 'posterior.png',
                      mass = 1000,
                      label_size = 0.05,
                      xlabel = "", ylabel = "",
                      xLegend = .65, yLegend = .30,
                      legendWidth = 0.30, legendHeight = 0.35,
                      fillStyle = 3395,
                      drawOption = 'APL3'):

    _legend = '[makePosteriorPlot()]:'

    gROOT.SetStyle('Plain')

    tf = TFile(infile, 'read')
    c1 = tf.Get('c1')

    _h = c1.GetListOfPrimitives().First()#Get('MCMCposterior_hist')
    _h2 = c1.GetListOfPrimitives().At(1)#Get('MCMCposterior_hist')

    _h.GetXaxis().SetRangeUser(0,0.5)
    _h2.GetXaxis().SetRangeUser(0,0.5)

    _h.SetFillStyle(1001)
    _h.SetFillColor(ROOT.kYellow-9)

    c = TCanvas("c", "c", 600, 400)
    c.SetFrameLineWidth(1)
    c.SetRightMargin(0.05)
    c.cd()

    c.SetFrameFillColor(19)
    #c.SetFrameBorderMode(0)
    c.SetTickx(1)
    c.SetTicky(1)
    
    ROOT.gPad.SetLogy(setlogscale)

    _h.Draw()
    _h2.Draw('same')

    cmsPrel(40.0, 0.89, 0.90)

    latex = TLatex()
    latex.SetNDC()
    #latex.SetTextAlign(31) # align right
    latex.SetTextAlign(5) # align left
    latex.SetTextSize(0.055)
    latex.DrawLatex(0.30, 0.70, "#font[62]{Combined #mu^{+}#mu^{-} and e^{+}e^{-}} channels")
    latex.DrawLatex(0.30, 0.60, "#font[62]{M=}"+'%g'%mass+"GeV")

    #c.Draw()

    # axes labels
    _h.GetXaxis().SetTitle('')
    _h.GetYaxis().SetTitle('')
    _h.GetXaxis().SetNdivisions(407)
    _h.GetYaxis().SetNdivisions(405)
    _h.GetYaxis().SetLabelSize(label_size)
    _h.GetXaxis().SetLabelSize(label_size)
    latex.SetTextSize(label_size)
    latex.SetTextAlign(31) # align right
    latex.DrawLatex(0.93,0.02, xlabel)
    latex.SetTextAngle(90)
    latex.DrawLatex(0.05,0.95, ylabel)


    c.Modified()
    #c.GetListOfPrimitives().Print()
    c.RedrawAxis()

    c.SaveAs(outfile)



# testing 

#data = Data(15.1)
#data.add_expected('comb', 'scan_expected/zpr_exp_nosyst_mcmc_comb_05nov2010v1.txt',
#                  'combined', 1.0)
#line_color = 1, line_style = 1, line_width = 1,
#marker_color = 1, marker_size = 1.2, marker_style=8,
#error_type = 'nrms'):
