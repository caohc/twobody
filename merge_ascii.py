#!/usr/bin/env python

################################################
#
# This script merges ASCII output of zprime
# batch jobs into one ASCII file, sorted
# by mass (first column) in the ascending order
#
#
# Gena Kukartsev, November 2010
#
################################################

import fileinput
import sys

def compare_columns(a, b):
    # sort on ascending index 0, then 1
    return cmp(a[0], b[0]) or cmp(a[1], b[1])
    #return cmp(b[0], a[0]) or cmp(a[1], b[1])

def compare_columns_second(a, b):
    # sort on ascending index 0, then 1
    return cmp(a[1], b[1]) or cmp(a[0], b[0])

def fix_lines():
    #print legend, 'about to fix missing newlines'
    for line in fileinput.input():
        _words = line.strip().split()
        
        if _words[0][0] == '#':
            continue

        _i = 0
        for w in _words:
            print w,
            _i +=1
            if _i == 6:
                print
                #print 'HUJ!!!!!!!!!!!'
                _i=0

legend = '[ASCII merge]:'

# if an output file is specified, write in the file
# otherwise to stdout
sort_by_second_column = False
if sys.argv[1] == '-f':
    if len(sys.argv)>2:
        ofile = open(sys.argv[2], "w")
        sys.stdout = ofile
        sys.argv.pop(1)
        sys.argv.pop(1)

elif sys.argv[1] == '--fix':
    sys.argv.pop(1)
    fix_lines()
    #print legend, 'done fixing missing newlines'
    sys.exit(0)

if sys.argv[1] == '-r':
    sys.argv.pop(1)
    sort_by_second_column = True

#observed number
item_to_find = None
if sys.argv[1] == '-o':
    item_to_find = sys.argv.pop(2)
    sys.argv.pop(1)
    sys.argv.pop(1)

_limits =[]

for line in fileinput.input():
    #print legend, line.strip()

    _words = line.strip().split()

    if _words[0][0] == '#':
        continue

    _number = []
    for word in _words:
        _number.append(float(word))
    
    _limits.append(_number)

if sort_by_second_column:
    _limits.sort(compare_columns_second)
else:
    _limits.sort(compare_columns)

counter = 0
for l in _limits:
    #print l[0],l[1]
    counter +=1
    for n in l:
        print n, '   ',

    #print '#',item_to_find, l[1]
    if item_to_find and float(item_to_find) < float(l[1]):
        print counter, '(found',item_to_find, ')'

    print

if item_to_find:
    print item_to_find
