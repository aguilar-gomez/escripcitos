#!/usr/bin/env python

import os, subprocess, sys

if len(sys.argv) != 4 or sys.argv[1] in ["-h","-help","--help"]:
    sys.stderr.write("Usage: "+os.path.basename(sys.argv[0])+" <in.prefix> <out.prefix> <filter.percentid>\n")
    sys.stderr.write("This script filters plot files based on percent identity.\n")
    sys.exit(1)

fplot = open(sys.argv[1]+'.fplot', 'r')
rplot = open(sys.argv[1]+'.rplot', 'r')
gp = open(sys.argv[1]+'.gp', 'r')
out_fplot = open(sys.argv[2]+'.fplot', 'w')
out_rplot = open(sys.argv[2]+'.rplot', 'w')
out_gp = open(sys.argv[2]+'.gp', 'w')

blank = 0

def rm_filter(infile,outfile):
    for line in infile:
        if not line.strip():
            if blank == 0 or blank == 1:
                blank += 1
                outfile.write(line)
        elif line.startswith('#') or line.startswith('0 0'):
            blank =	0
            outfile.write(line)
        else:
            if float(line.rstrip().split(' ')[2]) >= float(sys.argv[3]):
                blank = 0
                outfile.write(line)

rm_filter(fplot,out_fplot)
rm_filter(rplot,out_rplot)

for line in gp:
    if line.startswith('set output'):
        out_gp.write('set output "'+sys.argv[2]+'.png"\n')
    elif 'fplot' in line:
        out_gp.write(' "'+sys.argv[2]+'.fplot" title "FWD" w lp ls 1, \\'+'\n')
    elif 'rplot' in line:
        out_gp.write(' "'+sys.argv[2]+'.rplot" title "REV" w lp ls 2'+'\n')
    elif line.startswith('set style line'):
        out_gp.write(line.replace('pt 6 ps 1','pt 1 ps 0.1'))
    elif line.startswith('set grid'):
        out_gp.write('set grid\n')
    else:
        out_gp.write(line)

fplot.close()
rplot.close()
gp.close()
out_fplot.close()
out_rplot.close()
out_gp.close()

subprocess.check_call('gnuplot '+sys.argv[2]+'.gp',shell=True)
