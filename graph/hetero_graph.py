import sys
import numpy as np
import matplotlib
import pylab
import random
from optparse import OptionParser

T=5516

parser = OptionParser()

parser.add_option("-d",
                  "--data_file",
                  dest="data_file",
                  help="Data file")

(options, args) = parser.parse_args()

if not options.data_file:
    parser.error('Data file not given')

f = open(options.data_file,'r')

R={}

for l in f:
    if l[0] != '#':
        a = l.rstrip().split(',')
        n_a = a[0].split('_')

        mut_p = float(n_a[1])
        wt_p = float(n_a[3])

        coverage = mut_p + wt_p
        proportion = mut_p / coverage

        if not coverage in R:
            R[coverage] = {}

        if not proportion  in R[coverage]:
            R[coverage][proportion] = {}
            
        R[coverage][proportion]['l_1_1'] = [float(a[1]),float(a[2])]
        R[coverage][proportion]['l_1_100'] = [float(a[3]),float(a[4])]
        R[coverage][proportion]['g'] = [float(a[5]),float(a[6])]
        R[coverage][proportion]['d'] = [float(a[7]),float(a[8])]
 
f.close()


width = 0.30   
matplotlib.rcParams.update({'font.size': 12})
fig = matplotlib.pyplot.figure(figsize=(10,5),dpi=300)
fig.subplots_adjust(wspace=.05,left=.01,bottom=.01)
#matplotlib.pyplot.suptitle("5516 1KG Deletions",fontsize=12)
c=1
for coverage in sorted(R.keys()):
    l_1_1 = []
    l_1_100 = []
    g = []
    d = []
    x = []

    for proportion in sorted(R[coverage].keys()):
        #l_1_1.append(R[coverage][proportion]['l_1_1'][0])
        #l_1_100.append(R[coverage][proportion]['l_1_100'][0])
        #g.append(R[coverage][proportion]['g'][0])
        l_1_1.append( float(R[coverage][proportion]['l_1_1'][0])/float(T))
        l_1_100.append(float(R[coverage][proportion]['l_1_100'][0])/float(T))
        g.append(float(R[coverage][proportion]['g'][0])/float(T))
        d.append(float(R[coverage][proportion]['d'][0])/float(T))
        x.append(proportion)

    N = len(l_1_1)
    ind = np.arange(N)
    #ax = fig.add_subplot(4,2,c)
    ax = fig.add_subplot(2,4,c)
    #ax.set_ylim([0,6000])
    ax.set_ylim([0,1])
    ax.set_xlim([-0.5,4])
    l_1_1_p = ax.bar(ind, l_1_1, width, color='k')
    #l_1_100_p = ax.bar(ind+width, l_1_100, width, color='r')
    g_p = ax.bar(ind+width*1, g, width, color='r')
    d_p = ax.bar(ind+width*2, d, width, color='g')

    for i in range(2,10,2):
        matplotlib.pyplot.hlines(i/10.0,-0.5,5,colors='w')

    matplotlib.pyplot.title(str(int(coverage)) + 'x')
    matplotlib.pyplot.tick_params(axis='x',length=0)
    matplotlib.pyplot.tick_params(axis='y',length=0)

    if c == 1:
        ax.set_ylabel('Sensitivity')

    if c != 1:
        ax.set_yticklabels([]) 
        ax.spines['left'].set_visible(False)

    if c == 7:
        ax.set_xticks(ind+width*2)
        ax.set_xlabel('Sequencing Coverage')
        ax.set_xticklabels(x) 
    else:
        ax.set_xticklabels(['']*len(x)) 

    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    c+=1


for coverage in sorted(R.keys()):
    print coverage
    l_1_1 = []
    l_1_100 = []
    g = []
    d = []
    x = []
    for proportion in sorted(R[coverage].keys()):
        #l_1_1.append( float(R[coverage][proportion]['l_1_1'][1]))
        #l_1_100.append(float(R[coverage][proportion]['l_1_100'][1]))
        #g.append(float(R[coverage][proportion]['g'][1]))
        if float(R[coverage][proportion]['l_1_1'][0])+\
                float(R[coverage][proportion]['l_1_1'][1]) > 0:
            l_1_1.append( float(R[coverage][proportion]['l_1_1'][1])/\
                (float(R[coverage][proportion]['l_1_1'][0])+\
                float(R[coverage][proportion]['l_1_1'][1])))
        else:
            l_1_1.append(0)

        if float(R[coverage][proportion]['l_1_100'][0])+\
                float(R[coverage][proportion]['l_1_100'][1]) > 0:
            l_1_100.append(float(R[coverage][proportion]['l_1_100'][1])/\
                (float(R[coverage][proportion]['l_1_100'][0])+
                float(R[coverage][proportion]['l_1_100'][1])))
        else:
            l_1_100.append(0)


        if float(R[coverage][proportion]['g'][0])+\
                float(R[coverage][proportion]['g'][1]) > 0:
            g.append(float(R[coverage][proportion]['g'][1])/\
                (float(R[coverage][proportion]['g'][0])+
                float(R[coverage][proportion]['g'][1])))
        else:
            g.append(0)
        if float(R[coverage][proportion]['d'][0])+\
                float(R[coverage][proportion]['d'][1]) > 0:
            d.append(float(R[coverage][proportion]['d'][1])/\
                (float(R[coverage][proportion]['d'][0])+
                float(R[coverage][proportion]['d'][1])))
        else:
            d.append(0)

        x.append(proportion)

    N = len(l_1_1)
    ind = np.arange(N)
    ax = fig.add_subplot(2,4,c)
    #ax.set_ylim([0,6000])
    ax.set_ylim([0,1])
    ax.set_xlim([-0.5,4])
    l_1_1_p = ax.bar(ind, l_1_1, width, color='k')
    #l_1_100_p = ax.bar(ind+width, l_1_100, width, color='r')
    g_p = ax.bar(ind+width*1, g, width, color='r')
    d_p = ax.bar(ind+width*2, d, width, color='g')

    for i in range(2,10,2):
        matplotlib.pyplot.hlines(i/10.0,-0.5,5,colors='w')

    #matplotlib.pyplot.title(pretty_type[sv_type])
    matplotlib.pyplot.tick_params(axis='x',length=0)
    matplotlib.pyplot.tick_params(axis='y',length=0)

    if c == 5:
        ax.set_ylabel('Positive Predictive Value')

    if c == 5:
        #lg = ax.legend( (l_1_1_p[0], l_1_100_p[0], g_p[0],d_p[0]),\
        lg = ax.legend( (l_1_1_p[0],  g_p[0],d_p[0]),\
                ('lumpy', 'GASVPro','DELLY'),\
                prop={'size':10},\
                loc=2)

        lg.draw_frame(False)

    if c != 5:
        ax.set_yticklabels([]) 
        ax.spines['left'].set_visible(False)

    ax.set_xticks(ind+width*2)
    ax.set_xlabel('SV Allele Frequency')
    ax.set_xticklabels(x) 
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    #ax.set_tick_params(which='x',length=0)
    #fig.autofmt_xdate()

    #c+=2
    c+=1

matplotlib.pyplot.savefig('heterozygous.png',bbox_inches='tight')
matplotlib.pyplot.savefig('heterozygous.pdf',bbox_inches='tight')
matplotlib.pyplot.savefig('heterozygous.eps',bbox_inches='tight')
