import sys
import numpy as np
import matplotlib
import pylab
import random
from optparse import OptionParser


lumpy_color='#003366'
gasv_color='#0080ff'
delly_color='#99ccff'

matplotlib.rc('font',**{'family':'sans-serif', \
    'sans-serif':['Helvetica']})
matplotlib.rc('text', usetex=True)

T={ 'del':1000, \
    'dup':1000, \
    'inr':2000, \
    'inv':1000 }

pretty_type={ 'del':'Deletion', \
              'dup':'Duplication', \
              'inr':'Insertion', \
              'inv':'Inversion' }

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
    a = l.rstrip().split(',')
    n_a = a[0].split('.')

    sv_type = n_a[1]
    coverage = int(n_a[2])

    #coverage = mut_p + wt_p
    #proportion = mut_p / coverage

    if not sv_type in R:
        R[sv_type] = {}

    if not coverage  in R[sv_type]:
        R[sv_type][coverage] = {}
        
    R[sv_type][coverage]['l_1_1'] = [float(a[1]),float(a[2])]
    R[sv_type][coverage]['l_1_100'] = [float(a[3]),float(a[4])]
    R[sv_type][coverage]['g'] = [float(a[5]),float(a[6])]
    R[sv_type][coverage]['d'] = [float(a[7]),float(a[8])]
    
f.close()


width = 0.30   
matplotlib.rcParams.update({'font.size': 12})
fig = matplotlib.pyplot.figure(figsize=(10,5),dpi=300)
fig.subplots_adjust(wspace=.05,left=.01,bottom=.01)
#matplotlib.pyplot.suptitle("5516 1KG Deletions",fontsize=12)
c=1
for sv_type in sorted(R.keys()):
    l_1_1 = []
    l_1_100 = []
    g = []
    d = []
    x = []
    for coverage in sorted(R[sv_type].keys()):
        l_1_1.append( float(R[sv_type][coverage]['l_1_1'][0])/
                float(T[sv_type]))
        l_1_100.append(float(R[sv_type][coverage]['l_1_100'][0])/
                float(T[sv_type]))
        g.append(float(R[sv_type][coverage]['g'][0])/
                float(T[sv_type]))
        d.append(float(R[sv_type][coverage]['d'][0])/
                float(T[sv_type]))

        x.append(coverage)

    N = len(l_1_1)
    ind = np.arange(N)
    ax = fig.add_subplot(2,4,c)



    minorLocator   = matplotlib.ticker.MultipleLocator(0.1)
    ax.yaxis.set_minor_locator(minorLocator)
    ax.yaxis.grid(b=True, which='minor', color='k', linestyle='--')
    ax.yaxis.grid(b=True, which='major', color='k', linestyle='-')
    ax.set_axisbelow(True) 



    ax.set_ylim([0,1])
    ax.set_xlim([-0.5,5])
    l_1_1_p = ax.bar(ind, l_1_1, width, color=lumpy_color)
    #l_1_100_p = ax.bar(ind+width, l_1_100, width, color='r')
    g_p = ax.bar(ind+width*1, g, width, color=gasv_color)
    d_p = ax.bar(ind+width*2, d, width, color=delly_color)

    #for i in range(2,10,2):
        #matplotlib.pyplot.hlines(i/10.0,-0.5,5,colors='w')

    matplotlib.pyplot.title(pretty_type[sv_type])
    matplotlib.pyplot.tick_params(axis='x',length=0)
    matplotlib.pyplot.tick_params(axis='y',length=0)

    if c == 1:
        ax.set_ylabel('Sensitivity')
        ax.set_yticklabels([0,0.2,0.4,0.6,0.8,1.0]) 

    if c != 1:
        ax.set_yticklabels([]) 
        #ax.spines['left'].set_visible(False)


    ax.set_xticklabels(['']*len(x)) 

    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_visible(False)
    c+=1

for sv_type in sorted(R.keys()):
    print sv_type
    l_1_1 = []
    l_1_100 = []
    g = []
    d = []
    x = []
    for coverage in sorted(R[sv_type].keys()):
        if float(R[sv_type][coverage]['l_1_1'][0]) + \
            float(R[sv_type][coverage]['l_1_1'][1]) > 0 :

            l_1_1.append( float(R[sv_type][coverage]['l_1_1'][1])/ \
                (float(R[sv_type][coverage]['l_1_1'][0])+ \
                float(R[sv_type][coverage]['l_1_1'][1])) ) 
        else:
            l_1_1.append(0)

        if float(R[sv_type][coverage]['l_1_100'][0])+ \
            float(R[sv_type][coverage]['l_1_100'][1]) > 0:
            l_1_100.append(float(R[sv_type][coverage]['l_1_100'][1])/ \
                (float(R[sv_type][coverage]['l_1_100'][0])+ \
                float(R[sv_type][coverage]['l_1_100'][1])))
        else:
            l_1_100.append(0)
            
        if float(R[sv_type][coverage]['g'][0])+ \
            float(R[sv_type][coverage]['g'][1]) > 0:
            g.append(float(R[sv_type][coverage]['g'][1])/ \
                (float(R[sv_type][coverage]['g'][0])+ \
                float(R[sv_type][coverage]['g'][1])))
        else:
            g.append(0)

        if float(R[sv_type][coverage]['d'][0])+ \
            float(R[sv_type][coverage]['d'][1]) > 0:
            d.append(float(R[sv_type][coverage]['d'][1])/ \
                (float(R[sv_type][coverage]['d'][0])+ \
                float(R[sv_type][coverage]['d'][1])))
        else:
            d.append(0)

        x.append(coverage)

    print l_1_1
    print l_1_100
    print g

    N = len(l_1_1)
    ind = np.arange(N)
    ax = fig.add_subplot(2,4,c)
    #ax.set_ylim([0,6000])
    ax.set_ylim([0,0.5])
    ax.set_xlim([-0.5,5])

    #matplotlib.axis.Tick(color='g')

    matplotlib.pyplot.tick_params(axis='x',length=0)
    matplotlib.pyplot.tick_params(axis='y',length=0)

    minorLocator   = matplotlib.ticker.MultipleLocator(0.05)
    ax.yaxis.set_minor_locator(minorLocator)
    ax.yaxis.grid(b=True, which='minor', color='k', linestyle='--')
    ax.yaxis.grid(b=True, which='major', color='k', linestyle='-')



    print ax.yaxis
    ax.set_axisbelow(True) 
    #ax.yaxis.grid(b=True, which='minor', color='b', linestyle='-')

    l_1_1_p = ax.bar(ind, l_1_1, width, color=lumpy_color)
    #l_1_100_p = ax.bar(ind+width, l_1_100, width, color='r')
    g_p = ax.bar(ind+width*1, g, width, color=gasv_color)
    d_p = ax.bar(ind+width*2, d, width, color=delly_color)

    #for i in range(1,5,1):
        #matplotlib.pyplot.hlines(i/10.0,-0.5,5,colors='w')

    #matplotlib.pyplot.title(pretty_type[sv_type])

    if c == 5:
        ax.set_ylabel('False Discovery Rate')
        ax.set_xticklabels([0,0.1,0.2,0.3,0.4,0.5])

    if c == 5:
        #lg = ax.legend( (l_1_1_p[0], l_1_100_p[0], g_p[0],d_p[0]),\
        #lg = ax.legend( (l_1_1_p[0],  g_p[0],d_p[0]),\
                #('lumpy', 'GASVPro','DELLY'),\
                #prop={'size':10},\
                #loc=4)

        #lg.draw_frame(False)
        ax.set_yticklabels([0,0.1,0.2,0.3,0.4,0.5]) 

    if c != 5:
        ax.set_yticklabels([]) 
        ax.spines['left'].set_visible(False)

    ax.set_xticks(ind+width*1.5)
    ax.set_xlabel('Sequencing Coverage')
    pretty_x = []
    for p_x in x:
        pretty_x.append(str(p_x) + 'x')
    ax.set_xticklabels(pretty_x)

    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_visible(False)
    c+=1

l=fig.legend((l_1_1_p[0],  g_p[0],d_p[0]),\
                ('LUMPY', 'GASVPro','DELLY'),\
                prop={'size':10},\
                loc=10,ncol=3)
l.draw_frame(False)


matplotlib.pyplot.savefig('homozygous.png',bbox_inches='tight')
matplotlib.pyplot.savefig('homozygous.pdf',bbox_inches='tight')
matplotlib.pyplot.savefig('homozygous.eps',bbox_inches='tight')
