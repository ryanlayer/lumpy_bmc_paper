import sys
import numpy as np
import matplotlib
import pylab
import random
from optparse import OptionParser

lumpy_color_1='#003366'
lumpy_color_2='#0080ff'
lumpy_color_3='#99ccff'
delly_color_1='#660000'
delly_color_2='#ff0000'



matplotlib.rc('font',**{'family':'sans-serif', \
    'sans-serif':['Helvetica']})
matplotlib.rc('text', usetex=True)


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
            
        R[coverage][proportion]['l_pe'] = [float(a[1]),float(a[2])]
        R[coverage][proportion]['l_sr'] = [float(a[3]),float(a[4])]
        R[coverage][proportion]['l_pesr'] = [float(a[5]),float(a[6])]
        R[coverage][proportion]['d_pe'] = [float(a[7]),float(a[8])]
        R[coverage][proportion]['d_sr'] = [float(a[9]),float(a[10])]

f.close()


width = 0.18   
matplotlib.rcParams.update({'font.size': 12})
fig = matplotlib.pyplot.figure(figsize=(10,5),dpi=300)
fig.subplots_adjust(wspace=.05,left=.01,bottom=.01)
#matplotlib.pyplot.suptitle("5516 1KG Deletions",fontsize=12)
c=1
for coverage in sorted(R.keys()):
    l_pe = []
    l_sr = []
    l_pesr = []
    d_pe = []
    d_sr = []
 
    x = []

    for proportion in sorted(R[coverage].keys()):
        #l_1_1.append(R[coverage][proportion]['l_1_1'][0])
        #l_1_100.append(R[coverage][proportion]['l_1_100'][0])
        #g.append(R[coverage][proportion]['g'][0])
        #l_1_1.append( float(R[coverage][proportion]['l_1_1'][0])/float(T))
        #l_1_100.append(float(R[coverage][proportion]['l_1_100'][0])/float(T))
        #g.append(float(R[coverage][proportion]['g'][0])/float(T))
        #d.append(float(R[coverage][proportion]['d'][0])/float(T))

        l_pe.append( float(R[coverage][proportion]['l_pe'][0])/ float(T))
        l_sr.append( float(R[coverage][proportion]['l_sr'][0])/ float(T))
        l_pesr.append( float(R[coverage][proportion]['l_pesr'][0])/ float(T))
        d_pe.append( float(R[coverage][proportion]['d_pe'][0])/ float(T))
        d_sr.append( float(R[coverage][proportion]['d_sr'][0])/ float(T))




        x.append(proportion)

    N = len(l_pe)
    ind = np.arange(N)
    #ax = fig.add_subplot(4,2,c)
    ax = fig.add_subplot(2,4,c)
    #ax.set_ylim([0,6000])
    ax.set_ylim([0,1.0])
    ax.set_xlim([-0.5,4])

    minorLocator   = matplotlib.ticker.MultipleLocator(0.1)
    ax.yaxis.set_minor_locator(minorLocator)
    ax.yaxis.grid(b=True, which='minor', color='k', linestyle='--')
    ax.yaxis.grid(b=True, which='major', color='k', linestyle='-')
    ax.set_axisbelow(True) 


    l_pe_p = ax.bar(ind, l_pe, width, color=lumpy_color_1)
    l_sr_p = ax.bar(ind+width*1, l_sr, width, color=lumpy_color_2)
    l_pesr_p = ax.bar(ind+width*2, l_pesr, width, color=lumpy_color_3)
    d_pe_p = ax.bar(ind+width*3, d_pe, width, color=delly_color_1)
    d_sr_p = ax.bar(ind+width*4, d_sr, width, color=delly_color_2)




    #for i in range(2,10,2):
        #matplotlib.pyplot.hlines(i/10.0,-0.5,5,colors='w')

    matplotlib.pyplot.title(str(int(coverage)) + 'x')
    matplotlib.pyplot.tick_params(axis='x',length=0)
    matplotlib.pyplot.tick_params(axis='y',length=0)

    if c == 1:
        ax.set_ylabel('Sensitivity')
        ax.set_yticklabels([0,0.2,0.4,0.6,0.8,1.0]) 
        ax.spines['left'].set_visible(False)

    if c != 1:
        ax.set_yticklabels([]) 
        ax.spines['left'].set_visible(False)

    #if c == 1:
        #lg = ax.legend( (l_1_1_p[0], l_1_100_p[0], g_p[0],d_p[0]),\
        #lg = ax.legend( (l_1_1_p[0],  g_p[0],d_p[0]),\
                #('lumpy', 'GASVPro','DELLY'),\
                #prop={'size':10},\
                #loc=2)
        #lg.draw_frame(False)


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
    l_pe = []
    l_sr = []
    l_pesr = []
    d_pe = []
    d_sr = []
 
    x = []
    for proportion in sorted(R[coverage].keys()):
        if float(R[coverage][proportion]['l_pe'][0]) + \
            float(R[coverage][proportion]['l_pe'][1]) > 0 :

            l_pe.append( float(R[coverage][proportion]['l_pe'][1])/ \
                (float(R[coverage][proportion]['l_pe'][0])+ \
                float(R[coverage][proportion]['l_pe'][1])) ) 
        else:
            l_pe.append(0)

        if float(R[coverage][proportion]['l_sr'][0])+ \
            float(R[coverage][proportion]['l_sr'][1]) > 0:
            l_sr.append(float(R[coverage][proportion]['l_sr'][1])/ \
                (float(R[coverage][proportion]['l_sr'][0])+ \
                float(R[coverage][proportion]['l_sr'][1])))
        else:
            l_sr.append(0)
            
        if float(R[coverage][proportion]['l_pesr'][0])+ \
            float(R[coverage][proportion]['l_pesr'][1]) > 0:
            l_pesr.append(float(R[coverage][proportion]['l_pesr'][1])/ \
                (float(R[coverage][proportion]['l_pesr'][0])+ \
                float(R[coverage][proportion]['l_pesr'][1])))
        else:
            l_pesr.append(0)

        if float(R[coverage][proportion]['d_pe'][0])+ \
            float(R[coverage][proportion]['d_pe'][1]) > 0:
            d_pe.append(float(R[coverage][proportion]['d_pe'][1])/ \
                (float(R[coverage][proportion]['d_pe'][0])+ \
                float(R[coverage][proportion]['d_pe'][1])))
        else:
            d_pe.append(0)

        if float(R[coverage][proportion]['d_sr'][0])+ \
            float(R[coverage][proportion]['d_sr'][1]) > 0:
            d_sr.append(float(R[coverage][proportion]['d_sr'][1])/ \
                (float(R[coverage][proportion]['d_sr'][0])+ \
                float(R[coverage][proportion]['d_sr'][1])))
        else:
            d_sr.append(0)


        x.append(proportion)

    N = len(l_pe)
    ind = np.arange(N)
    ax = fig.add_subplot(2,4,c)
    #ax.set_ylim([0,6000])
    ax.set_ylim([0,0.5])
    ax.set_xlim([-0.5,4])

    minorLocator   = matplotlib.ticker.MultipleLocator(0.05)
    ax.yaxis.set_minor_locator(minorLocator)
    ax.yaxis.grid(b=True, which='minor', color='k', linestyle='--')
    ax.yaxis.grid(b=True, which='major', color='k', linestyle='-')
    ax.set_axisbelow(True) 

    l_pe_p = ax.bar(ind, l_pe, width, color=lumpy_color_1)
    l_sr_p = ax.bar(ind+width*1, l_sr, width, color=lumpy_color_2)
    l_pesr_p = ax.bar(ind+width*2, l_pesr, width, color=lumpy_color_3)
    d_pe_p = ax.bar(ind+width*3, d_pe, width, color=delly_color_1)
    d_sr_p = ax.bar(ind+width*4, d_sr, width, color=delly_color_2)

    #for i in range(1,5,1):
        #matplotlib.pyplot.hlines(i/10.0,-0.5,5,colors='w')

    #matplotlib.pyplot.title(pretty_type[sv_type])
    matplotlib.pyplot.tick_params(axis='x',length=0)
    matplotlib.pyplot.tick_params(axis='y',length=0)

    if c == 5:
        ax.set_ylabel('False Discovery Rate')
        ax.set_yticklabels([0,0.1,0.2,0.3,0.4,0.5]) 
        ax.spines['left'].set_visible(False)

#    if c == 5:
#        #lg = ax.legend( (l_1_1_p[0], l_1_100_p[0], g_p[0],d_p[0]),\
#        lg = ax.legend( (l_1_1_p[0],  g_p[0],d_p[0]),\
#                ('lumpy', 'GASVPro','DELLY'),\
#                prop={'size':10},\
#                loc=2)
#
#        lg.draw_frame(False)

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

l=fig.legend((l_pe_p[0],  l_sr_p[0], l_pesr_p[0], d_pe_p[0], d_sr_p[0]),\
                ('LUMPY pe', 'LUMPY sr','LUMPY pesr','DELLY pe','DELLY sr'),\
                prop={'size':10},\
                loc=10,ncol=5)
l.draw_frame(False)

matplotlib.pyplot.savefig('signals.png',bbox_inches='tight')
matplotlib.pyplot.savefig('signals.pdf',bbox_inches='tight')
matplotlib.pyplot.savefig('signals.eps',bbox_inches='tight')
