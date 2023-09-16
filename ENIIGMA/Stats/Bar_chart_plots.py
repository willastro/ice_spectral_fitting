import matplotlib.pyplot as plt
import glob
from pandas import DataFrame
import numpy as np
import os
import pandas as pd
import matplotlib
from matplotlib.ticker import MultipleLocator, FormatStrFormatter, ScalarFormatter, LogLocator
from matplotlib.backends.backend_pdf import PdfPages

plt.rcParams['axes.linewidth'] = 1.5

matplotlib.rcParams['font.size'] = 16


def word_count(filename, word):
    with open(filename, 'r') as f:
        return f.read().count(word)


def bar(dir=os.getcwd() + '/', sig_level=1.):
    """
	Create Pie charts.

	Parameters
	-------------

	sig_level : 'float'
		Confidence interval
		Default = 1

	"""

    print('====================================================')
    print('Creating bar chart...')
    print('====================================================')

    pathdir = dir + 'Workspace/Processing/Interp_proc/'
    store_f = pathdir + 'Best_comb.csv'

    tb = pd.read_csv(store_f, sep=',')
    ##
    Best = store_f
    df = pd.read_csv(Best,sep=',')
    n_genes = df.shape[1] - 3 # number of genes
    n_lines = df.shape[0]
    #n_solutions = n_lines/n_genes
    #print('NN:', n_genes)
    
    header = []
    for h in range(n_genes):
      header.append('w'+str(h+1))
    
    data = pd.read_csv(Best, sep=',', low_memory=True, usecols=header, nrows=1)
    cmin = data.T.values.tolist()
    #print('cmin is:', cmin)
    
    
    ##
    name = tb['name']
    chi = tb['best_chi']
    deltachi = chi - min(chi)

    sig1 = sig_level
    
    n_solutions = len(np.where( deltachi <= sig1)[0])/n_genes

    fs = open(pathdir + 'select1.txt', 'w')
    count = 0
    for t0 in range(0, n_lines, n_genes):
      for z in range(n_genes):
        #print(z,t0, 'w'+str(z+1))
        #rename = os.path.splitext(os.path.basename(name[count]))[0]
        #print(tb['w'+str(z+1)][t0])
        if deltachi[t0] <= sig1 and tb['w'+str(z+1)][t0] > 0.:
          rename = os.path.splitext(os.path.basename(name[count]))[0]
          rename = rename.replace('_', ':').replace('V3', '').replace('H2', 'H$_2$').replace('H3', 'H$_3$').replace('H4', 'H$_4$').replace('H5', 'H$_5$').replace('H6', 'H$_6$').replace('H7', 'H$_7$').replace('H8','H$_8$').replace('H9', 'H$_9$').replace('O2', 'O$_2$').replace('O3', 'O$_3$').replace('O4', 'O$_4$').replace('O5','O$_5$').replace('O6', 'O$_6$').replace('O7', 'O$_7$').replace('O8', 'O$_8$').replace('O9', 'O$_9$').replace('C2','C$_2$').replace('C3', 'C$_3$').replace('C4', 'C$_4$').replace('C5', 'C$_5$').replace('C6', 'C$_6$').replace('C7','C$_7$').replace('C8', 'C$_8$').replace('C9', 'C$_9$').replace('N2', 'N$_2$').replace('N3', 'N$_3$').replace('N4','N$_4$').replace('N5', 'N$_5$').replace('N6', 'N$_6$').replace('N7', 'N$_7$').replace('N8', 'N$_8$').replace('N9','N$_9$')
          #print('P2', tb['w'+str(z+1)][t0], rename)
          fs.write('{0:s}\n'.format(rename))
        count = count + 1
    fs.close()
    
    
    store_f1 = pathdir + 'select1.txt'

    tb = pd.read_csv(store_f1, sep='\s+', header=None)
    nselec = tb[0]

    org = list(set(nselec))

    forg = open(pathdir + 'select2.txt', 'w')
    for jo in range(len(org)):
        # print org[jo]
        forg.write('{0:s}\n'.format(org[jo]))
    forg.close()

    store_f2 = pathdir + 'select2.txt'
    tb = pd.read_csv(store_f2, sep='\s+', header=None)
    nselec = tb[0]

    list2 = list(range(len(nselec)))

    
    f = open(pathdir + 'frequency_list_short.txt', 'w')
    for jj in list2:
        # print nselec[jj]
        wc = word_count(store_f2, nselec[jj])
        f.write('{0:s} {1:d}\n'.format(nselec[jj], wc))
    f.close()

    f = open(pathdir + 'frequency_list.txt', 'w')
    for jj in list2:
        # print nselec[jj]
        wc = word_count(store_f1, nselec[jj])
        f.write('{0:s} {1:d}\n'.format(nselec[jj], wc))
    f.close()
    

    path = pathdir + 'frequency_list.txt'
    paths = pathdir + 'frequency_list_short.txt'
    t = pd.read_csv(path, sep='\s+', header=None)
    ts = pd.read_csv(paths, sep='\s+', header=None)
    name = t[0]
    freq = t[1] / ts[1]
    
    #print(freq, t[1], ts[1])
    
    Data1 = {'name': name, 'freq': freq}
    df1 = DataFrame(Data1, columns=['name', 'freq'])
    df1.sort_values(by=['freq'], inplace=True)
    df1.to_csv('frequency_list.csv', index=False)

    plt.figure(1, figsize=(15, 10))
    cmap = plt.get_cmap('CMRmap_r')
    colors = [cmap(i) for i in np.linspace(0.2, 1.0, len(name))]

    
    #Create a bar chart
    
    pie = (df1['freq'].to_numpy() / n_solutions)*100
    
    barWidth = 0.5
    
    # Set position of bar on X axis
    r1 = np.arange(len(pie))
    cmap = plt.get_cmap('inferno')
    colors = [cmap(i) for i in np.linspace(0.2, 1.0, len(pie))]
    
    fig = plt.figure(figsize=(10,10), dpi=300)
    ax2=fig.add_subplot(111)
    ax2.bar(r1, pie, color=colors, width=barWidth, edgecolor='black', zorder=5)
    pps = ax2.bar(r1, pie, color=colors, width=barWidth, edgecolor='black', zorder=5)
    for p in pps:
      height = p.get_height()
      ax2.text(x=p.get_x() + p.get_width() / 2, y=height+.90, s="{}".format(height), ha='center')
    
    ax2.set_ylabel('Recurrence (%)', fontsize=16)
    ax2.tick_params(direction='in', which='both')
    sp = df1['name']
    plt.xticks([r + 0 for r in range(len(pie))], sp, rotation=90, fontsize=12)
    ax2.tick_params(which='major', length=0, width=1, direction='in', labelsize=12)
    plt.tight_layout()
    plt.grid(b=True, which='major', linestyle=':')
    plt.tight_layout()
    plt.savefig(pathdir + 'Bar_chart.pdf', format='pdf', bbox_inches='tight', pad_inches=0.1)

    
    
