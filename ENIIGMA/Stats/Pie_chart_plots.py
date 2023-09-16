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


def pie(dir=os.getcwd() + '/', sig_level=1.):
    """
	Create Pie charts.

	Parameters
	-------------

	sig_level : 'float'
		Confidence interval
		Default = 1

	"""

    print('====================================================')
    print('Creating pie chart...')
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
    #print('NN:', n_genes, df.shape)
    
    
    ##
    name = tb['name']
    chi = tb['best_chi']
    deltachi = chi - min(chi)
    

    sig1 = sig_level
    
    n_solutions = len(np.where( deltachi <= sig1)[0])/n_genes
    #print('Nsolutions',n_solutions, len(np.where(deltachi <= sig1)[0]))
    

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
    #print('N', nselec)

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
    
    #print('check:', name, freq, t[1], ts[1])
    
    Data1 = {'name': name, 'freq': freq}
    df1 = DataFrame(Data1, columns=['name', 'freq'])
    df1.sort_values(by=['freq'], inplace=True)
    df1.to_csv('frequency_list.csv', index=False)

    plt.figure(1, figsize=(15, 10))
    cmap = plt.get_cmap('CMRmap_r')
    colors = [cmap(i) for i in np.linspace(0.2, 1.0, len(name))]

    """
    def absolute_value(val):
        # a  = numpy.round(val/100.*(df1['freq']/max(freq)).sum(), 0)
        a = (df1['freq'] / max(freq))[
            np.abs((df1['freq'] / max(freq)) - val / 100. * (df1['freq'] / max(freq)).sum()).idxmin()]
        a = a * 100
        b = ('% 1.1f' % a)
        return b
    """
    
    a = (df1['freq'].to_numpy() / n_solutions)
    values = a * 100
    print(df1['freq'].to_numpy(), n_solutions, values)
    
    def make_autopct(values):
      def my_autopct(pct):
        total = sum(values)
        val = int(round(pct*total/100.0))
        return '{v:d}'.format(v=val)
      return my_autopct

    
    #print('a:', a)
    #print(df1['name'], df1['freq'])
    
    # Create a pie chart
    plt.pie(
        # using data total)arrests
        df1['freq']/ n_solutions,
        # with the labels being officer names
        labels=df1['name'],
        # with no shadows
        shadow=False,
        # with colors
        colors=colors,
        # with one slide exploded out
        # explode=(0.6, 0.45, 0.3, 0.15, 0.05, 0.0, 0.0, 0.0, 0.0),
        # with the start angle at 90%
        startangle=90,
        # with the percent listed as a fraction
        autopct=make_autopct(values),
        textprops={'fontsize': 16},
    )

    # plt.suptitle('Recorrence in %', fontsize=16)
    # View the plot drop above
    plt.axis('equal')

    # View the plot
    plt.tight_layout()
    plt.savefig(pathdir + 'Pie_chart.pdf', format='pdf', bbox_inches='tight', pad_inches=0)
# plt.show()

# pathdir = '/Users/will_rocha_starplan/Downloads/Data_StarPlan_Project/Fitting/WLTS/WL17_5_8/Workspace/Processing/Interp_proc/'
# pie(pathdir+'Best_comb.csv', .1)
