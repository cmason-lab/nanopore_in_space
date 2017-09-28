import pandas as pd
import h5py
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import glob
import sys
from scipy.stats import ks_2samp
matplotlib.rcParams['svg.fonttype'] = 'none'

def location_by_version(vers):
    if vers > 50: # == 0.51.1.39:
        loc="Analyses/Basecall_1D_000/BaseCalled_template"
    else:
        loc="Analyses/Basecall_2D_000/BaseCalled_template"
    return loc

def print_attrs(name, obj):
    print name
    for key, val in obj.attrs.iteritems():
        print "    %s: %s" % (key, val)

indirs=glob.glob(sys.argv[1]+'/*') #'/zenodotus/masonlab/epitranscriptome_scratch/abm237/runs_to_date/nasa_06_2016/*')
names=[dir.split('/')[-1] for dir in indirs]
name=indirs[0].split('/')[-2]

stays = {nm:[] for nm in names}
skips = {nm:[] for nm in names}
for indir,nm in zip(indirs,names):
   try:
       stays[nm] = [float(x) for x in open(nm+'_stays_list.txt','r').read().split('\n')]
       skips[nm] = [float(x) for x in open(nm+'_skips_list.txt','r').read().split('\n')]
   except (IOError,ValueError) as e:
       num_reads = 0 
       loc = None
       version = None
       #print indir
       infi_list = glob.glob(indir+'/*.fast5')
       if len(infi_list) == 0:
          infi_list = glob.glob(indir+'/downloads/*.fast5')
       if len(infi_list) == 0:
          infi_list = glob.glob(indir+'/downloads/pass/*.fast5') + glob.glob(indir+'/downloads/fail/*.fast5')
       print nm, '# files = ',len(infi_list) #[:2]
       for infi in infi_list:
           try:
               with h5py.File(infi, "r") as hf:
                   #hf.visititems(print_attrs)
                   num_stays = hf['Analyses/Basecall_1D_000/Summary/basecall_1d_template'].attrs.get('num_stays')
                   num_skips = hf['Analyses/Basecall_1D_000/Summary/basecall_1d_template'].attrs.get('num_skips')
                   num_bases = hf['Analyses/Basecall_1D_000/Summary/basecall_1d_template'].attrs.get('sequence_length') 
                   stays[nm].append(num_stays*1./num_bases)
                   skips[nm].append(num_skips*1./num_bases)
                   num_reads += 1
                   if num_reads%500 == 0:
                      print num_reads
           except IOError:
              continue
           except KeyError:
              continue
       #try:
       #   exp_distribution[nm] = pd.read_csv(name+'_'+nm+'_dataframe.csv')
       #except IOError:
       #   pass
       print 'done directory'
        
#print stays
#print skips
fig = plt.figure(figsize=(10,6))

cols = sns.color_palette(['#56d6e1','#2b2381','#db1818']) #'#66B3BA','#C17767','#A4F9C8','#5C6F68','#93A3B1']) #'#ff8c00','#08901e','#ffa6b8','#56d6e1','#2b2381','#db1818']) 
sns.set_style('white')
sns.set_palette(cols)
for sname,sdict in zip(["stays","skips"],[stays,skips]):
   fig = plt.figure()
   maximum = 0
   for nm in sorted(names):
      if len(sdict[nm]) > 0:
         cut_list = [x for x in sdict[nm] if x < 0.5]
         if max(sdict[nm]) > maximum:
            maximum = max(cut_list)
         sns.distplot(cut_list,label=nm)
   print 'max',sname,'per read =',maximum
   plt.legend(loc=1)
   plt.xlim([0,round(maximum,2)+0.01])
   plt.ylabel('Read density')
   plt.xlabel('# '+sname+'/base')
   if len(names) == 2:
      k_s = ks_2samp(sdict[names[0]],sdict[names[1]])
      plt.title(sname+', K-S test: D='+str(k_s[0])+', p='+str(k_s[1]))
   sns.despine()
   for ext in ['png','pdf']:
      plt.savefig(name+'_'+sname+'.'+ext,dpi=900,bbox_inches='tight')
   for nm in sorted(names):
      with open(nm+'_'+sname+'_list.txt','w') as listfi:
         listfi.write('\n'.join([str(x) for x in sdict[nm]]))
