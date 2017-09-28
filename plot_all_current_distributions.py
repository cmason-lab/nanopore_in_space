import pandas as pd
import h5py
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import glob
import sys
matplotlib.rcParams['svg.fonttype'] = 'none'

def extract_kmer_events(loc,fi):
    location = loc+"/Events"
    try:
        #pA, time, kmer, steps advanced?
        #print fi[location].value[0]
        kmer_events = [(a[0],a[1],a[4],a[6]) for a in fi[location].value]
        #print kmer_events[0]
        return kmer_events
    except KeyError:
        return None

def extract_model(loc,fi):
    location = loc+"/Model"
    try: 
        kmer_list = [(a[0],a[1]) for a in fi[location].value]
        print kmer_list[0]
        model = {}
        for kmer in kmer_list:
            model[kmer[0]] = kmer[1]
        return model
    except KeyError:
        return None

def location_by_version(vers):
    if vers > 50: # == 0.51.1.39:
        loc="Analyses/Basecall_1D_000/BaseCalled_template"
    else:
        loc="Analyses/Basecall_2D_000/BaseCalled_template"
    return loc

def convert2df(kmer_dict):
    kmers = {'kmer':[],'mn':[],'sd':[]}
    for kmer in kmer_dict:
        kmers['kmer'].append(kmer)
        kmers['mn'].append(np.mean(kmer_dict[kmer]))
        kmers['sd'].append(np.std(kmer_dict[kmer]))
    kmer_df = pd.DataFrame(kmers)
    return kmer_df


def print_attrs(name, obj):
    print name
    for key, val in obj.attrs.iteritems():
        print "    %s: %s" % (key, val)

#Directory containing fast5 directories for each run
indirs=[x for sublist in [glob.glob(arg+'/*') for arg in sys.argv[1:]] for x in sublist]
print indirs

#BacterialMix1_minion_1/ version 0.48.2.14
#nasa_4_1_2016/downloads/ 0.51.1.39
#nasa_5_12_2016/downloads/
#Hillier_Bank_R9/downloads/pass/ 0.51.1.66
#Flight_data_2 0.50.2.15
#nasa_6_1_2016/downloads/combined
#nasa_6_2_2016_72ng/downloads/
#Hillier_Sed_R9/downloads/fail/
#nasa_5_13_2016/downloads/combined

names=[dir.split('/')[-1] for dir in indirs]
name="current_distributions"

kmer_dict = {}
model = {}
exp_distribution = {}
for indir,nm in zip(indirs,names):
    #try:
    #   model = pd.read_csv('model_dataframe.csv')
    #except:
    #   pass
    if not model: #xcept IOError:
       num_reads = 0 
       kmer_dict[nm] = {}
       loc = None
       version = None
       #print indir
       infi_list = glob.glob(indir+'/*.fast5')
       if len(infi_list) == 0:
          infi_list = glob.glob(indir+'/downloads/*.fast5')
       if len(infi_list) == 0:
          infi_list = glob.glob(indir+'/downloads/pass/*.fast5') + glob.glob(indir+'/downloads/fail/*.fast5')
       print '# files ',indir,':',len(infi_list) #[:2]
       for infi in infi_list:
           if model:
               break
           try:
               with h5py.File(infi, "r") as hf:
                   #hf.visititems(print_attrs)
                   #print hf.attrs.items()[0][1]
                   #print hf["UniqueGlobalKey/tracking_id"].attrs.get("version")
                   if not version:
                       version = hf["UniqueGlobalKey/tracking_id"].attrs.get("version")
                       pref,vers,pnt,suf = [int(x) for x in version.split()[0].split('.')]
                   if not loc:
                       loc = location_by_version(vers)
                   #print hf.keys()
                   #print hf["/Analyses/Basecall_1D_000/BaseCalled_template"].keys()
                   #print hf["Analyses/Basecall_1D_000/Configuration/general"].attrs.get("template_model")
                   if not model and suf < 66:
                       model = extract_model(loc,hf)
                   #read_key = (hf["Analyses/EventDetection_000/Reads"].keys())[0]
                   #print hf.attrs.items()[0][1]
                   #print hf["UniqueGlobalKey/tracking_id"].attrs.get("version")
                   
                   """events = extract_kmer_events(loc,hf)
                   num_reads += 1
                   if num_reads%1000 == 0:
                      print num_reads
                   if events:
                      for e in events:
                         if e[2] not in kmer_dict[nm]:
                            kmer_dict[nm][e[2]] = [e[0]]
                         else:
                            kmer_dict[nm][e[2]].append(e[0])"""
           except IOError:
              continue
           if model:
              break
    if model:
       try:
          exp_distribution[nm] = pd.read_csv(name+'_'+nm+'_dataframe.csv')
       except IOError:
          pass
    print 'done directory'

print kmer_dict.keys()
      
if kmer_dict and not exp_distribution:
   for nm in names:
      exp_distribution[nm] = convert2df(kmer_dict[nm]) #pd.DataFrame(kmers) #read_table("hmm_temp_events_stdv.csv", sep =';', header = None)

if model:
    model_distribution = convert2df(model)

if model:
    order_distribution = model_distribution
else:
    order_distribution = exp_distribution[names[-1]]
    
gap = int(len(order_distribution.kmer)/50)
order_distribution = order_distribution.sort_values(by='mn')
order_distribution['ind'] = np.arange(0,len(order_distribution.kmer))
len_inds = len(order_distribution.ind)
order_kmers = np.asarray(order_distribution['kmer'])
tick_nums = np.arange(0,len_inds,gap)
ticks = order_kmers[tick_nums]
order_distribution['data_set'] = ['model']*len(order_distribution.index)

if model:
    model_distribution = order_distribution
    
order_distribution.head()


for nm in names:
    exp_distribution[nm]['kmer'] = pd.Categorical(exp_distribution[nm]['kmer'], categories=list(order_kmers), ordered=True)
    exp_distribution[nm].sort_values(by='kmer',inplace=True)

    inds = [list(order_distribution[order_distribution.kmer == k].ind)[0] for k in exp_distribution[nm].kmer]
    exp_distribution[nm]['ind'] = inds
    print exp_distribution[nm].head()

fig = plt.figure(figsize=(10,6))
ax1 = fig.add_subplot(1, 1, 1)
print len(exp_distribution[nm].mn) 
print len(exp_distribution[nm].sd)

colpal = zip(['#2b2381','#56d6e1','#db1818','#ff8c00','#6B0465','#ffa6b8','#BCBF76'],['#ada8e7','#b5e5ea','#d89494','#fbc079','#8ac994','#eecad2','#C7C999'])
for col,nm in zip(colpal,sorted(names,reverse=True)):
    if len(names) < 4:
       ax1.errorbar(exp_distribution[nm].ind,exp_distribution[nm].mn,exp_distribution[nm].sd,color=col[0],fmt='.', markersize='2', linestyle='None',ecolor=col[1],capsize=4,elinewidth=0.5,label=nm,zorder=1)
    else:
       ax1.scatter(exp_distribution[nm].ind,exp_distribution[nm].mn,color=col[0],marker='.',s=4,label=nm,zorder=1)
if model:
    ax1.scatter(model_distribution.ind,model_distribution.mn,color='black',marker='.',s=4, label='model', zorder=2)
if model or len(names) > 1:
    handles, labels = ax1.get_legend_handles_labels()
    #handles = handles[len(names):]
    #labels = labels[len(names):]
    ax1.legend(handles, labels, markerscale=10,loc=2)
else:
    ax1.legend().set_visible(False)
plt.xticks(tick_nums,ticks,rotation=90)
plt.xlim([-10,tick_nums[-1]+10])
plt.ylabel('Current (pA)')
plt.xlabel('k-mer')
plt.savefig(name+'_all.eps',dpi=500,bbox_inches='tight')
plt.savefig(name+'_all.png',dpi=500,bbox_inches='tight')

for nm in exp_distribution:
    if len(glob.glob(name+'_'+nm+'_dataframe.csv')) == 0:
       exp_distribution[nm].to_csv(name+'_'+nm+'_dataframe.csv')
if model:
    if len(glob.glob('model_dataframe.csv')) == 0:
       model_distribution.to_csv('model_dataframe.csv')

