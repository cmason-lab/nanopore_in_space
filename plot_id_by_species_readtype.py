from __future__ import division
import re
import sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
import glob

samfis = glob.glob("*.sam")
reads = ['2d','template','complement']
cols = sns.color_palette(['#56d6e1','#2b2381','#ffa6b8','#db1818']) #,'#a34c7c'
sns.set_palette(cols)
sns.set_style('ticks')
fig = plt.plot()
#outfi = open('graphmap_identities.tab','w')
species_ids = {sam:{read:{x:[] for x in ['Enterobacteria_phage_lambda','Escherichia_coli','Mus_musculus']} for read in reads} for sam in samfis}
#match_lengths = {sam:{x:[] for x in ['Enterobacteria_phage_lambda','Escherichia_coli','Mus_musculus']} for sam in samfis}

for sam in samfis:

   infi = open(sam,'r').read().split('\n') #open(sys.argv[1],'r').read().split('\n')
   #identities = []
   #species_ids = {x:[] for x in ['Enterobacteria_phage_lambda','Escherichia_coli','Mus_musculus']}

   for line in infi:
      if len(line) > 0 and line[0] != '@' and line.split('\t')[1] != '4':
         readtype = line.split('\t')[0].split('_')[-1]
         len_read = len(line.split('\t')[9])
         md = line.split('MD:Z:')[1].split('\t')[0]
         num_matches = sum([int(x) for x in re.split('A|C|G|T|N|\^',md) if x != ''])
         len_seq = len([x for x in re.split('A|C|G|T|N|\^',md)])+num_matches-len(md.split('^'))
         percent_id = num_matches/len_seq
         #identities.append(percent_id)
         for species in species_ids[sam][readtype]:
            #print line.split(species)
            if len(line.split(species)) == 2:
               #print line.split(species)
               species_ids[sam][readtype][species].append(percent_id)
               break
         #break
         
   #sns.kdeplot(np.array(identities),label=sam.split('.')[0]+', mean = '+str(round(np.mean(identities),2)),legend=False)
   #plt.title(sys.argv[1].split('.')[0])
   #print sam+'\n\tmax identity =',np.max(identities),'\n\tmean identity = ', np.mean(identities)
   #for species in species_ids:
   #   print species, np.mean(species_ids[species])

   #outfi.write(sam.split('.')[0]+'\t'+str(np.max(identities))+'\t'+str(np.mean(identities)))
   #outfi.write(sys.argv[1].split('.')[0]+'\t'+str(np.max(identities))+'\t'+str(np.mean(identities)))

#outfi.close()

plt.legend(loc=2)
plt.xlabel("Identity")
plt.ylabel("Density") ## Reads")
plt.xlim([0,1])
sns.despine()
for ext in ['png','pdf']:
   plt.savefig('graphmap_identity.pdf',dpi=900,bbox_inches='tight')

cols = sns.color_palette(['#FF934F','#E86412','#B56E42','#89B6A5','#3A4C45','#6A706E','#DD8BBD','#D66DAC','#AD7A99']) #,'#a34c7c'
sns.set_palette(cols)
for i, sam in enumerate(samfis):
   fig = plt.figure()
   for readtype in species_ids[sam]:
      for species in species_ids[sam][readtype]:
         print species, round(np.median(species_ids[sam][readtype][species]),2)
         sns.distplot(species_ids[sam][readtype][species],label=species+' '+readtype+', median = '+str(round(np.median(species_ids[sam][readtype][species]),2))+', mean = '+str(round(np.mean(species_ids[sam][readtype][species]),2)))
   plt.title(sam.split('.')[0])
   plt.legend(loc=2,prop={'size':6})
   plt.xlabel("Identity")
   plt.ylabel("Density")
   plt.xlim([0,1])
   sns.despine()
   for ext in ['png','pdf']:
      plt.savefig(sam.split('.')[0]+'_id_by_species_readtype.'+ext,dpi=900,bbox_inches='tight')

