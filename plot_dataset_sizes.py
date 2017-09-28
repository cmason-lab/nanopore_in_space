import matplotlib.pyplot as plt
import seaborn as sns

fig = plt.figure(figsize=(5,5))
sns.set_style('white')
ground, iss = [],[]
for ds,ap in zip(open('dataset_sizes.txt','r'),open('active_pores.txt','r')):
   try:
      ground.append(float(ds.split('\t')[0])/float(ap.split('\t')[0]))
      iss.append(float(ds.split('\t')[1].strip())/float(ap.split('\t')[1].strip()))
   except ValueError:
      pass

plt.scatter(iss,ground,color='#296E26')
maxn = 100
plt.plot([0,maxn],[0,maxn],color='#55B196')
plt.ylim([0,maxn])
plt.xlim([0,maxn])
plt.ylabel('# ground reads/pore')
plt.xlabel('# flight reads/pore')
#plt.ylabel('# ground reads')
#plt.xlabel('# ISS reads')
plt.savefig('numbers_of_reads_per_pore.pdf',dpi=300,transparent=True,bbox_inches='tight')
