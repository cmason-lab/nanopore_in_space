import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
matplotlib.rcParams['svg.fonttype'] = 'none'

results = open('graphmap_species_counts.tab','r').read().split('\n')
species_list = results[0].split('\t')[1:]
results_dict = {index:[] for index in ['percent','species','sample']}
for result in results[1:]:
        sample = result.split('\t')[0]
        if len(sample.split('_2_temp')) > 1:
            sample = sample.split('_2_template')[0]
        elif len(sample.split('_72ng_temp')) > 1:
            sample = sample.split('_72ng_template')[0]
        else:
            sample = sample.split('_template')[0]
        if len(sample.split('nasa_')) > 1:
            sample = sample.split('nasa_')[1]  
        species_counts = [int(count) for count in result.split('\t')[1:6]]
        tot = sum(species_counts)
        sample = sample + ' (n=' +str(tot)+')'
        species_percents = [count*100./tot for count in species_counts]
        for count,species in zip(species_percents,species_list):
            results_dict['percent'].append(count)
            results_dict['species'].append(' '.join(species.split('_')))
            results_dict['sample'].append(sample)

results_df = pd.DataFrame(results_dict)
results_df

sns.set_style('ticks')
cols = sns.color_palette(['#56d6e1','#2b2381','#ffa6b8','#db1818']) #,'#a34c7c'
sns.set_palette(cols)
g = plt.figure(figsize=(3,3))
sns.barplot(x="species",y="percent",hue="sample",data=results_df,palette=cols)
plt.xticks(rotation=90)
#g.set_xticklabels(rotation=45)
#g.set_titles("{col_name}")
#plt.legend(title='',loc=2)
plt.ylim([0,80])
plt.ylabel("Percent of reads")
#g.set(xticks=[])
#g.add_legend()
#g.fig.subplots_adjust(wspace=.05, hspace=.05)
sns.despine()
for ext in ['png','pdf']:
   plt.savefig('graphmap_results.'+ext,format=ext,dpi=900,bbox_inches='tight',transparent=True)

