import sys
import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd
import matplotlib.patches as patches
import matplotlib.lines as mlines

qual_scores = {}
boxes = {}
whiskers = {}
for fi in sys.argv[1:]:
   qdist = open(fi,'r').read().split('\n')
   qual_scores[fi] = {score:0 for score in range(100)}
   boxes[fi] = {x:0 for x in [0.25,0.50,0.75]}
   whiskers[fi] = {'min':0,'max':0}
   #qual_scores[fi] = []
   for line in qdist:
      if line != '':
         score,freq = tuple([int(x) for x in line.split('\t')[1:3] if x!=''])
         #print score,freq
         #qual_scores[fi] = qual_scores[fi] + [score]*freq
         qual_scores[fi][score] += freq
   total = sum(qual_scores[fi].values())
   num_bases = 0 
   minimum = None
   maximum = 0
   for score in sorted(qual_scores[fi].keys()):
      num_bases += qual_scores[fi][score]
      for quant in boxes[fi]:
         if num_bases >= total*quant and boxes[fi][quant] == 0:
            boxes[fi][quant] = score
      if minimum == None and qual_scores[fi][score] > 0:
         minimum = score
      if qual_scores[fi][score] > 0 and score > maximum:
         maximum = score
   iqr = boxes[fi][0.75]-boxes[fi][0.25]
   whiskers[fi]['min'] = max(boxes[fi][0.25]-1.5*iqr,minimum)
   whiskers[fi]['max'] = min(boxes[fi][0.75]+1.5*iqr,maximum)

sns.set_style('white')
fig = plt.figure()
ax1 = fig.add_subplot(111)
c = '#086788'
lc = '#8EB9C8'
for i,fi in enumerate(sorted(qual_scores.keys())):
   print fi
   print boxes[fi][0.25], boxes[fi][0.5], boxes[fi][0.75], whiskers[fi]
   ax1.add_patch(
    patches.Rectangle(
        (i+0.7, boxes[fi][0.25]), 0.6, boxes[fi][0.75]-boxes[fi][0.25], #boxes[fi][0.75]/boxes[fi][0.75])
        color = lc 
    )
   )
   hmid = mlines.Line2D([i+0.7,i+1.3],[boxes[fi][0.5],boxes[fi][0.5]],color=c)
   down = mlines.Line2D([i+1,i+1],[whiskers[fi]['min'],boxes[fi][0.25]],color=c)
   hdown = mlines.Line2D([i+0.8,i+1.2],[whiskers[fi]['min'],whiskers[fi]['min']],color=c)
   up = mlines.Line2D([i+1,i+1],[boxes[fi][0.75],whiskers[fi]['max']],color=c)
   hup = mlines.Line2D([i+0.8,i+1.2],[whiskers[fi]['max'],whiskers[fi]['max']],color=c)
   ax1.add_line(hmid)
   ax1.add_line(down)
   ax1.add_line(hdown)
   ax1.add_line(up)
   ax1.add_line(hup)

ax1.set_xlim(0,len(qual_scores.keys())+1)
ax1.set_ylim(0,21)

ticks = ax1.get_xticks()
ax1.set_xticks(range(1,len(qual_scores)+1))
labels = sorted([x.split('_')[0] for x in qual_scores])
print labels
ax1.set_xticklabels(labels)

plt.show()
   
#qualDF = pd.DataFrame(qual_scores.items(),columns=["dataset","score"]) 
#fig = sns.boxplot(x="dataset",y="score",data=qualDF)
fig.savefig("qualscores.pdf",dpi=500,bbox_inches='tight')
