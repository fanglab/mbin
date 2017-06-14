import numpy as np
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import random

def scatterplot(results, labels, plot_fn, sizes, title):
	if len(labels)>0:
		label_set = list(set(labels))
		label_set.sort()
		colors    = plt.get_cmap('rainbow')(np.linspace(0, 1.0, len(label_set)))
		shapes    = ["o","s","^"]
		fig       = plt.figure(figsize=[15,12])
		ax        = fig.add_axes([0.05, 0.05, 0.8, 0.8])
		ax.axis("off")
		if len(sizes)>0:
			sizes[sizes<100000] = 100000
			scaled_sizes = sizes**1.5 / max(sizes**1.5) * 2000
		res = []
		for k,target_lab in enumerate(label_set):
			if target_lab=="unknown" or target_lab=="Unlabeled":
				target_lab = "Unlabeled"
			
			idxs  = [j for j,label in enumerate(labels) if label==target_lab]
			X     = results[idxs,0]
			Y     = results[idxs,1]
			color = cm.spectral(float(k) / len(label_set))
			for i,x in enumerate(X):
				res.append( (x, Y[i], target_lab, color, shapes[k%len(shapes)]) ) 
			if len(sizes)>0:
				scaled_sizes_idx = np.array(scaled_sizes)[idxs]
				if target_lab=="Unlabeled":
					ax.scatter(X, Y, edgecolors="r", label=target_lab, marker="+", facecolors="r", lw=3, alpha=0.3, s=scaled_sizes_idx)
				else:
					try:
						ax.scatter(X, Y, edgecolors=color, label=target_lab, marker=shapes[k%len(shapes)], facecolors="None", lw=3, alpha=0.7, s=scaled_sizes_idx)
					except:
						ax.scatter(X, Y, edgecolors=color, label=target_lab, marker=shapes[k%len(shapes)], facecolors="None", lw=3, alpha=0.7, s=scaled_sizes_idx)
			else:
				if target_lab=="Unlabeled":
					ax.scatter(X, Y, marker="+", s=15 , edgecolors="grey", label=target_lab, facecolors="grey")
				else:
					ax.scatter(X, Y, marker=shapes[k%len(shapes)], s=15 , edgecolors="None", label=target_lab, facecolors=color)
					# ax.scatter(X, Y, marker="o", s=15 , edgecolors="None", label=target_lab, facecolors="k", alpha=0.5)
		box = ax.get_position()
		ax.set_position([box.x0, box.y0, box.width * 0.7, box.height])
		ax.legend(loc='center left', bbox_to_anchor=(1, 0.5), prop={'size':14}, frameon=False, scatterpoints=1)
		ax.set_title(title)
		plt.savefig(plot_fn)
	else:
		fig       = plt.figure(figsize=[12,12])
		ax        = fig.add_subplot(111)
		if len(sizes)>0:
			sizes[sizes<100000] = 100000
			scaled_sizes = sizes**1.5 / max(sizes**1.5) * 2000
		X     = results[:,0]
		Y     = results[:,1]
		if len(sizes)>0:
			scaled_sizes_idx = np.array(scaled_sizes)[idxs]
			ax.scatter(X, Y, marker="o", lw=3, alpha=0.5, s=scaled_sizes_idx)
		else:
			ax.scatter(X, Y, marker="o", s=15 ,edgecolors="None")
		ax.set_title(title)
		plt.savefig(plot_fn)
