import os,sys
import logging
import optparse
import subprocess
import numpy as np
import shutil
import warnings
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import random
import bhtsne

class Features:
	def __init__( self, opts, mfeats_fn, ofeats_fn ):
		self.opts = opts
		self.mfn  = mfeats_fn
		self.ofn  = ofeats_fn

		if self.mfn.find("contig_methyl_features.txt")>-1:
			self.opts.seq_type = "contig"
		else:
			self.opts.seq_type = "read"

		self.parse_fields()

	def parse_fields( self ):
		"""
		Read in the fields contained in the output files from methylprofiles.
		"""
		m = np.loadtxt(self.mfn, dtype="str", skiprows=1)
		o = np.loadtxt(self.ofn, dtype="str", skiprows=1)
		
		# Optional flags
		m = self.length_filter(m)
		o = self.length_filter(o)
		if self.opts.n_seqs!=None:
			m = self.subsample_seqs( m )
			o = self.subsample_seqs( o )
		
		# Pull out values
		self.ids     = m[:,0].astype("str")
		self.lens    = m[:,1].astype("int")
		self.mscores = m[:,2:].astype("float")
		if self.opts.seq_type=="contig":
			self.covs     = o[:,2].astype("float")
			self.covcomps = o[:,2:].astype("float")
			self.comps    = o[:,3:].astype("float")
		else:
			self.comps    = o[:,2:].astype("float")

	def length_filter ( self, data ):
		"""
		Remove all sequences shorter than --l_min
		"""
		mk = data[:,1].astype("int")>self.opts.l_min

		return data[mk,:]

	def subsample_seqs( self, data ):
		"""
		If --n_seqs specified, downsample to desired number of sequences.
		"""
		mk = random.sample(range(data.shape[0]), self.opts.n_seqs)
		
		return data[mk,:]

	def tsne_reduce( self, data ):
		"""
		Use Barnes-Hut approximation of t-SNE to reduce dimensionality
		of features from N-dim to 2D for plotting.
		"""
		map2d = bhtsne.run_bh_tsne(data, no_dims=self.opts.n_dims, max_iter=self.opts.n_iters, use_pca=False, verbose=True)

		return map2d

	def PCA_reduce( self ):
		"""
		Use principal component analysis to reduce dimensionality
		of features from N-dim to 2D for plotting.
		"""
		# TO BE ADDED

	def plot( self, map2d, plot_fn ):
		"""
		Make scatter plot of 2-dimensional map.
		"""
		def scale_sizes(lengths):
			# Scale marker sizes for plotting
			max_size     = 2500
			min_size     = 5
			norm_lengths = (lengths-min(lengths)) / float((max(lengths)-min(lengths)))
			scaled_sizes = max_size * norm_lengths + min_size
			
			return scaled_sizes

		if self.opts.labels!=None:
			# Read labels into dictionary
			seq_labels_m       = np.loadtxt(self.opts.labels, dtype="str", delimiter="\t")
			seq_names          = seq_labels_m[:,0]
			seq_labs           = seq_labels_m[:,1]
			lab_d              = dict( zip(seq_names, seq_labs) )

			fig                = plt.figure(figsize=[15,12])
			ax                 = fig.add_axes([0.1, 0.1, 0.6, 0.6])
			# ax.axis("off")
			label_list         = list(set(lab_d.values()))
			label_list.sort()
			colors             = plt.get_cmap('spectral')(np.linspace(0.0, 0.95, len(label_list)))
			shapes             = ["o", "v", "^", "s", "D"]
			
			if self.opts.size_markers:
				scaled_sizes = scale_sizes(self.lens)

			res                = []
			for k,target_lab in enumerate(label_list):
				idxs  = [j for j,label in enumerate(seq_labs) if label==target_lab]
				X     = map2d[idxs,0]
				Y     = map2d[idxs,1]
				color = colors[k]
				if self.opts.size_markers:
					scaled_sizes_idx = scaled_sizes[idxs]
					for i,x in enumerate(X):
						res.append( (x, Y[i], target_lab, color, scaled_sizes_idx[i], shapes[k%len(shapes)]) )
				else:
					for i,x in enumerate(X):
						m_size = 15
						res.append( (x, Y[i], target_lab, color, m_size, shapes[k%len(shapes)]) )

			np.random.shuffle(res)
			plotted_labs = set()
			legend_plots = []
			legend_labs  = []
			for i,(x,y, target_lab, color, size, shape) in enumerate(res):
				if target_lab=="Unlabeled":
					shape = "*"
					color = "gray"
				plot = ax.scatter(x,y, marker=shape, s=size ,edgecolors=color, label=target_lab, facecolors="None", linewidth=1.2, alpha=0.8)
				if target_lab not in plotted_labs:
					plotted_labs.add(target_lab)
					legend_plots.append(plot)
					legend_labs.append(target_lab)
			box     = ax.get_position()
			ax.set_position([box.x0, box.y0, box.width * 0.7, box.height])
			leg_tup = zip(legend_labs,legend_plots)
			leg_tup.sort(key=lambda x: x[0])

			# Put "Unlabeled" at end of legend 
			unlabeled = [tup for tup in leg_tup if tup[0]=="Unlabeled"]
			if len(unlabeled)>0:
				unlabeled = unlabeled[0]
				leg_tup.remove(unlabeled)
				leg_tup.append(unlabeled)

			legend_labs, legend_plots = zip(*leg_tup)
			leg = ax.legend(legend_plots, legend_labs, loc='center left', prop={'size':24}, ncol=1, bbox_to_anchor=(1, 0.5), frameon=False, scatterpoints=1)
			for i in range(len(legend_labs)):
				leg.legendHandles[i]._sizes = [100]
			xmin  = min(map2d[:,0])
			xmax  = max(map2d[:,0])
			ymin  = min(map2d[:,1])
			ymax  = max(map2d[:,1])
			xspan = xmax - xmin
			yspan = ymax - ymin
			ax.set_xlim([xmin-(0.05*xspan), xmax+(0.05*xspan)])
			ax.set_ylim([ymin-(0.05*yspan), ymax+(0.05*yspan)])
			plt.savefig(plot_fn)
		else:
			fig = plt.figure(figsize=[12,12])
			ax  = fig.add_subplot(111)
			X   = map2d[:,0]
			Y   = map2d[:,1]
			if self.opts.size_markers:
				scaled_sizes = scale_sizes(self.lens)
				ax.scatter(X, Y, marker="o", lw=3, alpha=0.5, s=scaled_sizes)
			else:
				ax.scatter(X, Y, marker="o", s=15 ,edgecolors="None")
			xmin  = min(map2d[:,0])
			xmax  = max(map2d[:,0])
			ymin  = min(map2d[:,1])
			ymax  = max(map2d[:,1])
			xspan = xmax - xmin
			yspan = ymax - ymin
			ax.set_xlim([xmin-(0.05*xspan), xmax+(0.05*xspan)])
			ax.set_ylim([ymin-(0.05*yspan), ymax+(0.05*yspan)])
			logging.info("  * Saving %s map %s" % (self.opts.dim_reduce, plot_fn))
			plt.savefig(plot_fn)

	def digest_large_seqs( self ):
		"""
		Digests large sequences (probably just contigs) to give them 
		additional weight during t-SNE dimensionality reduction. Otherwise
		a large 3 Mb contig would have the same weight as a tiny 5 kb contig.
		"""
		max_size        = 50000
		self.lens_d     = []
		self.ids_d      = []
		self.mscores_d  = []
		self.comps_d    = []
		self.covs_d     = []
		self.covcomps_d = []
		for i,l in enumerate(self.lens):
		
			# Add original vals for all sequences
			self.lens_d.append(max_size)
			self.ids_d.append(self.ids[i])
			self.mscores_d.append(self.mscores[i,:])
			self.comps_d.append(self.comps[i,:])
			if self.opts.seq_type=="contig":
				self.covs_d.append(self.covs[i])
				self.covcomps_d.append(self.covcomps[i,:])
			
			# For large seqs, add values for digested seqs
			if l>max_size:
				j = 0
				k = 0
				while j < (l-max_size):
					j+=max_size
					self.lens_d.append(max_size)
					self.ids_d.append(self.ids[i]+".digest.%s" % k)
					self.mscores_d.append(self.mscores[i,:])
					self.comps_d.append(self.comps[i,:])
					if self.opts.seq_type=="contig":
						self.covs_d.append(self.covs[i])
						self.covcomps_d.append(self.covcomps[i,:])
					k+=1

		self.lens_d     = np.asarray(self.lens_d)
		self.ids_d      = np.asarray(self.ids_d)
		self.mscores_d  = np.asarray(self.mscores_d)
		self.comps_d    = np.asarray(self.comps_d)
		self.covs_d     = np.asarray(self.comps_d)
		self.covcomps_d = np.asarray(self.covs_d)

	def remove_digested_seqs( self, coords_d_2d ):
		"""
		After the digested sequences have been used to give
		large contigs more weight during t-SNE dimensionality
		reduction, remove those digested contigs so that only
		the original contigs are shown in the plot.
		"""
		mk              = np.ones(len(self.ids_d), dtype=bool) - 1
		digest_idx      = np.where(np.char.find(self.ids_d, '.digest.') < 0)
		mk[digest_idx] += 1
		mk              = mk.astype("bool")
		
		return coords_d_2d[mk]

def launch():
	opts,mfeats_fn,ofeats_fn = __parseArgs()
	__initLog(opts)

	feats = Features(opts, mfeats_fn, ofeats_fn)

	if feats.opts.dim_reduce=="bhtsne":
		# First digest any large contigs to give them
		# more weight during t-SNE clustering
		logging.info("Processing large contigs prior to BH-tSNE")
		feats.digest_large_seqs()
		logging.info("Using BH-tSNE to reduce methylation scores matrix from %s to %s" % (feats.mscores.shape[1], \
					 																	  feats.opts.n_dims))
		mscores2d_d  = feats.tsne_reduce(feats.mscores_d)
		mscores2d    = feats.remove_digested_seqs( mscores2d_d )
		mscores2d_fn = feats.mfn.replace(".txt", ".%s.txt" % feats.opts.dim_reduce)
		np.savetxt(mscores2d_fn, mscores2d, delimiter="\t", fmt="%.5f")
		plot_fn      = feats.mfn.replace(".txt", ".%s.png" % feats.opts.dim_reduce)
		if feats.opts.n_dims==2:
			feats.plot(mscores2d, plot_fn)
		print ""

		logging.info("Using BH-tSNE to reduce methylation scores matrix from %s to %s" % (feats.comps.shape[1], \
					 																	  feats.opts.n_dims))
		comps2d_d  = feats.tsne_reduce(feats.comps_d)
		comps2d    = feats.remove_digested_seqs( comps2d_d )
		comps2d_fn = feats.ofn.replace(".txt", ".comp.%s.txt" % feats.opts.dim_reduce)
		np.savetxt(comps2d_fn, comps2d, delimiter="\t", fmt="%.5f")
		plot_fn    = feats.ofn.replace(".txt", ".comp.%s.png" % feats.opts.dim_reduce)
		if feats.opts.n_dims==2:
			feats.plot(comps2d, plot_fn)
		print ""

		if feats.opts.seq_type=="contig":
			logging.info("Using BH-tSNE to reduce composition matrix (+ coverage column) from %s to %s" % (feats.covcomps.shape[1], \
						 																				   feats.opts.n_dims))
			covcomps2d_d  = feats.tsne_reduce(feats.covcomps_d)
			covcomps2d    = feats.remove_digested_seqs( covcomps2d_d )
			covcomps2d_fn = feats.ofn.replace(".txt", ".covcomp.%s.txt" % feats.opts.dim_reduce)
			np.savetxt(covcomps2d_fn, covcomps2d, delimiter="\t", fmt="%.5f")
			plot_fn       = feats.ofn.replace(".txt", ".covcomp.%s.png" % feats.opts.dim_reduce)
			if feats.opts.n_dims==2:
				feats.plot(covcomps2d, plot_fn)

	elif feats.opts.dim_reduce=="pca":
		raise Exception("PCA not yet supported. Sorry!")

	print >> sys.stderr, "mBin feature mapping has finished running. See log for details."

def __parseArgs():
	"""Handle command line argument parsing"""

	usage = """%prog [--help] [options] <SEQ>_methyl_features.txt <SEQ>_other_features.txt

	Examples:

	mapfeatures -i --labels=<LABELS.txt> --size_markers contig_methyl_features.txt contig_other_features.txt

	mapfeatures -i --labels=<LABELS.txt> --l_min=500 align_methyl_features.txt align_other_features.txt

	mapfeatures -i --l_min=10000 --n_seqs=1000 read_methyl_features.txt read_other_features.txt

	"""

	parser = optparse.OptionParser( usage=usage, description=__doc__ )

	parser.add_option( "-d", "--debug", action="store_true", help="Increase verbosity of logging" )

	parser.add_option( "-i", "--info", action="store_true", help="Add basic logging" )

	parser.add_option( "--logFile", type="str", help="Write logging to file [log.controls]" )
	
	parser.add_option( "--prefix", type="str", help="Prefix to use for output files [None]" )

	parser.add_option( "--size_markers", action="store_true", help="Adjust marker size in plot according to sequence length [False]" )
	
	parser.add_option( "--dim_reduce", type="str", help="Dimensionality reduction algorithm to apply (bhtsne or pca) [bhtsne]" )
	
	parser.add_option( "--labels", type="str", help="Tab-delimited file (no header) of sequence labels (seq_name\tlabel_name) [None]" )

	parser.add_option( "--l_min", type="int", help="Minimum read length to include for analysis [0]" )

	parser.add_option( "--n_seqs", type="int", help="Number of sequences to subsample [all]" )
	
	parser.add_option( "--n_dims", type="int", help="Number of dimensions to reduce to for visualization (only n_dims=2 will be plotted) [2]" )
	
	parser.add_option( "--n_iters", type="int", help="Number of iterations to use for BH-tSNE [500]" )
	
	parser.set_defaults( logFile="log.mapfeatures",          \
						 info=False,                         \
						 debug=False,                        \
						 prefix="",                          \
						 size_markers=False,                 \
						 dim_reduce="bhtsne",                \
						 labels=None,                        \
						 l_min=0,                            \
						 n_seqs=None,                        \
						 n_dims=2,                           \
						 n_iters=500,                        \
						 )

	opts, args = parser.parse_args( )

	mfeats_fn,ofeats_fn = __check_input( opts, args, parser )

	opts.comp_only     = False
	opts.sam           = None
	opts.skip_motifs   = None
	opts.bas_whitelist = None

	if opts.prefix != "":
		opts.prefix+="_"		
	
	return opts,mfeats_fn,ofeats_fn

def __initLog( opts ):
	"""Sets up logging based on command line arguments. Allows for three levels of logging:
	logging.error( ): always emitted
	logging.info( ) : emitted with --info or --debug
	logging.debug( ): only with --debug"""

	if os.path.exists(opts.logFile):
		os.remove(opts.logFile)

	logLevel = logging.DEBUG if opts.debug \
				else logging.INFO if opts.info \
				else logging.ERROR

	logger = logging.getLogger("")
	logger.setLevel(logLevel)
	
	# create file handler which logs even debug messages
	fh = logging.FileHandler(opts.logFile)
	fh.setLevel(logLevel)
	
	# create console handler with a higher log level
	ch = logging.StreamHandler()
	ch.setLevel(logLevel)
	
	# create formatter and add it to the handlers
	logFormat = "%(asctime)s [%(levelname)s] %(message)s"
	formatter = logging.Formatter(logFormat, "%Y-%m-%d %H:%M:%S")
	ch.setFormatter(formatter)
	fh.setFormatter(formatter)
	
	# add the handlers to logger
	logger.addHandler(ch)
	logger.addHandler(fh)

def __check_input( opts, args, parser ):
	"""
	Make sure the input is in the form of either a cmp.h5 file of aligned reads
	or a FOFN of unaligned bas.h5 files. Also make sure that a reference fasta 
	file is specified if 
	"""
	if len(args)!=2:
		print "ERROR -- expecting two arguments: \
				 (1) <SEQ>_methyl_features.txt output from methylprofiles containing methylation features for mapping \
				 (2) <SEQ>_other_features.txt output from methylprofiles containing alternative sequence features for mapping"

	mfeats_fn  = args[0]
	ofeats_fn  = args[1]
	feature_type = None

	if not os.path.exists(mfeats_fn):
		parser.error("Can't find file of sequence features (methylprofiles output) for mapping: %s" % mfeats_fn)

	if not os.path.exists(ofeats_fn):
		parser.error("Can't find file of sequence features (methylprofiles output) for mapping: %s" % ofeats_fn)

	return mfeats_fn, ofeats_fn

if __name__ == "__main__":
	main()