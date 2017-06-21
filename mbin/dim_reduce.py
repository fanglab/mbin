import os
import logging
import subprocess
import numpy as np
import shutil
import warnings
import plotting

class dim_reducer():

	def __init__(self, opts, fns, h5_files):
		self.opts     = opts
		self.fns      = fns
		self.h5_files = h5_files

	def run_OS_command( self, CMD ):
		p         = subprocess.Popen(CMD, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
		stdOutErr = p.communicate()
		sts       = p.returncode
		if sts != 0:
			logging.warning("Failed command: %s" % CMD)
		return sts, stdOutErr

	def run_tSNE(self, CMD ):
		sts        = 1
		perplexity = 30
		while sts!=0 and perplexity>1:
			CMD            = CMD.split(" -p ")[0]
			CMD           += " -p %s" % perplexity
			logging.info("Running %s" % CMD)
			sts, stdOutErr = self.run_OS_command( CMD )
			perplexity -= 1
			if sts!=0:
				logging.info("   ...trying again (p = %s)..." % perplexity)
		return sts

	def digest_large_contigs_for_tSNE( self ):
		names      = np.atleast_1d(np.loadtxt(self.fns["contig_names"],     dtype="str"))
		labels     = np.atleast_1d(np.loadtxt(self.fns["contig_labels"],    dtype="str"))
		covs       = np.atleast_1d(np.loadtxt(self.fns["contig_cov"],       dtype="float"))
		lengths    = np.atleast_1d(np.loadtxt(self.fns["contig_lengths"],   dtype="int"))
		SCp        = np.atleast_1d(np.loadtxt(self.fns["contig_SCp"],       dtype="float"))
		SCp_N      = np.atleast_1d(np.loadtxt(self.fns["contig_SCp_N"],     dtype="float"))
		comp       = np.atleast_1d(np.loadtxt(self.fns["contig_comp"],      dtype="float"))

		# if len(SCp.shape)==1:
		# 	SCp   = SCp.reshape(SCp.shape[0],1)
		# 	SCp_N = SCp_N.reshape(SCp_N.shape[0],1)

		# if len(comp.shape)==1:
		# 	comp   = comp.reshape(comp.shape[0],1)

		max_size          = self.opts.subcontig_size
		ref_comp_fn       = "contigs.frag_comp"
		ref_SCp_fn        = "contigs.frag_SCp"
		ref_SCp_N_fn      = "contigs.frag_SCp_N"
		ref_names_fn      = "contigs.frag_names"
		ref_labels_fn     = "contigs.frag_labels"
		ref_covs_fn       = "contigs.frag_covs"
		ref_lengths_fn    = "contigs.frag_lengths"
		f_ref_comp        = open(ref_comp_fn,    "w")
		f_ref_SCp         = open(ref_SCp_fn,     "w")
		f_ref_SCp_N       = open(ref_SCp_N_fn,   "w")
		f_ref_names       = open(ref_names_fn,   "w")
		f_ref_labels      = open(ref_labels_fn,  "w")
		f_ref_covs        = open(ref_covs_fn,    "w")
		f_ref_lengths     = open(ref_lengths_fn, "w")
		for i,length in enumerate(lengths):
			if length>max_size:
				name            = names[i]
				label           = labels[i]
				cov             = covs[i]
				# Copy the original contig SCp and composition vectors; apply to fragments
				if len(SCp.shape)==1:
					contig_SCp_str    = "\t".join(map(lambda x: str(x), SCp[:]))
					contig_SCp_N_str  = "\t".join(map(lambda x: str(x), SCp_N[:]))
					contig_comp_str   = "\t".join(map(lambda x: str(x), comp[:]))
				else: 
					contig_SCp_str    = "\t".join(map(lambda x: str(x), SCp[i,:]))
					contig_SCp_N_str  = "\t".join(map(lambda x: str(x), SCp_N[i,:]))
					contig_comp_str   = "\t".join(map(lambda x: str(x), comp[i,:]))
				
				for j,(refName,contig_seq) in enumerate(self.fasta_iter(self.opts.contigs)):
				
					if refName.find("|quiver")>-1:
						# SMRT assemblies add |quiver to contig names, but this
						# gets dropped from the contig names in the cmp.h5 file.
						refName = refName.replace("|quiver","")
			
					if cmph5_read.slugify(refName) == name:
						seq_chunks = list(self.chunks(contig_seq, max_size))
						# Omit the chunk containing the remainder sequence
						seq_chunks = seq_chunks[:-1]
						for seq_chunk in seq_chunks:
							contig_comp  = []
							contig_kmers = read_scanner.kmer_freq( "cmp", seq_chunk, 0, self.opts )
							for kmer,count in contig_kmers.iteritems():
								kmer_normed_comp = math.log(float(count) / sum(contig_kmers.values()))
								contig_comp.append(kmer_normed_comp)
							# contig_comp_str = "\t".join(map(lambda x: str(round(x,6)), contig_comp))
							f_ref_comp.write(   "%s\n" %  contig_comp_str)
							f_ref_SCp.write(    "%s\n" %  contig_SCp_str)
							f_ref_SCp_N.write(  "%s\n" %  contig_SCp_N_str)
							f_ref_names.write(  "%s\n" %  name)
							f_ref_labels.write( "%s\n" %  label)
							f_ref_covs.write(   "%s\n" %  cov)
							f_ref_lengths.write("%s\n" %  len(seq_chunk))
		f_ref_comp.close()
		f_ref_SCp.close()
		f_ref_SCp_N.close()
		f_ref_names.close()
		f_ref_labels.close()
		f_ref_covs.close()
		f_ref_lengths.close()

		def append_file( orig, new ):
			cat_CMD = "cat %s %s > tmp.appended" % (orig, new)
			sts, stdOutErr = self.run_OS_command( cat_CMD )
			os.rename("tmp.appended", orig)
			os.remove(new)
			
		append_file( self.fns["contig_comp"],    ref_comp_fn )
		append_file( self.fns["contig_SCp"],     ref_SCp_fn )
		append_file( self.fns["contig_SCp_N"],   ref_SCp_N_fn )
		append_file( self.fns["contig_names"],   ref_names_fn )
		append_file( self.fns["contig_labels"],  ref_labels_fn )
		append_file( self.fns["contig_cov"],     ref_covs_fn )
		append_file( self.fns["contig_lengths"], ref_lengths_fn )

	def write_contig_length_labels_file( self, cmph5_file ):
		contigs           = np.atleast_1d(np.loadtxt(self.fns["contig_names"], dtype="str"))
		f                 = open(self.fns["contig_lengths"], "w")
		for contig in contigs:
			if self.opts.sam!=None:
				if contig=="unknown":
					length = 1000
				else:
					length = self.contig_fasta_lens[contig]
			else:
				length = self.opts.cmph5_contig_lens[cmph5_file][contig]
			f.write("%s\n" % length)
		f.close()

	def tsne_reduce(self):
		"""

		"""
		if self.opts.h5_type=="cmp" or self.opts.sam!=None:
		
			self.write_contig_length_labels_file( self.h5_files[0] )

			tsne_sts = {}
			##########################################################
			# Reduce contigs composition dimensionality from ND --> 2D
			##########################################################
			logging.info("Reducing dimensionality of the CONTIGS composition matrix...")
			logging.info("    - digesting large contigs to provide more weight in tSNE")
			self.digest_large_contigs_for_tSNE()
			
			sne_CMD    = "python ~/gitRepo/eBinning/src/bhtsne.py -v -i %s -l %s -o %s -s contigs.comp.scatter.png -z %s -n contigs_composition_digest" %    \
																																	(self.fns["contig_comp"],    \
																																	 self.fns["contig_labels"],  \
																																	 self.fns["contig_comp_2D"], \
																																	 self.fns["contig_lengths"])
			tsne_sts["comp_2D"] = self.run_tSNE(sne_CMD)
			
			##################################################################################
			# Adding contigs coverage column to 2D composition matrix, reducing from 3D --> 2D
			##################################################################################
			logging.info("Adding contig coverage column to 2D composition matrix and reducing to 2D...")
			paste_CMD      = "paste %s %s > %s" % (self.fns["contig_comp_2D"], self.fns["contig_cov"], self.fns["contig_cov_comp_3D"])
			sts, stdOutErr = self.run_OS_command( paste_CMD )
			sne_CMD  = "python ~/gitRepo/eBinning/src/bhtsne.py -v -i %s -l %s -o %s -s contigs.cov_comp.scatter.png -z %s -n contigs_composition_coverage_digest" %    \
																																	(self.fns["contig_cov_comp_3D"], \
																																	 self.fns["contig_labels"],  \
																																	 self.fns["contig_cov_comp_2D"], \
																																	 self.fns["contig_lengths"])
			tsne_sts["cov_comp_2D"] = self.run_tSNE(sne_CMD)

			##################################################
			# Reduce contigs IPD dimensionality from ND --> 2D
			##################################################
			logging.info("Reducing dimensionality of the CONTIGS SCp matrix...")
			sne_CMD = "python ~/gitRepo/eBinning/src/bhtsne.py -v -i %s -l %s -o %s -s contigs.SCp.tSNE.scatter.png -z %s -n contigs_SCp_tSNE_digest" % \
																																	(self.fns["contig_SCp"],          \
																																	 self.fns["contig_labels"],       \
																																	 self.fns["contig_SCp_2D"],       \
																																	 self.fns["contig_lengths"])
			tsne_sts["methyl_2D"] = self.run_tSNE(sne_CMD)
			logging.info("Done.")

			#########################################################################
			# Apply Z-score transformation for combining composition and IPD matrices
			#########################################################################
			if (tsne_sts["methyl_2D"]==0 and tsne_sts["cov_comp_2D"]==0 and tsne_sts["comp_2D"]==0):
				logging.info("Z-score tranforming the 2D matrices...")
				for fn in [self.fns["contig_SCp_2D"], self.fns["contig_cov_comp_2D"], self.fns["contig_comp_2D"]]:
					if fn==self.fns["contig_SCp_2D"] and self.opts.h5_type=="cmp" and self.opts.comp_only:
						shutil.copy(self.fns["contig_SCp"], self.fns["contig_SCp_2D"])
						shutil.copy(self.fns["contig_SCp_2D"], self.fns["contig_SCp_2D"]+".zscores")
					else:
						m         = np.loadtxt(fn, dtype="float")
						m_std     = stats.mstats.zscore(m, axis=0)
						out_fn    = fn + ".zscores"
						np.savetxt(out_fn, m_std, fmt='%.3f', delimiter="\t")
				logging.info("Done.")

				paste_CMD      = "paste %s %s > %s" % (self.fns["contig_cov_comp_2D_z"], self.fns["contig_SCp_2D_z"], self.fns["contig_4D_z"])
				sts, stdOutErr = self.run_OS_command( paste_CMD )
				sne_CMD        = "python ~/gitRepo/eBinning/src/bhtsne.py -v -i %s -l %s -o %s -s contigs.combined.scatter.png -z %s -n contigs_combined" % \
																																		(self.fns["contig_4D_z"],         \
																																		 self.fns["contig_labels"],       \
																																		 self.fns["contig_combo_2D"],     \
																																		 self.fns["contig_lengths"])
				tsne_sts["combined_2D"] = self.run_tSNE(sne_CMD)

			########################################################################
			# Make both digested and undigested versions of the contigs output files
			########################################################################
			def drop_appended_lines( fn, n_to_keep ):
				head_CMD = "head -%s %s > tmp.head" % (n_to_keep, fn)
				sts, stdOutErr = self.run_OS_command( head_CMD )
				os.rename("tmp.head", fn)

			logging.info("...removing the fragmented contig seeds...")
			n_contigs = len(np.atleast_1d(np.loadtxt(self.fns["contig_names"], dtype="str")))

			shutil.copy( self.fns["contig_comp"],          self.fns["contig_comp"]+".digest" )
			shutil.copy( self.fns["contig_names"],         self.fns["contig_names"]+".digest" )
			shutil.copy( self.fns["contig_labels"],        self.fns["contig_labels"]+".digest" )
			shutil.copy( self.fns["contig_cov"],           self.fns["contig_cov"]+".digest" )
			shutil.copy( self.fns["contig_lengths"],       self.fns["contig_lengths"]+".digest" )
			shutil.copy( self.fns["contig_SCp"],           self.fns["contig_SCp"]+".digest" )
			shutil.copy( self.fns["contig_SCp_N"],         self.fns["contig_SCp_N"]+".digest" )
			drop_appended_lines( self.fns["contig_comp"],          n_contigs )
			drop_appended_lines( self.fns["contig_SCp"],           n_contigs )
			drop_appended_lines( self.fns["contig_SCp_N"],         n_contigs )
			drop_appended_lines( self.fns["contig_names"],         n_contigs )
			drop_appended_lines( self.fns["contig_labels"],        n_contigs )
			drop_appended_lines( self.fns["contig_cov"],           n_contigs )
			drop_appended_lines( self.fns["contig_lengths"],       n_contigs )

			if tsne_sts["comp_2D"]==0:
				shutil.copy( self.fns["contig_comp_2D"],       self.fns["contig_comp_2D"]+".digest" )
				shutil.copy( self.fns["contig_comp_2D_z"],     self.fns["contig_comp_2D_z"]+".digest" )
				shutil.copy( self.fns["contig_cov_comp_3D"],   self.fns["contig_cov_comp_3D"]+".digest" )
				drop_appended_lines( self.fns["contig_comp_2D"],       n_contigs )
				drop_appended_lines( self.fns["contig_comp_2D_z"],     n_contigs )
				drop_appended_lines( self.fns["contig_cov_comp_3D"],   n_contigs )
			
			if tsne_sts["cov_comp_2D"]==0:
				shutil.copy( self.fns["contig_cov_comp_2D"],   self.fns["contig_cov_comp_2D"]+".digest" )
				shutil.copy( self.fns["contig_cov_comp_2D_z"], self.fns["contig_cov_comp_2D_z"]+".digest" )
				drop_appended_lines( self.fns["contig_cov_comp_2D"],   n_contigs )
				drop_appended_lines( self.fns["contig_cov_comp_2D_z"], n_contigs )
			
			if tsne_sts["methyl_2D"]==0:
				shutil.copy( self.fns["contig_SCp_2D"],        self.fns["contig_SCp_2D"]+".digest" )
				shutil.copy( self.fns["contig_SCp_2D_z"],      self.fns["contig_SCp_2D_z"]+".digest" )
				drop_appended_lines( self.fns["contig_SCp_2D"],        n_contigs )
				drop_appended_lines( self.fns["contig_SCp_2D_z"],      n_contigs )

			if tsne_sts["cov_comp_2D"]==0 and tsne_sts["methyl_2D"]==0:
				shutil.copy( self.fns["contig_4D_z"],          self.fns["contig_4D_z"]+".digest" )
				drop_appended_lines( self.fns["contig_4D_z"],          n_contigs )
			
			if tsne_sts.get("combined_2D")==0:
				shutil.copy( self.fns["contig_combo_2D"],      self.fns["contig_combo_2D"]+".digest" )
				drop_appended_lines( self.fns["contig_combo_2D"],      n_contigs )			
			
			#########
			warnings.filterwarnings("ignore", category=DeprecationWarning, module="matplotlib")
			#########

			SCp = np.loadtxt(self.fns["contig_SCp"], dtype="float")
			if len(SCp.shape)!=1:
				# There is ore than one contig for plotting
				if tsne_sts["comp_2D"]==0:
					plotting.scatterplot(np.loadtxt(self.fns["contig_comp_2D"],     dtype="float"), np.loadtxt(self.fns["contig_labels"], dtype="str"), "contigs.comp.noDigest.scatter.png",     np.loadtxt(self.fns["contig_lengths"], dtype="int"), "contigs_comp_noDigest")
				if tsne_sts["cov_comp_2D"]==0:
					plotting.scatterplot(np.loadtxt(self.fns["contig_cov_comp_2D"], dtype="float"), np.loadtxt(self.fns["contig_labels"], dtype="str"), "contigs.cov_comp.noDigest.scatter.png", np.loadtxt(self.fns["contig_lengths"], dtype="int"), "contigs_cov_comp_noDigest")
				if tsne_sts["methyl_2D"]==0:
					plotting.scatterplot(np.loadtxt(self.fns["contig_SCp_2D"],      dtype="float"), np.loadtxt(self.fns["contig_labels"], dtype="str"), "contigs.SCp.tSNE.noDigest.scatter.png", np.loadtxt(self.fns["contig_lengths"], dtype="int"), "contigs_SCp_noDigest")
				if tsne_sts.get("combined_2D")==0:
					plotting.scatterplot(np.loadtxt(self.fns["contig_combo_2D"],    dtype="float"), np.loadtxt(self.fns["contig_labels"], dtype="str"), "contigs.combined.noDigest.scatter.png", np.loadtxt(self.fns["contig_lengths"], dtype="int"), "contigs_combined_noDigest")
			else:
				pass
		
		else:
			#######################################################
			# Reduce reads composition dimensionality from ND --> 2D
			#######################################################
			logging.info("Reducing dimensionality of the READS composition matrix...")
			sne_CMD = "python ~/gitRepo/eBinning/src/bhtsne.py -v -i %s -l %s -o %s -s reads.comp.scatter.png -n reads_composition" %                     \
																																	(self.fns["read_comp"],   \
																																	 self.fns["read_labels"], \
																																	 self.fns["read_comp_2D"])
			tsne_sts["comp_2D"] = self.run_tSNE(sne_CMD)

			###############################################
			# Reduce reads IPD dimensionality from ND --> 2D
			###############################################
			logging.info("Reducing dimensionality of the READS SMp matrix...")
			sne_CMD = "python ~/gitRepo/eBinning/src/bhtsne.py -v -i %s -l %s -o %s -s reads.SMp.tSNE.scatter.png -n reads_SMp_tSNE_digest" %   \
																																	(self.fns["read_SMp"],    \
																																	 self.fns["read_labels"], \
																																	 self.fns["read_SMp_2D"])
			tsne_sts["methyl_2D"] = self.run_tSNE(sne_CMD)
		