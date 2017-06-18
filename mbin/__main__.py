import os,sys
import optparse
import logging
from mbin import mbin_runner

def main(args=None):
	"""The main routine."""

	opts,args = __parseArgs()
	__initLog(opts)

	mBin = mbin_runner(opts,args)
	mBin.run()

	print >> sys.stderr, "mBin has finished running. See log for details."

def __parseArgs():
	"""Handle command line argument parsing"""

	usage = """%prog [--help] [options]

	EXAMPLES 

	For contigs:

	mbin -i --procs=4 --h5_type=cmp --h5=cmp.h5.fofn --contigs=polished_assembly.fasta --real_mix --control_dir=./control_data

	where cmp.h5.fofn has the format:
	aligned_read.cmp.h5 labelname



	For unaligned reads:

	mbin -i --procs=4 --h5_type=bas --h5=bas.h5.fofn --readlength_min=20000 --real_mix --control_dir=./control_data --motifs_file=motifs.txt

	where bas.h5.fofn has the format:
	smrtcell1.bas.h5 labelname
	smrtcell2.bas.h5 labelname
	smrtcell3.bas.h5 labelname
	smrtcell4.bas.h5 labelname

	"""

	parser = optparse.OptionParser( usage=usage, description=__doc__ )

	parser.add_option( "-d", "--debug", action="store_true", help="Increase verbosity of logging" )

	parser.add_option( "-i", "--info", action="store_true", help="Add basic logging" )

	parser.add_option( "--logFile", type="str", help="Write logging to file [log.out]" )

	parser.add_option( "--subreadlength_min", type="int", help="Minimum subread length to include for analysis [100]" )

	parser.add_option( "--readlength_min", type="int", help="Minimum read length to include for analysis [100]" )

	parser.add_option( "--readlength_max", type="int", help="Maximum read length to include for analysis [10000000]" )

	parser.add_option( "--minQV", type="float", help="If base has QV < minQV, do not include [0]" )

	parser.add_option( "--min_IPD", type="float", help="If motif IPD is < min_IPD, set to zero [1.0]" )

	parser.add_option( "--min_pct", type="float", help="Remove motifs if they are significant in < min_pct of all reads [5.0]" )

	parser.add_option( "--min_kmer", type="int", help="Minimum motif size to scan (contiguous motifs) [4]" )

	parser.add_option( "--max_kmer", type="int", help="Maximum motif size to scan (contiguous motifs) [6]" )

	parser.add_option( "--mod_bases", type="str", help="String containing bases to query for mods ['A']" )

	parser.add_option( "--comp_kmer", type="int", help="Kmer size to use for sequence composition measurements [5]" )

	parser.add_option( "--minAcc", type="float", help="Min subread accuracy of read [0.8]" )

	parser.add_option( "--minMapQV", type="float", help="Min mapping QV of aligned read [240]" )

	parser.add_option( "--minReadScore", type="float", help="Min read score of an unaligned read [0.0]" )

	parser.add_option( "--maxPausiness", type="float", help="Max pausiness value of an unaligned read [1000]" )

	parser.add_option( "--bipartite", action="store_true", help="Search bipartite motifs also [False]" )

	parser.add_option( "--comp_only", action="store_true", help="Only run composition binning (no eBinning) [False]" )

	parser.add_option( "--minContigLength", type="int", help="Min length of contig to consider [10000]" )

	parser.add_option( "--minContigForMotifs", type="int", help="Min length of contig to use when getting top motifs [50000]" )

	parser.add_option( "--get_cmph5_stats", action="store_true", help="Generate read stats and exit [False]" )

	parser.add_option( "--skip_motifs", type="str", help="File containing specific motifs to SKIP [None]" )

	parser.add_option( "--motifs_file", type="str", help="File containing specific motifs to include [None]" )

	parser.add_option( "--tmp", type="str", help="Directory where numerous temporary files will be written [tmp]" )

	parser.add_option( "--procs", type="int", help="Number of cores to use [4]" )

	parser.add_option( "--N_reads", type="int", help="Number of qualifying reads to include in analysis [1000000000]" )

	parser.add_option( "--h5", type="str", help="Path to the cmp.h5 file to analyze (for --h5_type=cmp) or a fofn with a bas.h5 file on each line (for --h5_type=bas) [None]" )

	parser.add_option( "--h5_type", type="str", help="cmp (aligned reads) or bas (unaligned) [cmp]" )
	
	parser.add_option( "--contigs", type="str", help="Fasta file containing entries for the assembled contigs [None]" )
	
	parser.add_option( "--bas_whitelist", type="str", help="File containing reads to use in the read-level binning [None]" )
	
	parser.add_option( "--subcontig_size", type="int", help="Size to decompose contigs into for tSNE weighting [50000]" )
	
	parser.add_option( "--min_motif_N", type="int", help="Min number of IPDs from a read/contig to keep for motif filtering [20]" )
	
	parser.add_option( "--min_motif_reads", type="int", help="Min number of reads with motif hits to keep for motif filtering [20]" )
	
	parser.add_option( "--N_motif_reads", type="int", help="Number of reads to scan for motif discovery procedure (will use N_reads \
															if smaller than N_motif_reads) [20000]" )
	
	parser.add_option( "--minMotifIPD", type="float", help="Min motif contig IPD for inclusion of motif in final set [1.7]" )
	
	parser.add_option( "--printIPDs", action="store_true", help="Print out IPD values for SCp ROC curve [False]" )
	
	parser.add_option( "--subtract_control", type="str", help="Subtract control IPDs in final calculations [True]" )
	
	parser.add_option( "--sam", type="str", help="Path to SAM file with CCS alignments for assigning read-level IPDs to contigs [None]" )
	
	parser.add_option( "--control_dir", type="str", help="Location to write control IPD data from WGA sequencing (or where to retrieve previously written control data) [None]" )
	
	parser.add_option( "--cross_cov_bins", type="str", help="Path to file containing binning results from CONCOCT. Will use to \
															 improve motif discovery. Only works with contig-level analysis \
															 (cmp.h5 input) inputs. File format should be '<contig_name>,<bin_id>' \
															 [None]" )
	
	parser.add_option( "--motif_discov_only", action="store_true", help="Quit pipeline after motif discovery [False]" )
	parser.set_defaults( logFile="log.out",                  \
						 info=False,                         \
						 debug=False,                        \
						 printIPDs=False,                    \
						 subtract_control=True,              \
						 subreadlength_min=100,              \
						 readlength_min=100,                 \
						 readlength_max=10000000,            \
						 minQV=0.0,                          \
						 min_IPD=1.0,                        \
						 min_pct=5.0,                        \
						 min_kmer=4,                         \
						 max_kmer=6,                         \
						 comp_kmer=5,                        \
						 mod_bases="A",                      \
						 minAcc=0.8,                         \
						 minMapQV=240,                       \
						 minReadScore=0.0,                   \
						 maxPausiness=1000,                  \
						 minContigLength=10000,              \
						 minContigForMotifs=25000,           \
						 minMotifIPD=1.7,                    \
						 min_motif_N=20,                     \
						 min_motif_reads=20,                 \
						 N_motif_reads=20000,                \
						 get_cmph5_stats=False,              \
						 bipartite=False,                    \
						 comp_only=False,                    \
						 skip_motifs=None,                   \
						 motifs_file=None,                   \
						 tmp="tmp",                          \
						 procs=4,                            \
						 N_reads=1000000000,                 \
						 h5=None,                            \
						 h5_type=None,                       \
						 bas_whitelist=None,                 \
						 contigs=None,                       \
						 sam=None,                           \
						 control_dir="Not_set",              \
						 cross_cov_bins=None,                \
						 motif_discov_only=False,            \
						 subcontig_size=50000)

	opts, args = parser.parse_args( )

	opts.bipart_config = [(3,4), (5,6), (3,4)]

	if (opts.h5_type!="cmp" and opts.h5_type!="bas"):
		parser.error("Please specify --h5_type to be either cmp (aligned reads) or bas (unaligned reads)!")

	__check_input( opts, args, parser )
	
	if opts.control_dir=="Not_set":
		parser.error("Please specify where to build control IPD data (or where an existing build exists)!")

	opts.control_dir = os.path.abspath(opts.control_dir)

	if opts.sam!=None and opts.motifs_file==None:
		parser.error("Use of SAM file for read<-->contig mapping only supported with input motifs file using --motifs_file!" )

	if opts.h5_type=="bas" and opts.cross_cov_bins!=None:
		parser.error("Use of the --cross_cov_bins option is not compatible with bas.h5 inputs!")

	
	#####################
	# Point to a whole-genome amplified or other 
	# methylation-free files to use as controls
	# 
	# wga_cmp_h5 = an aligned_reads.cmp.h5 file of WGA alignments
	# wga_bas_h5 = an unaligned file of WGA reads
	opts.wga_cmp_h5  = "UNSPECIFIED/aligned_reads.cmp.h5"
	opts.wga_bas_h5  = "UNSPECIFIED/m*****.1.bax.h5"
	#####################

	return opts,args

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
	formatter = logging.Formatter(logFormat)
	ch.setFormatter(formatter)
	fh.setFormatter(formatter)
	
	# add the handlers to logger
	logger.addHandler(ch)
	logger.addHandler(fh)

def __check_input( opts, args, parser ):
	input_fn = args[0]

	if (opts.h5_type=="cmp" and input_fn[-6:]!="cmp.h5"):
		parser.error("If a --h5_type=cmp is specified, please input a *.cmp.h5 file of aligned reads")
	else:
		opts.h5_files = [input_fn]

	if (opts.h5_type=="cmp" and opts.contigs==None):
		parser.error("Please specify the fasta file corresponding to the alignments in %s!" % input_fn)

	if opts.h5_type=="bas":
		input_lines = open(input_fn, "rb").read().split("\n")
		input_lines = [l for l in input_lines if l!=""]
		bas_fns     = [l for l in input_lines if l[-6:]=="bas.h5"]
		if len(input_lines)!=len(bas_fns):
			parser.error("If --h5_type=bas is specified, please input a FOFN with a list of *.bas.h5 files to analyze")
		else:
			opts.h5_files = input_lines

	if len(args) != 1:
		parser.error( "Expected 1 argument." )




if __name__ == "__main__":
	main()