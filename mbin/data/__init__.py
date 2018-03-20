from pkg_resources import Requirement, resource_filename

CMPH5 = {"aligned_reads.cmp.h5" : "pHel3.fa"}

BASH5 = {"m130522_092457_42208_c100497142550000001823078008081323_s1_p0.bas.h5" : \
				["m130522_092457_42208_c100497142550000001823078008081323_s1_p0.1.bax.h5", \
				 "m130522_092457_42208_c100497142550000001823078008081323_s1_p0.2.bax.h5", \
				 "m130522_092457_42208_c100497142550000001823078008081323_s1_p0.3.bax.h5"]
				}

MOTIFS_FN    = "motifs.txt"

CONTROL_DIR  = "controls"

class Opts:
	def __init__(self):
		self.sam=None
		self.logFile='log.out'
		self.info=False
		self.debug=False
		self.minAcc=0.8
		self.h5=None, 
		self.subcontig_size=50000
		self.maxPausiness=1000
		self.min_motif_N=20
		self.minContigForMotifs=25000
		self.readlength_max=10000000
		self.contigs='pHel3.fa'
		self.minMotifIPD=1.7
		self.minQV=0.0
		self.subreadlength_min=100
		self.N_reads=1000000000
		self.h5_type='cmp'
		self.readlength_min=100
		self.min_kmer=4
		self.max_kmer=6
		self.mod_bases='A'
		self.N_motif_reads=20000
		self.minContigLength=0
		self.min_motif_reads=20
		self.comp_kmer=5
		self.minMapQV=240
		self.minReadScore=0.0
		self.procs=1
		self.tmp='tmp'
		self.min_IPD=1.0
		self.subtract_control=True
		self.comp_only=False
		self.motif_discov_only=False
		self.cross_cov_bins=None
		self.bipartite=True
		self.skip_motifs=None
		self.printIPDs=False
		self.min_pct=5.0
		self.bas_whitelist=None
		self.get_cmph5_stats=False
		self.tsne=False
		self.bipart_config=[(3,4), (5,6), (3,4)]
		self.wga_cmp_h5=getCmpH5()["cmph5"]
		self.motifs_file=getMotifsFile()
		self.control_dir=getControl()

def _getAbsPath(fname):
	return resource_filename(Requirement.parse('mbin'),'mbin/data/%s' % fname)

def getCmpH5():
	"""
	Return full path to the cmp.h5 and associated reference
	"""
	return [{'cmph5' : _getAbsPath(cmph5), 'fasta':  _getAbsPath(fasta)} for cmph5,fasta in CMPH5.items()][0]

def getBasH5():
	"""
	Return full path to the bas.h5 and associated bax.h5 files
	"""
	return [{'bash5' : _getAbsPath(bash5), 'baxh5':  map(_getAbsPath, baxh5)} for bash5,baxh5 in BASH5.items()][0]

def getMotifsFile():
	"""
	Return full path to file containing motifs to check
	"""
	return _getAbsPath(MOTIFS_FN)

def getControl():
	"""
	Return location of control data to be used for tests
	"""
	return _getAbsPath(CONTROL_DIR)

def set_opts():
	"""
	Set some options for the test call to mbin
	"""
	opts = Opts()
	return opts