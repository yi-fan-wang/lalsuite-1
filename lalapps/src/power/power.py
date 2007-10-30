# $Id$
#
# This program is free software; you can redistribute it and/or modify it under
# the terms of the GNU General Public License as published by the Free Software
# Foundation; either version 2 of the License, or (at your option) any later
# version.
#
# This program is distributed in the hope that it will be useful, but WITHOUT
# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
# FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
# details.
#
# You should have received a copy of the GNU General Public License along with
# this program; if not, write to the Free Software Foundation, Inc., 51
# Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.


#
# =============================================================================
#
#                                   Preamble
#
# =============================================================================
#


"""
Excess power pipeline construction tools.
"""


import errno
import math
import os
import random
import sys
import time


from glue import segments
from glue import segmentsUtils
from glue import pipeline
from glue.lal import CacheEntry
from pylal.date import LIGOTimeGPS
from pylal import burstsearch


__author__ = "Duncan Brown <duncan@gravity.phys.uwm.edu>, Kipp Cannon <kipp@gravity.phys.uwm.edu>"
__date__ = "$Date$"[7:-2]
__version__ = "$Revision$"[11:-2]


#
# =============================================================================
#
#                                   Helpers
#
# =============================================================================
#


def get_universe(config_parser):
	return config_parser.get("condor", "universe")


def get_executable(config_parser, name):
	return config_parser.get("condor", name)


def get_out_dir(config_parser):
	return config_parser.get("pipeline", "out_dir")


def get_cache_dir(config_parser):
	return config_parser.get("pipeline", "cache_dir")


def make_dag_directories(config_parser):
	try:
		os.mkdir(get_cache_dir(config_parser))
	except OSError, e:
		if e.errno != errno.EEXIST:
			# OK if directory exists, otherwise report error
			raise e
	try:
		os.mkdir(get_out_dir(config_parser))
	except OSError, e:
		if e.errno != errno.EEXIST:
			# OK if directory exists, otherwise report error
			raise e


def get_files_per_binjfind(config_parser):
	return config_parser.getint("pipeline", "files_per_binjfind")


def get_files_per_bucluster(config_parser):
	return config_parser.getint("pipeline", "files_per_bucluster")


def get_files_per_bucut(config_parser):
	return config_parser.getint("pipeline", "files_per_bucut")


def get_timing_parameters(config_parser):
	# initialize data structure
	resample_rate = config_parser.getfloat("lalapps_power", "resample-rate")
	params = burstsearch.XLALEPGetTimingParameters(
		window_length = config_parser.getint("lalapps_power", "window-length"),
		max_tile_length = int(config_parser.getfloat("lalapps_power", "max-tile-duration") * resample_rate),
		tile_stride_fraction = config_parser.getfloat("lalapps_power", "tile-stride-fraction"),
		psd_length = config_parser.getint("lalapps_power", "psd-average-points")
	)
	params.resample_rate = resample_rate

	# retrieve additional parameters from config file
	params.filter_corruption = config_parser.getint("lalapps_power", "filter-corruption")
	params.max_tile_bandwidth = config_parser.getfloat("lalapps_power", "max-tile-bandwidth")

	# NOTE:  in previous code, psd_length, window_length, window_shift,
	# filter_corruption, and psd_overlap were all floats in units of
	# seconds
	return params


def make_cache_entry(input_cache, description, path):
	# summarize segment information
	seglists = segments.segmentlistdict()
	for c in input_cache:
		seglists |= c.to_segmentlistdict()

	# obtain instrument list
	instruments = seglists.keys()
	if None in instruments:
		instruments.remove(None)
	instruments.sort()

	# remove empty segment lists to allow extent_all() to work
	for instrument in seglists.keys():
		if not seglists[instrument]:
			del seglists[instrument]

	# make the URL
	if path:
		url = "file://localhost%s" % os.path.abspath(path)
	else:
		url = None

	# construct a cache entry from the instruments and
	# segments that remain
	return CacheEntry("+".join(instruments) or None, description, seglists.extent_all(), url)


def match_nodes_to_caches(nodes, caches, verbose = False):
	"""
	For each cache, get the set of the nodes whose output files it
	contains.  A node is allowed to provide more than one output file,
	and thus can be listed in more than one set.
	"""
	# can't use [set()] * len(caches) for the normal reason
	node_groups = [set() for cache in caches]
	n_found = 0
	for n, node in enumerate(nodes):
		if verbose and not (n % 10):
			print >>sys.stderr, "\t%.1f%%\r" % (100.0 * n / len(nodes)),
		found = False
		for output in node.get_output_cache():
			# find the caches in which this output has been
			# placed, and add the node to the matching group
			for i, cache in enumerate(caches):
				if output in cache:
					node_groups[i].add(node)
					found = True
		if found:
			n_found += 1
	if verbose:
		print >>sys.stderr, "\t100.0%"
	# count the number of nodes whose output files were not found in
	# any caches
	unused = len(nodes) - n_found

	# done
	return node_groups, unused


#
# =============================================================================
#
#                            DAG Node and Job Class
#
# =============================================================================
#


class BurstInjJob(pipeline.CondorDAGJob):
	"""
	A lalapps_binj job used by the power pipeline. The static options
	are read from the [lalapps_binj] section in the ini file. The
	stdout and stderr from the job are directed to the logs directory.
	The job runs in the universe specified in the ini file.  The path
	to the executable is determined from the ini file.
	"""
	def __init__(self, config_parser):
		"""
		config_parser = ConfigParser object
		"""
		pipeline.CondorDAGJob.__init__(self, get_universe(config_parser), get_executable(config_parser, "lalapps_binj"))

		# do this many injections between flow and fhigh inclusively
		self.injection_bands = config_parser.getint("pipeline", "injection_bands")

		self.add_ini_opts(config_parser, "lalapps_binj")

		self.set_stdout_file(os.path.join(get_out_dir(config_parser), "binj-$(macrochannelname)-$(macrogpsstarttime)-$(macrogpsendtime)-$(cluster)-$(process).out"))
		self.set_stderr_file(os.path.join(get_out_dir(config_parser), "binj-$(macrochannelname)-$(macrogpsstarttime)-$(macrogpsendtime)-$(cluster)-$(process).err"))
		self.set_sub_file("lalapps_binj.sub")


class BurstInjNode(pipeline.CondorDAGNode):
	def __init__(self, job):
		pipeline.CondorDAGNode.__init__(self, job)
		self.__usertag = None
		self.output_cache = []

	def set_user_tag(self, tag):
		self.__usertag = tag
		self.add_var_opt("user-tag", self.__usertag)

	def get_user_tag(self):
		if self.output_cache:
			raise AttributeError, "cannot change attributes after computing output cache"
		return self.__usertag

	def set_start(self, start):
		if self.output_cache:
			raise AttributeError, "cannot change attributes after computing output cache"
		self.add_var_opt("gps-start-time", start)

	def set_end(self, end):
		if self.output_cache:
			raise AttributeError, "cannot change attributes after computing output cache"
		self.add_var_opt("gps-end-time", end)

	def get_start(self):
		return self.get_opts().get("macrogpsstarttime", None)

	def get_end(self):
		return self.get_opts().get("macrogpsendtime", None)

	def get_output_cache(self):
		"""
		Returns a LAL cache of the output file names.  This must be
		kept synchronized with the name of the output file in
		binj.c.  Note in particular the calculation of the "start"
		and "duration" parts of the name.
		"""
		if not self.output_cache:
			if not self.get_start() or not self.get_end():
				raise ValueError, "start time or end time has not been set"
			seg = segments.segment(LIGOTimeGPS(self.get_start()), LIGOTimeGPS(self.get_end()))
			if self.__usertag:
				filename = "HL-INJECTIONS_%s-%d-%d.xml" % (self.__usertag, int(self.get_start()), int(self.get_end() - self.get_start()))
			else:
				filename = "HL-INJECTIONS-%d-%d.xml" % (int(self.get_start()), int(self.get_end() - self.get_start()))
			self.output_cache = [CacheEntry("H1+H2+L1", self.__usertag, seg, "file://localhost" + os.path.abspath(filename))]
		return self.output_cache

	def get_output_files(self):
		raise NotImplementedError

	def get_output(self):
		raise NotImplementedError


class PowerJob(pipeline.CondorDAGJob, pipeline.AnalysisJob):
	"""
	A lalapps_power job used by the power pipeline. The static options
	are read from the [lalapps_power] and [lalapps_power_<inst>]
	sections in the ini file. The stdout and stderr from the job are
	directed to the logs directory.  The job runs in the universe
	specified in the ini file. The path to the executable is determined
	from the ini file.
	"""
	def __init__(self, config_parser):
		"""
		config_parser = ConfigParser object
		"""
		pipeline.CondorDAGJob.__init__(self, get_universe(config_parser), get_executable(config_parser, "lalapps_power"))
		pipeline.AnalysisJob.__init__(self, config_parser)
		self.add_ini_opts(config_parser, "lalapps_power")
		self.set_stdout_file(os.path.join(get_out_dir(config_parser), "lalapps_power-$(cluster)-$(process).out"))
		self.set_stderr_file(os.path.join(get_out_dir(config_parser), "lalapps_power-$(cluster)-$(process).err"))
		self.set_sub_file("lalapps_power.sub")


class PowerNode(pipeline.AnalysisNode):
	def __init__(self, job):
		pipeline.CondorDAGNode.__init__(self, job)
		pipeline.AnalysisNode.__init__(self)
		self.__usertag = None
		self.output_cache = []

	def set_ifo(self, instrument):
		"""
		Load additional options from the per-instrument section in
		the config file.
		"""
		if self.output_cache:
			raise AttributeError, "cannot change attributes after computing output cache"
		pipeline.AnalysisNode.set_ifo(self, instrument)
		for optvalue in self.job()._AnalysisJob__cp.items("lalapps_power_%s" % instrument):
			self.add_var_arg("--%s %s" % optvalue)

	def set_user_tag(self, tag):
		if self.output_cache:
			raise AttributeError, "cannot change attributes after computing output cache"
		self.__usertag = tag
		self.add_var_opt("user-tag", self.__usertag)

	def get_user_tag(self):
		return self.__usertag

	def get_output_cache(self):
		"""
		Returns a LAL cache of the output file name.  Calling this
		method also induces the output name to get set, so it must
		be at least once.
		"""
		if not self.output_cache:
			self.output_cache = [CacheEntry(self.get_ifo(), self.__usertag, segments.segment(LIGOTimeGPS(self.get_start()), LIGOTimeGPS(self.get_end())), "file://localhost" + os.path.abspath(self.get_output()))]
		return self.output_cache

	def get_output_files(self):
		raise NotImplementedError

	def get_output(self):
		if self._AnalysisNode__output is None:
			if None in (self.get_start(), self.get_end(), self.get_ifo(), self.__usertag):
				raise ValueError, "start time, end time, ifo, or user tag has not been set"
			seg = segments.segment(LIGOTimeGPS(self.get_start()), LIGOTimeGPS(self.get_end()))
			self.set_output("%s-POWER_%s-%d-%d.xml.gz" % (self.get_ifo(), self.__usertag, int(self.get_start()), int(self.get_end()) - int(self.get_start())))
		return self._AnalysisNode__output

	def set_mdccache(self, file):
		"""
		Set the LAL frame cache to to use. The frame cache is
		passed to the job with the --frame-cache argument.  @param
		file: calibration file to use.
		"""
		self.add_var_opt("mdc-cache", file)
		self.add_input_file(file)

	def set_inspinj(self, file):
		"""
		Set the LAL frame cache to to use. The frame cache is
		passed to the job with the --frame-cache argument.  @param
		file: calibration file to use.
		"""
		self.add_var_opt("inspiralinjection-file", file)
		self.add_input_file(file)

	def set_siminj(self, file):
		"""
		Set the LAL frame cache to to use. The frame cache is
		passed to the job with the --frame-cache argument.  @param
		file: calibration file to use.
		"""
		self.add_var_opt("siminjection-file", file)
		self.add_input_file(file)


class LigolwAddNode(pipeline.LigolwAddNode):
	def __init__(self, *args):
		pipeline.LigolwAddNode.__init__(self, *args)
		self.input_cache = []
		self.output_cache = []

	def __update_output_cache(self):
		del self.output_cache[:]
		self.output_cache.append(make_cache_entry(self.input_cache, None, self._AnalysisNode__output))

	def set_name(self, *args):
		pipeline.LigolwAddNode.set_name(self, *args)
		self.cache_name = os.path.join(self._CondorDAGNode__job.cache_dir, "%s.cache" % self.get_name())
		self.add_var_opt("input-cache", self.cache_name)

	def add_input_cache(self, cache):
		self.input_cache.extend(cache)
		self.__update_output_cache()

	def set_output(self, path):
		pipeline.LigolwAddNode.set_output(self, path)
		self.__update_output_cache()

	def set_output_segment(self, seg):
		for c in self.output_cache:
			c.segment = seg

	def add_preserve_cache(self, cache):
		for c in cache:
			self.add_var_arg("--remove-input-except %s" % c.path())

	def get_input_cache(self):
		return self.input_cache

	def get_output_cache(self):
		return self.output_cache

	def write_input_files(self, *args):
		f = file(self.cache_name, "w")
		for c in self.input_cache:
			print >>f, str(c)
		pipeline.LigolwAddNode.write_input_files(self, *args)

	def get_output_files(self):
		raise NotImplementedError

	def get_output(self):
		raise NotImplementedError


class BucutJob(pipeline.CondorDAGJob):
	def __init__(self, config_parser):
		pipeline.CondorDAGJob.__init__(self, "vanilla", get_executable(config_parser, "ligolw_bucut"))
		self.set_sub_file("ligolw_bucut.sub")
		self.set_stdout_file(os.path.join(get_out_dir(config_parser), "ligolw_bucut-$(cluster)-$(process).out"))
		self.set_stderr_file(os.path.join(get_out_dir(config_parser), "ligolw_bucut-$(cluster)-$(process).err"))
		self.add_condor_cmd("getenv", "True")
		self.add_condor_cmd("Requirements", "Memory > 1100")
		self.add_ini_opts(config_parser, "ligolw_bucut")


class BucutNode(pipeline.CondorDAGNode):
	def __init__(self, *args):
		pipeline.CondorDAGNode.__init__(self, *args)
		self.input_cache = []
		self.output_cache = self.input_cache

	def add_input_cache(self, cache):
		self.input_cache.extend(cache)
		for c in cache:
			filename = c.path()
			pipeline.CondorDAGNode.add_file_arg(self, filename)
			self.add_output_file(filename)

	def add_file_arg(self, filename):
		raise NotImplementedError

	def get_input_cache(self):
		return self.input_cache

	def get_output_cache(self):
		return self.output_cache

	def get_output_files(self):
		raise NotImplementedError

	def get_output(self):
		raise NotImplementedError


class BuclusterJob(pipeline.CondorDAGJob):
	def __init__(self, config_parser):
		pipeline.CondorDAGJob.__init__(self, "vanilla", get_executable(config_parser, "ligolw_bucluster"))
		self.set_sub_file("ligolw_bucluster.sub")
		self.set_stdout_file(os.path.join(get_out_dir(config_parser), "ligolw_bucluster-$(cluster)-$(process).out"))
		self.set_stderr_file(os.path.join(get_out_dir(config_parser), "ligolw_bucluster-$(cluster)-$(process).err"))
		self.add_condor_cmd("getenv", "True")
		self.add_condor_cmd("Requirements", "Memory > 1100")
		self.add_ini_opts(config_parser, "ligolw_bucluster")


class BuclusterNode(pipeline.CondorDAGNode):
	def __init__(self, *args):
		pipeline.CondorDAGNode.__init__(self, *args)
		self.input_cache = []
		self.output_cache = self.input_cache

	def set_name(self, *args):
		pipeline.CondorDAGNode.set_name(self, *args)
		self.cache_name = os.path.join(self._CondorDAGNode__job.cache_dir, "%s.cache" % self.get_name())
		self.add_var_opt("input-cache", self.cache_name)

	def add_input_cache(self, cache):
		self.input_cache.extend(cache)

	def add_file_arg(self, filename):
		raise NotImplementedError

	def write_input_files(self, *args):
		pipeline.CondorDAGNode.write_input_files(self, *args)
		f = file(self.cache_name, "w")
		for c in self.input_cache:
			print >>f, str(c)

	def get_input_cache(self):
		return self.input_cache

	def get_output_cache(self):
		return self.output_cache

	def get_output_files(self):
		raise NotImplementedError

	def get_output(self):
		raise NotImplementedError


class BinjfindJob(pipeline.CondorDAGJob):
	def __init__(self, config_parser):
		pipeline.CondorDAGJob.__init__(self, "vanilla", get_executable(config_parser, "ligolw_binjfind"))
		self.set_sub_file("ligolw_binjfind.sub")
		self.set_stdout_file(os.path.join(get_out_dir(config_parser), "ligolw_binjfind-$(cluster)-$(process).out"))
		self.set_stderr_file(os.path.join(get_out_dir(config_parser), "ligolw_binjfind-$(cluster)-$(process).err"))
		self.add_condor_cmd("getenv", "True")
		self.add_ini_opts(config_parser, "ligolw_binjfind")


class BinjfindNode(pipeline.CondorDAGNode):
	def __init__(self, *args):
		pipeline.CondorDAGNode.__init__(self, *args)
		self.input_cache = []
		self.output_cache = self.input_cache

	def add_input_cache(self, cache):
		self.input_cache.extend(cache)
		for c in cache:
			filename = c.path()
			pipeline.CondorDAGNode.add_file_arg(self, filename)
			self.add_output_file(filename)

	def add_file_arg(self, filename):
		raise NotImplementedError

	def get_input_cache(self):
		return self.input_cache

	def get_output_cache(self):
		return self.output_cache

	def get_output_files(self):
		raise NotImplementedError

	def get_output(self):
		raise NotImplementedError


class BurcaJob(pipeline.CondorDAGJob):
	def __init__(self, config_parser):
		pipeline.CondorDAGJob.__init__(self, "vanilla", get_executable(config_parser, "ligolw_burca"))
		self.set_sub_file("ligolw_burca.sub")
		self.set_stdout_file(os.path.join(get_out_dir(config_parser), "ligolw_burca-$(cluster)-$(process).out"))
		self.set_stderr_file(os.path.join(get_out_dir(config_parser), "ligolw_burca-$(cluster)-$(process).err"))
		self.add_condor_cmd("getenv", "True")
		self.add_condor_cmd("Requirements", "Memory > 1100")
		self.add_ini_opts(config_parser, "ligolw_burca")


class Burca2Job(pipeline.CondorDAGJob):
	def __init__(self, config_parser):
		pipeline.CondorDAGJob.__init__(self, "vanilla", get_executable(config_parser, "ligolw_burca"))
		self.set_sub_file("ligolw_burca2.sub")
		self.set_stdout_file(os.path.join(get_out_dir(config_parser), "ligolw_burca2-$(cluster)-$(process).out"))
		self.set_stderr_file(os.path.join(get_out_dir(config_parser), "ligolw_burca2-$(cluster)-$(process).err"))
		self.add_condor_cmd("getenv", "True")
		self.add_ini_opts(config_parser, "ligolw_burca2")


class BurcaNode(pipeline.CondorDAGNode):
	def __init__(self, *args):
		pipeline.CondorDAGNode.__init__(self, *args)
		self.input_cache = []
		self.output_cache = self.input_cache

	def add_input_cache(self, cache):
		self.input_cache.extend(cache)
		for c in cache:
			filename = c.path()
			pipeline.CondorDAGNode.add_file_arg(self, filename)
			self.add_output_file(filename)

	def add_file_arg(self, filename):
		raise NotImplementedError

	def get_input_cache(self):
		return self.input_cache

	def get_output_cache(self):
		return self.output_cache

	def get_output_files(self):
		raise NotImplementedError

	def get_output(self):
		raise NotImplementedError


class SQLiteJob(pipeline.CondorDAGJob):
	def __init__(self, config_parser):
		pipeline.CondorDAGJob.__init__(self, "vanilla", get_executable(config_parser, "ligolw_sqlite"))
		self.set_sub_file("ligolw_sqlite.sub")
		self.set_stdout_file(os.path.join(get_out_dir(config_parser), "ligolw_sqlite-$(cluster)-$(process).out"))
		self.set_stderr_file(os.path.join(get_out_dir(config_parser), "ligolw_sqlite-$(cluster)-$(process).err"))
		self.add_condor_cmd("getenv", "True")
		self.add_ini_opts(config_parser, "ligolw_sqlite")


class SQLiteNode(pipeline.CondorDAGNode):
	def __init__(self, *args):
		pipeline.CondorDAGNode.__init__(self, *args)
		self.input_cache = []
		self.output_cache = []

	def add_input_cache(self, cache):
		if self.output_cache:
			raise AttributeError, "cannot change attributes after computing output cache"
		self.input_cache.extend(cache)
		for c in cache:
			filename = c.path()
			pipeline.CondorDAGNode.add_file_arg(self, filename)
			self.add_output_file(filename)

	def add_file_arg(self, filename):
		raise NotImplementedError

	def set_output(self, filename):
		if self.output_cache:
			raise AttributeError, "cannot change attributes after computing output cache"
		self.add_macro("macrodatabase", filename)

	def get_input_cache(self):
		return self.input_cache

	def get_output_cache(self):
		if not self.output_cache:
			self.output_cache = [make_cache_entry(self.input_cache, None, self.get_opts()["macrodatabase"])]
		return self.output_cache

	def get_output_files(self):
		raise NotImplementedError

	def get_output(self):
		raise NotImplementedError


class BurcaTailorJob(pipeline.CondorDAGJob):
	def __init__(self, config_parser):
		pipeline.CondorDAGJob.__init__(self, "vanilla", get_executable(config_parser, "ligolw_burca_tailor"))
		self.set_sub_file("ligolw_burca_tailor.sub")
		self.set_stdout_file(os.path.join(get_out_dir(config_parser), "ligolw_burca_tailor-$(cluster)-$(process).out"))
		self.set_stderr_file(os.path.join(get_out_dir(config_parser), "ligolw_burca_tailor-$(cluster)-$(process).err"))
		self.add_condor_cmd("getenv", "True")
		self.add_ini_opts(config_parser, "ligolw_burca_tailor")


class BurcaTailorNode(pipeline.CondorDAGNode):
	def __init__(self, *args):
		pipeline.CondorDAGNode.__init__(self, *args)
		self.input_cache = []
		self.output_cache = []

	def set_name(self, *args):
		pipeline.CondorDAGNode.set_name(self, *args)
		self.cache_name = os.path.join(self._CondorDAGNode__job.cache_dir, "%s.cache" % self.get_name())

	def add_input_cache(self, cache):
		if self.output_cache:
			raise AttributeError, "cannot change attributes after computing output cache"
		self.input_cache.extend(cache)
		for c in cache:
			filename = c.path()
			pipeline.CondorDAGNode.add_file_arg(self, filename)
		self.add_output_file(filename)

	def add_file_arg(self, filename):
		raise NotImplementedError

	def set_output(self, description):
		if self.output_cache:
			raise AttributeError, "cannot change attributes after computing output cache"
		cache_entry = make_cache_entry(self.input_cache, description, "")
		filename = "%s-%s-%d-%d.xml.gz" % (cache_entry.observatory, cache_entry.description, int(cache_entry.segment[0]), int(abs(cache_entry.segment)))
		self.add_var_opt("output", filename)
		cache_entry.url = "file://localhost" + os.path.abspath(filename)
		del self.output_cache[:]
		self.output_cache.append(cache_entry)
		return filename

	def get_input_cache(self):
		return  self.input_cache

	def get_output_cache(self):
		if not self.output_cache:
			raise AttributeError, "must call set_output(description) first"
		return self.output_cache

	def write_input_files(self, *args):
		# oh.  my.  god.  this is fscked.
		for arg in self.get_args():
			if "--add-from-cache" in arg:
				f = file(self.cache_name, "w")
				for c in self.input_cache:
					print >>f, str(c)
				pipeline.CondorDAGNode.write_input_files(self, *args)
				break

	def get_output_files(self):
		raise NotImplementedError

	def get_output(self):
		raise NotImplementedError


#
# =============================================================================
#
#                                DAG Job Types
#
# =============================================================================
#


#
# This is *SUCH* a hack I don't know where to begin.  Please, shoot me.
#


datafindjob = None
binjjob = None
powerjob = None
lladdjob = None
binjfindjob = None
bucutjob = None
buclusterjob = None
burcajob = None
burca2job = None
sqlitejob = None
burcatailorjob = None


def init_job_types(config_parser, types = ["datafind", "binj", "power", "lladd", "binjfind", "bucluster", "bucut", "burca", "burca2", "sqlite", "burcatailor"]):
	"""
	Construct definitions of the submit files.
	"""
	global datafindjob, binjjob, powerjob, lladdjob, binjfindjob, buclusterjob, llb2mjob, bucutjob, burcajob, burca2job, sqlitejob, burcatailorjob

	# LSCdataFind
	if "datafind" in types:
		datafindjob = pipeline.LSCDataFindJob(get_cache_dir(config_parser), get_out_dir(config_parser), config_parser)

	# lalapps_binj
	if "binj" in types:
		binjjob = BurstInjJob(config_parser)

	# lalapps_power
	if "power" in types:
		powerjob = PowerJob(config_parser)

	# ligolw_add
	if "lladd" in types:
		lladdjob = pipeline.LigolwAddJob(get_out_dir(config_parser), config_parser)
		lladdjob.cache_dir = get_cache_dir(config_parser)

	# ligolw_binjfind
	if "binjfind" in types:
		binjfindjob = BinjfindJob(config_parser)
		binjfindjob.files_per_binjfind = get_files_per_binjfind(config_parser)
		if binjfindjob.files_per_binjfind < 1:
			raise ValueError, "files_per_binjfind < 1"

	# ligolw_bucut
	if "bucut" in types:
		bucutjob = BucutJob(config_parser)
		bucutjob.files_per_bucut = get_files_per_bucut(config_parser)
		if bucutjob.files_per_bucut < 1:
			raise ValueError, "files_per_bucut < 1"

	# ligolw_bucluster
	if "bucluster" in types:
		buclusterjob = BuclusterJob(config_parser)
		buclusterjob.files_per_bucluster = get_files_per_bucluster(config_parser)
		if buclusterjob.files_per_bucluster < 1:
			raise ValueError, "files_per_bucluster < 1"
		buclusterjob.cache_dir = get_cache_dir(config_parser)

	# ligolw_burca
	if "burca" in types:
		burcajob = BurcaJob(config_parser)

	# ligolw_burca
	if "burca2" in types:
		burca2job = Burca2Job(config_parser)

	# ligolw_sqlite
	if "sqlite" in types:
		sqlitejob = SQLiteJob(config_parser)

	# ligolw_burca_tailor
	if "burcatailor" in types:
		burcatailorjob = BurcaTailorJob(config_parser)
		burcatailorjob.cache_dir = get_cache_dir(config_parser)


#
# =============================================================================
#
#                                 Segmentation
#
# =============================================================================
#


def psds_from_job_length(timing_params, t):
	"""
	Return the number of PSDs that can fit into a job of length t
	seconds.  In general, the return value is a non-integer.
	"""
	if t < 0:
		raise ValueError, t
	# convert to samples, and remove filter corruption
	t = t * timing_params.resample_rate - 2 * timing_params.filter_corruption
	if t < timing_params.psd_length:
		return 0
	return (t - timing_params.psd_length) / timing_params.psd_shift + 1


def job_length_from_psds(timing_params, psds):
	"""
	From the analysis parameters and a count of PSDs, return the length
	of the job in seconds.
	"""
	if psds < 1:
		raise ValueError, psds
	# number of samples
	result = (psds - 1) * timing_params.psd_shift + timing_params.psd_length
	# add filter corruption
	result += 2 * timing_params.filter_corruption
	# convert to seconds
	return result / timing_params.resample_rate


def split_segment(timing_params, segment, psds_per_job):
	"""
	Split the data segment into correctly-overlaping segments.  We try
	to have the numbers of PSDs in each segment be equal to
	psds_per_job, but with a short segment at the end if needed.
	"""
	# in seconds
	joblength = job_length_from_psds(timing_params, psds_per_job)
	# in samples
	joboverlap = 2 * timing_params.filter_corruption + (timing_params.psd_length - timing_params.psd_shift)
	# in seconds
	joboverlap /= timing_params.resample_rate

	segs = segments.segmentlist()
	t = segment[0]
	while t + joblength <= segment[1]:
		segs.append(segments.segment(t, t + joblength) & segment)
		t += joblength - joboverlap

	extra_psds = int(psds_from_job_length(timing_params, float(segment[1] - t)))
	if extra_psds:
		segs.append(segments.segment(t, t + job_length_from_psds(timing_params, extra_psds)))
	return segs


def segment_ok(timing_params, segment):
	"""
	Return True if the segment can be analyzed using lalapps_power.
	"""
	return psds_from_job_length(timing_params, float(abs(segment))) >= 1.0


#
# =============================================================================
#
#                            Single Node Fragments
#
# =============================================================================
#


datafind_pad = 512


def make_datafind_fragment(dag, instrument, seg):
	node = pipeline.LSCDataFindNode(datafindjob)
	node.set_name("LSCdataFind-%s-%s-%s" % (instrument, int(seg[0]), int(abs(seg))))
	node.set_start(seg[0] - datafind_pad)
	node.set_end(seg[1] + 1)
	# FIXME: argh, shoot me, I need the node to know what instrument
	# it's for, but can't call set_ifo() because that adds a
	# --channel-name command line argument (!?)
	node._AnalysisNode__ifo = instrument
	node.set_observatory(instrument[0])
	if node.get_type() == None:
		node.set_type(datafindjob.get_config_file().get("datafind", "type_%s" % instrument))
	node.set_retry(3)
	dag.add_node(node)
	return set([node])


def make_lladd_fragment(dag, parents, tag, preserve_cache = [], input_cache = None):
	node = LigolwAddNode(lladdjob)

	# collect input cache entries.  if input_cache is None, use the
	# parents' output caches as the source of input files
	parent_output_cache = []
	for parent in parents:
		node.add_parent(parent)
		parent_output_cache += parent.get_output_cache()
	if input_cache is not None:
		node.add_input_cache(input_cache)
	else:
		node.add_input_cache(parent.get_output_cache())

	# construct a name for the output file
	cache_entry = make_cache_entry(node.get_input_cache(), None, None)
	node.set_name("lladd_%s_%s_%d_%d" % (tag, cache_entry.observatory, int(cache_entry.segment[0]), int(abs(cache_entry.segment))))
	node.set_output("%s-%s-%d-%d.xml" % (cache_entry.observatory, tag, int(cache_entry.segment[0]), int(abs(cache_entry.segment))))
	node.set_output_segment(cache_entry.segment)

	node.add_preserve_cache(preserve_cache)
	node.set_retry(3)
	dag.add_node(node)
	# NOTE:  code that calls this generally requires a single node to
	# be returned.  if this function is ever modified to return more
	# than one node, check and fix all calling codes.
	return set([node])


def make_power_fragment(dag, parents, instrument, seg, tag, framecache, injargs = {}):
	node = PowerNode(powerjob)
	node.set_name("lalapps_power_%s_%s_%s_%s" % (instrument, tag, int(seg[0]), int(abs(seg))))
	map(node.add_parent, parents)
	node.set_cache(framecache)
	node.set_ifo(instrument)
	node.set_start(seg[0])
	node.set_end(seg[1])
	node.set_user_tag(tag)
	for arg, value in injargs.iteritems():
		# this is a hack, but I can't be bothered
		node.add_var_arg("--%s %s" % (arg, value))
	dag.add_node(node)
	return set([node])


def make_binj_fragment(dag, seg, tag, offset, flow, fhigh):
	# determine frequency ratio from number of injections across band
	# (take a bit off to allow fhigh to be picked up despite round-off
	# errors)
	fratio = 0.9999999 * (fhigh / flow) ** (1.0 / binjjob.injection_bands)

	# one injection every time-step / pi seconds, taking
	# injection_bands steps across the frequency range, this is how
	# often the sequence "repeats"
	period = binjjob.injection_bands * float(binjjob.get_opts()["time-step"]) / math.pi

	# adjust start time to be commensurate with injection period
	start = seg[0] - seg[0] % period + period * offset

	node = BurstInjNode(binjjob)
	node.set_start(start)
	node.set_end(seg[1])
	node.set_name("lalapps_binj_%d_%d" % (int(start), int(flow)))
	node.set_user_tag(tag)
	node.add_macro("macroflow", flow)
	node.add_macro("macrofhigh", fhigh)
	node.add_macro("macrofratio", fratio)
	node.add_macro("macroseed", int(time.time() + start))
	dag.add_node(node)
	return set([node])


def make_binjfind_fragment(dag, parents, tag, verbose = False):
	input_cache = [cache_entry for parent in parents for cache_entry in parent.get_output_cache()]
	input_cache.sort(lambda a, b: cmp(a.segment, b.segment))
	nodes = []
	while input_cache:
		node = BinjfindNode(binjfindjob)
		node.set_name("ligolw_binjfind_%s_%d" % (tag, len(nodes)))
		node.add_input_cache(input_cache[:binjfindjob.files_per_binjfind])
		del input_cache[:binjfindjob.files_per_binjfind]
		node.add_macro("macrocomment", tag)
		dag.add_node(node)
		nodes.append(node)
	parent_groups, unused = match_nodes_to_caches(parents, [node.get_input_cache() for node in nodes], verbose = verbose)
	if unused:
		# impossible
		raise Error
	for node, parents in zip(nodes, parent_groups):
		for parent in parents:
			node.add_parent(parent)
	return set(nodes)


def make_bucluster_fragment(dag, parents, tag, verbose = False):
	input_cache = [cache_entry for parent in parents for cache_entry in parent.get_output_cache()]
	input_cache.sort(lambda a, b: cmp(a.segment, b.segment))
	nodes = []
	while input_cache:
		node = BuclusterNode(buclusterjob)
		node.set_name("ligolw_bucluster_%s_%d" % (tag, len(nodes)))
		node.add_input_cache(input_cache[:buclusterjob.files_per_bucluster])
		del input_cache[:buclusterjob.files_per_bucluster]
		node.add_macro("macrocomment", tag)
		node.set_retry(3)
		dag.add_node(node)
		nodes.append(node)
	parent_groups, unused = match_nodes_to_caches(parents, [node.get_input_cache() for node in nodes], verbose = verbose)
	if unused:
		# impossible
		raise Error
	for node, parents in zip(nodes, parent_groups):
		for parent in parents:
			node.add_parent(parent)
	return set(nodes)


def make_bucut_fragment(dag, parents, tag, verbose = False):
	input_cache = [cache_entry for parent in parents for cache_entry in parent.get_output_cache()]
	input_cache.sort(lambda a, b: cmp(a.segment, b.segment))
	nodes = []
	while input_cache:
		node = BucutNode(bucutjob)
		node.set_name("ligolw_bucut_%s_%d" % (tag, len(nodes)))
		node.add_input_cache(input_cache[:bucutjob.files_per_bucut])
		del input_cache[:bucutjob.files_per_bucut]
		node.add_macro("macrocomment", tag)
		dag.add_node(node)
		nodes.append(node)
	parent_groups, unused = match_nodes_to_caches(parents, [node.get_input_cache() for node in nodes], verbose = verbose)
	if unused:
		# impossible
		raise Error
	for node, parents in zip(nodes, parent_groups):
		for parent in parents:
			node.add_parent(parent)
	return set(nodes)


def make_burca_fragment(dag, parents, tag):
	parents = list(parents)
	nodes = set()
	while parents:
		parent = parents.pop()
		input_cache = parent.get_output_cache()
		seg = make_cache_entry(input_cache, None, None).segment
		node = BurcaNode(burcajob)
		node.set_name("ligolw_burca_%s_%d_%d" % (tag, int(seg[0]), int(abs(seg))))
		node.add_parent(parent)
		node.add_input_cache(input_cache)
		node.add_macro("macrocomment", tag)
		dag.add_node(node)
		nodes.add(node)
	return nodes


def make_burca_tailor_fragment(dag, input_cache, seg, tag):
	input_cache = list(input_cache)
	input_cache.sort(reverse = True)
	nodes = []
	while input_cache:
		node = BurcaTailorNode(burcatailorjob)
		node.set_name("ligolw_burca_tailor_%s_%d_%d_%d" % (tag, int(seg[0]), int(abs(seg)), len(nodes)))
		node.add_input_cache(input_cache[-5:])
		del input_cache[-5:]
		node.set_output(tag)
		dag.add_node(node)
		nodes.append(node)
	node = BurcaTailorNode(burcatailorjob)
	node.set_name("ligolw_burca_tailor_%s_%d_%d" % (tag, int(seg[0]), int(abs(seg))))
	for parent in nodes:
		node.add_parent(parent)
		node.add_input_cache(parent.get_output_cache())
	del node.get_args()[:]
	node.add_var_arg("--add-from-cache %s" % node.cache_name)
	node.set_output(tag)
	node.set_post_script("/bin/rm -f %s" % " ".join([c.path() for c in node.get_input_cache()]))
	dag.add_node(node)
	return [node]


def make_burca2_fragment(dag, parents, input_cache, tag):
	# note that the likelihood data collection step immediately
	# preceding this is a choke point in the pipeline, where all jobs
	# must complete before anything else can begin, so there is no
	# advantage in grouping files in time order at this stage.
	# therefore, to load balance the jobs, we randomize the input
	# files.
	input_cache = list(input_cache)
	random.shuffle(input_cache)
	nodes = []
	# FIXME:  this whole function just isn't the way I want this to
	# work.  I want to provide DAG nodes, not input files...
	while input_cache:
		node = BurcaNode(burca2job)
		node.set_name("ligolw_burca2_%s_%d" % (tag, len(nodes)))
		node.add_macro("macrocomment", tag)
		node.add_input_cache(input_cache[-10:])
		del input_cache[-10:]
		for parent in parents:
			node.add_parent(parent)
		nodes.append(node)
		dag.add_node(node)
	return nodes


#
# =============================================================================
#
#                              LSCdataFind Stage
#
# =============================================================================
#


def make_datafind_stage(dag, seglists, verbose = False):
	if verbose:
		print >>sys.stderr, "building LSCdataFind jobs ..."

	#
	# Fill gaps smaller than the padding added to each datafind job.
	# Filling in the gaps ensures that exactly 1 datafind job is
	# suitable for each lalapps_power job, and also hugely reduces the
	# number of LSCdataFind nodes in the DAG.
	#

	filled = seglists.copy().protract(datafind_pad / 2).contract(datafind_pad / 2)

	#
	# Build the nodes.  Do this in time order to assist depth-first job
	# submission on clusters.
	#

	segs = [(seg, instrument) for instrument, seglist in filled.iteritems() for seg in seglist]
	segs.sort()

	nodes = set()
	for seg, instrument in segs:
		if verbose:
			print >>sys.stderr, "making datafind job for %s spanning %s" % (instrument, seg)
		new_nodes = make_datafind_fragment(dag, instrument, seg)
		nodes |= new_nodes

		# add a post script to check the file list
		required_segs_string = ",".join(segmentsUtils.to_range_strings(seglists[instrument] & segments.segmentlist([seg])))
		for node in new_nodes:
			node.set_post_script(datafindjob.get_config_file().get("condor", "LSCdataFindcheck") + " --dagman-return $RETURN --stat --gps-segment-list %s %s" % (required_segs_string, node.get_output()))

	return nodes


#
# =============================================================================
#
#         DAG Fragment Combining Multiple lalapps_binj With ligolw_add
#
# =============================================================================
#


def make_multibinj_fragment(dag, seg, tag):
	flow = float(powerjob.get_opts()["low-freq-cutoff"])
	fhigh = flow + float(powerjob.get_opts()["bandwidth"])

	nodes = make_binj_fragment(dag, seg, tag, 0.0, flow, fhigh)
	return make_lladd_fragment(dag, nodes, tag)


#
# =============================================================================
#
#            Analyze One Segment Using Multiple lalapps_power Jobs
#
# =============================================================================
#


#
# Without injections
#


def make_power_segment_fragment(dag, datafindnodes, instrument, segment, tag, timing_params, psds_per_job, verbose = False):
	"""
	Construct a DAG fragment for an entire segment, splitting the
	segment into multiple power jobs.
	"""
	# only one frame cache file can be provided as input
	# the unpacking indirectly tests that the file count is correct
	[framecache] = [node.get_output() for node in datafindnodes]
	seglist = split_segment(timing_params, segment, psds_per_job)
	if verbose:
		print >>sys.stderr, "Segment split: " + str(seglist)
	nodes = set()
	for seg in seglist:
		nodes |= make_power_fragment(dag, datafindnodes, instrument, seg, tag, framecache)
	return nodes


#
# With injections
#


def make_injection_segment_fragment(dag, datafindnodes, binjnodes, instrument, segment, tag, timing_params, psds_per_job, verbose = False):
	# only one frame cache file can be provided as input, and only one
	# injection description file can be provided as input
	# the unpacking indirectly tests that the file count is correct
	[framecache] = [node.get_output() for node in datafindnodes]
	[simfile] = [cache_entry.path() for node in binjnodes for cache_entry in node.get_output_cache()]
	seglist = split_segment(timing_params, segment, psds_per_job)
	if verbose:
		print >>sys.stderr, "Injections split: " + str(seglist)
	nodes = set()
	for seg in seglist:
		nodes |= make_power_fragment(dag, datafindnodes | binjnodes, instrument, seg, tag, framecache, injargs = {"burstinjection-file": simfile})
	return nodes


#
# =============================================================================
#
#        Analyze All Segments in a segmentlistdict Using lalapps_power
#
# =============================================================================
#


#
# Without injections
#


def make_single_instrument_stage(dag, datafinds, seglistdict, tag, timing_params, psds_per_job, verbose = False):
	nodes = []
	for instrument, seglist in seglistdict.iteritems():
		for seg in seglist:
			if verbose:
				print >>sys.stderr, "generating %s fragment %s" % (instrument, str(seg))

			# find the datafind job this power job is going to
			# need
			dfnodes = set([node for node in datafinds if (node.get_ifo() == instrument) and (seg in segments.segment(node.get_start(), node.get_end()))])
			if len(dfnodes) != 1:
				raise ValueError, "error, not exactly 1 datafind is suitable for power job at %s in %s" % (str(seg), instrument)

			# power jobs
			nodes += make_power_segment_fragment(dag, dfnodes, instrument, seg, tag, timing_params, psds_per_job, verbose = verbose)

	# done
	return nodes


#
# With injections
#


def make_single_instrument_injections_stage(dag, datafinds, binjnodes, seglistdict, tag, timing_params, psds_per_job, verbose = False):
	nodes = []
	for instrument, seglist in seglistdict.iteritems():
		for seg in seglist:
			if verbose:
				print >>sys.stderr, "generating %s fragment %s" % (instrument, str(seg))

			# find the datafind job this power job is going to
			# need
			dfnodes = set([node for node in datafinds if (node.get_ifo() == instrument) and (seg in segments.segment(node.get_start(), node.get_end()))])
			if len(dfnodes) != 1:
				raise ValueError, "error, not exactly 1 datafind is suitable for power job at %s in %s" % (str(seg), instrument)

			# power jobs
			nodes += make_injection_segment_fragment(dag, dfnodes, binjnodes, instrument, seg, tag, timing_params, psds_per_job, verbose = verbose)

	# done
	return nodes


#
# =============================================================================
#
#       The Coincidence Fragment:  bucluster + lladd + bucluster + burca
#
# =============================================================================
#


def make_pre_coinc_lladd(dag, parents, tag, input_cache):
	return make_lladd_fragment(dag, parents, tag, input_cache = input_cache, preserve_cache = input_cache)


#
# =============================================================================
#
#                             ligolw_sqlite stage
#
# =============================================================================
#


def make_sqlite_stage(dag, parents, tag, verbose = False):
	if verbose:
		print >>sys.stderr, "generating ligolw_sqlite stage for tag %s ..." % tag
	nodes = []
	for parent in parents:
		for cache_entry in parent.get_output_cache():
			node = SQLiteNode(sqlitejob)
			node.set_name("ligolw_sqlite_%s_%d" % (tag, len(nodes)))
			node.add_parent(parent)
			node.add_input_cache([cache_entry])
			node.set_output(cache_entry.path().replace(".xml.gz", ".sqlite"))
			node.set_retry(3)
			dag.add_node(node)
			nodes.append(node)
	return nodes
