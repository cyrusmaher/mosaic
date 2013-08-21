# The MIT License (MIT)
#
# Copyright (c) 2013 M. Cyrus Maher
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# all copies or substantial portions of the Software.
#
# The above copyright notice and this permission notice shall be included in
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
# THE SOFTWARE.

import sys
from commands import getoutput, getstatusoutput
from Bio import AlignIO
import tempfile
import numpy as np
import os
import random


def pairwisealign(seq1, seq2, **kwargs):
    """
    Globally align two sequences.

    :param seq1: Sequence 1
    :type seq1: str
    :param seq2: Sequence 2
    :type seq2: str

    :parm AA: True if protein sequences, false otherwise.
    :param gapopen: The cost for opening a gap.
    :param gapextend: The cost for extending a gap.

    :returns: Sequence 1 aligned, sequence 2 aligned
    :rtype: tuple
    """

    if kwargs['AA']:
        flag1 = '-sprotein1'
        flag2 = '-sprotein2'
        matrix = 'EBLOSUM62'
    else:
        flag1 = '-snucleotide1'
        flag2 = '-snucleotide2'
        matrix = 'EDNAFULL'

    outfile = tempfile.NamedTemporaryFile()

    callstr = 'stretcher -outfile=%s -asequence=asis:%s \
        -bsequence=asis:%s -gapopen=%s -gapextend=%s -aformat fasta \
        -datafile %s %s %s' % (outfile.name, seq1,
                               seq2, kwargs['gapopen'], kwargs['gapextend'],
                               matrix, flag1, flag2)

    status, output = getstatusoutput(callstr)

    if status == 0:
        result = AlignIO.parse(outfile, 'fasta')
        alignmentobj = result.next()
        outfile.close()
        getoutput('rm %s' % outfile.name)
    else:
        print output
        print "There was an error in pairwisealign.\
                Is stretcher installed? See above output and check out", \
            outfile.name
        print callstr
        outfile.close()
        raise SystemError

    return alignmentobj[0].seq._data, alignmentobj[1].seq._data


class MultiMSA:
    """
    Data object for containing sets of multiple sequence alignments.

    :param MSAlist1: A list of multiple sequence alignments
        (class Bio.Align.MultipleSeqAlignment, read in with Bio.AlignIO.read)
    :param MSAlist2: The corresponding AA or DNA alignments (optional) \n
    :param MSAspec1: A list of lists specifying the order of species in each
        MSAlist1 alignment(This or specfunc1 must be set)
    :param MSAspec2: Same as above for MSAlist2.
        If neither MSAspec2 or specfunc2
        are set, species order is assumed to be the same as MSAlist1
    :param specfunc1: A parsing function that returns species names from
        sequence labels (This or MSAspec1 must be set)
    :param specfunc2: Same as above for MSAlist2.
        If neither MSAspec2 or specfunc2 are set,
        species order is assumed to be the same as MSAlist1
    :param methodnames1: The names of the methods used corresponding
        to alignments in MSAlist1 (optional)
    :param methodnames2: The names of the methods used corresponding
        to alignments in MSAlist2 (optional)
    :param specorder: The desired order of species for downstream output.

    :returns: A MultiMSA object.

    """

    def __init__(self, MSAlist1, MSAlist2=None, MSAspec1=None, MSAspec2=None,
                 specfunc1=None, specfunc2=None, methodnames1=None,
                 methodnames2=None, specorder=None):
        assert MSAlist1, 'MSAlist1 is empty!!'
        assert MSAspec1 or specfunc1, \
            'Mosaic needs to know the order of species in MSA 1.' \
            'Please either supply a ordered list of species\
            or a function to obtain this information from record labels.'

        if MSAspec2:
            assert MSAspec2 or specfunc2 or specfunc1, \
                'Mosaic needs to know the order of species in the MSA 2.' \
                'Please either supply a ordered list of species\
                or a function to obtain this information from record labels.'

        self.MSAlist1 = MSAlist1
        self.MSAlist2 = MSAlist2
        self.catcount = 0
        self.defaultnames = True if methodnames1 is None else False

        self.methodnames1 = ['Met%s' % xx
                             for xx in range(1, len(MSAlist1) + 1)] \
            if methodnames1 is None else methodnames1

        self.methodnames2 = ['Met%s' % xx
                             for xx in range(1, len(MSAlist2) + 1)] \
            if ((methodnames2 is None)
                and (MSAlist2 is not None)) \
            else methodnames2

        self.nummethods = len(self.methodnames1)
        self.specorder = specorder

        if self.methodnames2:
            assert set(self.methodnames2) == set(self.methodnames1), \
                'Alignment sets must include the same set of methods.'

        if MSAspec1:
            self.MSAlist1_specs = MSAspec1
        else:
            self.MSAlist1_specs = [[specfunc1(record.id) for record in MSA]
                                   for MSA in MSAlist1]

        if MSAlist2:
            if MSAspec2:
                self.MSAlist2_specs = MSAspec2
            elif specfunc2:
                self.MSAlist2_specs = [[specfunc2(record.id)
                                        for record in MSA] for MSA in MSAlist2]
            elif specfunc1:
                self.MSAlist2_specs = [[specfunc1(record.id)
                                        for record in MSA] for MSA in MSAlist2]
        self.storeMSAdicts()

    def _storeMSAdict(self, methodnames, MSAlist_specs, MSAlist):
        metdict = {}
        for metcount, met in enumerate(methodnames):
            specs = MSAlist_specs[metcount]
            MSA = MSAlist[metcount]
            metdict[met] = {}
            for speccount, spec in enumerate(specs):
                metdict[met][spec] = MSA._records[speccount]
                metdict[met][spec].description = spec
        return metdict

    def storeMSAdicts(self):
        self.MSAdict1 = self._storeMSAdict(self.methodnames1,
                                           self.MSAlist1_specs, self.MSAlist1)
        if self.MSAlist2:
            self.MSAdict2 = self._storeMSAdict(self.methodnames2,
                                               self.MSAlist2_specs,
                                               self.MSAlist2)
        else:
            self.MSAdict2 = {}

    def _strn(self, thedict, startstr=''):
        for xx in thedict:
            specs = thedict[xx].keys()
            specs.sort()
            startstr += '\tFor %s, alignments of length %s returned for %s\n'\
                        % (xx, len(thedict[xx][specs[0]]), specs)
        return startstr

    def __repr__(self):
        str1 = self._strn(self.MSAdict1, startstr='Set 1:\n') if self.MSAdict1\
            else 'ERROR! No alignments included for set 1.'
        if self.MSAlist2:
            str2 = self._strn(self.MSAdict2, startstr='Set 2:\n') \
                if self.MSAdict2 \
                else 'No alignments included for set 2.'

            return '%s\n%s' % (str1, str2)
        return str1

    def __add__(self, other):
        for xx in other.methodnames1:
            assert xx not in self.methodnames1, \
                'ERROR: Catenated multiMSAs' \
                'should not have overlapping method names!'

        self.MSAdict1.update(self._storeMSAdict(other.methodnames1,
                                                other.MSAlist1_specs,
                                                other.MSAlist1))

        self.methodnames1 = self.methodnames1 + other.methodnames1

        self.MSAlist1_specs = self.MSAlist1_specs + other.MSAlist1_specs

        if self.MSAlist2 and other.MSAlist2:
            self.MSAdict2.update(self._storeMSAdict(other.methodnames2,
                                                    other.MSAlist2_specs,
                                                    other.MSAlist2))

            self.methodnames2 = self.methodnames2 + other.methodnames2
            self.MSAlist2_specs = self.MSAlist2_specs + other.MSAlist2_specs

        elif other.MSAlist2:
            self.MSAdict2 = other.MSAdict2
            self.methodnames2 = other.methodnames2
            self.MSAlist2_specs = other.MSAlist2_specs

        self.nummethods = len(self.methodnames1)

        return self


class Mosaic:
    """
    Allows for integration of sets of multiple sequence alignments


    :param multiMSA: a MultiMSA object
    :param ref: The species to use as the reference
        to anchor the sequence cluster.
    :param useonlyspec: If specified, use only the supplied subset of species.
    :param speccutoffs: Pass a dictionary of cutoffs. Corresponds to
                species-specific cutoffs.
    :param edgefunc: The function to calculate similarities between species.
            Can be 'perID', 'bitscore', or a
            user-specified function of the form: func(seq1, seq2, \*\*kwargs)
    :param optrule ('pairwise'): The rule for optimization. Can be 'pairwise',
            'toref', or a user-specified function of the form:
            ``func(edgeweightmatrix, \*\*kwargs)``
    :param ignoregaps: Ignore gaps in alignment scoring?
    :param AA: True if the primary MSA set is amino acid. False otherwise.

    :param customoptfunc: Custom optimization function
    :param scoremat: Custom scoring matrix for 'bitscore'
    :param stretcher_gapopen:
        Gap opening penalty for global pairwise sequence alignment.
    :param stretcher_gapextend:
        Gap extension penalty for global pairwise sequence alignment.

    :returns: Mosaic object

    """

    def __init__(self, multiMSA, ref, useonlyspec=None, speccutoffs=None,
                 edgefunc='perID', optrule='pairwise', ignoregaps=False,
                 customoptfunc=None, AA=True, scoremat=None,
                 stretcher_gapopen=8, stretcher_gapextend=1,
                 similaritythresh=-1e6):

        self.multiMSA = multiMSA
        self.speccutoffs = speccutoffs
        self.thresh = similaritythresh
        self.AA = AA
        self.optrule = optrule

        self.ref = ref
        self.ignoregaps = ignoregaps

        self.allspecs = set(self.getallspecs())

        if useonlyspec:
            ## All species that were present and requested.
            self.allspecs = self.allspecs & set(useonlyspec)

        ## Make sure reference species comes first in order.
        self.allspecs = (self.ref,) + tuple(self.allspecs - set([self.ref]))
        self.methods = self.multiMSA.methodnames1
        self.nummethods = self.multiMSA.nummethods

        self.specmethoddict = None
        self.customoptfunc = customoptfunc

        self.stretcher_gapopen = stretcher_gapopen
        self.stretcher_gapextend = stretcher_gapextend

        if edgefunc == 'perID':
            self.edgefunc = self.getperID
        elif edgefunc == 'bitscore':
            self.edgefunc = self.getbitscore
            import pandas as pd

            if scoremat:
                self.scoremat = pd.read_csv(scoremat,
                                            header=True,
                                            index_col=0, sep="\s+")
            elif AA:
                self.scoremat = pd.read_csv(
                    os.path.join(os.path.dirname(__file__),
                                 'BLOSUM62.txt'),
                    header=0, index_col=0, sep="\s+")
            else:
                self.scoremat = pd.read_csv(os.path.join(
                    os.path.dirname(__file__),
                    'EDNAFULL.txt'),
                                            header=0, index_col=0, sep="\s+")

        elif hasattr(edgefunc, '__call__'):
            self.edgefunc = edgefunc
        else:
            raise AssertionError('Unrecognized value for edge function. '
                                 'You supplied %s' % edgefunc)

        assert ref is not None, \
            "A reference species is required in order to anchor the cluster."

        if ignoregaps:
            assert ref, 'A reference must be defined in order to ignore gaps\
                        (multiple edges between nodes are not implemented).'

        assert optrule == 'pairwise' or optrule == 'toref', \
            'ERROR: unrecognized optimization rule.\
             Supported rules are "toref" and "pairwise"'
        if useonlyspec:
            try:
                len(useonlyspec)
            except:
                raise AssertionError('useonlyspec should be an iterable\
                 containing the desired species')

        self.calc_edgeweights()
        self.optimize_cluster()

    def getbitscore(self, seq1, seq2):
        """
        Calculate a bitscore between two (unaligned) sequences.

        :param seq1: Sequence 1
        :param seq2: Sequence 2

        :var self.AA: True if sequence alphabet is amino acids,
            false otherwise.
        :var self.stretcher_gapopen: The penalty for opening
            a gap in the alignment.
        :var self.stretcher_gapextend: The penalty for extending
            a gap in the alignment.
        :var self.scoremat: The score matrix to use for
            the scoring of the alignment.

        :return: A bit score for the aligned sequences.

        .. note::
            `pandas <http://pandas.pydata.org/>`_
            is required to manage scoring matrices.

        """
        seq1_aligned, seq2_aligned = pairwisealign(seq1, seq2,
                                           AA=self.AA,
                                           gapopen=self.stretcher_gapopen,
                                           gapextend=self.stretcher_gapextend)

        score = 0
        for xx, yy in zip(seq1_aligned, seq2_aligned):
            xx = '*' if xx == '-' else xx
            yy = '*' if yy == '-' else yy
            score += self.scoremat.ix[xx, yy]
        return 1.0 * score / len(seq1_aligned)

    def getperID(self, seq1, seq2):
        """
        Calculate percent identity between two (unaligned) sequences.

        :param seq1: Sequence 1
        :param seq2: Sequence 2

        :var self.AA: True if sequence alphabet is amino acids,
            false otherwise.
        :var self.ignoregaps: True if gaps in the first (reference)
            sequence are to be ignored.

        :return: A percent identity for the aligned sequences.
        """

        seq1_aligned, seq2_aligned = pairwisealign(seq1, seq2,
                                            AA=self.AA,
                                            gapopen=self.stretcher_gapopen,
                                            gapextend=self.stretcher_gapextend)

        seq1 = np.array([ss for ss in seq1_aligned])
        seq2 = np.array([ss for ss in seq2_aligned])

        ## You only want to count valid matches
        if self.AA:
            seq1_notmissing = np.logical_and((seq1 != 'U'), (seq1 != 'X'))
            seq2_notmissing = np.logical_and((seq2 != 'U'), (seq2 != 'X'))
            if self.ignoregaps:
                seq1_notmissing = np.logical_and((seq1 != '-'),
                                                 seq1_notmissing)
        else:
            seq1_notmissing = np.logical_and((seq1 != '?'), (seq1 != 'N'))
            seq2_notmissing = np.logical_and((seq2 != '?'), (seq2 != 'N'))
            if self.ignoregaps:
                seq1_notmissing = np.logical_and((seq1 != '-'),
                                                 seq1_notmissing)

        pairmatches = (seq1 == seq2)

        ## Get the percent of bases present in human
        ## that match in comparison species.
        qualitymatches = np.logical_and(pairmatches, seq1_notmissing)
        qualitycomparisons = np.logical_and(seq1_notmissing, seq2_notmissing)

        return 1.0 * qualitymatches.sum() / (qualitycomparisons.sum())

    def getallspecs(self):
        allspecs = set()
        for mm in self.multiMSA.MSAlist1_specs:
            allspecs = allspecs | set(mm)
        return list(allspecs)

    def indextomet_spec(self, xx):
        met = self.multiMSA.methodnames1[xx % self.nummethods]
        spec = self.allspecs[xx / self.nummethods]
        return met, spec

    def toindex(self, mnum, snum):
        return snum * self.nummethods + mnum

    def index_fromnames(self, method, spec):
        mnum = (np.array(self.methods) == method).nonzero()[0]
        snum = (np.array(self.allspecs) == spec).nonzero()[0]
        return self.toindex(mnum, snum)

    def indices_fromnames(self, method, spec):
        mnum = (np.array(self.methods) == method).nonzero()[0]
        snum = (np.array(self.allspecs) == spec).nonzero()[0]
        return mnum, snum

    def calc_sim_mat_pairwise(self):
        """
        Internal function: For pairwise optimization,
            calculate the similarity matrix that
            will define the cluster of sequences.

        :param self.multiMSA.MSAdict1: A dictionary of
            dictionaries of sequences.
        :param self.multiMSA.methodnames1:
            A list of the names of the methods producing each MSA.
        :param self.edgefunc:
            The function used to calculate
            the similarity between two sequences.

        :return: A (nmethods*nspec) x (nmethods*nspec) matrix of edgeweights,
            stored to ``self.edgeweights_pairwise``

        .. note::
            | Sequences are blocked by method
                (according to ``self.methodnames``).
            | These blocks are ordered by the specified species order
                (``self.allspecs``).

        """

        nn = (len(self.allspecs) - 1) * self.nummethods
        self.edgeweights_pairwise = np.zeros((nn, nn)) - np.inf

        xcount = 0
        for s_x in self.allspecs[1:]:
            for m_x in self.multiMSA.methodnames1:
                entry_x = self.multiMSA.MSAdict1[m_x]
                seq_x = entry_x[s_x] if s_x in entry_x else None

                ycount = 0
                for s_y in self.allspecs[1:]:
                    for m_y in self.multiMSA.methodnames1:
                        if (s_x == s_y):
                            ycount += 1
                            continue
                        elif ycount <= xcount:
                            ycount += 1
                            continue

                        entry_y = self.multiMSA.MSAdict1[m_y]
                        seq_y = entry_y[s_y] if s_y in entry_y else None

                        if seq_x and seq_y:
                            print "Calculating edge weights:", \
                                '(%s, %s vs. %s, %s)' % (m_x, s_x, m_y, s_y)

                            self.edgeweights_pairwise[xcount, ycount] = \
                                self.edgefunc(seq_x.seq._data, seq_y.seq._data)

                            print "\tValue is", \
                                self.edgeweights_pairwise[xcount, ycount]

                        ycount += 1
                xcount += 1

    def calc_sim_mat_toref(self):
        """
        Internal function: For "to reference" optimization,
            calculate the similarity vector
            that will relate each sequence to the reference.

        :param self.multiMSA.MSAdict1:
            A dictionary of dictionaries of sequences.
        :param self.multiMSA.methodnames1:
            A list of the names of the methods producing each MSA.
        :param self.edgefunc: The function used to calculate
            the similarity between two sequences.
        :return: A (nspec) x (nspec) matrix of edgeweights,
            stored to ``self.edgeweights_toref``
        .. note::
            This is the stage at which filtering takes place.
            Any sequence below the similarity cutoff is not assigned
            to the ``self.edgeweights_toref`` matrix.
        """

        self.edgeweights_toref = np.zeros((len(self.allspecs),
                                           self.nummethods)) - np.inf

        met1 = self.multiMSA.methodnames1[0]

        seq_x = self.multiMSA.MSAdict1[met1][self.ref]

        for scount, s_y in enumerate(self.allspecs):
            for mcount, m_y in enumerate(self.multiMSA.methodnames1):

                if s_y == self.ref:
                    continue

                entry_y = self.multiMSA.MSAdict1[m_y]
                seq_y = entry_y[s_y] if s_y in entry_y else None

                if seq_x and seq_y:
                    print "Calculating edge weights:", scount, mcount, \
                        '(%s, %s vs. %s, %s)' % (met1, self.ref, m_y, s_y)

                    score = self.edgefunc(seq_x.seq._data, seq_y.seq._data)

                    if self.speccutoffs and score < self.speccutoffs[s_y]:
                        print "Filtering sequence for", s_y, m_y
                        continue

                    self.edgeweights_toref[scount, mcount] = score

                    print "\tValue is", \
                        self.edgeweights_toref[scount, mcount]

        self.allspecs = [s_y for scount, s_y in enumerate(self.allspecs)
                         if s_y == self.ref or
                            np.any(
                                self.edgeweights_toref[scount, :]
                                > self.thresh)]

    def calc_edgeweights(self):
        self.calc_sim_mat_toref()
        if self.optrule == 'pairwise':
            self.calc_sim_mat_pairwise()

    def opt_cluster_toref(self):
        """
        Internal function: Takes sequences from each species with
            the highest similarity to the reference.
        """

        retvals = {}
        # Columns are methods, rows are species.
        # You want the indices of the row maxes.
        for cc in range(len(self.allspecs)):

            if np.any(self.edgeweights_toref[cc, :] > self.thresh):
                yindex = np.argmax(self.edgeweights_toref[cc, :])
                retvals[self.allspecs[cc]] = \
                    (self.multiMSA.methodnames1[yindex],
                     cc,
                     yindex
                    )
            else:
                retvals[self.allspecs[cc]] = None
        self.specmethoddict = retvals

    def gatherpairwisedifs(self, spec, method, index):
        mnum, snum = self.indices_fromnames(method, spec)

        dist = 0
        dd = self.edgeweights_toref[snum, mnum]

        if dd < self.thresh:
            print 'Why is distance to reference unknown for %s, %s?' \
                  % (spec, method)
            raise IndexError

        dist += dd

        for ss in self.specmethoddict:
            if ss != spec:
                i2 = self.specmethoddict[ss][1]
                if i2 < index:
                    dd = self.edgeweights_pairwise[i2, index]
                else:
                    dd = self.edgeweights_pairwise[index, i2]

                if dd < self.thresh:
                    print 'Why is distance to bestspec unknown?\
                        (Index is %s vs %s)' % (i2, index)
                    print 'This corresponds to %s, %s vs %s, %s' \
                          % (ss, self.specmethoddict[ss][0], spec, method)
                    print self.indextomet_spec(i2), self.indextomet_spec(index)
                    raise IndexError
                dist += dd
        return dist

    def initdict_pairwise(self):
        self.specmethoddict = {}
        ## Initialize return dictionary.
        for scount, ss in enumerate(self.allspecs[1:]):
            ## Make sure you initialize with sequences that exist!
            for mcount, mm in enumerate(self.methods):
                index = self.toindex(mcount, scount)

                ## Don't even try to include this species/method
                ## if it has no sequences present that are
                ## above the similarity cutoff
                if self.edgeweights_toref[scount + 1, mcount] > self.thresh:
                    self.specmethoddict[ss] = [mm, index, np.inf]
                    break

        ## Now go back through and fix pairwise dist numbers for each sequence
        for ss in self.specmethoddict:
            index = self.specmethoddict[ss][1]
            self.specmethoddict[ss][2] = \
                self.gatherpairwisedifs(ss, self.specmethoddict[ss][0],
                                        index)

    def optloop_pairwise(self):
        """
        Internal function: optimize cluster using pairwise similarities
            and Gibbs sampling.
        """

        print "optimizing pairwise"
        count = 0
        specnums = range(1, len(self.allspecs))
        metnums = range(len(self.methods))

        while True:
            print "Iteration", count + 1

            nochange = 0
            random.shuffle(specnums)  # pick a random species order
            for scount in specnums:
                ss = self.allspecs[scount]
                met_old = self.specmethoddict[ss][0]  # initialize variables
                met_new = met_old

                random.shuffle(metnums)  # pick a random method order
                for mcount in metnums:
                    mm = self.methods[mcount]

                    ## The second condition below
                    ## skips methods that have been filtered
                    if (mm != met_old and
                            (self.edgeweights_toref[scount, mcount]
                                 > self.thresh)):

                        index = self.toindex(mcount, scount - 1)

                        if np.all(self.edgeweights_pairwise[:, index]
                                < self.thresh):
                            continue

                        dist_new = self.gatherpairwisedifs(ss, mm, index)

                        if dist_new > self.specmethoddict[ss][2]:
                            met_new = mm
                            print "updating best seq for %s from %s to %s " \
                                  % (ss, met_old, met_new)
                            self.specmethoddict[ss] = (mm, index, dist_new)

                if met_new == met_old:
                    nochange += 1

            if nochange == len(specnums):
                break

            count += 1

        print "Optimization complete after %s iterations!!!" % (count + 1)

    def opt_cluster_pairwise(self):
        self.initdict_pairwise()
        self.optloop_pairwise()

    def optimize_cluster(self):
        """
        Internal function: optimize the sequence cluster.

        If ``self.optrule`` is 'pairwise', optimize cluster by picking the
            sequence for each species
            that optimizes the pairwise distance to current best sequences.
            This is repeated cyclically until convergence is reached.

        If ``self.optrule`` is 'toref', pick the sequence from each species
            that is most similar to the reference sequence

        If ``self.optrule`` is 'custom', apply the function defined
            in ``self.customoptfunc`` to the ``self.edgeweights_toref``
            and/or ``self.edgeweights_pairwise`` matrices

        """

        if self.optrule == 'pairwise':
            self.opt_cluster_pairwise()
        elif self.optrule == 'toref':
            self.opt_cluster_toref()
        elif self.optrule == 'custom':
            self.specmethoddict = self.customoptfunc()
        else:
            raise AssertionError(
                'Did not recognize cluster optimization rule: %s.'
                'Available options are: pairwise, toref, and custom.'
                % self.optrule)

    def _writefunc(self, filename, MSAdict, toloop, labelfunc,
                   inclspec, inclmet):

        print "Writing %s" % filename
        with open(filename, 'w') as f1:
            for ss in toloop:
                method = self.methods[0] if ss == self.ref \
                    else self.specmethoddict[ss][0]

                seq = MSAdict[method][ss]
                if labelfunc:
                    label = labelfunc(seq.name, ss, method)
                else:
                    label = '%s' % seq.name
                    if inclspec and inclmet:
                        label += ' (%s, %s)' % (ss, method)
                    elif inclspec:
                        label += ' (%s)' % (method)
                    elif inclmet:
                        label += ' (%s)' % (method)
                f1.write('>%s\n%s\n' % (label, seq.seq._data.replace('-', '')))

    def write_unaligned(self, filename1, filename2=None, inclspec=False,
                        inclmet=False, specorder=None, labelfunc=None):
        """
        Write (unaligned) optimal sequences to a file.

        :param filename1: The file to which to write the primary MSAlist
        :param filename2: The file to which to write the secondary MSAlist
        :param inclspec: Whether to include the species name
            in the sequence labels
        :param inclmet: Whether to include the method name
            in the sequence labels
        :param specorder: If specified, a different order
            for species in the output
        :param labelfunc: If specified, a function to output sequence labels.
            Should be of the form labelfunc(seq.name, species, method)
        """

        self.unalignedfile1 = filename1
        self.unalignedfile2 = filename2

        toloop = self.allspecs if specorder is None else \
            [xx for xx in specorder if xx in self.allspecs]

        if len(toloop) <= 1:
            print "No species to write!"
            self.unalignedfile1 = None
            self.unalignedfile2 = None
            return

        self._writefunc(filename1, self.multiMSA.MSAdict1, toloop,
                        labelfunc, inclspec, inclmet)

        if filename2:
            self._writefunc(filename2, self.multiMSA.MSAdict2, toloop,
                            labelfunc, inclspec, inclmet)

    def alignfunc(self, f_in, f_out, c=5, ir=500, **kwargs):
        """
        Create multiple sequence alignment from unaligned sequences

        :param f_in: The file of unaligned sequence.
        :param f_out: The desired output filename.
        :param ir: Specifies the -ir flag to msaprobs
        :param c: Specifies the -c flag to msaprobs

        .. note::
            This function requires
            `msaprobs <http://msaprobs.sourceforge.net/homepage.htm#latest>`_.
        """
        annotfile = '%s.annot' % f_out

        try:
            print "Attempting to align using MSAProbsCommandline."
            from Bio.Align.Applications import MSAProbsCommandline

            cline = MSAProbsCommandline(infile=f_in, outfile=f_out,
                                        annot=annotfile, **kwargs)
            cline()
        except ImportError as e1:
            print e1
            print "Trying another way..."

            callstr = 'msaprobs -annot %s -c %s -ir %s' % (annotfile, c, ir)

            tf = tempfile.NamedTemporaryFile(delete=False)
            s, o = getstatusoutput("cat %s | sed 's/\*/X/g' > %s"
                                   % (f_in, tf.name))

            if s == 0:
                s, o = getstatusoutput('%s %s > %s'
                                       % (callstr, tf.name, f_out))

            try:
                os.remove(tf.name)
            except Exception as e:
                print "Error removing temporary file"
                print e

            if s != 0:
                print "ERROR in aligning sequence"
                print o
                raise OSError

    def AAtoDNA(self, f_AA_aligned, f_DNA_unaligned, f_DNA_out):
        """
        Create a DNA multiple sequence alignment
            from an amino acid multiple sequence alignment

        :param f_AA_aligned:  File containing aligned amino acid sequences
        :param f_DNA_unaligned: The file containing unaligned DNA sequences
        :param f_DNA_out: The desired output DNA filename.

        .. note::
            This function requires
            `pal2nal <http://www.bork.embl.de/pal2nal/>`_.
        """
        callstr = 'perl pal2nal.pl %s %s -output fasta > %s' \
                  % (f_AA_aligned, f_DNA_unaligned, f_DNA_out)
        s, o = getstatusoutput(callstr)
        if s != 0:
            print "Error in pal2nal"
            print o
            raise OSError

    def align(self, filename1, filename2=None, AAtoDNA=True):
        """
        Align orthologous sequences.

        :param filename1: Output filename for primary alignments
        :param filename2: Output filename for secondary alignments
        :param AAtoDNA: Specifies that secondary sequences
            are DNA and should be aligned based on AA alignment.
        """

        if self.unalignedfile1:
            self.alignfunc(self.unalignedfile1, filename1)
            if filename2 and AAtoDNA:
                self.AAtoDNA(filename1, self.unalignedfile2, filename2)
            elif filename2:
                self.alignfunc(self.unalignedfile2, filename2)
        else:
            print "No file to align!"
