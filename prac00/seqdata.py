'''
Module that provides methods and classes for working with genome sequence data.
For instance,
- BED files
- 2bit genome sequence files

Examples:

>>> import seqdata as sqd
>>> hg19 = sqd.TwoBitFile('/Users/mikael/simhome/hg19.2bit')
>>> mybed1 = sqd.BedFile('/Users/mikael/simhome/mcf7/h3k4me3.bed', format = 'Peaks')
>>> print mybed1[0]
    ('chr1', 540640, 540790)
>>> for e in mybed1[0:5]:
...    print e.getloc(hg19)

    ggcgttttcctgtaaagttgggcacacgcttcccacatgactcagcaattgcacttctgggtatgtacccgagagaaacaaaagcttatgttcacacaaaaacctacaacgcaaatgcacaaacagctctatccaacaaccctggaagca
    ATATAGTAAAACCCAGCCCATGGCCCCTAACAGGGGCCCTCTCAGCCCTCCTAATGACCTCCGGCCTAGCCATGTGATTTCACTTCCACTCCACAACCCTCCTCATACTAGGCCTACTAACCAACACACTAACCATATACCAATGATGGC
    ...
>>> for key in hg19:
        print key

    chr19_gl000208_random
    chr8_gl000197_random
    chr6_apd_hap1
    chr13
    ...
>>> hg19['chrX']
    Out[1]: <seqdata.TwoBitSequence at 0x10cde7d90>
>>> len(hg19['chrX'])
    Out[1]: 155270560
>>> hg19['chrX'][1000000:1000060]
    Out[1]: 'AAAcagctacttggaaggctgaagcaggaggattgtttgagtctaggagtttgaggctgc'
'''

def overlap(chromLoc1, chromLoc2):
    """ Check if chromosome locations described by tuples
        (chrom, chromStart, chromEnd) overlap. 
        If so return the number of positions that overlap.
        Return 0 in case of NO overlap.
    """
    if chromLoc1[0] == chromLoc2[0]:
        halfWidth1 = (chromLoc1[2] - chromLoc1[1]) / 2
        halfWidth2 = (chromLoc2[2] - chromLoc2[1]) / 2
        minWidth = min(halfWidth1, halfWidth2)
        minWidth = max(minWidth, 1)
        maxWidth = max(halfWidth1, halfWidth2)
        maxWidth = max(maxWidth, 1)
        centre1 = chromLoc1[1] + halfWidth1
        centre2 = chromLoc2[1] + halfWidth2
        diffCentres = abs(centre1 - centre2)
        if diffCentres + minWidth < maxWidth: # one fragment encompasses the other
            return minWidth * 2
        else:
            return max(0, halfWidth1 + halfWidth2 - diffCentres)
    else:
        return 0

def distance(chromLoc1, chromLoc2, minimum = True):
    """ Check the distance between two locations described by tuples
        (chrom, chromStart, chromEnd).
        If chromLoc1 is BEFORE chromLoc2 then the distance is positive, else negative.
        If not on same chromosome return None.
        minimum: if True (default), then use minimum distance, if False, use centre to centre
    """ 
    if chromLoc1[0] == chromLoc2[0]:
        halfWidth1 = (chromLoc1[2] - chromLoc1[1]) / 2
        halfWidth2 = (chromLoc2[2] - chromLoc2[1]) / 2
        minWidth = min(halfWidth1, halfWidth2)
        minWidth = max(minWidth, 1)
        maxWidth = max(halfWidth1, halfWidth2)
        maxWidth = max(maxWidth, 1)
        centre1 = chromLoc1[1] + halfWidth1
        centre2 = chromLoc2[1] + halfWidth2
        diffCentres = abs(centre1 - centre2)
        if not minimum:
            return centre2 - centre1
        if diffCentres + minWidth < maxWidth: # one fragment encompasses the other
            return 0
        elif halfWidth1 + halfWidth2 - diffCentres > 0: # fragments overlap A-to-B or B-to-A
            return 0
        else:
            loc1_is_1st = chromLoc2[1] - chromLoc1[2]
            loc1_is_2nd = chromLoc1[1] - chromLoc2[2]
            if loc1_is_1st > loc1_is_2nd:
                return loc1_is_1st
            else:
                return -loc1_is_2nd
    else:
        return None
    
class BedEntry():
    
    def __init__(self, chrom, chromStart, chromEnd):
        self.chrom = chrom
        self.chromStart = chromStart
        self.chromEnd = chromEnd
        self.blockCount = None
        
    def addOption(self, 
                  name = None, 
                  score = None, 
                  strand = None, 
                  thickStart = None, 
                  thickEnd = None, 
                  itemRgb = None, 
                  blockCount = None, 
                  blockSizes = None, 
                  blockStarts = None,
                  signalValue = None,
                  pValue = None,
                  qValue = None,
                  peak = None,
                  tags = None,
                  summit = None,
                  fold = None,
                  fdr = None,
                  zscore = None,
                  bg = None):
        if name: self.name = name
        if score: self.score = score
        if strand: self.strand = strand
        if thickStart: self.thickStart = thickStart
        if thickEnd: self.thickEnd = thickEnd
        if itemRgb: self.itemRgb = [int(color) for color in itemRgb.split(',')]
        if blockCount: 
            self.blockCount = max(0, blockCount)
            if blockCount > 0: 
                self.blockSizes = [int(sizeword) for sizeword in blockSizes.split(',')]
                self.blockStarts = [int(startword) for startword in blockStarts.split(',')]
                if len(self.blockSizes) != blockCount or len(self.blockStarts) != blockCount:
                    raise RuntimeError('Blockcount is incorrect in BED entry \"%s\"' % str(self)) 
        if signalValue: self.signalValue = signalValue
        if pValue: self.pValue = pValue
        if qValue: self.qValue = qValue
        if peak: self.peak = peak
        if tags: self.tags = tags
        if summit: self.summit = summit
        if fold: self.fold = fold
        if fdr: self.fdr = fdr
        if bg: self.bg = bg
        if zscore: self.zscore = zscore
        
    def __str__(self):
        return str((self.chrom, self.chromStart, self.chromEnd))
    
    def __getitem__(self, i):
        if self.blockCount:
            return (self.chrom, self.blockStarts[i], self.blockStarts[i] + self.blockSizes[i]) 

    def __iter__(self):
        if self.blockCount:
            for i in range(self.blockCount):
                if self.blockSizes[i] > 0:
                    yield (self.chrom, self.blockStarts[i], self.blockStarts[i] + self.blockSizes[i])

    def __len__(self):
        return self.blockCount

    def loc(self, genome = None, fixedwidth = None, usesummit = False):
        """ Retrieve the genomic location for BED entry, or sequence if genome is provided
            genome: a dictionary with keys for sequence names, e.g. 'chr1', 'chrX', etc, and values with indexed/sliceable strings
            fixedwidth: the width of the location/sequence if the width in the BED entry is ignored, and only its centre is used
            usesummit: centre a fixedwidth window around an assigned "summit" 
        """
        if fixedwidth:
            if usesummit:
                diff = self.summit - fixedwidth / 2
            else:
                diff = (self.chromEnd - self.chromStart) / 2 - fixedwidth / 2
            if genome:
                return genome[self.chrom][max(0, self.chromStart + diff) : min(len(genome[self.chrom]), self.chromStart + diff + fixedwidth)]
            else:
                return (self.chrom, self.chromStart + diff, self.chromStart + diff + fixedwidth)
        if genome:
            return genome[self.chrom][self.chromStart:self.chromEnd]
        else:
            return (self.chrom, self.chromStart, self.chromEnd)

    def setwidth(self, fixedwidth = None, usesummit = False):
        if fixedwidth:
            if usesummit:
                diff = self.summit - fixedwidth / 2
            else:
                diff = (self.chromEnd - self.chromStart) / 2 - fixedwidth / 2
            self.chromStart += diff
            self.chromStart += diff + fixedwidth
        return (self.chrom, self.chromStart, self.chromEnd)

class BedFile():
    """ Read BED file.

        See http://genome.ucsc.edu/FAQ/FAQformat#format1

        The first three required BED fields are (part of all supported sub-formats):

        chrom - The name of the chromosome (e.g. chr3, chrY, chr2_random) or scaffold (e.g. scaffold10671).
        chromStart - The starting position of the feature in the chromosome or scaffold. The first base in a chromosome is numbered 0.
        chromEnd - The ending position of the feature in the chromosome or scaffold. The chromEnd base is not included in the display of the feature. For example, the first 100 bases of a chromosome are defined as chromStart=0, chromEnd=100, and span the bases numbered 0-99.

        The 9 additional optional BED fields are (part of sub-format "Optional"):
        
        name - Defines the name of the BED line. This label is displayed to the left of the BED line in the Genome Browser window when the track is open to full display mode or directly to the left of the item in pack mode.
        score - A score between 0 and 1000. If the track line useScore attribute is set to 1 for this annotation data set, the score value will determine the level of gray in which this feature is displayed (higher numbers = darker gray). This table shows the Genome Browser's translation of BED score values into shades of gray:
        shade                                             
        strand - Defines the strand - either '+' or '-'.
        thickStart - The starting position at which the feature is drawn thickly (for example, the start codon in gene displays).
        thickEnd - The ending position at which the feature is drawn thickly (for example, the stop codon in gene displays).
        itemRgb - An RGB value of the form R,G,B (e.g. 255,0,0). If the track line itemRgb attribute is set to "On", this RBG value will determine the display color of the data contained in this BED line. NOTE: It is recommended that a simple color scheme (eight colors or less) be used with this attribute to avoid overwhelming the color resources of the Genome Browser and your Internet browser.
        blockCount - The number of blocks (exons) in the BED line.
        blockSizes - A comma-separated list of the block sizes. The number of items in this list should correspond to blockCount.
        blockStarts - A comma-separated list of block starts. All of the blockStart positions should be calculated relative to chromStart. The number of items in this list should correspond to blockCount.

        ENCODE also defines broadpeaks and narrowpeaks format (part of our "Peaks" sub-format):
        
        name - Defines the name of the BED line. This label is displayed to the left of the BED line in the Genome Browser window when the track is open to full display mode or directly to the left of the item in pack mode.
        score - Indicates how dark the peak will be displayed in the browser (0-1000). If all scores were '0' when the data were submitted to the DCC, the DCC assigned scores 1-1000 based on signal value. Ideally the average signalValue per base spread is between 100-1000.
        strand - +/- to denote strand or orientation (whenever applicable). Use '.' if no orientation is assigned.
        signalValue - Measurement of overall (usually, average) enrichment for the region.
        pValue - Measurement of statistical significance (-log10). Use -1 if no pValue is assigned.
        qValue - Measurement of statistical significance using false discovery rate (-log10). Use -1 if no qValue is assigned.
        peak - Point-source called for this peak; 0-based offset from chromStart. Use -1 if no point-source called.      
        
        MACS also defines a "summit" peaks format (part of our "Summit" sub-format)
        It contains the peak summits locations for every peaks. The 5th column in this file is the .
        In addition to the required three, the following fields follow:
        length         [redundant, ignored]
        summit         summit height of fragment pileup
        tags
        pValue         [-10*log10(pvalue)]    
        fold           [enrichment]  
        FDR            [%; optional]

        "CCAT" BED-like file format: 
        chromosome, 
        peakcenter    [converted to summit], 
        regionstart, 
        regionend, 
        tags          [tagcount], 
        bg            [bgcount], 
        zscore, 
        fdr
        
    """
    
    def __init__(self, filename, format = 'Limited'):
        """ Read a BED file.
            format: specifies the format of the file, 
            "Limited", e.g.
                chr22 1000 5000
                chr22 2000 6000
            "Optional", e.g.
                track name=pairedReads description="Clone Paired Reads" useScore=1
                chr22 1000 5000 cloneA 960 + 1000 5000 0 2 567,488, 0,3512
                chr22 2000 6000 cloneB 900 - 2000 6000 0 2 433,399, 0,3601
                ...
                (also handles the Limited + score format)
            "Peaks", e.g.
                chr1    569780    569930    .    0    .    19    6.07811    -1    -1
                chr1    713300    713450    .    0    .    54    49.1167    -1    -1
            "Summit", e.g.
                # d = 130
                chr      start    end   length summit  tags -10*log10(pvalue)    fold_enrichment    FDR(%)
                chr1     8250     8671    422    286    46    145.84    11.68    0.51
                chr1    36382    36984    603    405    46    315.23    27.05    0.24
            "CCAT", e.g.
                chr8    94747805    94747070    94749250    525     3    21.519196    0.002000
                chr17   55277895    55277070    55279280    560    18    21.283333    0.002000
            "Cropped", e.g.
                chr1    851602    10
                chr1    921184    18
                chr1    931838    9
        """
        f = open(filename)
        row = 0
        acceptHeaderRows = 1
        headerRow = None
        self.rows = []
        self.format = 'Limited'
        self.indices = None
        for line in f:
            row += 1
            words = line.strip().split()
            if len(words) == 0:
                continue # ignore empty lines
            if words[0].strip().startswith('#'):
                continue # comment
            if words[0].strip().startswith('browser'):
                continue # ignore
            if words[0].strip().startswith('track'):
                continue # ignore
            try:
                chrom = words[0]
                if format.lower().startswith('ccat'):
                    chromStart = int(words[2])
                    chromEnd = int(words[3])
                else: # all other standard BED formats
                    chromStart = int(words[1])
                    chromEnd = int(words[2])
                entry = BedEntry(chrom, chromStart, chromEnd)
                if format.lower().startswith('opt'):
                    self.format = 'Optional'
                    if len(words) >= 12:
                        entry.addOption(name = words[3], score = float(words[4]), strand = words[5], thickStart = int(words[6]), thickEnd = int(words[7]), itemRgb = words[8], blockCount = int(words[9]), blockSizes = words[10], blockStarts = words[11])
                    elif len(words) >= 9:
                        entry.addOption(name = words[3], score = float(words[4]), strand = words[5], thickStart = int(words[6]), thickEnd = int(words[7]), itemRgb = words[8])
                    elif len(words) >= 6:
                        entry.addOption(name = words[3], score = float(words[4]), strand = words[5])
                    else:
                        entry.addOption(name = '.', score = int(words[3]), strand = '.')
                elif format.lower().startswith('peak'):
                    self.format = 'Peaks'
                    if len(words) >= 10: # narrowpeaks
                        entry.addOption(name = words[3], score = int(words[4]), strand = words[5], signalValue = float(words[6]), pValue = float(words[7]), qValue = float(words[8]), peak = int(words[9]))
                    else: # broadpeaks
                        entry.addOption(name = words[3], score = int(words[4]), strand = words[5], signalValue = float(words[6]), pValue = float(words[7]), qValue = float(words[8]))
                elif format.lower().startswith('summit'):
                    self.format = 'Summit'
                    if len(words) >= 9:
                        entry.addOption(summit = int(words[4]), tags = int(words[5]), pValue = float(words[6]), fold = float(words[7]), fdr = float(words[8]))
                    else:
                        entry.addOption(summit = int(words[4]), tags = int(words[5]), pValue = float(words[6]), fold = float(words[7]))
                elif format.lower().startswith('ccat'):
                    self.format = 'CCAT'
                    entry.addOption(summit = int(words[1]) - entry.chromStart, tags = int(words[4]), bg = int(words[5]), zscore = float(words[6]), fdr = float(words[7]), name = '.', score = int(words[4]), strand = '.')
                elif format.lower().startswith('crop'):
                    entry.addOption(score = int(words[2]), name = '.', strand = '.')
                    entry.chromEnd = entry.chromStart + 1
                self.rows.append(entry)
            except RuntimeError as e:
                if not acceptHeaderRows:
                    raise RuntimeError('Error in BED file at row %d (%s)' % (row, e.strerror))
                else:
                    headerRow = words
                    acceptHeaderRows -= 1 # count down the number of header rows that can occur
        f.close()
        self.indices = self._createIndices()

    def __iter__(self):
        return self.rows.__iter__()
    
    def __getslice__(self, i, j):
        return self.rows.__getslice__(i, j)
    
    def __getitem__(self, i):
        return self.rows[i]
    
    def __len__(self):
        return len(self.rows)

    def _createIndices(self):
        index_start = {}
        index_centre = {}
        index_end = {}
        for i in range(len(self.rows)):
            row = self.rows[i]
            if not index_start.has_key(row.chrom): # seeing chromosome entry first time
                index_start[row.chrom] = []
            if not index_centre.has_key(row.chrom): # seeing chromosome entry first time
                index_centre[row.chrom] = []
            if not index_end.has_key(row.chrom): # seeing chromosome entry first time
                index_end[row.chrom] = []
            index_start[row.chrom].append((row.chromStart, row.chromEnd - row.chromStart, i))
            index_centre[row.chrom].append((row.chromStart + (row.chromEnd - row.chromStart) / 2, (row.chromEnd - row.chromStart) / 2, i))
            index_end[row.chrom].append((row.chromEnd, row.chromEnd - row.chromStart, i))
        for chr in index_start:
            index_start[chr].sort()
            index_centre[chr].sort()
            index_end[chr].sort()
        return (index_start, index_centre, index_end)
        
    def __contains__(self, elem):
        """ Test for containment: does the specified elem overlap with at least one of the BED entries.
            The method performs a binary search. """
        try:
            if isinstance(elem, BedEntry):
                elem = elem.loc()
            entries = self.indices[0][elem[0]] # use the start index
            upper = len(entries)    # keep an upper boundary
            lower = 0               # and a lower boundary
            inspect = (upper - lower) / 2 # start by looking in the middle
            while True:
                entry = self.rows[entries[inspect][2]]
                d = distance(entry.loc(), elem, minimum = True)
                delta = 0
                if d == 0:
                    return True
                elif d > 0:
                    lower = inspect + 1
                    delta = (upper - inspect) / 2 # splitting in half, potential speed improvements with some heuristic?
                    inspect += delta
                else:
                    upper = inspect
                    delta = (inspect - lower + 1) / 2
                    inspect -= delta
                if delta == 0:
                    return False
        except KeyError:
            return False

    def match(self, elem, name):
        """ Test for containment: does the specified elem overlap with at least one of the BED entries 
            that has the nominated name (label).""" 
        try:
            if isinstance(elem, BedEntry):
                elem = elem.loc()
            entries = self.indices[0][elem[0]] # use the start index
            upper = len(entries)    # keep an upper boundary
            lower = 0               # and a lower boundary
            inspect = (upper - lower) / 2 # start by looking in the middle
            while True:
                entry = self.rows[entries[inspect][2]]
                d = distance(entry.loc(), elem, minimum = True)
                delta = 0
                if d == 0:
                    delta = 0
                    while d == 0:
                        if entry.name == name:
                            return True
                        delta += 1
                        entry = self.rows[entries[inspect + delta][2]]
                        d = distance(entry.loc(), elem, minimum = True)
                    delta = -1
                    entry = self.rows[entries[inspect + delta][2]]
                    d = distance(entry.loc(), elem, minimum = True)
                    while d == 0:
                        if entry.name == name:
                            return True
                        delta -= 1
                        entry = self.rows[entries[inspect + delta][2]]
                        d = distance(entry.loc(), elem, minimum = True)
                    return False
                elif d > 0:
                    lower = inspect + 1
                    delta = (upper - inspect) / 2 # splitting in half, potential speed improvements with some heuristic?
                    inspect += delta
                else:
                    upper = inspect
                    delta = (inspect - lower + 1) / 2
                    inspect -= delta
                if delta == 0:
                    return False
        except KeyError:
            return False
        
    def closest(self, myloc, minimum = True):
        """ Find the closest entry in the current BedFile to a given location.
            Return a tuple with the absolute distance and the entry that is closest.
            If several entries are closest, then any of the closest entries are returned. 
            If no location is found on the same chromosome, the tuple None, None is returned.
            minimum: if True, use minimum distance, if False, use centre to centre distance.
        """
        mindist = None
        minentry = None
        try:
            if isinstance(myloc, BedEntry):
                myloc = myloc.loc()
            if minimum:
                entries = self.indices[0][myloc[0]] # use start index
                upper = len(entries)    # keep an upper boundary
                lower = 0               # and a lower boundary
                inspect = (upper - lower) / 2 # start by looking in the middle
                delta = None
                while not delta == 0:
                    entry = self.rows[entries[inspect][2]]
                    d = distance(entry.loc(), myloc, minimum = True)
                    if mindist == None:
                        mindist = abs(d)
                        minentry = entry
                    elif abs(d) < mindist:
                        mindist = abs(d)
                        minentry = entry
                    if d == 0:
                        return (mindist, minentry)
                    elif d > 0:
                        lower = inspect + 1
                        delta = (upper - inspect) / 2 # splitting in half, potential speed improvements with some heuristic?
                        inspect += delta
                    else:
                        upper = inspect
                        delta = (inspect - lower + 1) / 2
                        inspect -= delta
                # we may have missed the closest, so need to look around this point
                for i_dn in range(inspect + 1, len(entries)): # Look downstream since   
                    entry = self.rows[entries[i_dn][2]]
                    d = distance(entry.loc(), myloc, minimum = True)
                    if abs(d) < mindist:
                        mindist = abs(d)
                        minentry = entry
                    elif abs(d) > mindist:
                        break
                # also need to investigate upstream, doing so by using end index
                entries = self.indices[2][myloc[0]] # use end index
                upper = len(entries)    # keep an upper boundary
                lower = 0               # and a lower boundary
                inspect = (upper - lower) / 2 # start by looking in the middle
                delta = None
                while not delta == 0:
                    entry = self.rows[entries[inspect][2]]
                    d = distance(entry.loc(), myloc, minimum = True)
                    if abs(d) < mindist:
                        mindist = abs(d)
                        minentry = entry
                    if d == 0:
                        return (mindist, minentry)
                    elif d > 0:
                        lower = inspect + 1
                        delta = (upper - inspect) / 2 # splitting in half, potential speed improvements with some heuristic?
                        inspect += delta
                    else:
                        upper = inspect
                        delta = (inspect - lower + 1) / 2
                        inspect -= delta
                # we may have missed the closest, so need to look around this point
                for i_up in range(inspect - 1, 0, -1): # Look upstream since   
                    entry = self.rows[entries[i_up][2]]
                    d = distance(entry.loc(), myloc, minimum = True)
                    if abs(d) < mindist:
                        mindist = abs(d)
                        minentry = entry
                    elif abs(d) > mindist:
                        break
                return (mindist, minentry)
            else: # minimum == False, i.e. use centre-to-centre distance
                entries = self.indices[1][myloc[0]] # use centre index
                upper = len(entries)    # keep an upper boundary
                lower = 0               # and a lower boundary
                inspect = (upper - lower) / 2 # start by looking in the middle
                delta = None
                while not delta == 0:
                    entry = self.rows[entries[inspect][2]]
                    d = distance(entry.loc(), myloc, minimum = False)
                    if mindist == None:
                        mindist = abs(d)
                        minentry = entry
                    elif abs(d) < mindist:
                        mindist = abs(d)
                        minentry = entry
                    if d == 0:
                        return (mindist, minentry)
                    elif d > 0:
                        lower = inspect + 1
                        delta = (upper - inspect) / 2 # splitting in half, potential speed improvements with some heuristic?
                        inspect += delta
                    else:
                        upper = inspect
                        delta = (inspect - lower + 1) / 2
                        inspect -= delta
                # at bottom of search
                return (mindist, minentry)
        except KeyError:
            return None

    def write(self, filename, format = 'BED6'):
        """ Save the data
            format - the format to use for WRITING, currently only BED6 ('Optional' 6-col format) is supported.
        """
        f = open(filename, 'w')
        for row in self.__iter__():
            if self.format == 'Peaks':
                f.write("%s %d %d %s %d %s" % (row.chrom, row.chromStart, row.chromEnd, row.name, row.signalValue, row.strand))
            elif self.format == 'Limited':
                f.write("%s %d %d" % (row.chrom, row.chromStart, row.chromEnd))
            else:
                f.write("%s %d %d %s %d %s" % (row.chrom, row.chromStart, row.chromEnd, row.name, row.score, row.strand))
            f.write("\n")
        f.close()
        
def writeBedFile(entries, filename, format = 'BED6'):
    """ Save the BED entries to a BED file.
        format - the format to use for WRITING, currently only BED6 ('Optional' 6-col format) is supported.
    """
    f = open(filename, 'w')
    for row in entries:
        if format == 'Peaks':
            f.write("%s %d %d %s %d %s" % (row.chrom, row.chromStart, row.chromEnd, row.name, row.signalValue, row.strand))
        elif format == 'Limited':
            f.write("%s %d %d" % (row.chrom, row.chromStart, row.chromEnd))
        else:
            f.write("%s %d %d %s %d %s" % (row.chrom, row.chromStart, row.chromEnd, row.name, row.score, row.strand))
        f.write("\n")
    f.close()

        
"""
This following code is a modified version of twobitreader (which is under Perl Artistic License 2.0).
As per license restrictions, the code below indicates what has been modified in relation to the 
standard version (retrieved from https://bitbucket.org/thesylex/twobitreader on the 16 May 2012). 
No warranty is provided, express or implied

Modifications to package:
- removed download.py and __main__ because they were not used and __main__ had errors.
- removed command-line interface because the BED file functionality is implemented more extensively elsewhere
"""
from array import array
from bisect import bisect_right
from errno import ENOENT, EACCES
from os import R_OK, access
try:
    from os import strerror
except ImportError:
    strerror = lambda x: 'strerror not supported'
from os.path import exists
from itertools import izip

def true_long_type():
    """
        OS X uses an 8-byte long, so make sure L (long) is the right size
        and switch to I (int) if needed
    """
    for type_ in ['L', 'I']:
        test_array = array(type_, [0])
        long_size = test_array.itemsize
        if long_size == 4: return type_
    raise ImportError("Couldn't determine a valid 4-byte long type to use \
as equivalent to LONG")
LONG = true_long_type()

def byte_to_bases(x):
    """convert one byte to the four bases it encodes"""
    c = (x >> 4) & 0xf
    f = x & 0xf
    cc = (c >> 2) & 0x3
    cf = c & 0x3
    fc = (f >> 2) & 0x3
    ff = f & 0x3
    return map(bits_to_base, (cc, cf, fc, ff))

def bits_to_base(x):
    """convert integer representation of two bits to correct base"""
    if x is 0: return 'T'
    if x is 1: return 'C'
    if x is 2: return 'A'
    if x is 3: return 'G'

def base_to_bin(x):
    """
        provided for user convenience
        convert a nucleotide to its bit representation 
    """
    if x == 'T': return '00'
    if x == 'C': return '01'
    if x == 'A': return '10'
    if x == 'G': return '11'

def create_byte_table():
    """create BYTE_TABLE"""
    d = {}
    for x in xrange(2**8):
        d[x] = byte_to_bases(x)
    return d

def split16(x):
    """
        split a 16-bit number into integer representation
        of its course and fine parts in binary representation
    """
    c = (x >> 8) & 0xff
    f = x & 0xff
    return c, f

def create_twobyte_table():
    """create TWOBYTE_TABLE"""
    d = {}
    for x in xrange(2**16):
        c, f = split16(x)
        d[x] = byte_to_bases(c) + byte_to_bases(f)
    return d

BYTE_TABLE = create_byte_table()
TWOBYTE_TABLE = create_twobyte_table()

def longs_to_char_array(longs, first_base_offset, last_base_offset, array_size):
    """
        takes in a iterable of longs and converts them to bases in a char array
        returns a ctypes string buffer
    """
    longs_len = len(longs)
    # dna = ctypes.create_string_buffer(array_size)
    dna = array('c', 'N' * longs_len)
    # translate from 32-bit blocks to bytes
    # this method ensures correct endianess (byteswap as neeed)
    bytes = array('B')
    bytes.fromstring(longs.tostring())
    # first block
    first_block = ''.join([''.join(BYTE_TABLE[bytes[x]]) for x in range(4)])
    i = 16 - first_base_offset
    if array_size < i: i = array_size
    dna[0:i] = array('c', first_block[first_base_offset:first_base_offset + i])
    if longs_len == 1: return dna
    # middle blocks (implicitly skipped if they don't exist)
    for byte in bytes[4:-4]:
        dna[i:i + 4] = array('c', BYTE_TABLE[byte])
        i += 4
    # last block
    last_block = array('c', ''.join([''.join(BYTE_TABLE[bytes[x]]) for x in range(-4,0)]))
    dna[i:i + last_base_offset] = last_block[0:last_base_offset]
    return dna

class TwoBitFile(dict):
    """
        python-level reader for .2bit files (i.e., from UCSC genome browser)
        (note: no writing support)
        
        TwoBitFile inherits from dict
        You may access sequences by name, e.g.
        >>> genome = TwoBitFile('hg18.2bit')
        >>> chr20 = genome['chr20']
        
        Sequences are returned as TwoBitSequence objects
        You may access intervals by slicing or using str() to dump the entire entry
        e.g.
        >>> chr20[100100:100200]
        'ttttcctctaagataatttttgccttaaatactattttgttcaatactaagaagtaagataacttccttttgttggtat
        ttgcatgttaagtttttttcc'
        >>> whole_chr20 = str(chr20)
        
        Fair warning: dumping the entire chromosome requires a lot of memory
        
        See TwoBitSequence for more info
    """
    
    def __init__(self, foo):
        super(TwoBitFile, self).__init__()
        if not exists(foo):
            raise IOError(ENOENT, strerror(ENOENT), foo)
        if not access(foo, R_OK):
            raise IOError(EACCES, strerror(EACCES), foo)
        self._filename = foo
        self._file_handle = open(foo, 'rb')
        self._load_header()
        self._load_index()
        for name, offset in self._offset_dict.iteritems():
            self[name] = TwoBitSequence(self._file_handle, offset,
                                        self._byteswapped)
        return        
        
    def _load_header(self):
        file_handle = self._file_handle
        header = array(LONG)
        header.fromfile(file_handle, 4)
        # check signature -- must be 0x1A412743
        # if not, swap bytes
        byteswapped = False
        (signature, version, sequence_count, reserved) = header
        if not signature == 0x1A412743:
            byteswapped = True
            header.byteswap()
            (signature2, version, sequence_count, reserved) = header
            if not signature2 == 0x1A412743:
                raise TwoBitFileError('Signature in header should be 0x1A412743'
                                    + ', instead found 0x%X' % signature)
        if not version == 0: 
            raise TwoBitFileError('File version in header should be 0.')
        if not reserved == 0:
            raise TwoBitFileError('Reserved field in header should be 0.')
        self._byteswapped = byteswapped
        self._sequence_count = sequence_count
        
    def _load_index(self):
        file_handle = self._file_handle
        byteswapped = self._byteswapped
        remaining = self._sequence_count
        sequence_offsets = []
        file_handle.seek(16)
        while True:
            if remaining == 0: break
            name_size = array('B')
            name_size.fromfile(file_handle, 1)
            if byteswapped: name_size.byteswap()
            name = array('c')
            if byteswapped: name.byteswap()
            name.fromfile(file_handle, name_size[0])
            offset = array(LONG)
            offset.fromfile(file_handle, 1)
            if byteswapped: offset.byteswap()
            sequence_offsets.append((name.tostring(), offset[0]))
            remaining -= 1
        self._sequence_offsets = sequence_offsets
        self._offset_dict = dict(sequence_offsets)

    def sequence_sizes(self):
        """returns a dictionary with the sizes of each sequence"""
        d = {}
        file_handle = self._file_handle
        byteswapped = self._byteswapped
        for name, offset in self._offset_dict.iteritems():
            file_handle.seek(offset)
            dna_size = array(LONG)
            dna_size.fromfile(file_handle, 1)
            if byteswapped: dna_size.byteswap()
            d[name] = dna_size[0]
        return d

class TwoBitSequence(object):
    """
        A TwoBitSequence object refers to an entry in a TwoBitFile
        
        You may access intervals by slicing or using str() to dump the entire entry
        e.g.
        >>> genome = TwoBitFile('hg18.2bit')
        >>> chr20 = genome['chr20']
        >>> chr20[100100:100200] # slicing returns a string
        'ttttcctctaagataatttttgccttaaatactattttgttcaatactaagaagtaagataacttccttttgttggtat
        ttgcatgttaagtttttttcc'
        >>> whole_chr20 = str(chr20) # get whole chr as string
        
        Fair warning: dumping the entire chromosome requires a lot of memory
         
        Note that we follow python/UCSC conventions:
        Coordinates are 0-based, end-open
        (Note: The UCSC web-based genome browser uses 1-based closed coordinates)
        If you attempt to access a slice past the end of the sequence,
        it will be truncated at the end.
        
        Your computer probably doesn't have enough memory to load a whole genome
        but if you want to string-ize your TwoBitFile, here's a recipe:
        
        x = TwoBitFile('my.2bit')
        d = x.dict()
        for k,v in d.iteritems(): d[k] = str(v)
    """
    def __init__(self, file_handle, offset, byteswapped=False):
        self._file_handle = file_handle
        self._original_offset = offset
        self._byteswapped = byteswapped
        file_handle.seek(offset)
        header = array(LONG)
        header.fromfile(file_handle, 2)
        if byteswapped: header.byteswap()
        dna_size, n_block_count = header
        self._dna_size = dna_size
        self._packed_dna_size = (dna_size + 15) / 16 # this is 32-bit fragments
        n_block_starts = array(LONG)
        n_block_sizes = array(LONG)
        n_block_starts.fromfile(file_handle, n_block_count)
        if byteswapped: n_block_starts.byteswap()
        n_block_sizes.fromfile(file_handle, n_block_count)
        if byteswapped: n_block_sizes.byteswap()
        self._n_block_starts = n_block_starts
        self._n_block_sizes= n_block_sizes
        mask_rawc = array(LONG)
        mask_rawc.fromfile(file_handle, 1)
        if byteswapped: mask_rawc.byteswap()
        mask_block_count = mask_rawc[0]
        mask_block_starts = array(LONG)
        mask_block_starts.fromfile(file_handle, mask_block_count)
        if byteswapped: mask_block_starts.byteswap()
        mask_block_sizes = array(LONG)
        mask_block_sizes.fromfile(file_handle, mask_block_count)
        if byteswapped: mask_block_sizes.byteswap()
        self._mask_block_starts = mask_block_starts
        self._mask_block_sizes = mask_block_sizes
        file_handle.read(4)
        self._offset = file_handle.tell()

    def __len__(self):
        return self._dna_size

    def __getslice__(self, min_, max_=None):
        return self.get_slice(min_, max_)

    def get_slice(self, min_, max_=None):
        """
        get_slice returns only a sub-sequence
        """
        # handle negative coordinates
        dna_size = self._dna_size
        if max_ < 0:
            if max_ < -dna_size: raise IndexError('index out of range')
            max_ = dna_size + 1 + max_
        if min_ < 0:
            if max_ < -dna_size: raise IndexError('index out of range')
            min_ = dna_size + 1 + min_
        # make sure there's a proper range
        if min_ > max_ and max_ is not None: return ''
        if max_ == 0: return ''
        # load all the data
        if max_ > dna_size: max_ = dna_size
        file_handle = self._file_handle
        byteswapped = self._byteswapped
        n_block_starts = self._n_block_starts
        n_block_sizes = self._n_block_sizes
        mask_block_starts = self._mask_block_starts
        mask_block_sizes = self._mask_block_sizes
        offset = self._offset
        packed_dna_size = self._packed_dna_size

        # region_size is how many bases the region is       
        if max_ is None: region_size = dna_size - min_
        else: region_size = max_ - min_
        
        # start_block, end_block are the first/last 32-bit blocks we need
        # note: end_block is not read
        # blocks start at 0
        start_block = min_ / 16
        end_block = max_ / 16
        # don't read past seq end
        if end_block >= packed_dna_size: end_block = packed_dna_size - 1
        # +1 we still need to read block
        blocks_to_read = end_block - start_block + 1
        
        # jump directly to desired file location
        local_offset = offset + start_block * 4
        file_handle.seek(local_offset)
        
        # note we won't actually read the last base
        # this is a python slice first_base_offset:16*blocks+last_base_offset
        first_base_offset = min_ % 16
        last_base_offset = max_ % 16
        
        fourbyte_dna = array(LONG)
        fourbyte_dna.fromfile(file_handle, blocks_to_read)
        if byteswapped: fourbyte_dna.byteswap()
        string_as_array = longs_to_char_array(fourbyte_dna, first_base_offset,
                                              last_base_offset, region_size)
        for start, size in izip(n_block_starts, n_block_sizes):
            end = start + size
            if end <= min_: continue
            if start > max_: break
            if start < min_: start = min_
            if end > max_: end = max_ 
            start -= min_
            end -= min_
            string_as_array[start:end] = array('c', 'N'*(end-start))
        lower = str.lower
        first_masked_region = max(0,
                                  bisect_right(mask_block_starts, min_) - 1)
        last_masked_region = min(len(mask_block_starts),
                                 1 + bisect_right(mask_block_starts, max_,
                                                  lo=first_masked_region))
        for start, size in izip(mask_block_starts[first_masked_region:last_masked_region],
                                mask_block_sizes[first_masked_region:last_masked_region]):
            end = start + size
            if end <= min_: continue
            if start > max_: break
            if start < min_: start = min_
            if end > max_: end = max_ 
            start -= min_
            end -= min_
            string_as_array[start:end] = array('c', lower(string_as_array[start:end].tostring()))
        if not len(string_as_array) == max_ - min_:
            raise RuntimeError, "Sequence was longer than it should be"
        return string_as_array.tostring()

    def __str__(self):
        """
        returns the entire chromosome
        """
        return self.__getslice__(0, None)
    
class TwoBitFileError(StandardError):
    """
    Base exception for TwoBit module
    """
    def __init__(self, msg):
        errtext = 'Invalid 2-bit file. ' + msg
        return super(TwoBitFileError, self).__init__(errtext)

def print_specification():
    """
    Prints the twoBit file format specification I got from the Internet.
    This is only here for reference
    """
    return """
        From http://www.its.caltech.edu/~alok/reviews/blatSpecs.html
        
        .2bit files
        
        A .2bit file can store multiple DNA sequence (up to 4 gig total) in a compact \
        randomly accessible format. The two bit files contain masking information as \
        well as the DNA itself. The file begins with a 16 byte header containing the \
        following fields:
        
        signature - the number 0x1A412743 in the architecture of the machine that \
        created the file.
        version - zero for now. Readers should abort if they see a version number \
        higher than 0.
        sequenceCount - the number of sequences in the file
        reserved - always zero for now.
        All fields are 32 bits unless noted. If the signature value is not as given, \
        the reader program should byte swap the signature and see if the swapped \
        version matches. If so all multiple-byte entities in the file will need to be \
        byte-swapped. This enables these binary files to be used unchanged on \
        different architectures.
        
        The header is followed by a file index. There is one entry in the index for \
        each sequence. Each index entry contains three fields:
        
        nameSize - a byte containing the length of the name field
        name - this contains the sequence name itself, and is variable length \
        depending on nameSize.
        offset - 32 bit offset of the sequence data relative to the start of the file
        
        The index is followed by the sequence records. These contain 9 fields:
        
        dnaSize - number of bases of DNA in the sequence.
        nBlockCount - the number of blocks of N's in the file (representing unknown \
        sequence).
        nBlockStarts - a starting position for each block of N's
        nBlockSizes - the size of each block of N's
        maskBlockCount - the number of masked (lower case) blocks
        maskBlockStarts - starting position for each masked block
        maskBlockSizes - the size of each masked block
        packedDna - the dna packed to two bits per base as so: 00 - T, 01 - C, 10 - A, \
        11 - G. The first base is in the most significant 2 bits byte, and the last \
        base in the least significant 2 bits, so that the sequence TCAG would be \
        represented as 00011011. The packedDna field will be padded with 0 bits as \
        
        necessary so that it takes an even multiple of 32 bit in the file, as this \
        improves i/o performance on some machines.
        .nib files
    """