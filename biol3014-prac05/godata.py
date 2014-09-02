'''
Created on Jul 12, 2012
Module for managing Gene Ontology data, in particular gene:terms 
annotations and term definitions
@author: mikael
'''

from struct import pack, unpack, calcsize, error
import operator
import time
import os
import stats

# Character codes used by binary format to identify ontology 
onto_codes = {
              'P': 'Biological process',
              'F': 'Molecular function',
              'C': 'Cellular component'}

# Labels for edges in the ontology graph, index is used in binary format
onto_rel = ['is_a', 'isect', 'part_of', 'has_part', 'regulates']

# Evidence codes assigned to annotations, an index is assigned when creating binary file and is stored in its header
evid_codes = { # Experimental Evidence Codes
    'EXP': 'Inferred from Experiment',
    'IDA': 'Inferred from Direct Assay',
    'IPI': 'Inferred from Physical Interaction',
    'IMP': 'Inferred from Mutant Phenotype',
    'IGI': 'Inferred from Genetic Interaction',
    'IEP': 'Inferred from Expression Pattern',
               #Computational Analysis Evidence Codes
    'ISS': 'Inferred from Sequence or Structural Similarity',
    'ISO': 'Inferred from Sequence Orthology',
    'ISA': 'Inferred from Sequence Alignment',
    'ISM': 'Inferred from Sequence Model',
    'IGC': 'Inferred from Genomic Context',
    'IBA': 'Inferred from Biological aspect of Ancestor',
    'IBD': 'Inferred from Biological aspect of Descendant',
    'IKR': 'Inferred from Key Residues',
    'IRD': 'Inferred from Rapid Divergence',
    'RCA': 'inferred from Reviewed Computational Analysis',
    'TAS': 'Traceable Author Statement',
    'NAS': 'Non-traceable Author Statement',
               #Curator Statement Evidence Codes
    'IC': 'Inferred by Curator',
    'ND': 'No biological Data available',
               #Automatically-assigned Evidence Codes
    'IEA': 'Inferred from Electronic Annotation',
               #Obsolete Evidence Codes
    'NR': 'Not Recorded'}

class BinGO():
    
    # Structures to hold all data relevant to session, all keys are "encoded"
    annots = {}     # annotations: annots[gene] = (taxa, terms[term] = (evid, T/F))
    termdefs = {}   # definitions: termdefs[term] = (onto, terms[term] = relation, name)
    # Codes for encoding and decoding
    gene_code = None
    term_code = None
    evid_code = None
    # indices
    annot_index = {}
    # Files
    f = None
    
    def __init__(self, filename, taxa = None):
        """ The binary file contains all the data and will initialise 
            gene annotations (annots) and term definitions (termdefs)
            and the encoding/decoding keys. """
        self.f = self._readBitFile(filename, taxa = taxa)
    
    def _decodeGeneIDs(self, gene_codes):
        if type(gene_codes) != list and type(gene_codes) != set and type(gene_codes) != tuple:
            gene_codes = [gene_codes]
        ids = []
        for i in gene_codes:
            s = decode(i, self.gene_code)
            ids.append(s)
        return ids
        
    def _encodeGeneIDs(self, gene_names):
        if type(gene_names) != list and type(gene_names) != set and type(gene_names) != tuple:
            gene_names = [gene_names]
        ids = []
        for i in gene_names:
            y = encode(i, self.gene_code)
            ids.append(y)
        return ids
    
    def _getGeneEntry(self, gene):
        peek = self.annot_index[gene]
        self.f.seek(peek, 0)
        buf = self.f.read(calcsize('IIH'))
        (gene_int, taxa_int, nterms) = unpack('IIH', buf)
        buf = self.f.read(nterms * calcsize('?BI'))
        terms_dict = {}
        for pos in range(0, len(buf) - 1, calcsize('?BI')):
            (qual_bool, evid_int, term_int) = unpack('?BI', buf[pos:pos+calcsize('?BI')])
            terms_dict[term_int] = (evid_int, qual_bool)
        return (taxa_int, terms_dict) 

    def _getSuperTerms(self, term, rel = None):
        """ Recursively compute the transitive closure. """
        found = set()
        try:
            (_, closure, _) = self.termdefs[term]
            for (t, r) in closure.items():
                if (not rel) or r == rel: 
                    found.add(t)
                    found.update(self._getSuperTerms(t, rel))
        except KeyError:
            print 'Could not find GO:%s' % (''.join(decode(term, self.term_code)))
        return found
    
    def _getChildTerms(self, term, rel = None):
        found = set()
        for (child, termdef) in self.termdefs.items():
            (_, parents_dict, _) = termdef
            try:
                myrel = parents_dict[term]
                if rel == myrel or not rel: found.add(child)
            except KeyError:
                pass
        return found
        
    def _getSpecificTerms(self, term, rel = None):
        direct = self._getChildTerms(term, rel)
        found = set()
        for t in direct:
            found.add(t)
            found.update(self._getSpecificTerms(t, rel))
        return found
    
    def getTerms(self, genes, evid = None, onto = None, include_more_general = True):
        """
        Retrieve all GO terms for a given set of genes (or single gene).
        The result is given as a map (key=gene name, value=list of unique terms) OR
        in the case of a single gene as a list of unique terms.
        If include_more_general is True (default) then transitively related terms are included
        """
        mymap = dict()
        # STEP 1: Find all terms directly associated with specified genes
        direct = set() # all GO terms (encoded)
        ids = self._encodeGeneIDs(genes)
        for i in ids:
            gene_name = ''.join(decode(i, self.gene_code))
            mymap[gene_name] = set()
            try:
                (taxa, terms) = self._getGeneEntry(i)
                for (term, evid_and_qual) in terms.items():
                    if evid_and_qual[1] and not evid: # if True and no evidence is specified
                        direct.add(term)
                        mymap[gene_name].add(term)
                    elif self.evid_code[evid_and_qual[0]] == evid:
                        direct.add(term)
                        mymap[gene_name].add(term)
            except KeyError:
                pass
                #print 'Failed to find annotations for gene %s' % gene_name
        if include_more_general:
            # STEP 2: Find the transitive closure of each term identified, store as a dictionary
            indirect = {}
            for t in direct:
                if not indirect.has_key(t):
                    indirect[t] = set(self._getSuperTerms(t))
        # STEP 3: compile and return results
        for gene in mymap:
            term_ids = mymap[gene]
            all_ids = set(term_ids)
            if include_more_general:
                for term_id in term_ids:
                    all_ids.update(indirect[term_id])
            mymap[gene] = set()
            for term_enc in all_ids:
                mymap[gene].add('GO:'+''.join(decode(term_enc, self.term_code)))
        return mymap
    
    def getAllGenes(self):
        names = []
        for g in self._decodeGeneIDs(self.annot_index.keys()):
            names.append(''.join(g))
        return names
    
    def getGenes(self, terms, evid = None, taxa = None, rel = None, include_more_specific = True):
        """ Retrieve all genes that are annotated with specified terms, and qualified by evidence, taxa etc. """
        """ TODO: Debug--suspect this implementation is incorrect. """
        term_ids = set()
        for t in terms:
            term_ids.add(encode(t[3:], self.term_code))
        # STEP 1 (optional): determine more specific terms to be included in query
        if include_more_specific:
            myterms = set()
            for t in term_ids:
                myterms.add(t)
                children = self._getSpecificTerms(t, rel)
                myterms.update(children)
            term_ids = myterms
        # STEP 2: identify genes with those terms
        found = {}
        for g in self.annot_index:
            gene_name = decode(g, self.gene_code)
            (mytaxa, tdict) = self._getGeneEntry(g)
            if not taxa or taxa == mytaxa:
                for annot_term in tdict.keys():
                    if tdict[annot_term] == evid:
                        if annot_term in terms:
                            try:
                                added = found[gene_name]
                                added.add(annot_term)
                            except KeyError:
                                found[gene_name] = set([annot_term])
        # STEP 3: compile and return results
        for gene in found:
            term_ids = found[gene]
            all_ids = set(term_ids)
            found[gene] = set()
            for term_enc in all_ids:
                found[gene].add('GO:'+''.join(decode(term_enc, self.term_code)))
        return found
    
    def getTermdef(self, term):
        term_id = encode(term[3:], self.term_code)
        try:
            (onto_ch, terms_dict, name_peek) = self.termdefs[term_id]
            self.f.seek(name_peek, 0)
            term_name = self.f.readline()
            return (onto_codes[onto_ch], terms_dict, term_name)
        except KeyError:
            return ('Unknown', 'Unknown', 'Unknown')
 
    def _readBitFile(self, filename, taxa, termnames = False):
        f = open(filename, 'r')
        # STEP 1: header info
        ngene_code = None
        nterm_code = None
        nevid_code = None
        ngene_cnt = 0
        nterm_cnt = 0
        nevid_cnt = 0
        header = True
        total_gene_cnt = None
        current_gene_cnt = 0
        current_terms_cnt = 0
        annot_offset = 0
        obo_offset = 0
        while f:
            if not ngene_code:
                line = f.readline()
                fields = line.split()
                total_gene_cnt = int(fields[0])
                total_terms_cnt = int(fields[1])
                ngene_code = int(fields[2])
                nterm_code = int(fields[3])
                nevid_code = int(fields[4])
                self.gene_code = ['' for _ in range(ngene_code)]
                self.term_code = ['' for _ in range(nterm_code)]
                self.evid_code = ['' for _ in range(nevid_code)]
            elif ngene_cnt < ngene_code:
                line = f.readline()
                self.gene_code[ngene_cnt] = line.strip()
                ngene_cnt += 1
            elif nterm_cnt < nterm_code:
                line = f.readline()
                self.term_code[nterm_cnt] = line.strip()
                nterm_cnt += 1
            elif nevid_cnt < nevid_code:
                line = f.readline()
                self.evid_code[nevid_cnt] = line.strip()
                nevid_cnt += 1
            else: # we're not in the header
                if header: offset = f.tell() 
                header = False
                try:
                    if current_gene_cnt < total_gene_cnt: # we are reading gene:terms annotations
                        peek = f.tell()
                        buf = f.read(calcsize('IIH'))
                        (gene_int, taxa_int, nterms) = unpack('IIH', buf)
                        current_gene_cnt += 1
                        if (not taxa) or (taxa_int == taxa or taxa_int in taxa):
                            self.annot_index[gene_int] = peek
                        bufsize = calcsize('?BI') 
                        f.read(nterms * bufsize)
                    elif current_terms_cnt < total_terms_cnt: # we are reading term definitions (term is_a term, term, term, ...) 
                        buf = f.read(calcsize('IcH'))
                        (term_int, onto_ch, nterms) = unpack('IcH', buf)
                        current_terms_cnt += 1
                        bufsize = calcsize('BI')
                        buf = f.read(nterms * bufsize)
                        terms_dict = {}
                        for pos in range(0, len(buf) - 1, bufsize):
                            (rel_ndx, sup_int) = unpack('BI', buf[pos:pos+bufsize])
                            terms_dict[sup_int] = rel_ndx
                        name_peek = f.tell()
                        f.readline() # skip putting name in memory, instead refer to the position in the file
                        self.termdefs[term_int] = (onto_ch, terms_dict, name_peek)
                    else:
                        buf = f.read(calcsize('II'))
                        (annot_offset, obo_offset) = unpack('II', buf)
                        break
                except error as inst:
                    print "Problem reading binary file: ", inst, "at gene ", current_gene_cnt, "at definition ", current_terms_cnt, "at", f.tell()
                    exit(3)
        print "Read %d genes and %d term definitions" % (current_gene_cnt, current_terms_cnt)
        print "Annotations start at", annot_offset, "\nDefinitions start at", obo_offset
        return f

    #FIXME: write code to perform test of taxa enrichment
    
    def getGOReport_byScore(self, gene_score_map, negatives_score_map = {}, include_more_general = True, descending_order = True):
        """ Generate a complete GO term report for a set of genes with associated scores. 
            Uses the Wilcoxon Ranksum test for each GO term to assign a p-value,
            indicating the enrichment of term to "top" genes in descending order by score (by default).
        """
        fg_map = self.getTerms(gene_score_map.keys(), include_more_general = include_more_general)
        fg_list = []
        for id in fg_map:
            for t in fg_map[id]:
                fg_list.append(t)
        term_set = set(fg_list)
        term_pval = {}
        if len(negatives_score_map) > 0:
            bg_map = self.getTerms(negatives_score_map.keys(), include_more_general = include_more_general)
        for t in term_set:
            pos = []
            neg = []
            for gene in gene_score_map: 
                annot = fg_map[gene]
                if not annot == None:
                    if t in annot:
                        pos.append(gene_score_map[gene])
                    else:
                        neg.append(gene_score_map[gene])
            if len(pos) > 0 and len(neg) > 0:
                if descending_order: 
                    p = stats.getRSpval(neg, pos)
                else:
                    p = stats.getRSpval(pos, neg)
            if len(negatives_score_map) > 0 and p <= 0.05:
                mpos = pos # scores of foreground genes with matching GO term
                mneg = [] # scores of background genes with matching GO terms 
                for gene in negatives_score_map: 
                    annot = bg_map[gene]
                    if not annot == None:
                        if t in annot:
                            mneg.append(negatives_score_map[gene])
                if len(mneg) > 0:
                    if descending_order: 
                        p2 = stats.getRSpval(mneg, mpos)
                    else:
                        p2 = stats.getRSpval(mpos, mneg)
                else:
                    p2 = 0.0
                term_pval[t] = (p, p2)
            else:
                term_pval[t] = (p, 1.0)
                
        sorted_pval = sorted(term_pval.items(), key=lambda v: v[1][0], reverse=False)

        ret = []
        for t in sorted_pval:
            defin = self.getTermdef(t[0])
            if defin == None:
                print 'Could not find definition of %s' % t[0]
            else:
                ret.append((t[0], t[1][0], t[1][1], defin[2].strip(), defin[0]))
        return ret
        
    def getGOReport(self, positives, background = None, taxa = None, include_more_general = True):
        """ Generate a complete GO term report for a set of genes (positives).
            Each GO term is also assigned an enrichment p-value (on basis of background, if provided).
            Returns a list of tuples (GO_Term_ID[str], Foreground_no[int], Term_description[str]) with no background, OR
            (GO_Term_ID[str], E-value[float], Foreground_no[int], Background_no[int], Term_description[str]). 
            E-value is a Bonferroni-corrected p-value.
            """
        pos = set(positives)
        fg_map = self.getTerms(pos, include_more_general = include_more_general)
        fg_list = []
        for id in fg_map:
            for t in fg_map[id]:
                fg_list.append(t)
        bg_map = {}
        bg_list = []
        neg = set()
        if background != None:
            neg = set(background).difference(pos)
            bg_map = self.getTerms(neg, include_more_general = include_more_general)
            for id in bg_map:
                for t in bg_map[id]:
                    bg_list.append(t)
        term_set = set(fg_list)
        term_cnt = {}

        nPos = len(pos)
        nNeg = len(neg)
        if background == None:
            for t in term_set:
                term_cnt[t] = fg_list.count(t)
            sorted_cnt = sorted(term_cnt.items(), key=lambda v: v[1], reverse=True)
        else: # a background is provided
            for t in term_set:
                fg_hit = fg_list.count(t)
                bg_hit = bg_list.count(t)
                fg_nohit = nPos - fg_hit
                bg_nohit = nNeg - bg_hit
                term_cnt[t] = (fg_hit, fg_hit + bg_hit, stats.getFETpval(fg_hit, bg_hit, fg_nohit, bg_nohit, False))
            sorted_cnt = sorted(term_cnt.items(), key=lambda v: v[1][2], reverse=False)

        ret = []
        for t in sorted_cnt:
            defin = self.getTermdef(t[0])
            if defin == None:
                print 'Could not find definition of %s' % t[0]
            else:
                if background != None:
                    ret.append((t[0], t[1][2] * len(term_set), t[1][0], t[1][0]+t[1][1], defin[2], defin[0]))
                else:
                    ret.append((t[0], t[1], defin[2], defin[0]))
        return ret

def encode(code_me, encode_strings):
    code = 0
    accum = 1
    try:
        for pos in range(len(code_me)):
            codelen = len(encode_strings[pos])
            for i in range(codelen):
                if encode_strings[pos][i] == code_me[pos]:
                    code += accum * i
                    accum *= codelen
                    break
    except IndexError as e:
        print e, code_me
    return code

def decode(code, encode_strings):
    npos    = len(encode_strings)
    accum   = [1 for _ in range(npos)]
    try:
        for pos in range(1, npos): accum[pos] = accum[pos - 1] * len(encode_strings[pos - 1]) 
        indices = [-1 for _ in range(npos)]
        for pos in range(npos - 1, -1, -1): # go backwards, start at last (most significant) position
            indices[pos] = code / accum[pos] 
            code -= accum[pos] * indices[pos] 
        string = [encode_strings[pos][indices[pos]] for pos in range(len(encode_strings))]
    except IndexError as e:
        print e, code
    return string

def _extractAnnotFields(line):
    fields = line.strip().split()
    offset = 0
    gene = fields[1]
    if len(gene) != 6:
        gene = None
    symb = fields[2]
    qual = (fields[3] != 'NOT')
    if not qual: offset += 1
    term = fields[3 + offset]
    if not term.startswith('GO:'):
        term = None
        for field in fields[4 + offset:]:
            offset += 1
            if field.startswith('GO:'):
                term = field
                break
    evid = fields[5 + offset]
    if not evid_codes.has_key(evid):
        evid = None
        for field in fields[6 + offset:]:
            offset += 1
            if evid_codes.has_key(field):
                evid = field
                break
    onto = fields[6 + offset]
    if not onto_codes.has_key(onto):
        onto = None
        for field in fields[7 + offset:]:
            offset += 1
            if onto_codes.has_key(field):
                onto = field
                break
    taxa = fields[9 + offset]
    if taxa.find('taxon:') == -1:
        taxa = None
        for field in fields[10 + offset:]:
            offset += 1
            if field.find('taxon:') > -1:
                taxa = field
                break
    if taxa != None:
        taxa_line = taxa.split(':')
        taxa = int(taxa_line[len(taxa_line) - 1]) # pick last taxon ID
    return (gene, symb, qual, term, evid, onto, taxa)

def readOBOFile(obofile):
    """
    http://www.geneontology.org/GO.format.obo-1_2.shtml
    """
    src = open(obofile, 'r')
    terms = {}
    in_term_def = False
    in_type_def = False
    for line in src:
        if in_term_def:
            if line.startswith('id: '):
                term_id = line[4:14]
                term_is = set()
            elif line.startswith('name: '):
                term_name = line[6:].strip()
            elif line.startswith('def: '):
                # Note this is a multi-line field, delimited by "'s
                pass
            elif line.startswith('namespace: '):
                if   line[11] == 'b': term_onto = 'P'
                elif line[11] == 'm': term_onto = 'F'
                elif line[11] == 'c': term_onto = 'C'
            elif line.startswith('is_a: '):
                term_is.add((line[6:16], 'is_a'))
            elif line.startswith('relationship: '):
                fields = line.split()
                term_is.add((fields[2], fields[1]))
            elif line.startswith('intersection_of: '):
                fields = line.split()
                if fields[1].startswith('GO:'):
                    term_is.add((fields[1], 'isect'))
                else:
                    term_is.add((fields[2], fields[1]))
            elif line.startswith('is_obsolete: '):
                in_term_def = False # ignore this entry
        if line.startswith('[Term]'):
            if in_term_def: # already defining one, stash it before moving on to the next...
                terms[term_id] = (term_name, term_onto, term_is)
            elif in_type_def:
                in_type_def = False
            in_term_def = True
        if line.startswith('[Typedef]'):
            if in_term_def: # already defining one, stash it before moving on to the next...
                in_term_def= False
            in_type_def = True
    if in_term_def: #  defining one, stash it
        terms[term_id] = (term_name, term_onto, term_is)
    return terms

def writeBitFile(annotFile, obofile, destFile, taxas = None):
    print "Started at", time.asctime()
    # open annotation file to analyse and index data 
    src = open(annotFile, 'r')
    gene_index = [{} for _ in range(6)]  # count different characters in different positions 
    term_index = [{} for _ in range(7)]  # count different characters in different positions
    evid_index = {}
    gene_cnt = 0
    cnt = 0
    prev_gene = None
    for line in src:
        cnt += 1
        #if cnt > 100000:
        #    break
        if line.startswith('!'):
            continue
        (gene, symb, qual, term, evid, onto, taxa) = _extractAnnotFields(line)
        if not (taxas and ((taxa == taxas) or (taxa in taxas))):  # The gene does NOT belong to a nominated taxon
            continue
        if not (gene == prev_gene): # not the same gene
            gene_cnt += 1
        try:
            evid_index[evid]
        except:
            evid_index[evid] = len(evid_index)
        pos = 0
        for ch in gene[0:6]:
            try:
                gene_index[pos][ch]
            except KeyError: # no match
                gene_index[pos][ch] = len(gene_index[pos])
            pos += 1
        pos = 0
        for ch in term[3:10]:
            try:
                term_index[pos][ch]
            except KeyError: # no match
                term_index[pos][ch] = len(term_index[pos])
            pos += 1
        prev_gene = gene
    src.close()
    print "Read annotations for %d genes" % gene_cnt

    gene_code = ['' for _ in range(6)]
    term_code = ['' for _ in range(7)]
    for d in range(len(gene_index)):
        arr = ['?' for _ in gene_index[d]]
        for (ch, index) in gene_index[d].items():
            arr[index] = ch
        gene_code[d] = ''.join(arr)
    for d in range(len(term_index)):
        arr = ['?' for _ in term_index[d]]
        for (ch, index) in term_index[d].iteritems():
            arr[index] = ch
        term_code[d] = ''.join(arr)
    evid_code = ['' for _ in range(len(evid_index))]
    for (e, ndx) in evid_index.items():
        evid_code[ndx] = e
    
    # Get GO definitions
    terms = readOBOFile(obofile)
    print "Read %d GO definitions" % len(terms)
    
    # re-open, now with the aim of copying info
    src = open(annotFile, 'r')
    dst = open(destFile, 'w')
    # STEP 1: header info
    dst.write("%d\t%d\t%d\t%d\t%d\n" % (gene_cnt, len(terms), len(gene_code), len(term_code), len(evid_index)))
    for code_str in gene_code:
        dst.write(code_str+"\n")
    for code_str in term_code:
        dst.write(code_str+"\n")
    for e_str in evid_code:
        dst.write(e_str+'\n')
    print "Wrote header %d\t%d\t%d\t%d\t%d, now at @%d" % (gene_cnt, len(terms), len(gene_code), len(term_code), len(evid_index), dst.tell())

    # STEP 2: write annotations
    annot_offset = dst.tell()
    prev_gene = None
    concat_terms = {}
    cnt = 0
    for line in src:
        cnt += 1
        #if cnt > 100000:
        #    break
        if line.startswith('!'):
            continue
        (gene, symb, qual, term, evid, onto, taxa) = _extractAnnotFields(line)
        if not (taxas and ((taxa == taxas) or (taxa in taxas))): # The gene does NOT belong to a nominated taxon
            continue
        if gene != prev_gene: # new gene is found
            if prev_gene != None:
                # write prev data
                s = pack('IIH', encode(prev_gene, gene_code), taxa, len(concat_terms))
                dst.write(s)
                for t in concat_terms:
                    (o, q, e) = concat_terms[t]
                    s = pack('?BI', q, evid_index[e], encode(t, term_code))
                    dst.write(s)
            # re-init
            prev_gene = gene
            concat_terms = {}
        concat_terms[term[3:]] = (onto, qual, evid) 
    if len(concat_terms) > 0:
        # write data in buffer
        s = pack('IIH', encode(prev_gene, gene_code), taxa, len(concat_terms))
        dst.write(s)
        for t in concat_terms:
            (o, q, e) = concat_terms[t]
            s = pack('?BI', q, evid_index[e], encode(t, term_code))
            dst.write(s)
    print "Wrote GO annotations, now at @%d" % dst.tell()

    # Next, the ontology definition...
    obo_offset = dst.tell() # remember the position where the OBO starts
    sorted_terms = sorted(terms.iteritems(), key=operator.itemgetter(0))
    for [t, _] in sorted_terms:
        (term_name, term_onto, term_is) = terms[t]
        s = pack('IcH', encode(t[3:], term_code), term_onto, len(term_is))
        dst.write(s)
        for (sup_term, sup_rel) in term_is:
            try:
                index = onto_rel.index(sup_rel)
            except ValueError:
                index = 9
            s = pack('BI', index, encode(sup_term[3:], term_code))
            dst.write(s)
        dst.write(term_name + '\n')
    print "Wrote %d GO definitions, now at @%d" % (len(sorted_terms), dst.tell())
    
    # Finally, write the offsets to quickly access annotations and definitions, resp
    dst.write(pack('II', annot_offset, obo_offset))
    # done, close
    dst.close()
    print "Completed at", time.asctime()

if __name__ == '__main__0':
    writeBitFile('/Users/mikael/simhome/TNR/GO/gene_association.tair', 
                 '/Users/mikael/simhome/TNR/GO/gene_ontology_ext.obo',
                 '/Users/mikael/simhome/TNR/GO/tair.bit')
    """
    Started at Tue Jul 17 09:24:25 2012
    Read annotations for 14973326 genes
    Read 35980 GO definitions
    Wrote header 14973326    35980    6    7    19, now at @302
    Wrote GO annotations, now at @773175002
    Wrote 35980 GO definitions, now at @775623250
    Completed at Tue Jul 17 10:59:05 2012
    """
    bgo = BinGO('/Users/mikael/simhome/TNR/GO/tair.bit')
    #bgo = BinGO('/Users/mikael/simhome/gene_association.bit', taxa = 39947)
    print "Done loading index with %d genes annotated" % len(bgo.annot_index)

if __name__ == '__main__':
    print os.getcwd()
    #writeBitFile('/Users/mikael/simhome/gene_association.goa_uniprot', 
    #             '/Users/mikael/simhome/gene_ontology_ext.obo',
    #             '/Users/mikael/simhome/gene_association_mammal.bit',
    ##             [9606,10090])
    """
    Started at Tue Jul 17 09:24:25 2012
    Read annotations for 14973326 genes
    Read 35980 GO definitions
    Wrote header 14973326    35980    6    7    19, now at @302
    Wrote GO annotations, now at @773175002
    Wrote 35980 GO definitions, now at @775623250
    Completed at Tue Jul 17 10:59:05 2012
    """
    bgo = BinGO('/Users/mikael/simhome/gene_association.bit')
    #bgo = BinGO('/Users/mikael/simhome/gene_association.bit', taxa = 39947)
    print "Done loading index with %d genes annotated" % len(bgo.annot_index)
    #pos = [id.strip() for id in open('/Users/mikael/simhome/nls/identifiers_streptophyta_pos.txt')]
    #neg = [id.strip() for id in open('/Users/mikael/simhome/nls/identifiers_streptophyta_neg.txt')]
    #rows = bgo.getGOReport(pos, background = neg, taxa = 39947, include_more_general = True)
    f_bg = open('/Users/mikael/simhome/homoaa/uniprotID/S288C-sgd.proteinID.list')
    f_fg = open('/Users/mikael/simhome/homoaa/uniprotID/S288C-sgd.homo')
    bg = [s.strip() for s in f_bg]
    print 'Background has %d proteins' % len(bg)
    fg = []
    for tract in f_fg:
        field = tract.split('\t')
        id = field[1]
        aa_cnt = int(field[5])
        tn_cnt = int(field[6])
        if aa_cnt >= 7:
            fg.append(id)
    print 'Foreground has %d proteins' % len(fg)
    f_bg.close()
    f_fg.close()
    rows = bgo.getGOReport(fg, bg, include_more_general = True)
    for row in rows[0:100]:
        if len(row) > 4:
            print "%s\t%4.2E\t%3d\t%6d\t%s (%s)" % (row[0], row[1], row[2], row[3], row[4].strip(), row[5])
        else:
            print "%s\t%3d\t%s (%s)" % (row[0], row[1], row[2].strip(), row[3])
