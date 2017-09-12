#!/usr/bin/env python2.7

'''
gkb2FiveColTbl
This script converts a GenBank file (.gbk or .gb) into a Sequin 5-column feature table (.tbl).
PACKAGE REQUIREMENTS:
    - BioPython
    - argparse
USAGE:
    python2 gbk2FiveColTbl.py -i gbk-file 2> stderr
INPUTS:
    A GenBank file.
OUTPUTS:
    infile_prefix.tbl: the five-column feature table
    infile_prefix.fsa: the DNA sequence in FASTA format
    infile_prefix.errors: an error file
ARGUMENTS:
    --i / -infile: gbk-file
'''

###############
# AUTHOR INFO #
###############
__author__ = 'Michael Gruenstaeudl <m.gruenstaeudl@fu-berlin.de>'
__copyright__ = 'Copyright (C) 2017 Michael Gruenstaeudl'
__info__ = 'gkb2FiveColTbl'
__version__ = '2017.09.08.1900'

########
# TODO #
########
''' foo bar '''

#############
# DEBUGGING #
#############
import pdb
#pdb.set_trace()

#####################
# IMPORT OPERATIONS #
#####################
from Bio import SeqIO

import argparse
import Bio
import collections
import os

####################
# GLOBAL VARIABLES #
####################
allowed_qualifiers = ['gene', 'product', 'transl_table', 'protein_id', 'exception', 'note', 'transl_except', 'trans_splicing'] # Do NOT change order of items in list! The list is position-sensitive! See section "ADD QUALIFIER INFO".

min_feature_len = 6
transl_table = 11
valid_start_codon = 'ATG'
valid_stop_codons = ['TAG', 'TAA', 'TGA']

###########
# CLASSES #
###########
class MyException(Exception):
    pass

###############
# DEFINITIONS #
###############
def open_unless_exists(*args, **kwargs):
    ''' This functions wraps open() to check if file exists. '''
    if os.path.exists(args[0]):
        exit('\nERROR: File `%s` already exists.' %(args[0]))
    f = open(*args, **kwargs)
    return f

def is_feature_transl_valid(seq, feature, error_fh):
    ''' Foo bar '''
    extracted_seq = feature.location.extract(seq)
    loc = feature.location
    feature_name = '%d %d %s' % (loc.nofuzzy_start+1, loc.nofuzzy_end, feature.type) ## Improvement necessary !

    if feature.type == 'CDS':
        try:
            validate_translation(extracted_seq, feature_name, test_internal_stops=True)
        except MyException as e:
            print >> error_fh, '\n%s\nLocation: %s\nFeatureType: %s\nQualifiers: %s\nSequence:\n%s\n' % (e, loc, feature.type, feature.qualifiers, str(extracted_seq))
        
    if feature.type == 'gene':
        #if 'trn' not in feature.qualifiers['gene'][0] and 'trans_splicing' not in feature.qualifiers:
        try:
            if 'trn' not in feature.qualifiers['gene'][0] and 'rrn' not in feature.qualifiers['gene'][0]:
                try:
                    validate_translation(extracted_seq, feature_name, test_internal_stops=False) # Do not test internal stop codons, because "gene" attributes do not have info on start and end position of introns stored in them; thus, they would falsely indicate internal stop codons when translated.
                except MyException as e:
                    print >> error_fh, '\n%s\nLocation: %s\nFeatureType: %s\nQualifiers: %s\nSequence:\n%s\n' % (e, loc, feature.type, feature.qualifiers, str(extracted_seq))
        except KeyError, e:
            print >> error_fh, '\nAttribute not found: %s\nLocation: %s\nFeatureType: %s\nQualifiers: %s\n' % (e, loc, feature.type, feature.qualifiers)

def validate_translation(extracted_seq, feature_name, test_internal_stops):
    ''' This function performs three checks on a coding region (CDS): 
    (i) tests if the feature starts with a valid start codon, 
    (ii) tests if the feature ends with a valid stop codon, and 
    (iii) tests if the sequence length is a multiple of three and that 
    there is only a single in frame stop codon (namely at the end).
    If any of these tests fail, an exception is raised immediately. 
    Step (iii) should not be tested on "gene" attributes, as such attributes do not have info on start and end position of introns stored in them; thus, they would falsely indicate internal stop codons when translated.'''
    
    # Testing if the feature starts with a valid start codon
    if extracted_seq.startswith(valid_start_codon): # Function "startswith()" returns True/False
        pass
    else:
        raise MyException('ERROR: Feature `%s` does not start with a valid start codon.' % (feature_name))

    # Testing if the feature ends with a valid stop codon
    if any([extracted_seq.endswith(stop_codon) for stop_codon in valid_stop_codons]): # Function "endswith()" returns True/False
        pass
    else:
        raise MyException('ERROR: Feature `%s` does not end with a valid stop codon.' % (feature_name))

    if test_internal_stops:
        # Testing if the sequence length is a multiple of three and that there is only a single in frame stop codon (namely at the end)
        try:
            result = extracted_seq.translate(table=transl_table, to_stop=False, cds=True)
        except Exception, e:
            transl_err = str(extracted_seq.translate(table=transl_table, to_stop=False, cds=False)) # enclosing "str()" is important !
            raise MyException('\nEXCEPTION MESSAGE: %s \nERROR: Feature `%s` has an issue: either length of sequence is not equal to a multiple of three or sequence has an internal stop codon.\nTranslation:\n%s' % (e, feature_name, transl_err))


def make_table(inFn, infmt='genbank'):
    
    # Define output file names
    file_stem = os.path.splitext(inFn)[0]
    fasta_outFn = file_stem + ".fasta"
    feature_outFn = file_stem + ".tbl"
    error_outFn = file_stem + ".errors"

    with open_unless_exists(feature_outFn, 'w') as feature_fh, \
    open_unless_exists(fasta_outFn, 'w') as fasta_fh, \
    open_unless_exists(error_outFn, 'w') as error_fh:
        
        # READ INDATA
        rec = SeqIO.read(inFn, infmt) # Compared to function "parse", "read" expects only a single record per gbk file
        
        # EXTRACT RECORD NAME
        try:
            record_name = rec.id.split('.')[0] or rec.annotations['accessions'][0] # NCBI number is preferential over taxon name (as extracted in except)
        except:
            try:
                record_name = rec.annotations['organism'].replace(' ', '_') or rec.name
            except:
                exit('\nERROR: Cannot extract a valid record name from file `%s`.' %(args[0]))
        
        # GROUP SEQFEATURES BY FEATURELOCATION
        featuredict_unsort = collections.OrderedDict()
        for f in rec.features:
            if isinstance(f, Bio.SeqFeature.SeqFeature) and f.location:
                loc_tag = '[%i:%i](%s)' % (f.location.start, f.location.end, f.location.strand)
                if loc_tag not in featuredict_unsort.keys():
                    featuredict_unsort[loc_tag] = []
                if loc_tag in featuredict_unsort.keys(): # Do not use else-statement, as updated featuredict_unsort would not be used
                    featuredict_unsort[loc_tag].append(f)
        
        # SORT EACH FEATURELOCATION SUCH THAT ATTRIBUTE VALUE "GENE" COMES FIRST
        featuredict_sorted = collections.OrderedDict()
        for k, v in featuredict_unsort.iteritems():
            if v[0].type: # See if first element in values list has attribute type
                sorted_v = []
                sorted_v.extend([e for e in v if e.type == 'gene'])
                sorted_v.extend([e for e in v if e.type != 'gene'])
                featuredict_sorted[k] = sorted_v
            if not v[0].type:
                featuredict_sorted[k] = v
        
        # GENERATE OUTFILE: FASTA FILE
        print >> fasta_fh, '> %s' % (record_name)
        print >> fasta_fh, '%s' % (rec.seq.lower())

        # GENERATE OUTFILE: FEATURE TABLE
        print >> feature_fh, '>Feature %s' % (record_name) # Write the organism name as the first line of the feature table
        
        for fts in featuredict_sorted.values():
            for f in fts:
                
                if isinstance(f.location, Bio.SeqFeature.FeatureLocation):
                    if f.location.strand == 1:
                        is_feature_transl_valid(rec.seq, f, error_fh)
                        print >> feature_fh, '%d\t%d\t%s' % (f.location.nofuzzy_start+1, f.location.nofuzzy_end, f.type)
                    if f.location.strand == -1:
                        is_feature_transl_valid(rec.seq, f, error_fh)
                        print >> feature_fh, '%d\t%d\t%s' % (f.location.nofuzzy_end, f.location.nofuzzy_start+1, f.type)
                    
                if isinstance(f.location, Bio.SeqFeature.CompoundLocation):
                    n_exons = len(f.location.parts)
                    
                    is_feature_transl_valid(rec.seq, f, error_fh) # validation of translation can only be performed on full feature, not on individual parts
                    if f.location.strand == 1:
                        for i in range(0,1):
                            p = f.location.parts[i]
                            print >> feature_fh, '%d\t%d\t%s' % (p.start+1, p.end, f.type)
                        for i in range(1,n_exons):
                            p = f.location.parts[i]
                            print >> feature_fh, '%d\t%d' % (p.start+1, p.end)
                    if f.location.strand == -1:
                        for i in range(0,1):
                            p = f.location.parts[i]
                            print >> feature_fh, '%d\t%d\t%s' % (p.end, p.start+1, f.type)
                        for i in range(1,n_exons):
                            p = f.location.parts[i]
                            print >> feature_fh, '%d\t%d' % (p.end, p.start+1)

                # INTEGRATE MIN FEATURE LENGTH
                #min_feature_len

                # ADD STRING IF QUALIFIERS INCOMPLETE
                if (f.type == 'CDS') and ('product' not in f.qualifiers):
                    f.qualifiers['product'] = ['hypothetical protein']

                # ADD QUALIFIER INFO
                for q in allowed_qualifiers:  # Loop over allowed_qualifiers, not over f.qualifiers, because only looping over the first ensures that the order of alloed_qualifiers is used (e.g., that "product" is first attribute for CDS).
                    if q in f.qualifiers.keys():
                        vals = f.qualifiers[q]
                        for v in vals:
                            print >> feature_fh, '%s\t%s\t%s\t%s\t%s' % ('', '', '', q, v)
                
            # ADD NEWLINE AFTER EVERY GENE
            print >> feature_fh, '\n'

    # INSTRUCTIONS FOR NECESSARY MANUAL STEPS
    print   '\nPERSONAL POST-ANALYSIS INSTRUCTIONS (2017.09.08):'\
            '\nDo the following steps to finish the feature table:'\
            '\n(a) Move one of the two copies of "rps12" (which is a transspliced gene) to the top of the table; the other copy shall be listed in its position.'\
            '\nConvert "trans-splicing" to "exception \t trans_splicing"'\
            '\n(c) Correct "repeat_region" and "misc_features"'\
            '\n(d) Remove "source" feature (if present).'\
            '\n--END OF RUN--\n'

############
# ARGPARSE #
############
if __name__ == '__main__':
    
    parser = argparse.ArgumentParser(description="  --  ".join([__author__, __copyright__, __info__, __version__]))
    # REQUIRED
    parser.add_argument('-i',
                        '--infile',
                        type = str,
                        help='absolute path to infile; infile in gbk format; Example: /path_to_input/test.gbk',
                        default='/home/username/Desktop/test.gbk',
                        required=True)
    # OPTIONAL
#    parser.add_argument('--topol',
#                        help='`circular` or `linear`', 
#                        default='linear',
#                        required=False)
    parser.add_argument('--version', 
                        help='Print version information and exit',
                        action='version',
                        version='%(prog)s ' + __version__)
    args = parser.parse_args()

########
# MAIN #
########
    make_table(args.infile)
    
#    try:
#        make_table(args.infile)
#    except:
#        print "ERROR: %s" %(args.infile)
