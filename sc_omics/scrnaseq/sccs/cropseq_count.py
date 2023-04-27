#!/usr/bin/env python3
'''
Copyright (c) 2019 Wei Li laboratory (weililab.org)
This code is free software; you can redistribute it and/or modify it
under the terms of the BSD License (see the file COPYING included with
the distribution).

Special thanks to Andrew Hill on part of the source codes.
'''

import sys
import re
import glob
import os
import itertools
import argparse
import logging

def crop_parseargs():
  """
  Parsing arguments
  """
  parser = argparse.ArgumentParser(description='cropseq-count: counting sgRNAs from CROPseq experiments.')
  
  parser.add_argument('-v', '--version',action='version',version='0.1')
  
  parser.add_argument('-n','--output-prefix',default='sample1',help='The prefix of the output file(s). Default sample1.')

  parser.add_argument('--lib-grna',required=True,help='A gRNA library file containing the list of sgRNA names, their sequences and associated genes, separated by tab.')
  parser.add_argument('--file-type',required=True,choices=['fastq','bam','paired-fastq'],default='fastq',help='The type of the file to search guide RNA sequence. ')
  parser.add_argument('--no-reverse-complement',action='store_true',help='Do not perform reverse complement search of sgRNAs.')
  parser.add_argument('-m','--max-mismatch',type=int,default=2,help='Maximum number of mismatches to be considered in sgRNA search. Default 2. Not recommended for values greater than 2. Decrease this value to speed up the search.')
  parser.add_argument('--anchor-before',default='GAAACACCG',help='Anchor sequence before the sgRNA. Default GAAACACCG (at the end of U6 promoter).')
  parser.add_argument('--anchor-after',default='GTTTTAGAG',help='Anchor sequence after the sgRNA. Default GTTTTAGAG.')
  parser.add_argument('--files',nargs='+',help='A list of fastq files, SAM/BAM files, or a wildcard of filenames.')
  


  # post-processing
  args=parser.parse_args()
    
  logging.basicConfig(level=10,
    format='%(levelname)-5s @ %(asctime)s: %(message)s ',
    datefmt='%a, %d %b %Y %H:%M:%S',
    # stream=sys.stderr,
    filename=args.output_prefix+'.log',
    filemode='w'
  )
  console = logging.StreamHandler()
  console.setLevel(logging.INFO)
  # set a format which is simpler for console use
  formatter = logging.Formatter('%(levelname)-5s @ %(asctime)s: %(message)s ','%a, %d %b %Y %H:%M:%S')
  #formatter.formatTime('%a, %d %b %Y %H:%M:%S')
  # tell the handler to use this format
  console.setFormatter(formatter)
  # add the handler to the root logger
  logging.getLogger('').addHandler(console)
  
  # add paramters
  logging.info('Parameters: '+' '.join(sys.argv))
  
  return args


# two versions of rev comp, depending on different versions of python
'''
Reverse complement
'''
if sys.version_info >= (3,1):
  trans_table=str.maketrans("ACGT","TGCA")
  def count_revcomp(x):
    return x.translate(trans_table)[::-1]
else:
  trans_table=string.maketrans("ACGT","TGCA")
  def count_revcomp(x):
    return x.translate(trans_table)[::-1]
#u6_before='GAAACACCG'
#u6_after='GTTTTAGAG'
# get grna_list
#grna_list_file='/groups/ligrp/weili/JinZhang/180912_RawDataForSingleCell/180912/metadata/grna_list.txt'

def gen_mismatches(sequence,num_mismatches):
  """
  generate mismatches of a certain sequence
  borrow from https://github.com/shendurelab/single-cell-ko-screens/blob/master/get_barcodes.py
  """
  letters = 'ACGT'
  mismatches=[]
  for locs in itertools.combinations(range(len(sequence)), num_mismatches):
      sequence_list = [[char] for char in sequence]
      for loc in locs:
          orig_char = sequence[loc]
          sequence_list[loc] = [l for l in letters if l != orig_char]

      for poss in itertools.product(*sequence_list):
          mismatches.append(''.join(poss))

  return mismatches

def gen_mismatch_dict(grna_table,num_mismatches):
  """
  Generate a dictionary of grna and its mismatches
  """

  grdict={}
  gmismatchdict={}
  g_orig_seq_cnt={}
  g_mis_seq_cnt={}
  for gr in grna_table:
    orig_seq=gr[1].upper()
    if orig_seq in g_orig_seq_cnt:
      logging.warning('Identical sequences '+orig_seq+ ' for '+gr[0]+' and '+g_orig_seq_cnt[orig_seq]+'. Skip '+gr[0]+'.')
      continue
    g_orig_seq_cnt[orig_seq]=gr[0]
    #grdict[gr[0]]=u6_before+gr[1]+u6_after
    grdict[gr[0]]=[orig_seq]
    nms=gen_mismatches(orig_seq,num_mismatches)
    gmismatchdict[gr[0]]=nms
    for sq in nms+[orig_seq]:
      g_mis_seq_cnt[sq]=(lambda s: g_mis_seq_cnt[s]+1 if s in g_mis_seq_cnt else 1)(sq)
  # remove those that overlap with each other
  n_skipped=0
  n_totalcnt=0
  for gr in grna_table:
    if gr[0] not in grdict:
      continue
    orig_seq=gr[1].upper()
    nms=gen_mismatches(orig_seq,num_mismatches)
    for sq in nms:
      if g_mis_seq_cnt[sq]>1:
        n_skipped=1
      else:
        gmismatchdict[gr[0]]+=[sq]
        n_totalcnt+=1
  logging.info('Total sgRNAs:'+str(len(grdict)))
  logging.info('Total sgRNAs with mismatches:'+str(n_totalcnt))
  logging.info(str(n_skipped)+ ' mismatched sequences are skipped due to overlap with other sequences.')
  
  
  return (grdict,gmismatchdict)

def process_input_files(args):
 
  in_fqfilelist=[]
  for sf in args.files:
    if os.path.isfile(sf):
      in_fqfilelist+=[sf]
    else:
      for sf2 in glob.glob(sf):
        if os.path.isfile(sf2):
          in_fqfilelist+=[sf2]
        else:
          logging.error(sf+' does not exist.')
          sys.exit(-1)
    
  return in_fqfilelist

def process_library_file(args): 
  grna_list_file=args.lib_grna
  grna=[line.strip().split() for line in open(grna_list_file)]
  ng=0
  for grl in grna:
    ng=ng+1
    if len(grl)<3:
      logging.error('Error in processing line'+str(ng)+':')
      logging.error('\t'.join(grl))
      logging.error('Library file must be at least three fields including: gRNA_ID, sequence, gene_ID.')
      sys.exit(-1)
  grna=grna[1:] # skip the header
  logging.info('Total gRNA:'+str(len(grna)))
  return grna

def search_sequence(line,grdict,gmismatchdict,args):
  """
  Perform sequence search between anchor_before and anchor_after
  Parameters
    line
      a sequence to be search
    grdict
    gmismatchdict
      a dictionary of {sgrna_name:[sequence]}
    args
      arguments
  Returns
    (found_rec, found_seq)
    found_rec: True/False whether a hit was found
    found_seq: the list of matched gRNA sequences. May contain more than 1 sequence 
	  
  """
  #line=count_revcomp(line)
  # u6 before or u6 after must be present
  u6_before=args.anchor_before
  u6_after=args.anchor_after
  if u6_before not in line and u6_after not in line:
    return (False,[])
  # first, count perfect matches
  found_rec=False
  found_seq=[]
  for (gn,gseqlist) in grdict.items():
    for gseq in gseqlist:
      if u6_before+gseq in line or gseq+u6_after in line:
        #gr_count[gn]+=1
        found_seq+=[gn]
        found_rec=True
        break
  # if perfect matches are not found, go to the next time-consuming step for mismatches 
  if found_rec == False:
    for (gn,gseqlist) in gmismatchdict.items():
      for gseq in gseqlist:
        if u6_before+gseq in line or gseq+u6_after in line:
          #gr_count[gn]+=1
          found_seq+=[gn]
          found_rec=True
          break
  return (found_rec,found_seq)
  
def process_pair_fastq_files(in_fqfile1,in_fqfile2,grdict,gmismatchdict,args):
    """
    Process a paired-end fastq files and search for possible hits.
    In 10X genomics platform, fq_file1 is usually the the file that contains molecular and cellular barcode, and fq_file2 is the actual RNA sequence
    """
    # use to determine cell barcode and umi barcode start and end; should place to argparse later
    cell_barcode_start=0
    cell_barcode_end=15
    umi_barcode_start=16
    umi_barcode_end=25
    
    nl=0
    #gr_count={s:0 for s in grdict.keys()}
    gr_count_universal={} # return {cell_bc:{sgRNA:{umi_bc:cnt}}} structure
    in_fqfile1_f=open(in_fqfile1)
    
    for line in open(in_fqfile2):
      nl+=1
      line_1st=in_fqfile1_f.readline().strip()
      if nl%4!=2:
        continue
      (hashit,found_seqlist)=search_sequence(line,grdict,gmismatchdict,args)
      if hashit==False and args.no_reverse_complement == False:
        line=count_revcomp(line)
        (hashit,found_seqlist)=search_sequence(line,grdict,gmismatchdict,args)
      if hashit:
        cell_bc=line_1st[cell_barcode_start:(cell_barcode_end+1)]
        umi_bc=line_1st[umi_barcode_start:(umi_barcode_end+1)]
        if cell_bc not in gr_count_universal:
          gr_count_universal[cell_bc]={}
        #if umi_bc not in gr_count_universal[cell_bc]:
        #  gr_count_universal[cell_bc][umi_bc]=0
        gr_count=gr_count_universal[cell_bc]
        for fs in found_seqlist:
          if fs not in gr_count:
            gr_count[fs]={}
          if umi_bc not in gr_count[fs]:
            gr_count[fs][umi_bc]=0
          gr_count[fs][umi_bc]+=1
    # end for
    nct=0
    for (cell_bc,c_dict) in gr_count_universal.items():
      for (sg,umicnt) in c_dict.items():
        ct=len(umicnt) # this is the number of umis
        if ct>0:
          print(cell_bc+'\t'+sg+'\t'+str(ct))
          nct+=ct
    logging.info(str(nct)+' umis found..')
    return gr_count_universal
 

def process_fastq_files(in_fqfile,grdict,gmismatchdict,args):
    """
    Process a single-end fastq file and search for possible hits
    """
    nl=0
    gr_count={s:0 for s in grdict.keys()}
    for line in open(in_fqfile):
      nl+=1
      if nl%4!=2:
        continue
      (hashit,found_seqlist)=search_sequence(line,grdict,gmismatchdict,args)
      for fs in found_seqlist:
        gr_count[fs]+=1
      if hashit==False and args.no_reverse_complement == False:
        line=count_revcomp(line)
        (hashit,found_seqlist)=search_sequence(line,grdict,gmismatchdict,args)
        for fs in found_seqlist:
          gr_count[fs]+=1
    
    nct=0
    for (sg,ct) in gr_count.items():
      if ct>0:
        print(sg+'\t'+str(ct))
        nct+=ct
    logging.info(str(nct)+' hits found..')
    return gr_count

def get_total_reads(bamfile):
    """
    Get the total number of reads from a BAM file.
    """
    import pysam
    read_counts = 0
    index_stats = pysam.idxstats(bamfile)

    # Versions of pysam differ in their output so check to make sure have a list
    if not isinstance(index_stats, list):
        index_stats = index_stats.split('\n')

    for line in index_stats:
        entries = line.strip().split('\t')

        if len(entries) == 4:
            read_counts += int(entries[2]) + int(entries[3])

    return read_counts


def process_bam_files(in_bamfile,grdict,gmismatchdict,args):
    """
    Processing bam files
    """
    import pysam
    total_reads = get_total_reads(in_bamfile)
    logging.info('Total reads:'+str(total_reads))

    read_number = 0
    gr_count_universal={} # return {cell_bc:{sgRNA:{umi_bc:cnt}}} structure

    for read in pysam.Samfile(in_bamfile):
        read_number += 1
        if read_number %1000000 ==0:
            logging.info('Processing '+str(int(read_number/1000000))+'/'+str(int(total_reads/1000000))+'M reads')
        # Ignore mapped reads (can't belong to guide transcripts...)
        if not read.is_unmapped:
            continue

        seq = read.seq.upper()
        tags = dict(read.tags)
        cell_bc = tags.get('CB', None)
        umi_bc = tags.get('UB', None)

        if not cell_bc or not umi_bc:
            continue
        # 
        # search sequences
        (hashit,found_seqlist)=search_sequence(seq,grdict,gmismatchdict,args)
        if hashit==False and args.no_reverse_complement == False:
          seq=count_revcomp(seq)
          (hashit,found_seqlist)=search_sequence(seq,grdict,gmismatchdict,args)
        if hashit:
          if cell_bc not in gr_count_universal:
            gr_count_universal[cell_bc]={}
          #if umi_bc not in gr_count_universal[cell_bc]:
          #  gr_count_universal[cell_bc][umi_bc]=0
          gr_count=gr_count_universal[cell_bc]
          for fs in found_seqlist:
            if fs not in gr_count:
              gr_count[fs]={}
            if umi_bc not in gr_count[fs]:
              gr_count[fs][umi_bc]=0
            gr_count[fs][umi_bc]+=1
    # end for
    nct=0
    for (cell_bc,c_dict) in gr_count_universal.items():
      for (sg,umicnt) in c_dict.items():
        ct=len(umicnt) # this is the number of umis
        if ct>0:
          print(cell_bc+'\t'+sg+'\t'+str(ct))
          nct+=ct
    logging.info(str(nct)+' UMIs found..')
    return gr_count_universal
 



 
def output_to_file(args,gr_count_dict,grdict,librecord=None):
  """
  Write count output to file
  """  
  out_file=args.output_prefix+'.output.txt'
  outff=open(out_file,'w')
  if args.file_type=='fastq':
    grdict_k=[x for x in grdict.keys()]
    print('\t'.join(['sample']+grdict_k),file=outff)
    for in_fqfile in gr_count_dict.keys():
      gr_count=gr_count_dict[in_fqfile]
      print('\t'.join([in_fqfile]+[str(gr_count[sg]) for sg in grdict_k]),file=outff)
  elif args.file_type=='paired-fastq' or args.file_type=='bam':
    # this would be a {cell_bc:{sgrna:{umi:count}}} structure
    if librecord!=None:
      sg_lib={}
      for sl in librecord:
        sg_lib[sl[0]]=(sl[1],sl[2])
    else:
      sg_lib={}
    print('cell\tbarcode\tsgrna\tgene\tread_count\tumi_count',file=outff)
    for (cell_bc, c_dict) in gr_count_dict.items():
      for (sg, sgcnt) in c_dict.items():
        umi_count=len(sgcnt)
        rcount=0
        for (umi, cnt) in sgcnt.items():
          rcount+=cnt
        if sg in sg_lib:
          (sgseq,geneid)=sg_lib[sg]
        else:
          sgseq='NA'
          geneid='NA'
        print('\t'.join([cell_bc,sg,sgseq,geneid,str(rcount),str(umi_count)]),file=outff)
  
  outff.close()


  
if __name__ == '__main__':
  args=crop_parseargs()
  in_fqfilelist=process_input_files(args)

  logging.info('Total number of files:'+str(len(in_fqfilelist)))
  logging.info(args.files)
  grna=process_library_file(args) 

  (grdict,gmismatchdict)=gen_mismatch_dict(grna,args.max_mismatch)



  gr_count_dict={}

  # process file
  if args.file_type=='fastq':
    for in_fqfile in in_fqfilelist:
      logging.info('Processing file '+in_fqfile)
      gr_count=process_fastq_files(in_fqfile,grdict,gmismatchdict,args)
      gr_count_dict[in_fqfile]=gr_count
  if args.file_type=='bam':
      if len(in_fqfilelist)!=1:
        logging.error('Only one bam file is accepted for bam file type.')
        sys.exit(-1)
      in_fqfile = in_fqfilelist[0] 
      logging.info('Processing file '+in_fqfile)
      gr_count_dict=process_bam_files(in_fqfile,grdict,gmismatchdict,args)
      #gr_count_dict[in_fqfile]=gr_count
  elif args.file_type=='paired-fastq':
    if len(in_fqfilelist)%2 != 0 or len(in_fqfilelist)<2:
      logging.error('Must provide two paired files for each sample.')
      sys.exit(-1)
    logging.info('First file:'+in_fqfilelist[0])
    logging.info('Second file:'+in_fqfilelist[1])
    gr_count_dict=process_pair_fastq_files(in_fqfilelist[0],in_fqfilelist[1],grdict,gmismatchdict,args)
   
  # output to file
  output_to_file(args,gr_count_dict,grdict,librecord=grna)
  
    
    
          




