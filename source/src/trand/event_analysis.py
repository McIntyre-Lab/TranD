#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
# Copyright © 2021 Oleksandr Moskalenko <om@rc.ufl.edu>
# Distributed under terms of the MIT license.
"""
TranD is a collection of tools to facilitate metrics of structural variation
for whole genome transcript annotation files (GTF) that pinpoint structural
variation to the nucleotide level.

TranD (Transcript Distances) can be used to calculate metrics of structural
variation within and between annotation files (GTF). Structural variation
reflects organismal complexity and three summary metrics for genome level
complexity are calculated for every gene in a GTF file: 1) the number of
transcripts per gene; 2) the number of exons per transcript; and 3) the number
of unique exons (exons with unique genomic coordinates) per gene. From each
these metrics distributions a summary statistics such as mean, median,
variance, inter-quartile range are calculated. With 1GTF file input, gene mode
can be used to generates these metrics for each gene and summary statistics and
distributions across genes. Distributions are visualized in a series of plots.
For 1 GTF and 2GTF a pairwise mode calculates distance metrics between 2
transcripts to the nucleotide. In 1 GTF this is all possible pairs within the
gene and in 2 GTF model this is all possible pairs among GTF files. The
distribution of these metrics across genes are visualized and summary
statistics for structural variations between pairs are calculated and reported.
Visualizations of the distributions of the frequency of intron retention,
alternative exon usage, donor/acceptor variation and 5', 3' variation in UTR
regions are provided as well as tabular formatted nucleotide level distances
for both 1GTF and 2 GTF.
"""

import itertools
import pandas as pd
from collections import namedtuple
from dataclasses import dataclass
from loguru import logger
from pybedtools import BedTool
from multiprocessing import Pool

# Import general trand functions
from .io import write_output
from .io import write_gtf
from .io import open_output_files
from .io import read_exon_data_from_file
from .bedtools import prep_bed_for_ea

# Import transcript distance functions and variables
from . import transcript_distance as TD
from . import minimum_distance as MD

# Import plotting functions
from . import plot_two_gtf_pairwise as P2GP
from . import plot_one_gtf_pairwise as P1GP
from . import plot_one_gtf_gene as P1GG

# Import complexity calculation function
from . import calculate_complexity as COMP

# Import consolidation function
from . import consolidation as CONSOL

# CONFIGURATION
jct_df_cols = ['gene_id', 'transcript_id', 'coords']
er_df_cols = ['gene_id', 'er_id', 'er_chr', 'er_start', 'er_end', 'er_strand', 'er_exon_ids',
              'er_transcript_ids', 'gene_transcript_ids', 'exons_per_er', 'transcripts_per_er',
              'transcripts_per_gene',  'er_ir_flag', 'er_annotation_frequency']
ea_df_cols = ['gene_id', 'transcript_1', 'transcript_2', 'transcript_id', 'ef_id', 'ef_chr',
              'ef_start', 'ef_end', 'ef_strand', 'ef_ir_flag', 'er_id', 'er_chr', 'er_start',
              'er_end', 'er_strand']
ef_df_cols = ['gene_id', 'er_id', 'ef_id', 'ef_chr', 'ef_start', 'ef_end', 'ef_strand',
              'ef_exon_ids', 'ef_transcript_ids', 'exons_per_ef', 'transcripts_per_ef',
              'ef_ir_flag', 'ea_annotation_frequency']
ir_df_cols = ['er_transcript_ids']
ue_df_cols = ['gene_id', 'num_uniq_exon']

# Data structures for ERs and EFs. Use a mutable dataclass, so we could add transcript names and
# counts during EA
@dataclass
class ER:
    gene_id: str
    id: str
    chrom: str
    start: int
    end: int
    strand: str
    ex_ids: list
    tx_ids: list
    gene_tx_ids: list
    ex_num: int = 0
    tx_num: int = 0
    gene_tx_num: int = 0
    ir_flag: int = 0
    er_freq: str = ''


@dataclass
class EF:
    gene_id: str
    er_id: str
    ef_id: str
    chrom: str
    start: int
    end: int
    strand: str
    ex_ids: list
    tx_ids: list
    ex_num: int = 0
    tx_num: int = 0
    ir_flag: int = 0
    ef_freq: str = ''


# TODO refactor module level variables used for first-pass parallelization
ea_list = []
jct_list = []
td_list = []


def chunks(lst, n):
    # Split list into n chunks and return list of lists
    splitSize = (len(lst)//n) + ((len(lst) % n) > 0)
    list_of_lists = []
    for element in range(0, len(lst), splitSize):
        list_of_lists.append(lst[element:element + splitSize])
    return list_of_lists


def callback_results(results):
    # Callback function to append result to list of results
    ea_data_cat, jct_data_cat, td_data_cat = results
    ea_list.append(ea_data_cat)
    jct_list.append(jct_data_cat)
    td_list.append(td_data_cat)


def write_td_data(data, out_fhs):
    """Write results of a transcript distance comparison to output files."""
    logger.debug("Writing TD data")


def do_td(data):
    """Perform Transcript Distance analysis on a pair of transcripts"""
    raise NotImplementedError
    # tx_names = data['transcript_id'].unique()
    # logger.debug("TD on {}", tx_names)
    td_data = ""
    return td_data


def format_ea_identical(gene, tx_names, data):
    """Format exons as EFs for identical transcripts"""
    # Start ef_id numbering from 1
    # ef_id = 1
    for feature in data:
        pass
    output = data
    return(output)


def create_junction_catalog(gene, tx, tx_data):
    """Create a junction catalog for a transcript"""
    junctions = []
    id = 0
    for e in tx_data:
        id += 1
        if id == 1:
            left_end = int(e.end) - 10
            continue
        right_start = int(e.start + 10)
        jct_coords = f"{e.chrom}:{left_end}:{right_start}:{e.strand}"
        junctions.append([gene, tx, jct_coords])
        left_end = int(e.end) - 10
    return junctions


def create_junction_catalog_str(gene, tx, tx_data):
    """Create a junction catalog for a transcript"""
    junctions = []
    id = 0
    tx_data_df = pd.DataFrame(tx_data, columns=['chrom', 'start', 'end', 'name',
                              'score', 'strand'])
    tx_data_df['start'] = tx_data_df['start'].astype(int)
    tx_data_df['end'] = tx_data_df['end'].astype(int)
    tx_data_df = tx_data_df.sort_values('start')
    for index, e in tx_data_df.iterrows():
        id += 1
        if id == 1:
            left_end = int(e.end) - 10
            continue
        right_start = int(e.start + 10)
        jct_coords = f"{e.chrom}:{left_end}:{right_start}:{e.strand}"
        junctions.append([gene, tx, jct_coords])
        left_end = int(e.end) - 10
    return junctions


def get_intron_retention_efs(ers_bed, efs_bed, common_efs):
    """
    Produce a list of exon fragments that participate in intron retention events.
    In essence, retained introns are ER-internal transcript-specific exonic fragments.
    So, they cannot be on the ER borders and cannot be shared between two transcripts.
    For a tiny speedup check that we have more than two EFs in an ER.
    """
    ir_efs = []
    for er in ers_bed:
        er_bed_str = f"{er.chrom}\t{er.start}\t{er.end}\t{er.name}\t{er.score}\t{er.name}"
        er_bed = BedTool(er_bed_str, from_string=True)
        er_efs = efs_bed.intersect(er_bed, wb=True)
        if len(er_efs) >= 3:
            for ef in er_efs:
                # Intron cannot be on the ER border by definition
                if ef.start != er.start and ef.end != er.end:
                    if ef.name not in common_efs:
                        ir_efs.append(ef.name)
    return ir_efs


def remove_ir_transcripts(tx_data, introns):
    """Remove transcripts that contain Intron Retention events i.e. an
    intron is encompassed by an exon
    """
    to_remove = []
    for e_tx in tx_data:
        REMOVED = False
        e_ints = [(int(x.start), int(x.end)) for x in tx_data[e_tx]]
        e_ints.sort(key=lambda x: x[0])
        for i in introns:
            if REMOVED:
                break
            for e in e_ints:
                # This is the exon this intron may overlap with
                if i[0] >= e[0] and i[0] <= e[1]:
                    # Intron is fully encompassed within the exon
                    if i[1] <= e[1]:
                        logger.debug("IR detected, removing {}", e_tx)
                        # logger.debug("Exon: {}:{}, intron: {}:{}", e[0], e[1], i[0], i[1])
                        to_remove.append(e_tx)
                        REMOVED = True
                        break
                    else:
                        # Stop processing this transcript's exons, move on to next intron
                        continue
    for ir_tx in to_remove:
        del tx_data[ir_tx]
    return tx_data


def list_ir_exon_transcript(tx_data, introns):
    """List exons that contain Intron Retention events i.e. an
    intron is encompassed by an exon
    """
    ir_exon_list = []
    ir_transcript_list = []
    for e_tx in tx_data:
        e_ints = [(int(x.start), int(x.end), x.name) for x in tx_data[e_tx]]
        e_ints.sort(key=lambda x: x[0])
        for i in introns:
            for e in e_ints:
                # Check if exon has already been added to list
                if e[2] in ir_exon_list:
                    continue
                else:
                    # This is the exon this intron may overlap with
                    if i[0] >= e[0] and i[0] <= e[1]:
                        # Intron is fully encompassed within the exon
                        if i[1] <= e[1]:
                            logger.debug("IR detected, {}", e)
                            ir_exon_list.append(e[2])
                            if e_tx not in ir_transcript_list:
                                ir_transcript_list.append(e_tx)
    return ir_exon_list, ir_transcript_list


def do_ea_gene(data, keep_ir):
    """
    Event Analysis on a full gene
    """
    gene_id = data['gene_id']
    logger.debug("Gene EA on {}", gene_id)
    # logger.debug(data)
    tx_names = data['transcript_list']
    logger.debug("Transcripts in {}: {}", gene_id, tx_names)
    if len(tx_names) == 1:
        logger.info("Gene {} has a single transcript.", gene_id)
        er_data, ef_data = single_transcript_ea(gene_id, tx_names[0], data[tx_names[0]])
        junction_data = create_junction_catalog_str(gene_id, tx_names[0], data[tx_names[0]])
        junction_df = pd.DataFrame(junction_data, columns=jct_df_cols)
        ir_transcripts = []
    else:
        tx_data = {}
        tx_coords = {}
        all_introns = []
        for tx_name in tx_names:
            tx_bed_str = data[tx_name]
            starts = [int(x[1]) for x in tx_bed_str]
            ends = [int(x[2]) for x in tx_bed_str]
            chrom = data[tx_name][0][0]
            strand = data[tx_name][0][5]
            feature = f"{tx_name}_transcript"
            left_edge = min(starts)
            right_edge = max(ends)
            tx_exons = BedTool(tx_bed_str).saveas()
            tx_data[tx_name] = tx_exons
            tx_coords[tx_name] = [left_edge, right_edge]
            tx_full = BedTool(f"{chrom} {left_edge} {right_edge} {feature} 0 {strand}",
                              from_string=True)
            tx_introns = tx_full.subtract(tx_exons)
            i_coords = [(int(x.start), int(x.end)) for x in tx_introns]
            all_introns.extend(i_coords)
        intron_data = list(set(all_introns))
        intron_data.sort(key=lambda x: x[0])
        # Removal of IR containing transcripts
        if not keep_ir:
            tx_data = remove_ir_transcripts(tx_data, intron_data)
            old_tx_names = tx_names
            tx_names = tx_data.keys()
            tx_coords_delete = [k for k in old_tx_names if k not in tx_names]
            for k in tx_coords_delete:
                del tx_coords[k]
            ir_exons = []
            ir_transcripts = []
        else:
            ir_exons, ir_transcripts = list_ir_exon_transcript(tx_data, intron_data)
            logger.info("Not attempting to remove IR-containing transcripts")
        # Junction Catalog creation
        junction_data = []
        for tx_name in tx_names:
            tx_jc_data = create_junction_catalog(gene_id, tx_name, tx_data[tx_name])
            junction_data.extend(tx_jc_data)
            for jc in tx_jc_data:
                tx_coords[tx_name].append(jc[2])
        junction_df = pd.DataFrame(junction_data, columns=jct_df_cols)
        # Event Analysis
        if not tx_data:
            logger.warning("Missing transcript data for {} gene, skipping", gene_id)
            return None, None, None
        er_data, ef_data = ea_analysis(gene_id, tx_data, tx_coords, ir_exons)
    er_df = pd.DataFrame(er_data, columns=er_df_cols)
    ef_df = pd.DataFrame(ef_data, columns=ef_df_cols)
    return er_df, ef_df, junction_df, ir_transcripts


def single_transcript_ea(gene_id, tx_name, tx_data):
    """
    Generate ERs (Exonic Regions) and EFs (Exonic Fragments) for single transcript genes.

    Required outputs:
    gene_id
    gene_transcript_ids
    transcripts_per_gene
    er_id
    er_chr
    er_start
    er_end
    er_strand
    er_flag_ir ( zero if IRs were removed)
    exons_per_er
    er_exon_ids = Piped ("|") list of exon IDs that are contained within the exon region
    transcripts_per_er
    er_transcript_ids = Piped ("|") list of transcript IDs the exon region is present in
    er_annotation_frequency = “unique”, “common”, “constitutive” (one, many, all) per txs
    ef_id
    ef_chr
    ef_start
    ef_end
    ef_strand
    ef_flag_ir (zero if IRs were removed)
    ef_exon_ids
    exons_per_ef
    transcripts_per_ef
    ef_transcript_ids = Piped ("|") list of transcript IDs the exon fragment is present in
    ef_annotation_frequency = “unique”, “common”, “constitutive” (one, many, all) per txs
    """
    er_out_data = []
    ef_out_data = []
    er_id = 1
    ef_id = 1
    for i in tx_data:
        er_name = f"{gene_id}:ER{er_id}"
        ef_name = f"{gene_id}:ER{er_id}:EF{ef_id}"
        chrom = i[0]
        start = i[1]
        end = i[2]
        strand = i[5]
        ex_id = i[3]
        ex_num = 1
        tx_num = 1
        gene_tx_num = 1
        ir_flag = 0
        er_freq = "unique"
        ef_freq = "unique"
        er_out_data.append([gene_id, er_name, chrom, start, end, strand, ex_id, tx_name, tx_name,
                            ex_num, tx_num, gene_tx_num, ir_flag, er_freq])
        ef_out_data.append([gene_id, er_name, ef_name, chrom, start, end, strand, ex_id, tx_name,
                            ex_num, tx_num, ir_flag, ef_freq])
        er_id += 1
        ef_id += 1
    return er_out_data, ef_out_data


def do_ea_pair(data):
    """
    Event Analysis on a pair of transcripts
    """
    gene_id = data['gene_id']
    tx_names = data['transcript_list']
    tx1_name, tx2_name = tx_names[0], tx_names[1]
    tx1_bed_str = data[tx1_name]
    tx2_bed_str = data[tx2_name]
    tx1_bed_df = pd.DataFrame(tx1_bed_str, columns=['chrom', 'start', 'end', 'name', 'score',
                                                    'strand'])
    tx1_bed_df['start'] = tx1_bed_df['start'].astype(int)
    tx1_bed_df['end'] = tx1_bed_df['end'].astype(int)
    tx1_bed_df = tx1_bed_df.sort_values('start')
    tx2_bed_df = pd.DataFrame(tx2_bed_str, columns=['chrom', 'start', 'end', 'name', 'score',
                                                    'strand'])
    tx2_bed_df['start'] = tx2_bed_df['start'].astype(int)
    tx2_bed_df['end'] = tx2_bed_df['end'].astype(int)
    tx2_bed_df = tx2_bed_df.sort_values('start')
    logger.debug("Comparing {} vs {}", tx1_name, tx2_name)
    junction_data1 = create_junction_catalog_str(gene_id, tx1_name, tx1_bed_str)
    junction_data2 = create_junction_catalog_str(gene_id, tx2_name, tx2_bed_str)
    junction_str1 = "|".join(pd.DataFrame(junction_data1, columns=['gene_id', 'transcript_id',
                                                                   'junction_id'])['junction_id'])
    junction_str2 = "|".join(pd.DataFrame(junction_data2, columns=['gene_id', 'transcript_id',
                                                                   'junction_id'])['junction_id'])
    junction_data = junction_data1 + junction_data2
    # Check if transcripts are full-splice match (all junctions are the same)
    same_start = False
    same_end = False
    is_fsm = False
    if junction_str1 == junction_str2:
        tx1_min = tx1_bed_df['start'].min()
        tx2_min = tx2_bed_df['start'].min()
        tx1_max = tx1_bed_df['end'].max()
        tx2_max = tx2_bed_df['end'].max()
        # Check for non-overlapping mono-exon transcript pairs
        if tx1_min < tx2_min and tx1_max <= tx2_min:
            # T1 completely upstream of T2
            tx1_bed = BedTool(tx1_bed_str).saveas()
            tx2_bed = BedTool(tx2_bed_str).saveas()
            ea_data = er_ea_analysis(tx1_bed, tx2_bed, tx1_name, tx2_name, gene_id)
        elif tx2_min < tx1_min and tx2_max <= tx1_min:
            # T2 completely upstream of T1
            tx1_bed = BedTool(tx1_bed_str).saveas()
            tx2_bed = BedTool(tx2_bed_str).saveas()
            ea_data = er_ea_analysis(tx1_bed, tx2_bed, tx1_name, tx2_name, gene_id)
        else:
            # Some overlap present between transcripts
            is_fsm = True
            # Check if start coordinates of the transcripts the same
            if tx1_bed_df['start'].min() == tx2_bed_df['start'].min():
                same_start = True
            # Check if end coordinates of transcripts are also the same
            if tx1_bed_df['end'].max() == tx2_bed_df['end'].max():
                same_end = True
            tx_names_str = " / ".join(tx_names)
            if same_start and same_end:
                # Transcripts are identical
                logger.info("Exons are identical for {}.", tx_names_str)
                ea_data = format_fsm_pair_ea(tx1_bed_df, tx2_bed_df, tx1_name, tx2_name, gene_id,
                                             side_diff="none")
            elif same_start and not same_end:
                # Start coordinates and junctions are the same
                logger.info("Junctions and start coordinates are identical for {}.", tx_names_str)
                ea_data = format_fsm_pair_ea(tx1_bed_df, tx2_bed_df, tx1_name, tx2_name, gene_id,
                                             side_diff="end")
            elif same_end and not same_start:
                # End coordinates and junctions are the same
                logger.info("Junctions and end coordinates are identical for {}.", tx_names_str)
                ea_data = format_fsm_pair_ea(tx1_bed_df, tx2_bed_df, tx1_name, tx2_name, gene_id,
                                             side_diff="start")
            else:
                # Junctions are the same and ends are both different
                logger.info("Junctions are identical for {}.", " / ".join(tx_names))
                ea_data = format_fsm_pair_ea(tx1_bed_df, tx2_bed_df, tx1_name, tx2_name, gene_id,
                                             side_diff="both")
    else:
        # Junctions are not identical (there is some alternate donor/acceptor/exon)
        tx1_bed = BedTool(tx1_bed_str).saveas()
        tx2_bed = BedTool(tx2_bed_str).saveas()
        ea_data = er_ea_analysis(tx1_bed, tx2_bed, tx1_name, tx2_name, gene_id)
    out_df = pd.DataFrame(ea_data, columns=ea_df_cols)
    junction_df = pd.DataFrame(junction_data, columns=jct_df_cols)
    td_df = pd.DataFrame(TD.calculate_distance(out_df, junction_df, gene_id, tx1_name, tx2_name,
                                               fsm=is_fsm)).T
    return out_df, junction_df, td_df


def er_ea_analysis(tx1_bed, tx2_bed, tx1_name, tx2_name, gene_id):
    """
    Generate ERs (Exonic Regions) and EFs (Exonic Fragments) and analyze events in a pair of
    transcripts.
    """
    ea_data = []
    strand = list(set([x.strand for x in tx1_bed]))[0]
    # logger.debug("TX1: {}: \n{}", tx1_name, tx1_bed)
    # logger.debug("TX2: {}: \n{}", tx2_name, tx2_bed)
    raw_ers_bed = tx1_bed.cat(tx2_bed, postmerge=True)
    raw_ers_list = [str(x).split() for x in raw_ers_bed]
    er_id = 1
    ers_list = []
    for i in raw_ers_list:
        i.extend([f"{gene_id}:ER{er_id}", '0', strand])
        ers_list.append("\t".join(i))
        er_id += 1
    ers_str = "\n".join(ers_list)
    # logger.debug("ERs for {} and {}: \n{}", tx1_name, tx2_name, ers_str)
    # logger.debug("ERs str:\n{}", ers_str)
    ers_bed = BedTool(ers_str, from_string=True)
    er_data = {}
    er_datum = namedtuple('ER', ['chrom', 'start', 'end', 'strand'])
    for er in ers_bed:
        er_data[er.name] = er_datum(er.chrom, er.start, er.end, er.strand)
    ef_tx1 = tx1_bed.subtract(tx2_bed)
    ef_tx2 = tx2_bed.subtract(tx1_bed)
    ef_both = tx1_bed.intersect(tx2_bed)
    efs_raw = ef_tx1.cat(ef_tx2, postmerge=False).cat(ef_both, postmerge=False).sort()
    # logger.debug("Raw EFs: \n{}", efs_raw)
    efs_er = ers_bed.intersect(efs_raw).sort()
    # logger.debug("EFs in ERs: \n{}", efs_er)
    efs_list = []
    ef_id = 1
    er_name = ''
    for i in efs_er:
        if i.name != er_name:
            ef_id = 1
        er_name = i.name
        ef_name = f"{er_name}:EF{ef_id}"
        ef_str = f"{i.chrom}\t{i.start}\t{i.end}\t{ef_name}\t {i.score}\t{i.strand}"
        efs_list.append(ef_str)
        ef_id += 1
    efs_str = "\n".join(efs_list)
    # logger.debug("Final EFs: \n{}", efs_str)
    efs_bed = BedTool(efs_str, from_string=True)
    common_efs_set = set([x.name for x in efs_bed.intersect(ef_both)])
    common_efs = list(common_efs_set)
    tx1_efs = list(set([x.name for x in efs_bed.intersect(tx1_bed)]).difference(common_efs_set))
    tx2_efs = list(set([x.name for x in efs_bed.intersect(tx2_bed)]).difference(common_efs_set))
    ir_efs = get_intron_retention_efs(ers_bed, efs_bed, common_efs)

    for ef in efs_bed:
        ir_flag = '0'
        ef_name = ef.name
        er_name = ":".join(ef_name.split(':')[:-1])
        if ef_name in ir_efs:
            ir_flag = '1'
        if ef_name in common_efs:
            tx_list = f"{tx1_name}|{tx2_name}"
        elif ef_name in tx1_efs:
            tx_list = f"{tx1_name}"
        elif ef_name in tx2_efs:
            tx_list = f"{tx2_name}"
        else:
            logger.error(f"Cannot find the fragment {ef_name} in any fragment lists, skipping")
        ea_datum = [gene_id, tx1_name, tx2_name, tx_list, ef_name, ef.chrom, ef.start, ef.end,
                    ef.strand, ir_flag, er_name, er_data[er_name].chrom, er_data[er_name].start,
                    er_data[er_name].end, er_data[er_name].strand]
        ea_data.append(ea_datum)
    logger.debug("Final EA data: \n{}", ea_data)
    return ea_data


def format_identical_pair_ea(tx1_bed, tx2_bed, tx1_name, tx2_name, gene_id):
    ea_data = []
    er_id, ef_id = 1, 1
    for i in tx1_bed:
        er_name = f"{gene_id}:ER{er_id}"
        ef_name = f"{gene_id}:ER{er_id}:EF{ef_id}"
        ea_data.append([gene_id, tx1_name, tx2_name, f"{tx1_name}|{tx2_name}", ef_name,
                        i.chrom, i.start, i.end, i.strand, 0, er_name, i.chrom, i.start,
                        i.end, i.strand])
        er_id += 1
    ea_df = pd.DataFrame(ea_data, columns=ea_df_cols)
    return ea_df


# !!! Below function potentially can replace format_identical_pair_ea (side = "none", or anything
# not "start", "end", or "both")
def format_fsm_pair_ea(tx1_bed_df, tx2_bed_df, tx1_name, tx2_name, gene_id, side_diff):
    ea_data = []
    er_id = 1
    max_er_id = len(tx1_bed_df)
    for index, i in tx1_bed_df.iterrows():
        ef_id = 1
        # Check first ER
        if er_id == 1:
            # Monoexon - Make EF and ER for monoexon transcript pairs (must both be monoexon since
            # junctions are identical)
            if er_id == max_er_id:
                # Check if which ends are different
                if side_diff == "start":
                    tx1_start = tx1_bed_df['start'].min()
                    tx2_start = tx2_bed_df['start'].min()
                    min_start = min(tx1_start, tx2_start)
                    max_start = max(tx1_start, tx2_start)
                    er_name = f"{gene_id}:ER{er_id}"
                    ef_name = f"{gene_id}:ER{er_id}:EF{ef_id}"
                    if min_start == tx1_start:
                        tx_name = tx1_name
                    else:
                        tx_name = tx2_name
                    ea_data.append([gene_id, tx1_name, tx2_name, f"{tx_name}", ef_name,
                                    i.chrom, min_start, max_start, i.strand, 0, er_name, i.chrom,
                                    min_start, i.end, i.strand])
                    ef_id += 1
                    ef_name = f"{gene_id}:ER{er_id}:EF{ef_id}"
                    ea_data.append([gene_id, tx1_name, tx2_name, f"{tx1_name}|{tx2_name}", ef_name,
                                    i.chrom, max_start, i.end, i.strand, 0, er_name, i.chrom,
                                    min_start, i.end, i.strand])
                    ef_id += 1
                elif side_diff == "end":
                    tx1_end = tx1_bed_df['end'].max()
                    tx2_end = tx2_bed_df['end'].max()
                    min_end = min(tx1_end, tx2_end)
                    max_end = max(tx1_end, tx2_end)
                    er_name = f"{gene_id}:ER{er_id}"
                    ef_name = f"{gene_id}:ER{er_id}:EF{ef_id}"
                    ea_data.append([gene_id, tx1_name, tx2_name, f"{tx1_name}|{tx2_name}", ef_name,
                                    i.chrom, i.start, min_end, i.strand, 0, er_name, i.chrom,
                                    i.start, max_end, i.strand])
                    ef_id += 1
                    ef_name = f"{gene_id}:ER{er_id}:EF{ef_id}"
                    if max_end == tx1_end:
                        tx_name = tx1_name
                    else:
                        tx_name = tx2_name
                    ea_data.append([gene_id, tx1_name, tx2_name, f"{tx_name}", ef_name,
                                    i.chrom, min_end, max_end, i.strand, 0, er_name, i.chrom,
                                    i.start, max_end, i.strand])
                elif side_diff == "both":
                    tx1_start = tx1_bed_df['start'].min()
                    tx2_start = tx2_bed_df['start'].min()
                    min_start = min(tx1_start, tx2_start)
                    max_start = max(tx1_start, tx2_start)
                    tx1_end = tx1_bed_df['end'].max()
                    tx2_end = tx2_bed_df['end'].max()
                    min_end = min(tx1_end, tx2_end)
                    max_end = max(tx1_end, tx2_end)
                    er_name = f"{gene_id}:ER{er_id}"
                    ef_name = f"{gene_id}:ER{er_id}:EF{ef_id}"
                    if min_start == tx1_start:
                        tx_name = tx1_name
                    else:
                        tx_name = tx2_name
                    ea_data.append([gene_id, tx1_name, tx2_name, f"{tx_name}", ef_name, i.chrom,
                                    min_start, max_start, i.strand, 0, er_name, i.chrom, min_start,
                                    max_end, i.strand])
                    ef_id += 1
                    ef_name = f"{gene_id}:ER{er_id}:EF{ef_id}"
                    ea_data.append([gene_id, tx1_name, tx2_name, f"{tx1_name}|{tx2_name}", ef_name,
                                    i.chrom, max_start, min_end, i.strand, 0, er_name, i.chrom,
                                    min_start, max_end, i.strand])
                    ef_id += 1
                    if max_end == tx1_end:
                        tx_name = tx1_name
                    else:
                        tx_name = tx2_name
                    ef_name = f"{gene_id}:ER{er_id}:EF{ef_id}"
                    ea_data.append([gene_id, tx1_name, tx2_name, f"{tx_name}", ef_name, i.chrom,
                                    min_end, max_end, i.strand, 0, er_name, i.chrom, min_start,
                                    max_end, i.strand])
                else:
                    # No difference on ends so monoexon are identical
                    er_name = f"{gene_id}:ER{er_id}"
                    ef_name = f"{gene_id}:ER{er_id}:EF{ef_id}"
                    ea_data.append([gene_id, tx1_name, tx2_name, f"{tx1_name}|{tx2_name}", ef_name,
                                    i.chrom, i.start, i.end, i.strand, 0, er_name, i.chrom, i.start,
                                    i.end, i.strand])
            # Multiexon - Make EF for first ER in multiexon tx pairs that have different starts
            elif (er_id != max_er_id) and (side_diff == "start" or side_diff == "both"):
                tx1_start = tx1_bed_df['start'].min()
                tx2_start = tx2_bed_df['start'].min()
                min_start = min(tx1_start, tx2_start)
                max_start = max(tx1_start, tx2_start)
                er_name = f"{gene_id}:ER{er_id}"
                ef_name = f"{gene_id}:ER{er_id}:EF{ef_id}"
                if min_start == tx1_start:
                    tx_name = tx1_name
                else:
                    tx_name = tx2_name
                ea_data.append([gene_id, tx1_name, tx2_name, f"{tx_name}", ef_name, i.chrom,
                                min_start, max_start, i.strand, 0, er_name, i.chrom, min_start,
                                i.end, i.strand])
                ef_id += 1
                ef_name = f"{gene_id}:ER{er_id}:EF{ef_id}"
                ea_data.append([gene_id, tx1_name, tx2_name, f"{tx1_name}|{tx2_name}", ef_name,
                                i.chrom, max_start, i.end, i.strand, 0, er_name, i.chrom, min_start,
                                i.end, i.strand])
                ef_id += 1
            # Make single EF in ER in multiexon transcript pairs that have no difference in starts
            else:
                er_name = f"{gene_id}:ER{er_id}"
                ef_name = f"{gene_id}:ER{er_id}:EF{ef_id}"
                ea_data.append([gene_id, tx1_name, tx2_name, f"{tx1_name}|{tx2_name}", ef_name,
                                i.chrom, i.start, i.end, i.strand, 0, er_name, i.chrom, i.start,
                                i.end, i.strand])
        # Make one EF for middle ER (since all junctions are shared, all internal ER are identical)
        if er_id != 1 and er_id != max_er_id:
            er_name = f"{gene_id}:ER{er_id}"
            ef_name = f"{gene_id}:ER{er_id}:EF{ef_id}"
            ea_data.append([gene_id, tx1_name, tx2_name, f"{tx1_name}|{tx2_name}", ef_name,
                            i.chrom, i.start, i.end, i.strand, 0, er_name, i.chrom, i.start,
                            i.end, i.strand])
        # Check last ER (for multiexon only)
        if er_id == max_er_id and er_id != 1:
            # Make EF for last ER in multiexon transcript pairs that have different ends
            if (side_diff == "end" or side_diff == "both"):
                tx1_end = tx1_bed_df['end'].max()
                tx2_end = tx2_bed_df['end'].max()
                min_end = min(tx1_end, tx2_end)
                max_end = max(tx1_end, tx2_end)
                er_name = f"{gene_id}:ER{er_id}"
                ef_name = f"{gene_id}:ER{er_id}:EF{ef_id}"
                ea_data.append([gene_id, tx1_name, tx2_name, f"{tx1_name}|{tx2_name}", ef_name,
                                i.chrom, i.start, min_end, i.strand, 0, er_name, i.chrom, i.start,
                                max_end, i.strand])
                ef_id += 1
                ef_name = f"{gene_id}:ER{er_id}:EF{ef_id}"
                if max_end == tx1_end:
                    tx_name = tx1_name
                else:
                    tx_name = tx2_name
                ea_data.append([gene_id, tx1_name, tx2_name, f"{tx_name}", ef_name,
                                i.chrom, min_end, max_end, i.strand, 0, er_name, i.chrom, i.start,
                                max_end, i.strand])
                ef_id += 1
            # Make single EF in ER in multiexon transcript pairs that have no difference in ends
            else:
                er_name = f"{gene_id}:ER{er_id}"
                ef_name = f"{gene_id}:ER{er_id}:EF{ef_id}"
                ea_data.append([gene_id, tx1_name, tx2_name, f"{tx1_name}|{tx2_name}", ef_name,
                                i.chrom, i.start, i.end, i.strand, 0, er_name, i.chrom, i.start,
                                i.end, i.strand])

        er_id += 1
    ea_df = pd.DataFrame(ea_data, columns=ea_df_cols)
    return ea_df


def ea_analysis(gene_id, tx_data, tx_coords, ir_exons):
    """
    Generate ERs (Exonic Regions) and EFs (Exonic Fragments) and analyze events in transcripts.
    Generalize to do both full-gene and transcript pairs to replace er_ea_analysis.

    Required outputs:
    gene_id
    gene_transcript_ids
    transcripts_per_gene
    er_id
    er_chr
    er_start
    er_end
    er_strand
    er_flag_ir ( zero if IRs were removed)
    exons_per_er
    er_exon_ids = Piped ("|") list of exon IDs that are contained within the exon region
    transcripts_per_er
    er_transcript_ids = Piped ("|") list of transcript IDs the exon region is present in
    er_annotation_frequency = “unique”, “common”, “constitutive” (one, many, all) per txs
    ef_id
    ef_chr
    ef_start
    ef_end
    ef_strand
    ef_flag_ir (zero if IRs were removed)
    ef_exon_ids
    exons_per_ef
    transcripts_per_ef
    ef_transcript_ids = Piped ("|") list of transcript IDs the exon fragment is present in
    ef_annotation_frequency = “unique”, “common”, “constitutive” (one, many, all) per txs
    """
    # Use the junction catalog and transcripts start/end to check for identical transcripts
    logger.debug("Performing EA analysis for {} gene", gene_id)
    logger.debug("TX Names: {}", list(tx_data.keys()))
    # logger.debug("TX Coords: {}", tx_coords)
    gene_transcript_ids = list(tx_data.keys())
    total_number_of_transcripts = len(gene_transcript_ids)
    tx_ident_check = {}
    for i in tx_coords:
        tx_structure_str_list = [str(x) for x in tx_coords[i]]
        tx_structure = ",".join(tx_structure_str_list)
        if tx_structure in tx_ident_check:
            tx_ident_check[tx_structure].append(i)
        else:
            tx_ident_check[tx_structure] = [i]
    identical_tx_list = [x for x in tx_ident_check.values() if len(x) > 1]
    identical_transcripts = {x[0]: x[1:] for x in identical_tx_list}
    if identical_transcripts:
        logger.debug("Identical records: {}", identical_transcripts)
    # Stash names for EA output, but remove duplicate data
    duplicate_txs = [x[1:] for x in identical_tx_list]
    if duplicate_txs:
        logger.debug("Duplicate records: {}", duplicate_txs)
    for dupe_list in duplicate_txs:
        for item in dupe_list:
            del tx_data[item]
    logger.debug("Unique TX Data: {}", list(tx_data.keys()))
    # Combine transcripts and generate ERs
    strands = []
    all_tx = None
    for i in tx_data:
        tx = tx_data[i]
        if not all_tx:
            all_tx = tx
        else:
            all_tx = all_tx.cat(tx, postmerge=True)
        strand = list(set([x.strand for x in tx]))[0]
        strands.append(strand)
    strand = list(set(strands))[0]
    # Is this a valid check?
    if len(strand) > 1:
        raise ValueError("Multiple strands detected when there should be one")
    logger.debug("Strand: '{}'", strand)
    raw_ers_list = [str(x).split() for x in all_tx]
    er_id = 1
    ers_list = []
    for i in raw_ers_list:
        i.extend([f"{gene_id}:ER{er_id}", '0', strand])
        ers_list.append("\t".join(i))
        er_id += 1
    ers_str = "\n".join(ers_list)
    ers_bed = BedTool(ers_str, from_string=True)
    er_data = {}
    # Coordinates of ER ends to check when creating EFs
    er_ends = []
    er_range_sets = {}
    for er in ers_bed:
        er_data[er.name] = ER(gene_id=gene_id, id=er.name, chrom=er.chrom, start=er.start,
                              end=er.end, strand=er.strand, ex_ids=[], tx_ids=[],
                              gene_tx_ids=gene_transcript_ids)
        er_ends.append(er.end)
        er_range_sets[er.name] = set(range(er.start, er.end + 1))
    # Generate EFs
    # Catalog starting and ending points (0-coordinate as we go from left side only)
    sts_ends = []
    tx_ranges = {}
    for i in tx_data:
        tx = tx_data[i]
        for ex in tx:
            # logger.debug(f"{ex.name}: {ex.start}/{ex.end}")
            # There are many ways to check if exon start is in ER, pick one at random
            if i in tx_ranges:
                tx_ranges[ex.name].append(ex.start, ex.end, i)
            else:
                tx_ranges[ex.name] = (ex.start, ex.end, i)
            if ex.start not in sts_ends:
                sts_ends.append(ex.start)
            if (ex.end) not in sts_ends:
                sts_ends.append(ex.end)
    sts_ends.sort()
    tx_range_sets = {}
    tx_exon_map = {}
    for i in tx_ranges:
        tx = tx_ranges[i]
        tx_range_set = set(range(tx[0], tx[1] + 1))
        tx_range_sets[i] = tx_range_set
        tx_exon_map[i] = tx[2]
    # Define EFs
    ef_list = []
    ef_start = sts_ends[0]
    ef_end = None
    for i in sts_ends[1:]:
        # Start of an ER
        if not ef_start and not ef_end:
            ef_start = i
            continue
        ef_end = i
        ef_list.append((ef_start, ef_end))
        if ef_end in er_ends:
            # Skip the intron
            ef_start = None
            ef_end = None
            continue
        else:
            ef_start = ef_end
    # Process all ERs and EFs vs transcript exons/transcripts to generate requested EA output
    # Sets are in 'er_range_sets' and 'tx_range_sets'
    for er_id in er_data:
        for ex_id in tx_exon_map:
            exon_data = tx_ranges[ex_id]
            # Check if exon start is in the ER
            if exon_data[0] in er_range_sets[er_id]:
                # Exon is in the ER, record exon, tx, and duplicates if any
                tx_id = exon_data[2]
                if ex_id not in er_data[er_id].ex_ids:
                    er_data[er_id].ex_ids.append(ex_id)
                if tx_id not in er_data[er_id].tx_ids:
                    er_data[er_id].tx_ids.append(tx_id)
                # Add duplicate txs and their exons if exist
                if tx_id in identical_transcripts:
                    duplicate_tx_list = identical_transcripts[tx_id]
                    for duplicate_tx in duplicate_tx_list:
                        if duplicate_tx not in er_data[er_id].tx_ids:
                            er_data[er_id].tx_ids.append(duplicate_tx)
                        ex_number = ex_id.split('_exon_')[1]
                        duplicate_ex = f"{duplicate_tx}_exon_{ex_number}"
                        if duplicate_ex not in er_data[er_id].ex_ids:
                            er_data[er_id].ex_ids.append(ex_id)
    for er in er_data:
        ex_num = len(er_data[er].ex_ids)
        tx_num = len(er_data[er].tx_ids)
        er_data[er].ex_num = ex_num
        er_data[er].tx_num = tx_num
        er_data[er].gene_tx_num = total_number_of_transcripts
        if tx_num == 1:
            er_data[er].er_freq = 'unique'
        elif tx_num == total_number_of_transcripts:
            er_data[er].er_freq = 'constitutive'
        else:
            er_data[er].er_freq = 'common'

    # Process EFs vs Exons and Transcripts
    ef_data = {}
    ef_start_set = set([x[0] for x in ef_list])
    ef_start_end_map = {x[0]: x[1] for x in ef_list}
    er_ef_map = {}
    # Match agains ERs to get id and common data
    for er in er_range_sets:
        er_set = er_range_sets[er]
        er_ef_starts = list(er_set.intersection(ef_start_set))
        for ef_start in er_ef_starts:
            er_ef_map[er] = er_ef_starts
    for er_id in er_ef_map:
        ef_id_ord = 1
        for ef_start in er_ef_map[er_id]:
            ef_id = f"{er_id}:EF{ef_id_ord}"
            ef_id_ord += 1
            ef_data[ef_id] = EF(gene_id=gene_id, er_id=er_id, ef_id=ef_id,
                                chrom=er_data[er_id].chrom, start=ef_start,
                                end=ef_start_end_map[ef_start], strand=er_data[er_id].strand,
                                ex_ids=[], tx_ids=[])
    # Match against exons
    for ex_id in tx_ranges:
        ex_set = set(range(tx_ranges[ex_id][0], tx_ranges[ex_id][1] + 1))
        for ef_id in ef_data:
            ef_start = ef_data[ef_id].start
            if ef_start in ex_set:
                # EF is a part of an exon - are there any duplicate exons and transcripts?
                tx_id = tx_exon_map[ex_id]
                if ex_id not in ef_data[ef_id].ex_ids:
                    ef_data[ef_id].ex_ids.append(ex_id)
                if tx_id not in ef_data[ef_id].tx_ids:
                    ef_data[ef_id].tx_ids.append(tx_id)
                # Add duplicate txs and their exons if exist
                if tx_id in identical_transcripts:
                    duplicate_tx_ids = identical_transcripts[tx_id]
                    ex_number = ex_id.split('_exon_')[1]
                    for dupe_tx_id in duplicate_tx_ids:
                        if dupe_tx_id not in ef_data[ef_id].tx_ids:
                            ef_data[ef_id].tx_ids.append(dupe_tx_id)
                        dupe_ex_id = f"{dupe_tx_id}_exon_{ex_number}"
                        if dupe_ex_id not in ef_data[ef_id].ex_ids:
                            ef_data[ef_id].ex_ids.append(dupe_ex_id)
    for ef in ef_data:
        ex_num = len(ef_data[ef].ex_ids)
        tx_num = len(ef_data[ef].tx_ids)
        ef_data[ef].ex_num = ex_num
        ef_data[ef].tx_num = tx_num
        if tx_num == 1:
            ef_data[ef].ef_freq = 'unique'
        elif tx_num == total_number_of_transcripts:
            ef_data[ef].ef_freq = 'constitutive'
        else:
            ef_data[ef].ef_freq = 'common'
    er_out_data = []
    for er in er_data:
        i = er_data[er]
        # Flag exon regions that contain IR exons
        if len(ir_exons) != 0 and len(set(i.ex_ids).intersection(set(ir_exons))) > 0:
            i.ir_flag = 1
        else:
            i.ir_flag = 0
        er_out_data.append([i.gene_id, i.id, i.chrom, i.start, i.end, i.strand, "|".join(i.ex_ids),
                            "|".join(i.tx_ids), "|".join(i.gene_tx_ids), i.ex_num, i.tx_num,
                            i.gene_tx_num, i.ir_flag, i.er_freq])
    ef_out_data = []
    for ef in ef_data:
        i = ef_data[ef]
        # Flag exon fragments that contain IR exons
        if len(ir_exons) != 0 and len(set(i.ex_ids).intersection(set(ir_exons))) > 0:
            i.ir_flag = 1
        else:
            i.ir_flag = 0
        ef_out_data.append([i.gene_id, i.er_id, i.ef_id, i.chrom, i.start, i.end, i.strand,
                            "|".join(i.ex_ids), "|".join(i.tx_ids), i.ex_num, i.tx_num, i.ir_flag,
                            i.ef_freq])
    return er_out_data, ef_out_data


def do_ea(tx_data, mode='pairwise', keep_ir=False):
    """Perform event analysis using pybedtools on bed string data"""
    try:
        bed_data = prep_bed_for_ea(tx_data)
    except ValueError:
        raise
    if mode == 'pairwise':
        ea_results, jct_catalog, transcript_distance = do_ea_pair(bed_data)
        return ea_results, jct_catalog, transcript_distance
    # full-gene, more transcripts than two
    elif mode == 'gene':
        # full-gene, all transcripts
        er_df, ef_df, jct_catalog, ir_transcripts = do_ea_gene(bed_data, keep_ir)
        return er_df, ef_df, jct_catalog, ir_transcripts
    else:
        raise ValueError("Wrong EA mode")


def ea_pairwise(data):
    "Do EA (Event Analysis) for a pair transcripts."
    transcripts = data.groupby("transcript_id")
    transcript_groups = transcripts.groups
    tx_data = {}
    for transcript in transcript_groups:
        transcript_df = data[data['transcript_id'] == transcript]
        tx_data[transcript] = transcript_df
    transcript_pairs = list(itertools.combinations(tx_data.keys(), 2))
    ea_df = pd.DataFrame(columns=ea_df_cols)
    jct_df = pd.DataFrame(columns=jct_df_cols)
    td_df = pd.DataFrame(columns=TD.td_df_cols)
    for tx_pair in transcript_pairs:
        # try:
        tx_df_1 = tx_data[tx_pair[0]]
        tx_df_2 = tx_data[tx_pair[1]]
        tx_pair_data = pd.concat([tx_df_1, tx_df_2])
        # try:
        ea_data, jct_data, td_data = do_ea(tx_pair_data, mode='pairwise')
        ea_df = ea_df.append(ea_data)
        jct_df = jct_df.append(jct_data)
        td_df = td_df.append(td_data)
        # except ValueError as e:
        #     logger.error(e)
        #     exit(1)
        #     #     # DEBUG, exit for now, just write to rejects later
        #     #     # continue
    return ea_df, jct_df, td_df


def process_single_file(infile, ea_mode, keep_ir, outdir, outfiles, complexity_only, skip_plots,
                        skip_interm, consolidate, consol_prefix, consol_outfiles):
    """Compare all transcript pairs in a single GTF file."""
    logger.info("Input file: {}", infile)
    if skip_interm:
        del(outfiles['ea_fh'])
        del(outfiles['er_fh'])
        del(outfiles['ef_fh'])
        del(outfiles['jc_fh'])
        del(outfiles['ir_fh'])
        del(outfiles['ue_fh'])
    else:
        if ea_mode == 'gene':
            del(outfiles['ea_fh'])
        else:
            del(outfiles['er_fh'])
            del(outfiles['ef_fh'])
            del(outfiles['ir_fh'])
            del(outfiles['ue_fh'])

    data = read_exon_data_from_file(infile)
    genes = data.groupby("gene_id")
    transcripts = data.groupby(["gene_id", "transcript_id"])
    logger.info("Found {} genes and {} transcripts", len(genes), len(transcripts))

    # If requested, consolidate transcripts with identical junctions
    #   (remove 5'/3' variation in redundantly spliced transcripts)
    if consolidate:
        data, genes = CONSOL.consolidate_transcripts(data,  outdir, consol_prefix, consol_outfiles,
                                                     genes, skip_interm,)
    # Output complexity measures using GTF data
    uniq_ex = COMP.calculate_complexity(outdir, data, skip_plots)

    # If requested, skip all other functions
    if complexity_only:
        logger.info("Complexity only option selected. Skipping all other functions.")
        return

    # Full processing
    out_fhs = open_output_files(outdir, outfiles)
    logger.debug("Output files: {}".format(outfiles))
    # Write out csv file headers
    if not skip_interm:
        if ea_mode == 'gene':
            out_fhs['er_fh'].write_text(",".join(er_df_cols) + '\n')
            out_fhs['ef_fh'].write_text(",".join(ef_df_cols) + '\n')
            out_fhs['ir_fh'].write_text(",".join(ir_df_cols) + '\n')
            out_fhs['ue_fh'].write_text(",".join(ue_df_cols) + '\n')
        else:
            out_fhs['ea_fh'].write_text(",".join(ea_df_cols) + '\n')
            out_fhs['td_fh'].write_text(",".join(TD.td_df_cols) + '\n')
        out_fhs['jc_fh'].write_text(",".join(jct_df_cols) + '\n')
    # Initialize concatenated pairwise transcript distance dataframe
    if ea_mode == "pairwise":
        td_data_cat = pd.DataFrame()
    else:
        er_data_cat = pd.DataFrame(columns=er_df_cols)
        ef_data_cat = pd.DataFrame(columns=ef_df_cols)
        ir_data_cat = pd.DataFrame(columns=['er_transcript_ids'])
    for gene in genes.groups:
        gene_df = data[data['gene_id'] == gene]
        transcripts = gene_df.groupby("transcript_id")
        transcript_groups = transcripts.groups
        number_of_transcripts = len(transcript_groups)
        if ea_mode == 'gene':
            try:
                er_data, ef_data, jct_data, ir_transcripts = do_ea(gene_df, mode='gene',
                                                                   keep_ir=keep_ir)
                if er_data is None:
                    continue
                if not skip_interm:
                    write_output(er_data, out_fhs, 'er_fh')
                    write_output(ef_data, out_fhs, 'ef_fh')
                    write_output(jct_data, out_fhs, 'jc_fh')
                er_data_cat = pd.concat([er_data_cat, er_data], ignore_index=True)
                ef_data_cat = pd.concat([ef_data_cat, ef_data], ignore_index=True)
                if len(ir_transcripts) > 0:
                    ir_data = pd.DataFrame(ir_transcripts, columns=['er_transcript_ids'])
                    ir_data_cat = pd.concat([ir_data_cat, ir_data], ignore_index=True)
            except ValueError as e:
                logger.error(e)
                continue
        else:
            if number_of_transcripts == 1:
                logger.info("Gene {} has a single transcript. Skipping", gene)
                continue
            ea_data, jct_data, td_data = ea_pairwise(gene_df)
            td_data_cat = pd.concat([td_data_cat, td_data], ignore_index=True)
            if not skip_interm:
                write_output(ea_data, out_fhs, 'ea_fh')
                write_output(jct_data, out_fhs, 'jc_fh')
                write_output(td_data, out_fhs, 'td_fh')
    # Output additional intermediate files for gene mode
    if ea_mode == 'gene' and not skip_interm:
        write_output(ir_data_cat, out_fhs, 'ir_fh')
        write_output(ue_data_cat, out_fhs, 'ue_fh')
    if not skip_plots:
        if ea_mode == 'pairwise':
            P1GP.plot_one_gtf_pairwise(outdir, td_data_cat)
        elif ea_mode == 'gene':
            P1GG.plot_one_gtf_gene(er_data_cat, ef_data_cat, ir_data_cat, uniq_ex, outdir)


def ea_two_files(f1_data, f2_data, gene_id, name1, name2):
    "Do EA (Event Analysis) for pairs of transcripts from two files for a gene."
    f1_transcripts = list(set(f1_data['transcript_id']))
    f2_transcripts = list(set(f2_data['transcript_id']))
    transcript_combos = list(itertools.product(f1_transcripts, f2_transcripts))
    logger.debug("Transcript combinations to process for {}: \n{}", gene_id, transcript_combos)
    ea_df = pd.DataFrame(columns=ea_df_cols)
    jct_df = pd.DataFrame(columns=jct_df_cols)
    td_df = pd.DataFrame(columns=TD.td_df_cols)
    for pair in transcript_combos:
        # try:
        tx_df_1 = f1_data[f1_data['transcript_id'] == pair[0]]
        tx_df_1_s = tx_df_1.assign(transcript_id=lambda x: x.transcript_id + '_' + name1)
        tx_df_2 = f2_data[f2_data['transcript_id'] == pair[1]]
        tx_df_2_s = tx_df_2.assign(transcript_id=lambda x: x.transcript_id + '_' + name2)
        tx_pair_data = pd.concat([tx_df_1_s, tx_df_2_s])
        try:
            ea_data, jct_data, td_data = do_ea(tx_pair_data, mode='pairwise')
            ea_df = ea_df.append(ea_data)
            jct_df = jct_df.append(jct_data)
            td_df = td_df.append(td_data)
        except ValueError:
            raise
        # except ValueError as e:
        #     logger.error(e)
        #     exit(1)
        #     #     # DEBUG, exit for now, just write to rejects later
        #     #     # continue
    return ea_df, jct_df, td_df


def process_two_files(infiles, outdir, outfiles, cpu_cores, out_pairs, complexity_only, skip_plots,
                      skip_interm, name1, name2):
    """Compare transcript pairs between two GTF files."""
    logger.debug("Input files: {}", infiles)
    if out_pairs != 'all':
        # Do not create full transcript distance output
        del(outfiles['td_fh'])
    else:
        # Do not create subset minimum distance output
        del(outfiles['md_fh'])
    # Do not make gene mode files (remove from outfiles)
    del(outfiles['er_fh'])
    del(outfiles['ef_fh'])
    del(outfiles['ir_fh'])
    del(outfiles['ue_fh'])

    infile_1 = infiles[0]
    infile_2 = infiles[1]
    in_f1 = read_exon_data_from_file(infile_1)
    in_f2 = read_exon_data_from_file(infile_2)

    # Calculate complexity of individual transcriptomes
    COMP.calculate_complexity(outdir, in_f1, skip_plots, name1)
    COMP.calculate_complexity(outdir, in_f2, skip_plots, name2)

    # Complexity only processing
    if complexity_only:
        logger.info("Complexity only option selected. Skipping all other functions.")
        return

    # Full processing
    if skip_interm:
        # Skip junction and ER/EF files (intermediate files)
        del(outfiles['ea_fh'])
        del(outfiles['jc_fh'])
        del(outfiles['gtf1_fh'])
        del(outfiles['gtf2_fh'])
        out_fhs = open_output_files(outdir, outfiles)
    else:
        out_fhs = open_output_files(outdir, outfiles)
        out_fhs['ea_fh'].write_text(",".join(ea_df_cols) + '\n')
        out_fhs['jc_fh'].write_text(",".join(jct_df_cols) + '\n')
    logger.debug("Output files: {}".format(outfiles))
    if out_pairs == 'all':
        out_fhs['td_fh'].write_text(",".join(MD.get_md_cols(name1, name2)) + '\n')
    else:
        out_fhs['md_fh'].write_text(",".join(MD.get_md_cols(name1, name2)) + '\n')
    f1_gene_names = set(in_f1['gene_id'])
    f2_gene_names = set(in_f2['gene_id'])
    only_f1_genes = f1_gene_names.difference(f2_gene_names)
    only_f2_genes = f2_gene_names.difference(f1_gene_names)
    odd_genes = only_f1_genes.union(only_f2_genes)
    f1_odds = in_f1[in_f1['gene_id'].isin(only_f1_genes)].copy()
    f2_odds = in_f2[in_f2['gene_id'].isin(only_f2_genes)].copy()
    if not skip_interm:
        write_gtf(f1_odds, out_fhs, 'gtf1_fh')
        write_gtf(f2_odds, out_fhs, 'gtf2_fh')
    common_genes = f1_gene_names.difference(odd_genes)
    valid_f1 = in_f1[in_f1['gene_id'].isin(common_genes)]
    valid_f2 = in_f2[in_f2['gene_id'].isin(common_genes)]
    f1_genes = valid_f1.groupby("gene_id")
    f1_transcripts = valid_f1.groupby(["gene_id", "transcript_id"])
    f2_genes = valid_f2.groupby("gene_id")
    f2_transcripts = valid_f2.groupby(["gene_id", "transcript_id"])
    logger.info("Found {} genes and {} transcripts in {} file", len(f1_genes),
                len(f1_transcripts), infile_1)
    logger.info("Found {} genes and {} transcripts in {} file", len(f2_genes),
                len(f2_transcripts), infile_2)
    gene_list = list(set(valid_f1['gene_id']))
    logger.debug("Genes to process: \n{}", gene_list)
    # Serial processing
    if cpu_cores == 1:
        ea_data, jct_data, td_data = loop_over_genes(gene_list, valid_f1, valid_f2, out_fhs,
                                                     cpu_cores, name1, name2)
        if not skip_interm:
            write_output(ea_data, out_fhs, 'ea_fh')
            write_output(jct_data, out_fhs, 'jc_fh')
        # Identify minimum pairs using transcript distances
        md_data = MD.identify_min_pair(td_data, out_pairs, name1, name2)
        if out_pairs != 'all':
            write_output(md_data, out_fhs, 'md_fh')
        else:
            write_output(md_data, out_fhs, 'td_fh')
        # Generate 2 GTF pairwise plots
        if not skip_plots:
            P2GP.plot_two_gtf_pairwise(outdir, md_data, f1_odds, f2_odds, name1=name1,
                                       name2=name2)
    # Parallel processing
    else:
        # Get lists for each process based on cpu_cores value
        geneLists = chunks(gene_list, cpu_cores)
        # Generate multiprocess Pool with specified number of cpus
        #     to loop through genes and calculate distances
        pool = Pool(cpu_cores)
        for genes in geneLists:
            subset_f1 = valid_f1[valid_f1['gene_id'].isin(genes)]
            subset_f2 = valid_f2[valid_f2['gene_id'].isin(genes)]
            pool.apply_async(loop_over_genes, args=(genes, subset_f1, subset_f2, out_fhs, cpu_cores,
                                                    name1, name2), callback=callback_results)
        pool.close()
        pool.join()
        ea_cat = pd.concat(ea_list)
        jct_cat = pd.concat(jct_list)
        td_cat = pd.concat(td_list)
        if not skip_interm:
            write_output(ea_cat, out_fhs, 'ea_fh')
            write_output(jct_cat, out_fhs, 'jc_fh')
        # Identify minimum pairs using transcript distances
        md_data = MD.identify_min_pair(td_cat, out_pairs, name1, name2)
        if out_pairs != 'all':
            write_output(md_data, out_fhs, 'md_fh')
        else:
            write_output(md_data, out_fhs, 'td_fh')
        # Generate 2 GTF pairwise plots
        if not skip_plots:
            P2GP.plot_two_gtf_pairwise(outdir, md_data, f1_odds, f2_odds, name1=name1,
                                       name2=name2)


def loop_over_genes(gene_list, valid_f1, valid_f2, out_fhs, cpu_cores, name1, name2):
    ea_data_list = []
    jct_data_list = []
    td_data_list = []
    # Loop over genes in provided gene list
    for gene in gene_list:
        f1_data = valid_f1[valid_f1['gene_id'] == gene]
        f2_data = valid_f2[valid_f2['gene_id'] == gene]
        try:
            ea_data, jct_data, td_data = ea_two_files(f1_data, f2_data, gene, name1, name2)
        except ValueError as e:
            logger.error(e)
            continue
        # Append output to list
        ea_data_list.append(ea_data)
        jct_data_list.append(jct_data)
        td_data_list.append(td_data)
    # Concatenate outputs
    ea_data_cat = pd.concat(ea_data_list)
    jct_data_cat = pd.concat(jct_data_list)
    td_data_cat = pd.concat(td_data_list)
    return [ea_data_cat, jct_data_cat, td_data_cat]
