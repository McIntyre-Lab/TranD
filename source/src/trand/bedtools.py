#! /usr/bin/env python3
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Copyright © 2021 Oleksandr Moskalenko <om@rc.ufl.edu>
#
# Distributed under terms of the MIT license.

"""
Accessing pybedtools/bedtools functionality for the first-pass approach to event analysis in paired
transcript analysis. We changed to direct interval analysis for full-gene event analysis. Paired
transcript analsysis is being refactored into direct interval analysis for performance and to
streamline both types of analyses.
"""

from loguru import logger
from pybedtools import BedTool


def prep_bed_for_ea(data):
    """Event analysis"""
    exons = {}
    gene_id = data['gene_id'].unique()[0]
    tx_names = list(data['transcript_id'].unique())
    try:
        verify_same_strand(tx_names, data)
        verify_same_reference(tx_names, data)
    except ValueError:
        raise
    exons = {'gene_id': gene_id}
    exons['transcript_list'] = tx_names
    for tx in tx_names:
        exon_id = 1
        tx_data = data[data['transcript_id'] == tx]
        tx_bed = []
        for index, row in tx_data.iterrows():
            # Score is optional, use zero to keep the strand data in the next column
            tx_bed.append((row.seqname, str(int(row.start)-1), str(int(row.end)),
                           f"{row.transcript_id}_exon_{str(exon_id)}", '0', row.strand))
            exon_id += 1
        exons[tx] = tx_bed
    return exons


def verify_same_reference(tx_names, data):
    """Verify that exons are on the same contig/chrom"""
    locs = data['seqname'].unique()
    if len(locs) > 1:
        raise ValueError(f"Multiple contig/chromn locations found for {tx_names}, skipping.")


def verify_same_strand(tx_names, data):
    """Verify exons are on the same strand"""
    strands = data['strand'].unique()
    if len(strands) > 1:
        raise ValueError(f"Multiple strands found in {tx_names}, skipping.")


def convert_tx_bed_to_coords(data, mode="pairwise", keep_ir=True):
    """
    Data pre-processing for all EA analyses. Convert bed into coordinate data.
    """
    # Multiple transcripts
    tx_names = data['transcript_list']
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
    return tx_data, tx_coords, intron_data


def remove_ir_transcripts(tx_data, introns, tx_names, tx_coords):
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
    old_tx_names = tx_names
    tx_names = tx_data.keys()
    tx_coords_delete = [k for k in old_tx_names if k not in tx_names]
    for k in tx_coords_delete:
        del tx_coords[k]
    return tx_data, tx_names, tx_coords
