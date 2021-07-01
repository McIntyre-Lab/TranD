#! /usr/bin/env python3
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Copyright © 2021 Oleksandr Moskalenko <om@rc.ufl.edu>
#
# Distributed under terms of the MIT license.

"""
Visualize fragments and events in a terminal
"""

def print_region(data):
    """
    Print a scaled region to screen with labels and dashes for visualization.
    """
    label = ""
    region = ""
    r_end = 0
    for i in (sorted(data.keys())):
        start = i
        name = data[start][0]
        end = data[start][1]
        ilen = end - start - 1
        if ilen < 0:
            ilen = 0
        gap = start - r_end - 1
        if gap <= 0:
            gap = 0
        region += gap * " " + '|'
        label += (gap + 1) * " "
        region += (ilen * "-") + '|'
        r_end = end + 1
        l_pad = round((end - start) / 2 - len(name) / 2) - 1
        if l_pad < 0:
            l_pad = 0
        label += l_pad * " "
        label += name
        # print(f"ER: {name}, {start}, {end}, {ilen}, {r_end}, {gap}, {l_pad}")
        # print(label)
        # print(region)
        l_fill = len(region) - len(label)
        label += l_fill * " "
    print(label)
    print(region)


def plot_ers(ers, tx1, tx2):
    """
    Plot ERs visually when debugging.
    """
    start, end, er_dict, tx1_dict, tx2_dict = [], [], {}, {}, {}
    for er in ers:
        start.append(er.start)
        end.append(er.end)
    g_start = min(start)
    g_end = max(end)
    er_range = g_end - g_start + 1
    scale = 200.0 / er_range
    print(tx1_dict)
    for i in ers:
        i_name = i.name.split(':')[1]
        er_dict[round((i.start - g_start) * scale)] = [i_name, round((i.end - g_start) * scale)]
    for i in tx1:
        tx1_name = i.name.split("_")[0]
        i_name = "E" + i.name.split('_')[-1]
        tx1_dict[round((i.start - g_start) * scale)] = [i_name, round((i.end - g_start) * scale)]
    for i in tx2:
        tx2_name = i.name.split("_")[0]
        i_name = "E" + i.name.split('_')[-1]
        tx2_dict[round((i.start - g_start) * scale)] = [i_name, round((i.end - g_start) * scale)]
    gene_name = next(iter(er_dict.values()))[0].split(':')[0]
    print(f"Exonic Regions for {gene_name} gene")
    print('|' + '-' * 198 + '|')
    print_region(er_dict)
    print('|' + '-' * 98 + '|')
    print(f"Transcript 1: {tx1_name}")
    print_region(tx1_dict)
    print('|' + '-' * 98 + '|')
    print(f"Transcript 2: {tx2_name}")
    print_region(tx2_dict)
    print('|' + '-' * 98 + '|')



