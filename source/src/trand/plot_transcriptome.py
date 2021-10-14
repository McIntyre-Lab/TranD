#!/usr/bin/env python
# coding: utf-8

import pandas as pd
import numpy as np
import os
import matplotlib.pyplot as plt
import matplotlib
import seaborn as sns
import math


def calculate_donor_ratio_newINF(unit, feature, count_col, data):
    df = data.groupby([unit, feature])[count_col].count().reset_index().rename(columns = {count_col : "num_col_id"})
    df["flag_donor_acceptor"] = np.where(df["num_col_id"] > 1, 1, 0)
    df2 = df.groupby(unit).agg({feature : "count", "flag_donor_acceptor" : 'sum'})
    res = df2["flag_donor_acceptor"]/df2[feature]
    return res

def uniqueList(test_list):
    res = []
    for i in test_list:
        if i not in res:
            res.append(i)
    return res

def calculate_number(unit, feature, data, calLen = True, numbers = True):
    nr = data.shape[0]
    temp = {}
    dic_len = {}
    if calLen:
        startcol = feature.replace("_id", "_start")
        stopcol = feature.replace("_id", "_end")
    for i in range(0, nr):
        uni = data[unit][i]
        fea = data[feature][i].split("|")
        if calLen:
            length = data[stopcol][i] - data[startcol][i]
        else:
            length = 0
        if uni not in temp:
            dic_len[uni] = length
            temp[uni] = fea
        else:
            dic_len[uni] += length
            for j in fea:
                if j not in temp[uni]:
                    temp[uni].append(j)

    for key in temp:
        temp[key] = uniqueList(temp[key])                

    if numbers:             
        for key in temp:
            temp[key] = len(temp[key])

    return pd.Series(temp), pd.Series(dic_len)


def calculate_IR_ratio(data_er, irs):
    dat, length = calculate_number("gene_id", "er_transcript_ids", data_er, calLen=False, numbers = False)
    dat = pd.DataFrame(dat)
    dat = dat.reset_index()
    dat.columns = ["gene_id", "er_transcript_ids"]
    dat = dat.explode("er_transcript_ids")
    
    mergedDF = pd.merge(dat, irs, how = 'outer', on = 'er_transcript_ids', indicator = 'merge_check')
    mergedDF['flag_IR_transcript'] = np.where(mergedDF['merge_check'] == 'both', 1, 0)
    mergedDF['flag_all_transcript'] = 1
    num_transcript = mergedDF.groupby("gene_id")['flag_all_transcript'].sum()
    num_IR_transcript = mergedDF.groupby("gene_id")['flag_IR_transcript'].sum()
    df_ir_ratio = pd.DataFrame({"gene_id": num_transcript.index,
                                "num_transcript": num_transcript,
                                "num_IR_transcript": num_IR_transcript,
                                "prop_IR_transcript": num_IR_transcript/num_transcript}).reset_index(drop=True)
    return df_ir_ratio


def calculate_length(dat, feature):
    start = feature + '_start'
    stop = feature + '_end'
    ids = feature + "_id"
    res = dat[stop] - dat[start]
    res.index = dat[ids]
    return res

def generate_bins(data, nb = 10, merge_tail = 99.5, ceil = True, hardcut=False):
    bins = [min(data)]
    if not hardcut:
        cut = np.percentile(data, merge_tail)
        itv = (cut - min(data))/(nb-1)
    else:
        cut = 100
        if max(data) > cut:
            itv = (cut - min(data))/(nb-1)
        else:
            cut = np.percentile(data, merge_tail)
            itv = (cut - min(data))/(nb-1)
    
    for i in range(1, nb):
        bins.append(min(data) + i*itv)
    
    if not hardcut:
        bins = bins + [max(data)+0.1]
    else:
        if max(data) > cut:
            bins = bins + [max(data)+0.1]
        else:
            bins = bins + [100]

    if ceil:
        bins = [math.ceil(x) for x in bins]
    return bins


def set_boxAttri(box, color_list):
    for patch, filers, whis, med, color in zip(box['boxes'], box['fliers'], box['whiskers'], box['medians'], color_list):
        patch.set_color(color)
        filers.set_markeredgecolor(color)
        whis.set_color(color)
        med.set_color('black')
    return

def calculate_varFeatureRatio(varList, allList):
    df1 = allList.to_frame()
    df1.columns = ['n']
    df1.insert(0, 'name', allList.index)
    df2 = varList.to_frame()
    df2.columns = ['n']
    df2.insert(0, 'name', varList.index)
    df_merge = df1.merge(df2, how='left', on = 'name', suffixes=('_totalFeature', '_varFeature')).fillna(0)
    df_merge['ratio'] = df_merge['n_varFeature'] / df_merge['n_totalFeature']
    return df_merge

def hist_plot(dat, bins, colors, text=True, omitText = []):
    x_coords = []
    y_coords = []
    first_height = len(dat[dat.values < bins[1]])
    #second_height = len(dat[(dat.values >= bins[1]) & (dat.values < bins[2])])
    #ylim = max(first_height, second_height) * 1.1
    xlim = len(bins) * 1.1
    shrink = (len(bins)-1)/10
    plt.xlim(1, xlim)
    plt.plot(xlim, first_height)
    plt.grid(False)
    plt.xticks([])
    omitIter = 0
    for i in range(0, len(bins)-1):
        x = i+2
        y = len(dat[(dat.values < bins[i+1]) & (dat.values >= bins[i])])
        col = colors[i]
        plt.vlines(x, 0, y, color = col, linewidth=44/shrink)
        if text and (omitIter not in omitText):
            plt.text(x, y, y, fontsize = 10/shrink, horizontalalignment='center',verticalalignment='bottom')
        x_coords.append(x)
        y_coords.append(y)
        omitIter+=1
    plt.ylim(0, max(y_coords)*1.1)
    return x_coords, y_coords



def counts_projectToColor(lowCol, highCol, counts, firstBinColor = False):
    sumC = sum(counts)
    if sumC == 0:
        return tuple(map(lambda x: x/255 ,list(lowCol)))
    if firstBinColor:
        return tuple(map(lambda x: x/255, list(lowCol)))
    ratio = [i/sumC for i in counts]
    heat_colors = [((lowCol[0] + (highCol[0] - lowCol[0]) * k)/255, (lowCol[1] + (highCol[1] - lowCol[1]) * k)/255, (lowCol[2] + (highCol[2] - lowCol[2]) * k)/255) for k in ratio]
    return heat_colors



def heatmap_ratio_visualize(data, ax, heat_nrows, tick_pos, lowCol=(0, 255, 0), highCol=(255, 0, 0), kb = 10):
    heat_bins = generate_bins([0,1], nb = heat_nrows, merge_tail=100, ceil=False)
    heat_x_coords = heat_bins[0:heat_nrows]
    heat_x_coords = [e+0.05 for e in heat_x_coords]
    heat_starts = tick_pos[0:(kb+1)]
    heat_ends = tick_pos[1:(kb+2)]
    res = []
    for b in range(0, kb+1):
        ind = pd.Series(data[b])
        counts = []
        for i in range(0, len(heat_bins)-1):
            each = ind[(ind.values < heat_bins[i+1]) & (ind.values >= heat_bins[i])]
            counts.append(each.count())
        if b == 0:
            firstBinColor = True
        else:
            firstBinColor = False
        heat_colors = counts_projectToColor(lowCol, highCol, counts, firstBinColor = firstBinColor)
        ax.hlines(heat_x_coords, heat_starts[b], heat_ends[b], color = heat_colors, linewidth = 65/heat_nrows)
        res.append(counts)
    return res


def find_break_lim(y_coords):
    n = len(y_coords)
    sorted_ys = np.flip(np.sort(y_coords))
    sorted_idx = np.flip(np.argsort(y_coords))
    differences = sorted_ys[0:(n-1)] - sorted_ys[1:n]
    max_idx = np.argmax(differences)
    break_y_up = y_coords[sorted_idx[max_idx]]
    break_y_down = y_coords[sorted_idx[max_idx+1]]
    return break_y_up, break_y_down

def extract_genes_from_bins(geneList, data, results, types):
    if types == 'feature':
        each = data[data.index.isin(geneList)].tolist()
        each = [x for x in each if str(x) != 'nan']
    elif types == "ratio":
        each = data[data['name'].isin(geneList)]['ratio'].tolist()
    results.append(each)
    

def visualize_combined(list_main, list_e, list_eu, df_ratio_IR, f_donor_ratio, df_ratio_e, df_ratio_f_lensum, out, kbins=4):
    dat = list_main
    dat2 = list_e
    dat3 = list_eu

    xlabel = 'Number of transcripts per gene'
    ylabel = 'Number of genes'
    ylabel2 = 'Exon regions\nper gene'
    ylabel3 = 'Unique exons\nper gene'
    ylabel4 = 'Proportion of transcripts\nwith intron retention'
    ylabel5 = 'Proportion of exon regions\nwith alt. don/acc'
    ylabel6 = 'Proportion of variable\nexon region'
    ylabel7 = 'Proportion of variable \n nucleotides'
    titile = ''
    #lowColor = (0, 255, 0)
    lowColor = (224, 224, 224)
    highColor = (255, 0, 0)

    font2 = {'family' : 'serif',
             'weight' : 'normal',
             'size'   : 13,
    }
    font3 = {'family' : 'serif',
             'weight' : 'normal',
             'size'   : 8,
    }
    if max(dat) > 16:
        lastbin = max(dat)
    else:
        lastbin = 17
    bins = [1, 2, 4, 8, 16, lastbin]

    ticklinewidth = 0.5
    colors = sns.color_palette("hls", kbins+1)
    xlim = (kbins + 2) * 1.1
    shrink = (len(bins)-1)/10
    
    sns.set_style("whitegrid")
    fig = plt.figure(figsize=(13,10), facecolor='white',edgecolor='black')
    grid = plt.GridSpec(9, 1, wspace=0, hspace=0.15)
    y_coords = []
    x_coords = []
    ###############################
    # plot section 1
    ax1 = fig.add_subplot(grid[0,0])
    ax1.patch.set_facecolor('white')
    ax1.axes.get_yaxis().set_visible(True)
    ax1.spines['bottom'].set_visible(False)
    ax1.set_title(titile, fontsize=14)

    x_coords, y_coords = hist_plot(dat, bins, colors, text=False)

    break_y_up, break_y_down = find_break_lim(y_coords)
    
    ax1.set_ylim(break_y_up * 0.8, max(y_coords) * 1.2)
    for i in range(0, len(bins)-1):
        if y_coords[i] >= break_y_up:
            ax1.text(x_coords[i], y_coords[i], y_coords[i], fontsize=10/shrink,
                     horizontalalignment='center', verticalalignment='bottom')

    tick_pos = [1] + x_coords
    tick_pos = [i+0.5 for i in tick_pos]
    ax1.tick_params(axis='y', which='major', pad=-5, labelsize=6, direction='in')
    for tick in ax1.yaxis.get_majorticklabels():
        tick.set_horizontalalignment("left")

    ##################################
    # plot section 2
    ax2 = fig.add_subplot(grid[1:3,0])
    ax2.patch.set_facecolor('white')
    ax2.spines['top'].set_visible(False)
    ax2.spines['bottom'].set_visible(False)
#    ax2.set_ylabel(ylabel,font2)
    ax2.set_ylabel(ylabel, font3, labelpad = 50, rotation=20, verticalalignment='top')
    #ax2.plot(x_coords, y_coords, color = 'steelblue', lw = 1) # fit a line
    
    hist_plot(dat, bins, colors, text=False)
    ax2.set_ylim(0, break_y_down * 1.3)
    for i in range(0, len(bins)-1):
        if y_coords[i] <= break_y_down:
            ax2.text(x_coords[i], y_coords[i], y_coords[i], fontsize = 10/shrink, horizontalalignment='center',verticalalignment='bottom')

    ax2.tick_params(axis='y', which='major', pad=-5, labelsize=6, direction='in')
    for tick in ax2.yaxis.get_majorticklabels():
        tick.set_horizontalalignment("left")
    ################################
    # extract genes based on bins
    feature2s = []
    feature3s = []
    ratio_IRs = []
    ratio_donors = []
    ratio_es = []
    ratio_flensum = []
    for i in range(0, len(bins)-1):
        genes = dat[(dat.values < bins[i+1]) & (dat.values >= bins[i])].index.tolist()        
        extract_genes_from_bins(genes, dat2, feature2s, 'feature')
        extract_genes_from_bins(genes, dat3, feature3s, 'feature')
        extract_genes_from_bins(genes, df_ratio_IR, ratio_IRs, 'ratio')
        extract_genes_from_bins(genes, f_donor_ratio, ratio_donors, 'ratio')
        extract_genes_from_bins(genes, df_ratio_e, ratio_es, 'ratio')
        extract_genes_from_bins(genes, df_ratio_f_lensum, ratio_flensum, 'ratio')
        
    ################################
    # plot section 3
    ax3 = fig.add_subplot(grid[3,0])
    ax3.patch.set_facecolor('white')
    ax3.spines['top'].set_visible(False)
    ax3.spines['bottom'].set_visible(False)
    ax3.grid(False)
    ax3.set_ylim(0, max(dat2))
    ax3.set_ylabel(ylabel2, font3, labelpad = 50, rotation=20, verticalalignment='top')
    box = ax3.boxplot(feature2s, positions = x_coords, notch = False, patch_artist=True)
    set_boxAttri(box, colors)
    
    y_arrow = max(dat2)
    ax3.set_xlim(1, xlim)
    ax3.set_xticklabels([])
    ax3.vlines(tick_pos, 0, max(dat2), color = 'grey', linewidth=ticklinewidth)
    ax3.set_yticks([y_arrow * 0.1, round(y_arrow/2), y_arrow*0.95])
    ax3.set_yticklabels([0, round(y_arrow/2), y_arrow], ha = 'left')
    ax3.tick_params(axis='y', which='major', pad=-5, labelsize=6, direction='in')

    for i in range(0, len(bins)-1): 
        ax3.arrow(x_coords[i], y_arrow, 0, -y_arrow/20, overhang=y_arrow/20, head_width=0.2, head_length=1, width = 1, shape="full",fc=colors[i], ec=colors[i],alpha=0.9)


    ##################################
    # plot section 4
    ax4 = fig.add_subplot(grid[4,0])
    ax4.patch.set_facecolor('white')
    ax4.spines['top'].set_visible(False)
    ax4.spines['bottom'].set_visible(False)
    ax4.grid(False)
    ax4.set_ylim(0, max(dat3))
    ax4.set_ylabel(ylabel3, font3, labelpad = 50, rotation=20, verticalalignment='top')
    box2 = ax4.boxplot(feature3s, positions = x_coords, notch = False, patch_artist=True)
    set_boxAttri(box2, colors)

    y_arrow = max(dat3)
    ax4.set_xticklabels([])
    ax4.set_xlim(1, xlim)
    ax4.set_yticks([y_arrow * 0.1, round(y_arrow/2), y_arrow*0.95])
    ax4.set_yticklabels([0, round(y_arrow/2), y_arrow], ha = 'left')
    ax4.tick_params(axis='y', which='major', pad=-5, labelsize=6, direction='in')
    
    ax4.vlines(tick_pos, 0, max(dat3), color = 'grey', linewidth=ticklinewidth)
    for i in range(0, len(bins)-1): 
        ax4.arrow(x_coords[i], y_arrow, 0, -y_arrow/20, overhang=y_arrow/20, head_width=0.2, head_length=1, width = 1, shape="full",fc=colors[i], ec=colors[i], alpha=0.9)


    
    ###########################################
    # set heatmap resolution
    heat_nrows = 20
    heat_ylim = 1.1
    
    ####### plot section 5
    
    for hm, ratio, lab in zip([5, 6, 7], [ratio_IRs, ratio_donors, ratio_es], [ylabel4, ylabel5, ylabel6]):
        
        ax = fig.add_subplot(grid[hm,0])
        ax.patch.set_facecolor('white')
        ax.spines['top'].set_visible(False)
        ax.spines['bottom'].set_visible(False)
        ax.grid(False)
        ax.set_ylim(0, heat_ylim)
        
        heatmap_ratio_visualize(ratio, ax, heat_nrows, tick_pos, kb = kbins, lowCol = lowColor, highCol = highColor)
        ax.set_ylabel(lab, font3, labelpad = 60, rotation=20, verticalalignment='top')
        ax.set_yticks([0.1, 0.5, 0.95])
        ax.set_yticklabels([0, 0.5, 1], ha = 'left')
        ax.tick_params(axis='y', which='major', pad=-5, labelsize=6, direction='in')
        
        ax.set_xticklabels([])
        ax.set_xlim(1, xlim)
        ax.vlines(tick_pos, 0, heat_ylim-0.03, color = '#4682b4', linewidth=ticklinewidth)
        ax.hlines(0, tick_pos[0], tick_pos[-1], color = '#D3D3D3', linewidth=ticklinewidth*3)

    
     ##################################
    ####### plot section 6
    ax6 = fig.add_subplot(grid[8,0])
    ax6.patch.set_facecolor('white')
    ax6.spines['top'].set_visible(False)
    ax6.grid(False)
    ax6.set_ylim(0, heat_ylim)
    
    heatmap_ratio_visualize(ratio_flensum, ax6, heat_nrows, tick_pos, kb = kbins, lowCol = lowColor, highCol = highColor)
    
    ax6.set_xlim(1, xlim)
    ax6.set_xlabel(xlabel, font2, labelpad = 10) 
    ax6.set_xticks(tick_pos)
    ax6.set_xticklabels([str(round(i, 2)) for i in bins])
    
    ax6.set_ylabel(ylabel7, font3, labelpad = 50, rotation=20, verticalalignment='top')
    ax6.set_yticks([0.1, 0.5, 0.95])
    ax6.set_yticklabels([0, 0.5, 1], ha = 'left')
    ax6.tick_params(axis='y', which='major', pad=-5, labelsize=6, direction='in', length=1)
    ax6.vlines(tick_pos, 0, heat_ylim-0.03, color = '#4682b4', linewidth=ticklinewidth)
    
#    outfile = out + "/transcriptome_summary_plot.pdf"
    outfile2 = out + "/transcriptome_summary_plot.png" 
#    plt.savefig(outfile, dpi=600, format='pdf')
    plt.savefig(outfile2, dpi=600, format='png')
    return

def plot_transcriptome(er_data, ef_data, ir_data, uniqex_data, outdir):
    er_data['er_start'] = er_data['er_start'].astype(int)
    er_data['er_end'] = er_data['er_end'].astype(int)
    ef_data['ef_start'] = ef_data['ef_start'].astype(int)
    ef_data['ef_end'] = ef_data['ef_end'].astype(int)
    matplotlib.rcParams['pdf.fonttype'] = 42
    matplotlib.rcParams['ps.fonttype'] = 42
    ###############
    ####### process gene transcript input using er input

    num_list_main, g_lensum = calculate_number('gene_id', 'er_transcript_ids', er_data, calLen = False)

    ###############
    ###### process exon region input
    er_data_var = er_data[er_data['er_annotation_frequency'] != 'constitutive'].reset_index().drop('index', axis=1)

    num_list_e, e_lensum = calculate_number('gene_id', 'er_id', er_data)
    num_list_e_varXcript, e_lensum_varXcript = calculate_number('gene_id', 'er_id', er_data_var)

    df_ratio_e = calculate_varFeatureRatio(num_list_e_varXcript, num_list_e)
    ################
    ###### process fragment input
    ef_data_var = ef_data[ef_data['ea_annotation_frequency'] != 'constitutive'].reset_index().drop('index', axis=1)


    # calculate donor/acceptor ratio for each gene
    f_donor_ratio = pd.DataFrame()
    seri_fs = calculate_donor_ratio_newINF('gene_id', 'er_id', 'ef_id', ef_data)
    f_donor_ratio['name'] = seri_fs.index
    f_donor_ratio['ratio'] = seri_fs.values

    num_list_f, f_lensum = calculate_number('gene_id', 'ef_id', ef_data)
    num_list_f_varXcript, f_lensum_varXcript = calculate_number('gene_id', 'ef_id', ef_data_var)
    del(ef_data)
    df_ratio_f_lensum = calculate_varFeatureRatio(f_lensum_varXcript, f_lensum)
    ##################

    num_list_ue = uniqex_data.set_index('gene_id')['num_uniq_exon']

    df_ratio_IR = calculate_IR_ratio(er_data, ir_data)
    df_ratio_IR.columns = ['name', 'num_transcript', 'num_IR_transcript', 'ratio']
    del(er_data)

    visualize_combined(
        num_list_main,
        num_list_e,
        num_list_ue,
        df_ratio_IR,
        f_donor_ratio,
        df_ratio_e,
        df_ratio_f_lensum,
        out = os.path.normpath(outdir)
    )






