#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
#  care_graph_step_06.py
#
#  Copyright 2019 Carlos Eduardo Sequeiros Borja <casebor@gmail.com>
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with this program; if not, write to the Free Software
#  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
#  MA 02110-1301, USA.
#
#

import os
import json
import logging
import collections
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Circle, ConnectionPatch

from libs import utils
from libs import paths
from libs.configuration import initialize_configuration, TASK, ENV

logger = logging.getLogger(__name__)
logging.getLogger('matplotlib.font_manager').disabled = True


def _get_atom_color(atom_name):
    if "H" in atom_name:
        return "cyan"
    elif "Cl" in atom_name:
        return "lime"
    elif "I" in atom_name:
        return "purple"
    elif "Br" in atom_name:
        return "darkred"
    elif "X" in atom_name:
        return "green"
    elif "C" in atom_name:
        return "gray"
    elif "O" in atom_name:
        return "red"
    elif "N" in atom_name:
        return "blue"
    elif "S" in atom_name:
        return "orange"
    else:
        return "gray"


def _func(x, data):
    return str(int(x))


def _get_colors(labels):
    colors = {"4K2A":"tab:blue", "3A2M":"tab:orange", "2QVB":"tab:green",
              "1CQW":"tab:red", "1MJ5":"tab:purple", "4C6H":"tab:olive",
              "4H77":"tab:pink", "4E46":"tab:cyan", "4HZG":"tab:gray",
              "conformer_03":"tab:blue"}
    out = []
    for l in labels:
        out.append(colors.get(l, "tab:gray"))
    return out


def main():
    utils.init_logging()
    initialize_configuration(__file__)
    utils.set_logging_level(ENV["console_log_level"])
    paths.set_task_name(TASK["name"])
    out_path = os.path.join(paths.task_root(), "graphs", "filters")
    os.makedirs(out_path, exist_ok=True)
    plt.rcParams.update({'font.size': 6})
    logging.getLogger('matplotlib.font_manager').disabled = True

    clusters = collections.defaultdict(set)
    for clu_file in os.listdir(paths.substrates_clusters()):
        cluster = []
        with open(os.path.join(paths.substrates_clusters(), clu_file), 'r') as clu_stream:
            for line in clu_stream:
                cluster.extend(json.loads(line))
        for c in cluster:
            protein = c["name"].split('-')[0]
            fragment = "{}-{}".format(c["name"].split('_')[-1], c["model_index"])
            clusters[protein].add(fragment)
    bar_graph_data = [None] * len(clusters)
    pos = 0
    for key, value in clusters.items():
        print("Number of fragments for {}: {}".format(key, len(value)))
        bar_graph_data[pos] = (key, len(value))
        pos += 1
    bgd_labels = [bgd[0] for bgd in bar_graph_data]
    bgd_data = [bgd[1] for bgd in bar_graph_data]
    bgd_pos = [i for i in range(len(bgd_labels))]
    bgd_colors = _get_colors(bgd_labels)

    for f in range(len(os.listdir(paths.filters()))):
        try:
            file_fil = str(f).zfill(3) + ".json"
            clu_file = os.path.join(paths.substrates_clusters(), file_fil + "l")
            clu_inpath = os.path.join(paths.substrates_clusters(), clu_file)
            cluster = []
            with open(clu_inpath, 'r') as clu_stream:
                for line in clu_stream:
                    cluster.extend(json.loads(line))
            clu_dict = collections.defaultdict(int)
            for c in cluster:
                label = c["ref"]["mapping"].split(',')[-1]
                label = label.split('-')[0]
                clu_dict[label] += 1
            clu_labels = []
            clu_values = []
            for key, value in clu_dict.items():
                clu_labels.append(key)
                clu_values.append(value)
            if os.path.isfile(os.path.join(paths.filters(), file_fil)):
                inpath = os.path.join(paths.filters(), file_fil)
                with open(inpath, 'r') as instream:
                    built_filters = json.load(instream)
                cols = 1
                qtty = len(built_filters["filters"]) + 2
                while qtty - 3 > 0:
                    cols += 1
                    qtty -= 3
                fig, ax = plt.subplots(cols, 3)
                for i, c_filter in enumerate(built_filters["filters"]):
                    role = c_filter["keys"][0][2]
                    atom1 = c_filter["keys"][0][3]
                    atom2 = c_filter["keys"][0][4]
                    conns = []
                    if c_filter["prop"] == "distances":
                        filter_value = c_filter["intervals"][0][1]
                        conns.append(ConnectionPatch(xyA=(2, 0), xyB=(3, 0),
                                                     coordsA="data", coordsB="data",
                                                     arrowstyle='-'))
                        conns[0].set_linewidth(4)
                    else:
                        filter_value = c_filter["intervals"][0][0]
                        conns.append(ConnectionPatch(xyA=(1, 0), xyB=(4, 0),
                                                     coordsA="data", coordsB="data",
                                                     arrowstyle='->'))
                        conns.append(ConnectionPatch(xyA=(1, 0), xyB=(4, 1),
                                                     coordsA="data", coordsB="data",
                                                     arrowstyle='->'))
                        conns.append(ConnectionPatch(xyA=(1, 0), xyB=(4, -1),
                                                     coordsA="data", coordsB="data",
                                                     arrowstyle='->'))
                        conns[0].set_linewidth(1)
                        conns[1].set_linewidth(1)
                        conns[2].set_linewidth(1)
                    at1 = Circle((1, 0), 1, label=atom1, color=_get_atom_color(atom1))
                    at2 = Circle((4, 0), 1, label=atom2, color=_get_atom_color(atom2))
                    f_title = "{}\n{} from {} to {}".format(role, c_filter["prop"],
                                                            atom1, atom2)
                    row = i//3
                    col = int(i)
                    while col >= 3:
                        col -= 3
                    if c_filter["prop"] == "distances":
                        ax[row, col].text(1.3, 1, str(filter_value) + r'$\AA$', fontsize=12)
                    else:
                        ax[row, col].text(1.3, 1, str(filter_value), fontsize=12)
                    ax[row, col].text(0.7, 1.6, "{}\n{}".format(c_filter["keys"][0][1], atom1), fontsize=8)
                    ax[row, col].text(3.7, 1.4, atom2, fontsize=8)
                    ax[row, col].set_xlim(left=-1, right=6)
                    ax[row, col].set_ylim(bottom=-2, top=3)
                    ax[row, col].set_title(f_title)
                    ax[row, col].add_patch(at1)
                    ax[row, col].add_patch(at2)
                    for c in conns:
                        ax[row, col].add_artist(c)
                    ax[row, col].tick_params(axis='x', which="both", bottom=False,
                                             top=False, labelbottom=False)
                    ax[row, col].tick_params(axis='y', which="both", left=False,
                                             right=False, labelleft=False)
                pie_colors = _get_colors(clu_labels)
                ax[-1, -1].pie(clu_values, labels=clu_labels, colors=pie_colors,
                               autopct=lambda x: _func(x, clu_values), textprops={'fontsize':6})
                ax[-1, -2].bar(bgd_pos, bgd_data, color=bgd_colors)
                ax[-1, -2].set_xticks(bgd_pos)
                ax[-1, -2].set_xticklabels(bgd_labels, fontdict={"fontsize":4}, rotation=90)
                fig.set_size_inches(8.4, 2.4*cols)
                fig_path = os.path.join(out_path, "filter_{}.png".format(str(f).zfill(3)))
                plt.savefig(fig_path, dpi=300, format="png")
                plt.cla()
                plt.close()
        except IndexError:
            print("Probably there is an error in the file filter_{}.png\n"
                  "Please check the file.".format(str(f).zfill(3)))
            fig, ax = plt.subplots(cols, 3)
            for i, c_filter in enumerate(built_filters["filters"]):
                role = c_filter["keys"][0][2]
                atom1 = c_filter["keys"][0][3]
                atom2 = c_filter["keys"][0][4]
                conns = []
                if c_filter["prop"] == "distances":
                    filter_value = c_filter["intervals"][0][1]
                    conns.append(ConnectionPatch(xyA=(2, 0), xyB=(3, 0),
                                                 coordsA="data", coordsB="data",
                                                 arrowstyle='-'))
                    conns[0].set_linewidth(4)
                else:
                    filter_value = c_filter["intervals"][0][0]
                    conns.append(ConnectionPatch(xyA=(1, 0), xyB=(4, 0),
                                                 coordsA="data", coordsB="data",
                                                 arrowstyle='->'))
                    conns.append(ConnectionPatch(xyA=(1, 0), xyB=(4, 1),
                                                 coordsA="data", coordsB="data",
                                                 arrowstyle='->'))
                    conns.append(ConnectionPatch(xyA=(1, 0), xyB=(4, -1),
                                                 coordsA="data", coordsB="data",
                                                 arrowstyle='->'))
                    conns[0].set_linewidth(1)
                    conns[1].set_linewidth(1)
                    conns[2].set_linewidth(1)
                at1 = Circle((1, 0), 1, label=atom1, color=_get_atom_color(atom1))
                at2 = Circle((4, 0), 1, label=atom2, color=_get_atom_color(atom2))
                f_title = "{}\n{} from {} to {}".format(role, c_filter["prop"],
                                                        atom1, atom2)
                if c_filter["prop"] == "distances":
                    ax[i].text(2, 1, str(filter_value) + r'$\AA$')
                else:
                    ax[i].text(2, 1, str(filter_value))
                ax[i].text(0.7, 1.4, "{}\n{}".format(c_filter["keys"][0][1], atom1))
                ax[i].text(3.7, 1.4, atom2)
                ax[i].set_xlim(left=0, right=5)
                ax[i].set_ylim(bottom=-2, top=2)
                ax[i].set_title(f_title)
                ax[i].add_patch(at1)
                ax[i].add_patch(at2)
                for c in conns:
                    ax[i].add_artist(c)
                ax[i].tick_params(axis='x', which="both", bottom=False,
                                  top=False, labelbottom=False)
                ax[i].tick_params(axis='y', which="both", left=False,
                                  right=False, labelleft=False)
            pie_colors = _get_colors(clu_labels)
            ax[-1].pie(clu_values, labels=clu_labels, colors=pie_colors,
                       autopct=lambda x: _func(x, clu_values), textprops={'fontsize': 6})
            ax[-2].bar(bgd_pos, bgd_data, color=bgd_colors)
            ax[-2].set_xticks(bgd_pos)
            ax[-2].set_xticklabels(bgd_labels, fontdict={"fontsize":4}, rotation=90)
            fig.set_size_inches(8.4, 2.4)
            fig_path = os.path.join(out_path, "filter_{}.png".format(str(f).zfill(3)))
            plt.savefig(fig_path, dpi=300, format="png")
            plt.cla()
            plt.close()


if __name__ == "__main__":
    main()
