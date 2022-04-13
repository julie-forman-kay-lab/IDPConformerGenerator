"""
Different plotting functions for analysis of generated ensembles.

Inspired by https://github.com/joaomcteixeira/taurenmd/pull/60

Parameters described in each function concern data representation
and are considered of highest importance because their incorrect
use can mislead data analysis and consequent conclusions.

Plot style parameters concernning only plot style, i.e., colors,
shapes, fonts, etc... and which do not distort the actual data,
are not listed in the paremeter list bellow. We hope these
parameter names are self-explanatory and are listed in the function
definition.
"""
import numpy as np
from matplotlib import pyplot as plt

def get_xtick_lbl(xticks):
    """
    Mainly for residues, reworks the first x-tick label to
    residue 1 since it's default labeled as 0.
    
    Parameters
    ----------
    xticks : ndarray
        Array of evenly spaced values depends on the x-ticks.
    
    Returns
    -------
    list
        List of renamed xtick-labels
    """
    xticks[0] = 1
    return list(xticks)
    

def plot_torsions(
        residues,
        angles,
        n_conf,
        *,
        title=None,
        xlabel="Residues",
        ylabel="Torsions (deg)",
        xlabel_fs=20,
        ylabel_fs=20,
        xticks=None,
        yticks=None,
        xticks_labels=None,
        yticks_labels=None,
        increment=10,
        xticks_fs=15,
        yticks_fs=15,
        fig_size=(10, 6),
        filename='plot_torsions.png',
        dpi=300,
        ):
    """
    Plot all torsion angle distributions as a scatter plot.
    
    Parameters
    ----------
    residues : integer
        Total number of residues of a protein in the ensemble.
    angles : np.ndarray, shape=(angles, len(residues))
        Container of the Y axis data.
    n_conf : integer
        Number of conformers we're processing.
    filename : str, optional
        The file name with which the plot figure will be saved
        in disk. Defaults to plot_torsions.png
        You can change the file type by specifying its extention in
        the file name.
    fig_size : tuple of float or int
        The size ratio of the subplot in the figure.
    """

    plt.figure(figsize=fig_size)
    for i in range(1, len(residues)):
        res_ang = [i]*n_conf
        plt.scatter(res_ang, angles, s=10, facecolors='none', edgecolors='k')
    if xticks == None:
        xticks = np.arange(0,len(residues), increment)
        xticks_labels = get_xtick_lbl(xticks)
    plt.xticks(ticks= xticks, labels=xticks_labels, fontsize=xticks_fs)
    plt.yticks(fontsize = yticks_fs)
    plt.ylabel(ylabel, fontsize=ylabel_fs)
    plt.xlabel(xlabel, fontsize=xlabel_fs)
    plt.title(title)
    plt.savefig(filename, dpi=dpi, transparent=False, bbox_inches='tight')
    plt.close("all")

    return

def plot_fracSS_DSSP(
        residues,
        frac_dssp,
        *,
        title=None,
        xlabel="Residues",
        ylabel="Fraction Secondary Structure",
        xlabel_fs=20,
        ylabel_fs=20,
        xticks=None,
        yticks=None,
        xticks_labels=None,
        yticks_labels=None,
        increment=10,
        colors=('b', 'g', 'r', 'c', 'm', 'y', 'k'),
        xticks_fs=15,
        yticks_fs=15,
        fig_size=(6, 4),
        filename='plot_fracSS_6AIO.png',
        dpi=300,
        ):
    """
    Plot all torsion angle distributions as a scatter plot.
    
    Parameters
    ----------
    residues : integer
        Total number of residues of a protein in the ensemble.
    frac_DSSP : dictionary
        Contains fraction of each secondary structure code per residue.
    n_conf : integer
        Number of conformers we're processing.
    filename : str, optional
        The file name with which the plot figure will be saved
        in disk. Defaults to plot_torsions.png
        You can change the file type by specifying its extention in
        the file name.
    fig_size : tuple of float or int
        The size ratio of the subplot in the figure.
    """
    
    plt.figure(figsize=fig_size)
    for i in len(colors):
        for ss in frac_dssp:
            plt.plot(residues, frac_dssp[ss], label=f'DSSP {ss}', color=colors[i])
    plt.legend(bbox_to_anchor=(1.05,1), loc='upper left', borderaxespad=0)
    ax = plt.gca()
    if xticks == None:
        xticks = np.arange(0, residues, increment)
        xticks_labels = get_xtick_lbl(xticks)
    ax.set_xticks(xticks)
    ax.set_xticklabels(labels=xticks_labels)
    ax.tick_params(axis='x', labelsize=xticks_fs)
    ax.tick_params(axis='y', labelsize=yticks_fs)
    plt.title(title)
    plt.xlabel(xlabel, fontsize=xlabel_fs)
    plt.ylabel(ylabel, fontsize=ylabel_fs)
    plt.savefig(filename, dpi=dpi, transparent=False, bbox_inches='tight')
    plt.close("all")

    return