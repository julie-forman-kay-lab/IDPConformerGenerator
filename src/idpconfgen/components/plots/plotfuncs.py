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
    """# noqa: D205, D400
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
        degrees,
        n_conf,
        *,
        angtype=None,
        title=None,
        xlabel="Residues",
        colors='k',
        ylabel=None,
        title_fs=20,
        xlabel_fs=15,
        ylabel_fs=15,
        xticks=None,
        yticks=None,
        xticks_labels=None,
        yticks_labels=None,
        increment=None,
        xticks_fs=15,
        yticks_fs=15,
        fig_size=(10, 6),
        filename='plot_torsions.png',
        dpi=300,
        ):
    """# noqa: D205, D400
    Plot all torsion angle distributions as a scatter plot.
    Defaults to Phi torsion angles.

    Parameters
    ----------
    residues : integer
        Total number of residues of a protein in the ensemble.

    angles : np.ndarray, shape=(n_confs, residues)
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

    Returns
    -------
    errmsg : list
        List of errors to print to the log while plotting.
    """
    errmsg = []

    plt.figure(figsize=fig_size)
    for i in range(1, residues):
        res_ang = [i] * n_conf
        plt.scatter(res_ang,
                    angles[:, i - 1],
                    s=10,
                    facecolors='none',
                    edgecolors=colors
                    )

    if increment is None:
        increment = int(residues / fig_size[0])
        if increment == 0:
            increment = 1

    if xticks is None:
        xticks = np.arange(0, residues, increment)
        xticks_labels = get_xtick_lbl(xticks)
    plt.xticks(ticks=xticks, labels=xticks_labels, fontsize=xticks_fs)
    plt.yticks(fontsize=yticks_fs)

    if degrees: yunits = "deg"  # noqa: E701
    else: yunits = "rad"        # noqa: E701
    if not ylabel: ylabel = f'{angtype} ({yunits})'  # noqa: E701

    plt.ylabel(ylabel, fontsize=ylabel_fs)
    plt.xlabel(xlabel, fontsize=xlabel_fs)
    plt.title(title, fontsize=title_fs)
    try:
        plt.savefig(filename, dpi=dpi, transparent=False, bbox_inches='tight')
    except Exception as e:
        errmsg.append("There was an issue with your figure saving parameters:")
        errmsg.append(f"{e}")
        errmsg.append("Saving figure with 300 dpi in the working directory as `plot_torsions.png`...")  # noqa: E501
        plt.savefig('plot_torsions.png',
                    dpi=300,
                    transparent=False,
                    bbox_inches='tight'
                    )
    plt.close("all")

    return errmsg


def plot_fracSS(
        residues,
        frac_ss,
        *,
        sstype="Generic",
        title=None,
        xlabel="Residues",
        ylabel="Fraction Secondary Structure",
        title_fs=20,
        xlabel_fs=15,
        ylabel_fs=15,
        xticks=None,
        yticks=None,
        xticks_labels=None,
        yticks_labels=None,
        increment=None,
        colors=['#D55E00', '#0072B2', 'k', 'g', 'r', 'c', 'm', 'y', 'b'],  # noqa: B006, E501
        xticks_fs=15,
        yticks_fs=15,
        fig_size=(6, 4),
        filename='plot_fracSS.png',
        dpi=300,
        ):
    """# noqa: D205, D400
    Plot the fractional secondary structure information as a line graph.
    Made to handel DSSP as well as Ramachandran.

    Parameters
    ----------
    residues : integer
        Total number of residues of a protein in the ensemble.

    frac_ss : dictionary
        Contains fraction of each secondary structure per residue.

    n_conf : integer
        Number of conformers we're processing.

    filename : str, optional
        The file name with which the plot figure will be saved
        in disk. Defaults to plot_torsions.png
        You can change the file type by specifying its extention in
        the file name.

    fig_size : tuple of float or int
        The size ratio of the subplot in the figure.

    Returns
    -------
    errmsg : list
        List of errors to print to the log while plotting.
    """
    errmsg = []

    aa = [x + 1 for x in range(residues)]
    plt.figure(figsize=fig_size)

    if len(colors) < len(frac_ss):
        errmsg.append("Number of colors is less than number of secondary structures.")  # noqa: E501
        errmsg.append("Reverting to default colorset...")
        colors = ['#D55E00', '#0072B2', 'k', 'g', 'r', 'c', 'm', 'y', 'b']

    clr = 0
    for ss in frac_ss:
        plt.plot(aa, frac_ss[ss], label=f'{sstype} {ss}', color=colors[clr])
        clr += 1
    plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0)
    ax = plt.gca()

    if increment is None:
        increment = int(residues / fig_size[0])
        if increment == 0:
            increment = 1

    if xticks is None:
        xticks = np.arange(0, residues, increment)
        xticks_labels = get_xtick_lbl(xticks)
    ax.set_xticks(xticks)
    ax.set_xticklabels(labels=xticks_labels)
    ax.tick_params(axis='x', labelsize=xticks_fs)
    ax.tick_params(axis='y', labelsize=yticks_fs)
    plt.title(title, fontsize=title_fs)
    plt.xlabel(xlabel, fontsize=xlabel_fs)
    plt.ylabel(ylabel, fontsize=ylabel_fs)
    try:
        plt.savefig(filename, dpi=dpi, transparent=False, bbox_inches='tight')
    except Exception as e:
        errmsg.append("There was an issue with your figure saving parameters:")
        errmsg.append(f"{e}")
        errmsg.append("Saving figure with 300 dpi in the working directory as `plot_fracSS.png`...")  # noqa: E501
        plt.savefig('plot_fracSS.png',
                    dpi=300,
                    transparent=False,
                    bbox_inches='tight'
                    )
    plt.close("all")

    return errmsg


def plot_bend_angles(
        angles,
        names,
        *,
        title=None,
        xlabel="Backbone angles",
        colors=['#D55E00', '#0072B2', 'k', 'g'],  # noqa: B006
        ylabel="Bend Angle (Radians)",
        title_fs=20,
        xlabel_fs=15,
        ylabel_fs=15,
        xticks=None,
        yticks=None,
        xticks_labels=None,
        yticks_labels=None,
        increment=None,
        xticks_fs=15,
        yticks_fs=15,
        fig_size=(10, 6),
        filename='plot_bgeo.png',
        dpi=300,
        ):
    """
    Plot backbone bond angle distribution as a boxplot.

    Parameters
    ----------
    names : string
        Names of the bend types

    angles : np.ndarray, shape=(names, n_trimers)
        Container of the Y axis data.

    filename : str, optional
        The file name with which the plot figure will be saved
        in disk. Defaults to plot_torsions.png
        You can change the file type by specifying its extention in
        the file name.

    Returns
    -------
    errmsg : list
        List of errors to print to the log while plotting.
    """
    errmsg = []

    plt.figure(figsize=fig_size)
    plt.boxplot(angles, labels=names)
    plt.ylabel(ylabel, fontsize=ylabel_fs)
    plt.xlabel(xlabel, fontsize=xlabel_fs)
    plt.title(title, fontsize=title_fs)
    try:
        plt.savefig(filename, dpi=dpi, transparent=False, bbox_inches='tight')
    except Exception as e:
        errmsg.append("There was an issue with your figure saving parameters:")
        errmsg.append(f"{e}")
        errmsg.append("Saving figure with 300 dpi in the working directory as `plot_bgeo.png`...")  # noqa: E501
        plt.savefig('plot_bgeo.png',
                    dpi=300,
                    transparent=False,
                    bbox_inches='tight'
                    )

    plt.close("all")

    return errmsg
