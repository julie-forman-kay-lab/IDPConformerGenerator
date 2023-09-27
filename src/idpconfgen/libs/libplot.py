"""Plot."""
import numpy as np
from matplotlib import pyplot as plt


def plot_distribution(
        data,
        subplots=1,
        title='Box Plot',
        ylabel='ylabel',
        xlabel='xlabel',
        vert=False,
        usermedians=None,
        labels=None,
        figname='boxplot.pdf',
        figdpi=300,
        ):
    """."""
    fig, axs = plt.subplots(ncols=subplots, figsize=[10 * subplots, data.shape[1] * 0.2])  # noqa: E501
    try:
        axs = axs.ravel()
    except AttributeError:
        axs = [axs]  # needed to facility the for loop bellow
    # plt.subplots_adjust(left=0.1, bottom=0.03, right=0.98, top=0.98)

    for i, ax in zip(range(subplots), axs):
        ax.set_title(title[i], fontsize=15, weight='bold')
        ax.invert_yaxis()
        ax.boxplot(
            data[:, i::subplots],
            vert=vert,
            usermedians=usermedians,
            labels=labels[i::subplots],
            )
        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)
        ax.tick_params(axis='both', labelsize=13)

    plt.tight_layout(h_pad=0.01, w_pad=0.01)
    fig.savefig(figname, dpi=figdpi)


def plot_distribution_list(
        data,
        subplots=1,
        title='Box Plot',
        ylabel='ylabel',
        xlabel='xlabel',
        vert=False,
        usermedians=None,
        usermean=None,
        userstd=None,
        labels=None,
        figname='boxplot.pdf',
        figdpi=300,
        ):
    """."""
    h = len(labels[0])
    w = 20
    fig, axs = plt.subplots(ncols=subplots, figsize=[w, h])
    try:
        axs = axs.ravel()
    except AttributeError:
        axs = [axs]  # needed to facility the for loop bellow

    for i, ax in zip(range(subplots), axs):
        ax.set_title(title[i], fontsize=15, weight='bold')
        ax.invert_yaxis()
        ax.boxplot(
            data[i],
            vert=vert,
            usermedians=usermedians,
            labels=labels[i],
            zorder=2,
            )
        ax.set_xlabel(xlabel, weight='bold', fontsize=13)
        ax.set_ylabel(ylabel, weight='bold', fontsize=13)
        ax.tick_params(axis='both', labelsize=13)

        if usermean:
            ax.axvline(usermean[i], color='grey', zorder=0)

        if userstd:
            M = usermean[i]
            ax.axvline(M + userstd[i], color='lightgrey', zorder=0)
            ax.axvline(M - userstd[i], color='lightgrey', zorder=0)
            ax.axvline(M + userstd[i] * 2, color='lightgrey', zorder=0)
            ax.axvline(M - userstd[i] * 2, color='lightgrey', zorder=0)
            ax.axvline(M + userstd[i] * 3, color='lightgrey', zorder=0)
            ax.axvline(M - userstd[i] * 3, color='lightgrey', zorder=0)

    plt.tight_layout(h_pad=0.01, w_pad=0.04)
    fig.savefig(figname, dpi=figdpi)


def plot_contacts_matrix(matrix, sequence, output, dpi=300):
    """
    Plot matrix heatmap from `contact_matrix` function.

    Parameters
    ----------
    matrix : np.ndarray
        Probability matrix of contacts
    
    sequence : str or list
        Single sequence for intra- or two sequences for inter-
    
    output : str
        Path/filename to output the plot
    
    dpi : int
        dpi of plot to save
    """
    lw = 10
    plt.figure(figsize=(lw, lw))
    im = plt.imshow(matrix, cmap='plasma', interpolation='nearest')
    plt.title('Contacts Frequency Heatmap', fontsize=18)
    plt.colorbar().set_label(label="Frequency", size=18)
    im.figure.axes[1].tick_params(axis="y", labelsize=18)
    
    if type(sequence) is list:
        seq1 = sequence[0]
        seq2 = sequence[1]
        
        len_s1 = len(seq1)
        len_s2 = len(seq2)
        if len_s1 > 100:
            fs1 = len_s1 / 30
        else:
            fs1 = len_s1 / lw
        if len_s2 > 100:
            fs2 = len_s2 / 30
        else:
            fs2 = len_s2 / lw
        plt.xticks(np.arange(len_s1), seq1, fontsize=fs1)
        plt.yticks(np.arange(len_s2), seq2, fontsize=fs2)
    else:
        len_seq = len(sequence)
        if len_seq > 100:
            fs = len_seq / 30
        else:
            fs = len_seq / lw
        plt.xticks(np.arange(len_seq), sequence, fontsize=fs)
        plt.yticks(np.arange(len_seq), sequence, fontsize=fs)
    
    plt.tight_layout(h_pad=0.01, w_pad=0.04)
    plt.savefig(output, dpi=dpi)
