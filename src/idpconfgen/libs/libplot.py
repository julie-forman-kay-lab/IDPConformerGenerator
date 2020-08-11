"""
Plots :o
"""

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
    fig, axs = plt.subplots(ncols=subplots, figsize=[10 * subplots, data.shape[1] * 0.2])
    try:
        axs = axs.ravel()
    except AttributeError:
        axs = [axs]  # needed to facility the for loop bellow
    #plt.subplots_adjust(left=0.1, bottom=0.03, right=0.98, top=0.98)


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
