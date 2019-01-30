import numpy as np
import copy



# Plot HiC matrix in a way we can store it in subplots
def plotHiC(matrix1, bad_color=None, axe=None, transform=np.log2, 
            rescale_zeros=True, title = None, **kwargs):
    """
    Plot HiC matrix

    :param matrix: list of lists with values to be plotted
    :param None bad_color: plots NaNs in a given color
    :param kwargs: extra parameters for the imshow function of matplotlib

    :returns: Nothing but pain and despair
    """

    matrix = copy.deepcopy(matrix1)

    if bad_color is not None:
        kwargs['cmap'] = plt.get_cmap(kwargs.get('cmap', None))
        kwargs['cmap'].set_bad(bad_color, 1.)

    if not isinstance(matrix, (np.ndarray, np.generic)):
        matrix = np.asarray(matrix)

    # remove zeroes from the matrix in order to avoid -inf with log transform
    if rescale_zeros:
        try:
            mini = min(matrix[np.nonzero(matrix)]) / 2.
        except ValueError:
            mini = 0.
        matrix[matrix==0] = mini

    matrix = np.ma.masked_where(np.isnan(matrix), transform(matrix))
    if axe == None:
        plt.imshow(matrix, interpolation='None', origin='lower', **kwargs)

        plt.xlim(0 - 0.5, len(matrix[0]) - 0.5)
        plt.ylim(-0.5, len(matrix) - 0.5)
        
        if title != None:
            plt.title(title)
    else:
        axe.imshow(matrix, interpolation='None', origin='lower', **kwargs)

        axe.set_xlim(0 - 0.5, len(matrix[0]) - 0.5)
        axe.set_ylim(-0.5, len(matrix) - 0.5)

        if title != None:
            axe.set_title(title)
