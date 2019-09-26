import mpl_toolkits.mplot3d as plt3d
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import pyplot as plt
import matplotlib.backends.backend_pdf
import matplotlib.cm as cm
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
import seaborn as sns


def plot3Dmatrix(longi, focusMultiGroups2, regionStartBin, filterPerc=False, title='',
                 saveFig=False, rangeVal=[False, False], toZscore=False, track=False,
                trackAlpha=0.2, elevAzim=[None, None], fig=False, show=True):
    
    '''
    Function to plot in a 3D cube al 3-wise interactions in the genome. Just one 6th of
        the matrix will have data to avoid redundancy
    
    :param longi: x, y or z axis length (should be the same)
    :param focusMultiGroups2: Dictionary with (x,y,z) list with 3D coordinates as key
        and interaction frequency value as value
    :param regionStartBin: Bin position at which our matrix starts
    :param False filterPerc: Float between 0 and 100 inclusive to choose the percentile 
        bellow which wont show any interaction data
    :param '' title: String containing the title for the plot
    :param False saveFig: Set to True if want to store data. Before calling the function
        use pdf = matplotlib.backends.backend_pdf.PdfPages(PATH)
        pdf.close() to save the file
    :param [False, False] rangeVal: Set minimum and maximum values for the colorMap
    :param False toZscore: Transform values to ZScore
    :param False track: List with (location, color) lists to be marked as a plane. Location
        should be writen in bin position integers. E.j. [(23, 'red'), (45, 'blue')]. Used to
        check the 3-wise interactions in relation with a bin of interest
    :param 0.2 trackAlpha: floating value indicating alpha applyed to track
    :param [False, False] elevAzim: elevation and azimut for the orientation of the 3D plot
    :param False fig: plt.figure() object
    :param True show: Whether to show or not the plot by plt.show()
    '''

    # Store maxVal and minVal
    maxVal = rangeVal[1]
    minVal = rangeVal[0]
    # Copy dictionary to be safe from modifications
    focusMultiGroups = copy.copy(focusMultiGroups2)
    # get keys from dictionary
    keys = sorted(focusMultiGroups.keys())

    # If we want to convert to Zscore values
    if toZscore == True:
        values = [focusMultiGroups[k] for k in keys]

        values = stats.zscore(values)
        for nk, k in enumerate(keys):
            focusMultiGroups[k] = values[nk]


    # H is a HiC matrix
    w = longi

    # Build a w x w x w sized matrix
    T = np.empty((w,w,w))
    T[:] = np.nan
    #T = np.zeros([w, w, w])



    if filterPerc == False:
        for i, j, k in keys:
            # Just show i <= j <= k
            T[i - regionStartBin][j - regionStartBin][k - regionStartBin] = focusMultiGroups[(i,j,k)]

            # test show all permutations
            #for p in itertools.permutations((i,j,k), 3):
            #    T[p[0] - regionStartBin][p[1] - regionStartBin][p[2] - regionStartBin] = focusMultiGroups[(i,j,k)]
    else:
        thress = np.percentile(focusMultiGroups.values(), filterPerc)
        for i, j, k in keys:
            # Just show i <= j <= k
            if focusMultiGroups[(i,j,k)] >= thress:
                T[i - regionStartBin][j - regionStartBin][k - regionStartBin] = focusMultiGroups[(i,j,k)]

            # test show all permutations
            #for p in itertools.permutations((i,j,k), 3):
            #    T[p[0] - regionStartBin][p[1] - regionStartBin][p[2] - regionStartBin] = focusMultiGroups[(i,j,k)]


    #np.save("data/Triangular_3Depint_50kb",T)

    # build variable with all matrix row positions
    t = np.arange(0, w)
    # transform tensor into SPARSE format
    # Build list with x,y,x coordinates and value
    #(x, y, z, value)
    Sp = np.zeros([w**3, 4])
    for i in t:
        for j in t:
            for k in t:
                s = (w**2)*i + w*j + k
                Sp[s,0] = i
                Sp[s,1] = j
                Sp[s,2] = k
                Sp[s,3] = T[i,j,k]

    #np.save("data/Triangular_3Depint_50kb_sparse.npy", Sp)


    #Sp[:,3] = 100*Sp[:,3]

    if fig == False:
        fig = plt.figure(figsize=(12,10))
        
    ax = fig.add_subplot(111, projection='3d')

    #colors = cm.hsv(the_fourth_dimension/max(the_fourth_dimension))



    # Store variables
    xs = Sp[:,0]
    ys = Sp[:,1]
    zs = Sp[:,2]
    xyz = Sp[:,3]

    ## colorbar
    # Set range for colors in cmap
    colmap = cm.ScalarMappable(cmap=cm.viridis)
    colmap.set_array(xyz)
    
    # check if there are negative numbers (important for the vlim
    #and cmap). cmap values will be set to the real minimum, but the 
    #coloring and alpha will be made just with adjusted positive numbers
    negative = False
    if np.nanmin(xyz) < 0:
        negative = True

    if maxVal == False:
        if minVal == False:
            maxVal = np.nanmax(xyz)
            if negative == False:
                colmap.set_clim(0, maxVal)
            else:
                minVal1 = np.nanmin(xyz)
                minVal = abs(minVal1)
                colmap.set_clim(minVal1, maxVal)
                # Need to deal with negative numbers for alpha
                minValR = abs(np.nanmin(xyz))
                maxValZS = maxVal + minVal
        else:
            maxVal = np.nanmax(xyz)
            if negative == False:
                colmap.set_clim(minVal, maxVal)
            else:
                colmap.set_clim(minVal, maxVal)
                # Need to deal with negative numbers for alpha
                minValR = abs(np.nanmin(xyz))
                maxValZS = maxVal + minVal
    else:
        if minVal == False:
            if negative == False:
                colmap.set_clim(0, maxVal)
            else:
                minVal1 = np.nanmin(xyz)
                minVal = abs(minVal1)
                colmap.set_clim(minVal1, maxVal)
                # Need to deal with negative numbers for alpha
                minValR = abs(np.nanmin(xyz))
                maxValZS = maxVal + minValR
        else:
            if negative == False:
                colmap.set_clim(minVal, maxVal)
            else:
                # Need to deal with negative numbers for alpha
                minValR = abs(np.nanmin(xyz))
                colmap.set_clim(minVal, maxVal)
                maxValZS = maxVal + minVal


    # Trick to modify alpha depending on value
    #positively proportional
    noVals = tuple([0] * 4)
    # set min(i/float(maxVal), 1) for cases when up limit is below the maximum value
    if toZscore == False:
        colors = [((colmap.to_rgba(i)[0], colmap.to_rgba(i)[1], 
                  colmap.to_rgba(i)[2], min(i / float(maxVal), 1)) if (i == i and i != 0)
                 else noVals) for i in xyz]
    else:
        # Take into account negative numbers of Zscore with the alpha
        colors = [((colmap.to_rgba(i)[0], colmap.to_rgba(i)[1], 
                  colmap.to_rgba(i)[2], min((i + minValR) / float(maxValZS), 1)) if (i == i and i != 0)
                 else noVals) for i in xyz]

    xs, ys, zs, colors = zip(*[(x, y, z, c) for x, y, z, s, c in zip(xs, ys, zs, xyz, colors) 
                             if s==s])

    ax.scatter(xs, ys, zs, marker='s', s=50, c=colors, depthshade=False, cmap='viridis',
              vmin=0, vmax=maxVal)

    #ax.scatter(xs, ys, zs, marker='s', s=50, c=Sp[:,3], depthshade=False, cmap='viridis',
    #          vmin=0, vmax=np.nanmax(Sp[:,3]))

    # Add colorbar
    #from mpl_toolkits.axes_grid1 import make_axes_locatable
    #divider = make_axes_locatable(ax)
    # create an axes on the right side of ax. The width of cax will be 5%
    # of ax and the padding between cax and ax will be fixed at 0.05 inch.
    #cax = divider.append_axes("right", size="5%", pad=0.05)
    #cb = plt.colorbar(colmap, cax=cax)

    
    cb = plt.colorbar(colmap,fraction=0.046, pad=0.04, ax=ax)

    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.set_zlabel('z')
    # set axis margins
    ax.set_xlim3d(0,w-1)
    ax.set_ylim3d(0,w-1)
    ax.set_zlim3d(0,w-1)

    # set title
    fig.suptitle(title)

    # Add lines to positions of interest
    if track != False:
        for tr, color in track:
            # Add small triangle
            verts = [np.array([0,0,tr]), np.array([0,tr,tr]), 
                     np.array([tr,tr,tr])]
            collection = Poly3DCollection([verts], alpha=trackAlpha, linewidths=1,
                                         edgecolors='grey')
            collection.set_facecolor(color)
            ax.add_collection3d(collection)

            # Add big triangle
            verts = [np.array([tr,tr,tr]), np.array([tr,tr,longi]), 
                     np.array([tr,longi,longi])]
            collection = Poly3DCollection([verts], alpha=trackAlpha, linewidths=1,
                                         edgecolors='grey')
            collection.set_facecolor(color)
            ax.add_collection3d(collection)

            # Add rectangle
            verts = [np.array([0,tr,tr]), np.array([tr,tr,tr]), np.array([tr,tr,longi]), 
                     np.array([0,tr,longi])]
            collection = Poly3DCollection([verts], alpha=trackAlpha, linewidths=1,
                                         edgecolors='grey')
            collection.set_facecolor(color)
            ax.add_collection3d(collection)

                        
            # Fill all y and x postion for all z positions
            #for i in range(0, li):
            #    xs = (li, li)
            #    ys = (li, longi)
            #    zs = (i, longi)
            #    for a in itertools.permutations([xs, ys, zs], 3):
            #        line_ = plt3d.art3d.Line3D(xs, ys, zs)
            #        line_.set_color(color)
            #        line_.set_linestyle('--')
            #        ax.add_line(line_)
            
    # orient the 3D matrix
    # Set the elevation and azimuth of the axes
    #'elev' stores the elevation angle in the z plane.
    #'azim' stores the azimuth angle in the x,y plane.
    elev = elevAzim[0]
    azim = elevAzim[1]
    ax.view_init(elev, azim)

    if saveFig == True:
        pdf.savefig(fig , bbox_inches='tight')

    if show == True:
        plt.show()
    #else:
    #    return colmap, ax
        
    
# Python program to check  
# if two lists have at-least  
# one element common 
# using set and property 
  
def common_member(a, b): 
    a_set = set(a) 
    b_set = set(b) 
    if (a_set & b_set): 
        return True 
    else: 
        return False


def plotRidgePlot(df, longi, locusCh, viewPoint, viewPointReal, regionStartBin,
                title='', heigh=10, wide=20, ysize=8):
    
    # change plotting style
    oldStyle = sns.axes_style()
    sns.set(style="white", rc={"axes.facecolor": (0, 0, 0, 0)})
    
    
    # calculate height
    height = heigh / float(longi)  # looking for a height of 10
    wide = wide / height  # looking for a wide of 20

    # Initialize the FacetGrid object
    #pal = sns.cubehelix_palette(longi, rot=-.25, light=.7)
    pal = [[66/255., 109/255., 185/255.]] * longi
    g = sns.FacetGrid(df, row="g", hue="g", aspect=wide, height=height, palette=pal)

    ## Draw the densities in a few steps
    # tracks
    #g.map(sns.kdeplot, "x", clip_on=False, shade=True, alpha=0.8, lw=1.5, bw=.2, gridsize=100)
    g.map(plt.plot, locusCh, alpha=0.8, lw=1.5)

    # draw a point to know where is the second view
    g.map(plt.plot, locusCh, alpha=0.8, lw=1.5)


    # fill plot background
    for nax, ax in enumerate(g.axes.flat):
        ax.fill_between(ax.lines[0].get_xdata().astype(int),
                        ax.lines[0].get_ydata(0),
                        facecolor=pal[nax], alpha=0.6)

    #g = g.map(plt.fill_between, locusCh, 'g', alpha=1)
    # white border for the tracs
    #g.map(sns.kdeplot, "x", clip_on=False, color="w", lw=2, bw=.2)
    # baseline zero for the tracks
    #g.map(plt.axhline, y=0, lw=2, clip_on=False)

    # Define and use a simple function to label the plot in axes coordinates
    def label(x, color, label):
        # change color when viewPoint
        if int(label) in viewPointReal:
            color = 'red'
        ax = plt.gca()
        #ax.set_ylim(min(x), max(x))
        up = height * 0.25
        left = -1/wide
        ax.text(left, up, '{:,}'.format(int(label)), fontweight="bold", color=color, 
                ha="left", va="center", transform=ax.transAxes, fontsize=ysize)

    g.map(label, locusCh)


    # Remove axes details that don't play will with overlap
    g.set_titles("")
    g.set(yticks=[])
    g.despine(bottom=True, left=True)

    # add red line in viewPoint and dot for second bin position
    for vi in viewPoint:
        for nax, ax in enumerate(g.axes.ravel()):
            ylim = ax.get_ylim()
            ax.plot((vi-regionStartBin, vi-regionStartBin), 
                    (0, ax.get_ylim()[1]), 
                    ls='-', color='red')
            # plot dot for position
            ax.plot(nax, 0, color='black', marker='o')
            ax.set_ylim(ylim)



    # change ticks to genomic coordinates
    # iterate over axes of FacetGrid
    #for ax in g.axes.flat:
    #    labels = ax.get_xticklabels() # get y labels for x axis
    #    nlabels = len(labels)
    #    adjust = nlabels = 10.0
    #    for i,l in enumerate(labels):
    #        if i != 0:
    #            labels[i] = int(labels[i].get_position()[0])
    #    ax.set_xticklabels(labels, rotation=30,
    #                      fontweight='bold') # set new labels

    for ax in g.axes.flat:
        labels0 = ax.get_xticklabels() # get y labels for x axis
        labels = list(df['g'])
        nlabels = len(labels)
        adjust = nlabels = 10.0
        newLabels = [0] * int(nlabels)
        for ni,i in enumerate(labels0):
            if ni != 0:
                pos = int(i.get_position()[0])
                newLabels[ni] = '{:,}'.format(labels[pos])
        ax.set_xticklabels(newLabels, rotation=30,
                          fontweight='bold') # set new labels

    # Set the subplots to overlap
    g.fig.subplots_adjust(hspace=-0.8)

    # add title
    g.fig.suptitle(title)
    
    # go back to old ploting style
    sns.set(oldStyle)

    return g.fig
