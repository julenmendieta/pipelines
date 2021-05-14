import mpl_toolkits.mplot3d as plt3d
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import pyplot as plt
import matplotlib.backends.backend_pdf
import matplotlib.cm as cm
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
import seaborn as sns
import matplotlib.ticker as plticker


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

def plotRidgePlot(df, longi, locusCh, viewPointReal,
                  resol, signifBinedM, positionsToMark={},
                  trackMatrx=[],
                title='', height=10, wide=20, ysize=8,
                 nxlabel = 10, ymax=False, cRanges=False):
    
    '''
    Function to plot 3-wise interaction data in a sea-plot or ridge plot form
        :param df: Pandas dataframe with two columns: column 'g' indicates 
            starting point (in genomic coordinates) of the bin of intrest.
            Second column (named by the chromosome of the region) indicates 
            frequency value at which we found the interactions.
        :param longi: Number of bins our region has
        :param viewPointReal: list with real coordinates (indicating starting
            point of bin) at wich we had the viewPoint or capture
        :param resol: resolution we are working with (in bp)
        :param signifBinedM: List ot lists indicating Zscore or significance values 
            for each interaction (to paint the interaction frequencies by these values)
        :param {} positionsToMark: Dictonary with bin postions to mark as key and label
            as value
        :param [] trackMatrx: List of lists with same data as df, but more handy to locate
            the significant points from signifBinedM
        :param '' title: Title for the plot
        :param 10 height: Integer indicating plot height
        :param 20 wide: Integer indicating plot width
        :param 8 ysize: Integer indicating size of column labels
        :param 10 nxlabel: Integer indicating number of x labels to be shown
        :param False ymax: Float indicating the upper limit of each line in the plot.
            If not set maximum value in dataset will be used, and shown in the title.
        :param False cRanges: Add dictionary with list of [starting, ending) point
            as key and color as value, to color by your desired range. Ej:
            {(0.0000001,0.5):'#fee5d9', (0.5,1):'#fcae91'}


    '''

    if cRanges == False:
        cRanges = {(0.0000001,0.5):'#fee5d9', (0.5,1):'#fcae91', (1,1.5):'#fb6a4a', 
                   (1.5,2):'#de2d26', (2,100):'#a50f15',
                  (-0.5,-0.0000001):'#eff3ff', (-1, -0.5):'#bdd7e7', (-1.5, -1):'#6baed6', 
                   (-2, -1.5):'#3182bd', (-100, -2):'#08519c'}


    overlap = 0.8
    regionStart1 = min(df['g'])
    posMarkReal = [(p * resol) + regionStart1 for p in positionsToMark]
    # change plotting style
    oldStyle = sns.axes_style()
    sns.set(style="white", rc={"axes.facecolor": (0, 0, 0, 0)})


    # calculate height
    height = height / float(longi)  # looking for a height of 10
    wide = wide / height  # looking for a wide of 20

    # Initialize the FacetGrid object
    #pal = sns.cubehelix_palette(longi, rot=-.25, light=.7)
    #pal = [[66/255., 109/255., 185/255.]] * longi  # blue
    pal = [[189/255., 189/255., 189/255.]] * longi  # grey
    g = sns.FacetGrid(df, row="g", hue="g", aspect=wide, height=height, palette=pal)

    ## Draw the densities in a few steps
    # tracks
    #g.map(sns.kdeplot, "x", clip_on=False, shade=True, alpha=0.8, lw=1.5, bw=.2, gridsize=100)
    g.map(plt.plot, locusCh, alpha=0.8, lw=1.5)

    # draw a point to know where is the second view
    #g.map(plt.plot, locusCh, alpha=0.8, lw=1.5)


    ## fill plot background
    if signifBinedM != False:
        for nax, ax in enumerate(g.axes.flat):
            ax.fill_between(ax.lines[0].get_xdata().astype(int),
                            ax.lines[0].get_ydata(0),
                            facecolor=pal[nax], alpha=0.6)
        for cr in cRanges:

            # This version would be to fill colors by significance
            #for nax, ax in enumerate(g.axes.flat):
            #    # select fill color
            #    toColor = [False for i in ax.lines[0].get_ydata(0)]
            #    for nx, x in enumerate(signifBinedM[nax]):
            #        if  cr[0] <= x < cr[1]:
            #            toColor[nx] = True

            #    ax.fill_between(ax.lines[0].get_xdata().astype(int),
            #                    ax.lines[0].get_ydata(0),
            #                    where=toColor,
            #                    facecolor=cRanges[cr], alpha=0.6)

            # this version to add dots by significance
            for nax, ax in enumerate(g.axes.ravel()):
                for nx, x in enumerate(signifBinedM[nax]):
                    if  cr[0] <= x < cr[1]:
                        xx = trackMatrx[nax][nx]
                        #xx = df.loc[df['g'] == (nx * resol) + regionStart].iloc[(nax + 1) / len(signifBinedM)][locusCh]
                        # plot dot for position
                        ax.plot(nx, xx, color=cRanges[cr], marker='o', markersize=6)



    else:
        for nax, ax in enumerate(g.axes.flat):
            ax.fill_between(ax.lines[0].get_xdata().astype(int),
                            ax.lines[0].get_ydata(0),
                            facecolor=pal[nax], alpha=0.6)


    #g = g.map(plt.fill_between, locusCh, 'g', alpha=1)
    # white border for the tracs
    #g.map(sns.kdeplot, "x", clip_on=False, color="w", lw=2, bw=.2)
    # baseline zero for the tracks
    #g.map(plt.axhline, y=0, lw=2, clip_on=False)

    # Define and use a simple function to label the plot in y axis coordinates
    def label(x, color, label):
        # change color when we have a label
        if int(label) in posMarkReal:
            color = 'green'
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


    # Remove axes details that don't play with overlap
    g.set_titles("")
    g.set(yticks=[])
    g.despine(bottom=True, left=True)

    # add red line in viewPoint and dot for second bin position
    for vi in viewPointReal:
        for nax, ax in enumerate(g.axes.ravel()):
            ylim = ax.get_ylim()
            #ax.axvline(x=(vi-regionStart1) / resol )
            if nax == 0:
                ax.plot(((vi-regionStart1) / resol, (vi-regionStart1) / resol), 
                    (0, max(trackMatrx[nax])), 
                    ls='-', color='red')
            else:
                ax.plot(((vi-regionStart1) / resol, (vi-regionStart1) / resol), 
                        (0, ax.get_ylim()[1] * (1 - overlap)), 
                        ls='-', color='red')
            # plot dot for position
            ax.plot(nax, 0, color='#bdbdbd', marker='o')
            if ymax == False:
                ymax = ylim[1]

            ax.set_ylim(ylim[0], ymax)

    # mark positions of interest
    if positionsToMark != False:
        for p in positionsToMark:
            for nax, ax in enumerate(g.axes.ravel()):
                ylim = ax.get_ylim()
                if nax == 0:
                    ax.plot((p, p), 
                        (0, max(trackMatrx[nax])), 
                        ls='--', color='green')
                else:
                    ax.plot((p, p), 
                        (0, ax.get_ylim()[1] * (1 - overlap)), 
                        ls='--', color='green')

                ax.set_ylim(ylim[0], ymax)

        print [(p, positionsToMark[p]) for p in sorted(positionsToMark)]



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

    loc = plticker.MultipleLocator(base=float(longi/nxlabel)) # this locator puts ticks at regular intervals
    for ax in g.axes.flat:
        ax.xaxis.set_major_locator(loc)
        labels0 = ax.get_xticks() # get x labels for x axis
        if len(labels0) != 0:
            nlabels = len(labels0) 
            newLabels = [0] * int(nlabels)
            for ni,i in enumerate(labels0):
                if ni != 0:
                    newLabels[ni] = '{:,}'.format(int(min(df['g']) + (i * resol)))
                else:
                    newLabels[ni] = 0
            ax.set_xticklabels(newLabels, rotation=30,
                              fontweight='bold') # set new labels
            xlim = ax.get_xlim()
            ax.set_xlim(xlim[0] ,longi - 1)

    # Set the subplots to overlap
    g.fig.subplots_adjust(hspace=-0.8)

    # add title
    g.fig.suptitle(title + ' maxVal=%s' %round(ymax, 4))

    # go back to old ploting style
    sns.set(oldStyle)


    return g.fig
