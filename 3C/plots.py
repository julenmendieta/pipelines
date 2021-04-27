import numpy as np
import copy
from matplotlib import pyplot as plt
import matplotlib.patches as patch
from matplotlib.colors import LinearSegmentedColormap
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.patches as mpatches

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
        im = plt.imshow(matrix, interpolation='None', origin='lower', **kwargs)

        plt.xlim(0 - 0.5, len(matrix[0]) - 0.5)
        plt.ylim(-0.5, len(matrix) - 0.5)
        
        if title != None:
            plt.title(title)
    else:
        im = axe.imshow(matrix, interpolation='None', origin='lower', **kwargs)

        axe.set_xlim(0 - 0.5, len(matrix[0]) - 0.5)
        axe.set_ylim(-0.5, len(matrix) - 0.5)

        if title != None:
            axe.set_title(title)

    return im

################## Functions for the stats analysis
def running_mean(xbin,x, N):
    cumsum = np.cumsum(np.insert(x, 0, 0))
    y=(cumsum[N:] - cumsum[:-N]) / N
    for ll in range(N/2):
    	print ll,len(xbin)-ll-1,xbin[-1],xbin[-2]
    	print xbin[ll]
    	print xbin[len(xbin)-ll-1]

    	y=np.insert( y,ll,  np.nan)
    	y=np.append( y,  np.nan)
    print len(xbin)
    print len(y)
    return zip(xbin,y)

#def GetXY(fileINname,diff=0):
#	list_dat_contac=[]
#	fileIN=open(fileINname,'r')
#	for l in fileIN.readlines()[1:]:
#	#mod=int(l.split()[0])
#	#print l.split()
#		part=int(l.split()[0])+diff
#		score=float(l.split()[1])
#		dev=float(l.split()[2])
#		list_dat_contac.append([part,score,score-dev,score+dev,dev])
#	xco,yco,eMco,ePco,dev=zip(*list_dat_contac)
#	return xco,yco,eMco,ePco,dev

def GetXY(fileINname,diff=0):
	list_dat_contac=[]
	fileIN=open(fileINname,'r')
	for l in fileIN.readlines()[1:]:
	#mod=int(l.split()[0])
	#print l.split()
		part=int(l.split()[0])+diff
		score=float(l.split()[1])
		list_dat_contac.append([part,score])
	xco,yco=zip(*list_dat_contac)
	return xco,yco

def GetXYA(fileINname,diff=0):
	list_dat_contac=[]
	fileIN=open(fileINname,'r')
	for l in fileIN.readlines()[1:]:
	#mod=int(l.split()[0])
	#print l.split()
		part=int(l.split()[0])+diff
		score=float(l.split()[1])*100.
		dev=float(l.split()[2])
		list_dat_contac.append([part,score,score-dev,score+dev,dev])
	xco,yco,eMco,ePco,dev=zip(*list_dat_contac)
	return xco,yco,eMco,ePco,dev

def GetXY_smooth2(fileINname,diff=0):
	list_dat_contac=[]
	fileIN=open(fileINname,'r')
	for l in fileIN.readlines()[1:]:
	#mod=int(l.split()[0])
	#print l.split()
		part=int(l.split()[0])+diff
		score=float(l.split()[9])
		dev=float(l.split()[10])
		list_dat_contac.append([part,score,score-dev,score+dev,dev])
	xco,yco,eMco,ePco,dev=zip(*list_dat_contac)
	return xco,yco,eMco,ePco,dev

def GetXY_smooth3(fileINname,diff=0):
	list_dat_contac=[]
	fileIN=open(fileINname,'r')
	for l in fileIN.readlines()[1:]:
	#mod=int(l.split()[0])
	#print l.split()
		part=int(l.split()[0])+diff
		score=float(l.split()[3])
		dev=float(l.split()[4])
		list_dat_contac.append([part,score,score-dev,score+dev,dev])
	xco,yco,eMco,ePco,dev=zip(*list_dat_contac)
	return xco,yco,eMco,ePco,dev

def GetXY_smooth4(fileINname,diff=0):
	list_dat_contac=[]
	fileIN=open(fileINname,'r')
	for l in fileIN.readlines()[1:]:
	#mod=int(l.split()[0])
	#print l.split()
		part=int(l.split()[0])+diff
		score=float(l.split()[5])
		dev=float(l.split()[6])
		list_dat_contac.append([part,score,score-dev,score+dev,dev])
	xco,yco,eMco,ePco,dev=zip(*list_dat_contac)
	return xco,yco,eMco,ePco,dev


# Plot stats for one model object
def plotearStats(indir, clust, marksAll, marksBar, regionStart, regionEnd, resol, rotation=90,
                 title='', titleSize = 20, figsize = (20, 10), diff=0, propor=[4], radius="?",
                 fromOri = 0, fromEnd = 0):
    fontSize = 15
    # propor is the proportion of size for normal plots against grey bar info plots
    # marksAll is to add the marks in all the plots, and MarksBar is to add marks just in the
    #top bar (is a list of lists)
    xa,ya,eMa,ePa,da=GetXYA(indir + 'accessibility.cluster%s.dat' %clust,diff=diff)
    xc,yc=GetXY(indir + 'consistency.cluster%s.dat' %clust,diff=diff)
    xi,yi,eMi,ePi,di=GetXY_smooth4(indir + 'interactionscluster%s.dat' %clust,diff=diff)

    nplot = len(marksBar) + 3 #we alway have 3 plots at least
    #f, (ax3, ax1, ax2, ax5) = plt.subplots(4, sharex=True, sharey=False,figsize=figsize)
    f, allAxis = plt.subplots(nplot, sharex=True, sharey=False,figsize=figsize,
                              gridspec_kw = {'height_ratios':[1] *  ((nplot) - 3) + propor * 3})
    # change tic size
    #f.tick_params(axis='both', which='major', labelsize=10)
    #f.tick_params(axis='both', which='minor', labelsize=8)
    #gs1 = gridspec.GridSpec(figsize[0], figsize[1])
    #gs1.update(wspace=0.025, hspace=0.05)

    allAxis[-3].plot(xc, yc,'k-')
    allAxis[-3].set_title('Consistency per particle', fontsize=fontSize)
    #allAxis[-3].axhline(y=np.mean(yc),c="blue")#,linewidth=0.2)
    allAxis[-3].set_ylabel("Consistency %",rotation=rotation, fontsize=fontSize)
    allAxis[-3].set_ylim(0, 100)

    allAxis[-2].plot(xa, ya,'k-')
    allAxis[-2].set_title('Accessibility per particle', fontsize=fontSize)
    allAxis[-2].axhline(y=np.mean(np.nan_to_num(ya)),c="blue")#,linewidth=0.2)
    allAxis[-2].set_ylabel("Accessibility",rotation=rotation,multialignment='center',
                          fontsize=fontSize)
    allAxis[-2].fill_between(xa,eMa,ePa,alpha=0.5, edgecolor='#bebebe', facecolor='#bebebe' )
    allAxis[-2].set_ylim(0, 100)

    allAxis[-1].plot(xi, yi,'k-')
    allAxis[-1].set_title('Contacts within sphere per particle', fontsize=fontSize)
    allAxis[-1].set_ylabel("Number of contact",rotation=rotation, multialignment='center',
                          fontsize=fontSize)
    allAxis[-1].fill_between(xi,eMi,ePi,alpha=0.5, edgecolor='#bebebe', facecolor='#bebebe' )
    allAxis[-1].axhline(y=np.mean(np.nan_to_num(yi)),c="blue")#,linewidth=0.2)

    # Turn to grey the info rectangles
    for axi in allAxis[:-3]:
        rect2 = patch.Rectangle((np.min(xa),0), np.max(xa)-np.min(xa), 1, color='lightgrey')
        axi.add_patch(rect2)

    # Add marks in all the plots
    for markStart, markEnd, color in marksAll:
        # Now we add rectnagles in ohur possitions of interest
        # The variables that matter are markStart and markEnd
        high=float(round(np.nanmax(ya)-np.nanmin(ya))+1000)
        min_s=float(round(np.nanmin(ya))-100)
        rect2 = patch.Rectangle((markStart,min_s), markEnd - markStart, high, color=color,alpha=0.3)
        allAxis[-2].add_patch(rect2)

        high=float(round(np.nanmax(yc)-np.nanmin(yc))+1000)
        min_s=float(round(np.nanmin(yc))-100)
        #rect2 = patch.Rectangle((np.min(b[k][0]),0), np.max(b[k][0])-np.min(b[k][0]), 1, color='#ff7f24')
        rect2 = patch.Rectangle((markStart,min_s), markEnd - markStart, high, color=color,alpha=0.3)
        allAxis[-3].add_patch(rect2)

        high=float(round(np.nanmax(yi)-np.nanmin(yi))+1000)
        min_s=float(round(np.nanmin(yi))-100)
        #rect2 = patch.Rectangle((np.min(b[k][0]),0), np.max(b[k][0])-np.min(b[k][0]), 1, color='#ff7f24')
        rect2 = patch.Rectangle((markStart,min_s), markEnd - markStart, high, color=color,alpha=0.3)
        allAxis[-1].add_patch(rect2)

        #rect2 = patch.Rectangle((markStart,0), markEnd - markStart, 1, color=color)
        #ax3.add_patch(rect2)

    # Add marks just in the top bar
    for naxi, marks in enumerate(marksBar):
        for markStart, markEnd, color in marks:
            rect2 = patch.Rectangle((markStart,0), markEnd - markStart, 1, color=color)
            allAxis[naxi].add_patch(rect2)

        #allAxis[naxi].set_ylim([0, 8])
        allAxis[naxi].axis('off')
        #allAxis[naxi].set_aspect('equal')

    f.subplots_adjust(hspace=0.3)

    # we want reference genome positions and enough labels to see start and end
    mini = 10
    # get the number of divisions where labels would be equally spread
    start = 0
    end = ((regionEnd - regionStart) / resol) + 1
    plt.xlim(start, end)
    #start, end = allAxis[-3].get_xlim()
    #print start, end
    binrange = range(int(start), int(end) + 1)
    posrange = range(regionStart, regionEnd + resol, resol)
    #for i in range(8, 22):
    #    divisor1 = len(binrange) / float(i)
    #    if mini >= divisor1/int(divisor1):
    #        mini = min(mini, divisor1/int(divisor1))
    #        divisor = int(divisor1)
    divisor = 9
    # use range to select the values in bin (value used at the time to plo) and
    #genomic position (value that we actually want to see)

    binx = [binrange[i] for i in range(0, len(binrange), divisor)]
    posx = ['{:,}'.format(posrange[i]) for i in range(0, len(posrange), divisor)]

    # check if last position is too close
    #if binx[-1] - binrange[-1] < divisor/1.5:
    #    binx = binx[:-1] + [binrange[-1]]
    #    posx = posx[:-1] + [posrange[-1]]
    #else:
    #    binx = binx + [binrange[-1]]
    #    posx = posx + [posrange[-1]]

    plt.xticks(binx,
               posx,
               rotation=45)
    plt.xlabel('Genomic coordinates', fontsize=fontSize)

    if title != '':
        f.suptitle(title, size=titleSize)

    plt.xlim(np.min(xa) + fromOri, np.max(xa) - fromEnd)

    plt.show()
    return f


# Example of call
# pdf = matplotlib.backends.backend_pdf.PdfPages(outpath + 'statsMultiplotsAll.pdf')
# # Beggin to plot
# for fi in modFiles:
#     flag = fi.split('/')[-1].split('_')[:2]
#     # Get region index
#     nregion = int(fi.split('/')[-2][3:])
#     # get regions start and end position (and chromosome)
#     reg = regiones['reg' + str(nregion)]
#     chrom = reg[0]
#     regionStart = reg[1]
#     regionEnd = reg[2]
    
#     region = fi.split('/')[-2]
#     fromOri = (originalRegionPos[region][1]/resol) - (reg[1] / resol)
#     fromEnd = (reg[2] / resol) - (originalRegionPos[region][2] / resol)

#     colors = ["red", "blue", "black", 'green']
#     marksAll = []
#     marksBar2 = []
#     # esto todavia no funcion, tengo q ver como poner el titulo en toda la figura con subplots
#     title = 'Model statistics in region %s of %s' %(flag[1][3:], flag[0])
    
#     # Add enhancers in the marking step
#     for en in enhAll[region]:
#         marksAll.append((en, en, colors[0]))

#     # If we want just promoters alone in their bins
#     #for i in names1prom:
#     for pro in promAll[region]:
#         marksAll.append((pro, pro, colors[1]))
        
#     # add captures in the marking step
#     for cap in catpAll[region]:
#         marksBar2.append((cap[0], cap[1], colors[2]))
        
#     # add interesting points in the marking step
#     marksBar3 = []
#     for po in interestAll[region]:
#         marksBar3.append((po, po, colors[3]))
    
#     inpath = outFile + '/%s%s_Stats/' %(flag[0], flag[1])
#     print flag
#     if len(marksBar3) != 0:
#         plotearStats(inpath, 1, marksAll, [marksBar2, marksBar3, marksAll], regionStart, regionEnd, resol, title=title, 
#                  figsize = (20, 10), diff=0, radius=50, fromOri=fromOri, fromEnd=fromEnd)
#     else:
#         plotearStats(inpath, 1, marksAll, [marksBar2, marksAll], regionStart, regionEnd, resol, title=title, 
#                      figsize = (20, 10), diff=0, radius=50, fromOri=fromOri, fromEnd=fromEnd)
#     print "#" * 80
# pdf.close()



# Plot compared stats plot
def plotearStatsCompare(values, indir, clust, marksAll, marksBar, regionStart, regionEnd, resol, rotation=90,
                 title='', titleSize = 20, figsize = (20, 10), diff=0, propor=[4], radius="?",
                 fromOri = 0, fromEnd = 0, ylim1 = (-30, 30), ylim2 = (-4, 4)):
    fontSize = 15
    # propor is the proportion of size for normal plots against grey bar info plots
    # marksAll is to add the marks in all the plots, and MarksBar is to add marks just in the
    #top bar (is a list of lists)
    xa,ya,xc,yc,xi,yi = values

    nplot = len(marksBar) + 3 #we alway have 3 plots at least
    #f, (ax3, ax1, ax2, ax5) = plt.subplots(4, sharex=True, sharey=False,figsize=figsize)
    f, allAxis = plt.subplots(nplot, sharex=True, sharey=False,figsize=figsize,
                              gridspec_kw = {'height_ratios':[1] *  ((nplot) - 3) + propor * 3})
    #gs1 = gridspec.GridSpec(figsize[0], figsize[1])
    #gs1.update(wspace=0.025, hspace=0.05)

    allAxis[-3].plot(xc, yc,'k-')
    allAxis[-3].set_title('Consistency per particle', fontsize=fontSize)
    allAxis[-3].axhline(y=np.mean(yc),c="blue")#,linewidth=0.2)
    allAxis[-3].set_ylabel("Consistency %",rotation=rotation, fontsize=fontSize)
    allAxis[-3].set_ylim(ylim1)

    allAxis[-2].plot(xa, ya,'k-')
    allAxis[-2].set_title('Accessibility per particle', fontsize=fontSize)
    allAxis[-2].axhline(y=np.mean(np.nan_to_num(ya)),c="blue")#,linewidth=0.2)
    allAxis[-2].set_ylabel("Accessibility %s" %radius,rotation=rotation,multialignment='center',
                          fontsize=fontSize)
    allAxis[-2].set_ylim(ylim1)

    allAxis[-1].plot(xi, yi,'k-')
    allAxis[-1].set_title('Contacts within sphere per particle', fontsize=fontSize)
    allAxis[-1].set_ylabel("Number of contact",rotation=rotation, multialignment='center',
                          fontsize=fontSize)
    allAxis[-1].axhline(y=np.mean(np.nan_to_num(yi)),c="blue")#,linewidth=0.2)
    allAxis[-1].set_ylim(ylim2)

    # Turn to grey the info rectangles
    for axi in allAxis[:-3]:
        rect2 = patch.Rectangle((np.min(xa),0), np.max(xa)-np.min(xa), 1, color='lightgrey')
        axi.add_patch(rect2)

    # Add marks in all the plots
    for markStart, markEnd, color in marksAll:
        # Now we add rectnagles in ohur possitions of interest
        # The variables that matter are markStart and markEnd
        high=float(round(np.nanmax(ya)-np.nanmin(ya))+1000)
        min_s=float(round(np.nanmin(ya))-100)
        rect2 = patch.Rectangle((markStart,min_s), markEnd - markStart, high, color=color,alpha=0.3)
        allAxis[-2].add_patch(rect2)

        high=float(round(np.nanmax(yc)-np.nanmin(yc))+1000)
        min_s=float(round(np.nanmin(yc))-100)
        #rect2 = patch.Rectangle((np.min(b[k][0]),0), np.max(b[k][0])-np.min(b[k][0]), 1, color='#ff7f24')
        rect2 = patch.Rectangle((markStart,min_s), markEnd - markStart, high, color=color,alpha=0.3)
        allAxis[-3].add_patch(rect2)

        high=float(round(np.nanmax(yi)-np.nanmin(yi))+1000)
        min_s=float(round(np.nanmin(yi))-100)
        #rect2 = patch.Rectangle((np.min(b[k][0]),0), np.max(b[k][0])-np.min(b[k][0]), 1, color='#ff7f24')
        rect2 = patch.Rectangle((markStart,min_s), markEnd - markStart, high, color=color,alpha=0.3)
        allAxis[-1].add_patch(rect2)

        #rect2 = patch.Rectangle((markStart,0), markEnd - markStart, 1, color=color)
        #ax3.add_patch(rect2)

    # Add marks just in the top bar
    for naxi, marks in enumerate(marksBar):
        for markStart, markEnd, color in marks:
            rect2 = patch.Rectangle((markStart,0), markEnd - markStart, 1, color=color)
            allAxis[naxi].add_patch(rect2)

        #allAxis[naxi].set_ylim([0, 8])
        allAxis[naxi].axis('off')
        #allAxis[naxi].set_aspect('equal')

    f.subplots_adjust(hspace=0.3)

    # we want reference genome positions and enough labels to see start and end
    mini = 10
    # get the number of divisions where labels would be equally spread
    start = 0
    end = ((regionEnd - regionStart) / resol) + 1
    plt.xlim(start, end)
    #start, end = allAxis[-3].get_xlim()
    #print start, end
    binrange = range(int(start), int(end) + 1)
    posrange = range(regionStart, regionEnd + resol, resol)
    #for i in range(8, 22):
    #    divisor1 = len(binrange) / float(i)
    #    if mini >= divisor1/int(divisor1):
    #        mini = min(mini, divisor1/int(divisor1))
    #        divisor = int(divisor1)
    divisor = 9
    # use range to select the values in bin (value used at the time to plo) and
    #genomic position (value that we actually want to see)

    binx = [binrange[i] for i in range(0, len(binrange), divisor)]
    posx = [posrange[i] for i in range(0, len(posrange), divisor)]

    # check if last position is too close
    #if binx[-1] - binrange[-1] < divisor/1.5:
    #    binx = binx[:-1] + [binrange[-1]]
    #    posx = posx[:-1] + [posrange[-1]]
    #else:
    #    binx = binx + [binrange[-1]]
    #    posx = posx + [posrange[-1]]

    plt.xticks(binx,
               posx,
               rotation=45)

    if title != '':
        f.suptitle(title, size=titleSize)

    # We add std lines
    for n, values in enumerate([yi, ya, yc]):
        values = np.nan_to_num(values)
        # get SD
        std = np.std(values)
        rect2 = patch.Rectangle((0, np.mean(values) - std), # x pos to start and y pos to start
                                len(values),# x pos to end
                                (np.mean(values) + std) - (np.mean(values) - std),  # y length till end
                                color='grey',alpha=0.3)
        allAxis[-(n+1)].add_patch(rect2)
        #print n+1, std

    plt.xlim(np.min(xa) + fromOri, np.max(xa) - fromEnd)

    plt.show()
    return f

# Example of call
# pdf = matplotlib.backends.backend_pdf.PdfPages(outpath + 'statsMultiplotsCellComparison.pdf')
# modFilesShort =modFiles
# for i in range(0, len(modFilesShort), 2):
#     allValues = {}
#     fi1 = modFilesShort[i]
#     flag1 = fi1.split('/')[-1].split('_')[:2]

#     colors = ["red", "blue", "black", 'green', 'yellow']
#     inpath2 = outFile + '/%s%s_Stats/' %(flag1[0], flag1[1])

#     clust = 1
#     diff=0
#     xa0,ya0,eMa0,ePa0,da0=GetXYA(inpath2 + 'accessibility.cluster%s.dat' %clust,diff=diff)
#     xc0,yc0=GetXY(inpath2 + 'consistency.cluster%s.dat' %clust,diff=diff)
#     xi0,yi0,eMi0,ePi0,di0=GetXY_smooth4(inpath2 + 'interactionscluster%s.dat' %clust,diff=diff)



#     # Beggin to plot
#     fi2 = modFilesShort[i + 1]
#     flag2 = fi2.split('/')[-1].split('_')[:2]
#     allValues[flag2[1]] = {}
#     # Get region index
#     region = flag2[1]
#     # get regions start and end position (and chromosome)
#     reg = regiones[region]
#     chrom = reg[0]
#     regionStart = reg[1]
#     regionEnd = reg[2]

#     fromOri = (originalRegionPos[region][1]/resol) - (reg[1] / resol)
#     fromEnd = (reg[2] / resol) - (originalRegionPos[region][2] / resol)

#     marksAll = []
#     marksBar2 = []
#     # esto todavia no funcion, tengo q ver como poner el titulo en toda la figura con subplots
#     title = 'Model statistics in %s, subtracting TB statistics to PI ' %(flag1[1])

#     # Add enhancers in the marking step
#     for en in enhAll[region]:
#         marksAll.append((en, en, colors[0]))

#     # If we want just promoters alone in their bins
#     #for i in names1prom:
#     for pro in promAll[region]:
#         marksAll.append((pro, pro, colors[1]))


#     # add captures in the marking step
#     for cap in catpAll[region]:
#         marksBar2.append((cap[0], cap[1], colors[2]))

#     # add interesting points in the marking step
#     marksBar3 = []
#     for po in interestAll[region]:
#         marksBar3.append((po, po, colors[3]))

#     inpath3 = outFile + '/%s%s_Stats/' %(flag2[0], flag2[1])

#     # Get values and remove to non CRISPRed values
#     xa,ya,eMa,ePa,da=GetXYA(inpath3 + 'accessibility.cluster%s.dat' %clust,diff=diff)
#     xc,yc=GetXY(inpath3 + 'consistency.cluster%s.dat' %clust,diff=diff)
#     xi,yi,eMi,ePi,di=GetXY_smooth4(inpath3 + 'interactionscluster%s.dat' %clust,diff=diff)

#     print len(ya)
#     print len(ya0)

#     if len(ya) == len(ya0):
#         ya = [ya0[ii] - ya[ii] for ii in range(len(ya))]
#         yc = [yc0[ii] - yc[ii] for ii in range(len(yc))]
#         yi = [yi0[ii] - yi[ii] for ii in range(len(yi))]
#     else:
#         ya = 0
#         yc = 0
#         yi = 0



#     allValues[region]['acces'] = ya
#     allValues[region]['consi'] = yc
#     allValues[region]['inter'] = yi

#     print flag1
#     plotearStatsCompare([xa, ya, xc, yc, xi, yi], outpath, 1, marksAll, [marksBar2, marksBar3, marksAll],
#                         regionStart, regionEnd, resol, title=title, figsize = (20, 10), diff=0, radius=50,
#                         fromOri=fromOri, fromEnd=fromEnd, ylim1 = (-40, 40), ylim2 = (-6, 6))

#     print "#" * 80
# pdf.close()


# Stats plot comparing two model objects
def convertRGB(value):
    return [i/255.0 for i in value]
def plotearStatsBoth(indirs, clust, marksAll, marksBar, regionStart, regionEnd, resol, rotation=90,
                 title='', titleSize = 20, figsize = (20, 10), diff=0, propor=[4], radius="?",
                 fromOri = 0, fromEnd = 0, colors=False):
    fontSize = 15
    # propor is the proportion of size for normal plots against grey bar info plots
    # marksAll is to add the marks in all the plots, and MarksBar is to add marks just in the
    #top bar (is a list of lists)
    xas = []
    yas = []
    eMas = []
    ePas = []
    das = []
    xcs = []
    ycs = []
    xis = []
    yis = []
    eMis = []
    ePis = []
    dis = []
    for ni, indi in enumerate(indirs):
        xa,ya,eMa,ePa,da=GetXYA(indirs[ni] + 'accessibility.cluster%s.dat' %clust,diff=diff)
        xc,yc=GetXY(indirs[ni] + 'consistency.cluster%s.dat' %clust,diff=diff)
        xi,yi,eMi,ePi,di=GetXY_smooth4(indirs[ni] + 'interactionscluster%s.dat' %clust,diff=diff)
        xas += [copy.copy(xa)]
        yas += [copy.copy(ya)]
        eMas += [copy.copy(eMa)]
        ePas += [copy.copy(ePa)]
        das += [copy.copy(da)]
        xcs += [copy.copy(xc)]
        ycs += [copy.copy(yc)]
        xis += [copy.copy(xi)]
        yis += [copy.copy(yi)]
        eMis += [copy.copy(eMi)]
        ePis += [copy.copy(ePi)]
        dis += [copy.copy(di)]

    nplot = len(marksBar) + 3 #we alway have 3 plots at least
    #f, (ax3, ax1, ax2, ax5) = plt.subplots(4, sharex=True, sharey=False,figsize=figsize)
    f, allAxis = plt.subplots(nplot, sharex=True, sharey=False,figsize=figsize,
                              gridspec_kw = {'height_ratios':[1] *  ((nplot) - 3) + propor * 3})
    # change tic size
    #f.tick_params(axis='both', which='major', labelsize=10)
    #f.tick_params(axis='both', which='minor', labelsize=8)
    #gs1 = gridspec.GridSpec(figsize[0], figsize[1])
    #gs1.update(wspace=0.025, hspace=0.05)
    if colors == False:
        colors = [(120,94,240), (220,38,127), (254,97,0)]
    colors2 = []
    for c in colors:
        colors2.append(convertRGB(c))
    for ni, indi in enumerate(indirs):
        allAxis[-3].plot(xcs[ni], ycs[ni],'-', color = colors2[ni], lw=1.75)
    allAxis[-3].set_title('Consistency per particle', fontsize=fontSize)
    #allAxis[-3].axhline(y=np.mean(yc),c="blue")#,linewidth=0.2)
    allAxis[-3].set_ylabel("Consistency %",rotation=rotation, fontsize=fontSize)
    allAxis[-3].set_ylim(0, 100)

    for ni, indi in enumerate(indirs):
        allAxis[-2].plot(xas[ni], yas[ni],'-', color = colors2[ni], lw=1.75)
    allAxis[-2].set_title('Accessibility per particle', fontsize=fontSize)
    #allAxis[-2].axhline(y=np.mean(np.nan_to_num(ya)),c="blue")#,linewidth=0.2)
    allAxis[-2].set_ylabel("Accessibility" ,rotation=rotation,multialignment='center',
                          fontsize=fontSize)
    #allAxis[-2].fill_between(xa,eMa,ePa,alpha=0.5, edgecolor='#bebebe', facecolor='#bebebe' )
    allAxis[-2].set_ylim(0, 100)

    for ni, indi in enumerate(indirs):
        allAxis[-1].plot(xis[ni], yis[ni],'-', color = colors2[ni], lw=1.75)
    allAxis[-1].set_title('Contacts within sphere per particle', fontsize=fontSize)
    allAxis[-1].set_ylabel("Number of contact",rotation=rotation, multialignment='center',
                          fontsize=fontSize)
    #allAxis[-1].fill_between(xi,eMi,ePi,alpha=0.5, edgecolor='#bebebe', facecolor='#bebebe' )
    #allAxis[-1].axhline(y=np.mean(np.nan_to_num(yi)),c="blue")#,linewidth=0.2)

    # Turn to grey the info rectangles
    for axi in allAxis[:-3]:
        rect2 = patch.Rectangle((np.min(xa),0), np.max(xa)-np.min(xa), 1, color='lightgrey')
        axi.add_patch(rect2)

    # Add marks in all the plots
    for markStart, markEnd, color in marksAll:
        # Now we add rectnagles in ohur possitions of interest
        # The variables that matter are markStart and markEnd
        high=float(round(np.nanmax(ya)-np.nanmin(ya))+1000)
        min_s=float(round(np.nanmin(ya))-100)
        rect2 = patch.Rectangle((markStart,min_s), markEnd - markStart, high, color=color,alpha=0.3)
        allAxis[-2].add_patch(rect2)

        high=float(round(np.nanmax(yc)-np.nanmin(yc))+1000)
        min_s=float(round(np.nanmin(yc))-100)
        #rect2 = patch.Rectangle((np.min(b[k][0]),0), np.max(b[k][0])-np.min(b[k][0]), 1, color='#ff7f24')
        rect2 = patch.Rectangle((markStart,min_s), markEnd - markStart, high, color=color,alpha=0.3)
        allAxis[-3].add_patch(rect2)

        high=float(round(np.nanmax(yi)-np.nanmin(yi))+1000)
        min_s=float(round(np.nanmin(yi))-100)
        #rect2 = patch.Rectangle((np.min(b[k][0]),0), np.max(b[k][0])-np.min(b[k][0]), 1, color='#ff7f24')
        rect2 = patch.Rectangle((markStart,min_s), markEnd - markStart, high, color=color,alpha=0.3)
        allAxis[-1].add_patch(rect2)

        #rect2 = patch.Rectangle((markStart,0), markEnd - markStart, 1, color=color)
        #ax3.add_patch(rect2)

    # Add marks just in the top bar
    for naxi, marks in enumerate(marksBar):
        for markStart, markEnd, color in marks:
            rect2 = patch.Rectangle((markStart,0), markEnd - markStart, 1, color=color)
            allAxis[naxi].add_patch(rect2)

        #allAxis[naxi].set_ylim([0, 8])
        allAxis[naxi].axis('off')
        #allAxis[naxi].set_aspect('equal')

    f.subplots_adjust(hspace=0.3)

    # we want reference genome positions and enough labels to see start and end
    mini = 10
    # get the number of divisions where labels would be equally spread
    start = 0
    end = ((regionEnd - regionStart) / resol) + 1
    plt.xlim(start, end)
    #start, end = allAxis[-3].get_xlim()
    #print start, end
    binrange = range(int(start), int(end) + 1)
    posrange = range(regionStart, regionEnd + resol, resol)
    #for i in range(8, 22):
    #    divisor1 = len(binrange) / float(i)
    #    if mini >= divisor1/int(divisor1):
    #        mini = min(mini, divisor1/int(divisor1))
    #        divisor = int(divisor1)
    divisor = 9
    # use range to select the values in bin (value used at the time to plo) and
    #genomic position (value that we actually want to see)

    binx = [binrange[i] for i in range(0, len(binrange), divisor)]
    posx = ['{:,}'.format(posrange[i]) for i in range(0, len(posrange), divisor)]

    # check if last position is too close
    #if binx[-1] - binrange[-1] < divisor/1.5:
    #    binx = binx[:-1] + [binrange[-1]]
    #    posx = posx[:-1] + [posrange[-1]]
    #else:
    #    binx = binx + [binrange[-1]]
    #    posx = posx + [posrange[-1]]

    plt.xticks(binx,
               posx,
               rotation=45)
    plt.xlabel('Genomic coordinates', fontsize=fontSize)

    if title != '':
        f.suptitle(title, size=titleSize)

    plt.xlim(np.min(xa) + fromOri, np.max(xa) - fromEnd)

    plt.show()
    return f


# Example of run
# pdf = matplotlib.backends.backend_pdf.PdfPages(outpath + 'statsMultiplotsBothCell_reg18.pdf')
# modFilesShort = modFiles
# for i in range(0, len(modFilesShort), 2):
#     colors = ["red", "blue", "black", 'green', 'yellow']
#     # Beggin to plot
#     fi1 = modFilesShort[i]
#     flag1 = fi1.split('/')[-1].split('_')[:2]

#     fi2 = modFilesShort[i + 1]
#     flag2 = fi2.split('/')[-1].split('_')[:2]
#     # Get region index
#     region = flag2[1]
#     # get regions start and end position (and chromosome)
#     reg = regiones[region]
#     chrom = reg[0]
#     regionStart = reg[1]
#     regionEnd = reg[2]

#     fromOri = (originalRegionPos[region][1]/resol) - (reg[1] / resol)
#     fromEnd = (reg[2] / resol) - (originalRegionPos[region][2] / resol)

#     marksAll = []
#     marksBar2 = []
#     # esto todavia no funcion, tengo q ver como poner el titulo en toda la figura con subplots
#     title = 'Model statistics in %s, subtracting TB statistics to PI ' %(flag1[1])

#     # Add enhancers in the marking step
#     for en in enhAll[region]:
#         marksAll.append((en, en, colors[0]))

#     # If we want just promoters alone in their bins
#     #for i in names1prom:
#     for pro in promAll[region]:
#         marksAll.append((pro, pro, colors[1]))


#     # add captures in the marking step
#     for cap in catpAll[region]:
#         marksBar2.append((cap[0], cap[1], colors[2]))

#     # add interesting points in the marking step
#     marksBar3 = []
#     for po in interestAll[region]:
#         marksBar3.append((po, po, colors[3]))


#     inpath2 = outFile + '/%s%s_Stats/' %(flag1[0], flag1[1])
#     inpath3 = outFile + '/%s%s_Stats/' %(flag2[0], flag2[1])


#     print flag1

#     plotearStatsBoth([inpath2, inpath3], 1, [], [marksAll], regionStart, regionEnd, resol, title=title,
#              figsize = (20, 10), diff=0, radius=50, fromOri=fromOri, fromEnd=fromEnd)


#     print "#" * 80
# pdf.close()


## Function to plot differential contact matrices with marqued points of interest
def plotDiffMtrx(mtComp, colorate=[], arrows=[], title='',
                 vRange = [-0.3, 0.3],
                 whiteLim=0.3,
                figsize=(10,10), titleAdj=1):

    '''
    Function to plot an interaction matrix (normaly used for differential ones)
        while showing some positions of interest
    :param mtComp: Interaction matrix in format of list of lists
    :param [] colorate: List with nested lists with position indexes to color the side bars.
        Coloring orther is blue, orange, green, red, purple, brown, pink, grey, greenish,
        lightblue, +white
    :param [] arrows: List with position indexes to add arrows in the sides of the matrix
    :param '' title: String with plot title
    :param [-0.3, 0.3] vRange: minimum and maximum values for the main plot colorbar
    :param 0.3 whiteLim: value till which we color in white the matrix values
    '''
    ## Prepare color map for side bars
    barCmap = plt.cm.get_cmap('tab10')
    colors = barCmap(np.arange(barCmap.N))
    colors = colors.tolist()
    colors.append([1, 1, 1, 1])
    barCmap = LinearSegmentedColormap.from_list('mycmap', colors)
    whitePos = len(colors)
    # Order od tab10 colors is:
    # blue, orange, green, red, purple, brown, pink, grey, greenish,
    #lightblue, +white

    ## Prepare colorbar for matrix
    cmap=plt.get_cmap('bwr')
    # get positions to turn white
    posi = int(round(cmap.N * whiteLim / 2))
    # get color change position
    midPos = cmap.N / 2
    # obtain color list
    ccolors = cmap(np.arange(cmap.N))
    ccolors = ccolors.tolist()

    # make range to select colors from beggining
    rangebeg = midPos / float(midPos - posi)
    rangebeg = [int(a) for a in np.arange(0, midPos, rangebeg)]

    # make range to select colors from end
    rangeend = midPos / float(midPos - posi)
    rangeend = [int(a) for a in np.arange(midPos, cmap.N, rangeend)]

    # join with whites in the middle
    colors =  [ccolors[c] for c in rangebeg] + [[1, 1, 1, 1]] * posi + [[1, 1, 1, 1]] * posi + \
    [ccolors[c] for c in rangeend]

    # create colormap
    cmap = LinearSegmentedColormap.from_list('mycmap', colors)

    ## Prepare sidebars content
    if len(colorate) <= len(colors) and len(colorate) != 0:
        # now we create a vector of mtComp size and add values to color
        binVector = [0] * len(mtComp)
        for b in range(len(binVector)):
            # go through each list in colorate in order
            added = False
            for nbb, bb in enumerate(colorate):
                if b in bb:
                    binVector[b] = nbb
                    added=True
            # if no mark in there leave white
            if added == False:
                binVector[b] = whitePos

    else:
        print 'To many things to colorate, not enough colors'

    # Add arrows
    arrowVect = [float('nan')] * len(mtComp)
    for a in arrows:
        arrowVect[a] = 1

    # create space for all plots
    fig, ax = plt.subplots(1,1, figsize=figsize)


    ### Main matrix plot
    mat = ax.imshow(mtComp, interpolation='nearest', origin='lower', vmin= vRange[0], vmax=vRange[1],
                    cmap=cmap, )


    ### Down right (colorbar)
    divider = make_axes_locatable(ax)
    cax1 = divider.append_axes("right", size="5%", pad=0.05)
    plt.colorbar(mat, cax=cax1)

    ### Color marks in the upper part
    barwide = max(2, len(mtComp) * 0.02)
    if len(colorate) != 0:
        ax2 = divider.append_axes("top", size="5%", pad=0.05)
        # turn it into a matrix
        mtBinVect = []
        for i in range(len(binVector)):
            mtBinVect.append(binVector)
        ax2.imshow(mtBinVect, interpolation='nearest', origin='lower',
                   extent=[0,len(binVector),0,barwide], cmap=barCmap,
                  vmin=0, vmax=whitePos)
        ax2.get_yaxis().set_visible(False)
        ax2.get_xaxis().set_visible(False)
        ax2.set_xlim(0, len(binVector))

    ### Arrows in the upper part
    ax3 = divider.append_axes("top", size="5%", pad=0.05)
    ax3.get_yaxis().set_visible(False)
    ax3.get_xaxis().set_visible(False)
    arrow = u'$\u2193$'
    ax3.plot(range(len(binVector)), arrowVect, linestyle='none', marker=arrow, markersize=20)
    ax3.set_xlim(0, len(binVector))

    ax3.axis('off')


    ### Color marks in the left part
    if len(colorate) != 0:
        ax4 = divider.append_axes("left", size="5%", pad=0.05)
        ax4.set_xticklabels([])
        ax4.imshow(np.matrix.transpose(np.array(mtBinVect)), interpolation='nearest',
                   origin='lower', extent=[0,barwide,0, len(binVector)], cmap=barCmap,
                  vmin=0, vmax=whitePos)
        ax4.get_xaxis().set_visible(False)
        ax4.set_ylim(0, len(binVector))

        # If there is color marks we hide y tics from main plot
        ax.set_yticklabels([])


    ### Arrows in the left part
    ax5 = divider.append_axes("left", size="5%", pad=0.05)
    ax5.set_xlim(0, len(binVector))
    #ax6.get_yaxis().set_visible(False)
    ax5.get_xaxis().set_visible(False)
    arrow = u'$\u2192$'
    ax5.plot(arrowVect, range(len(binVector)), linestyle='none', marker=arrow, markersize=20)
    ax5.set_ylim(0, len(binVector))
    ax5.set_xlim(-5, 8)
    ax5.axis('off')



    fig.suptitle(title,size=18)
    fig.subplots_adjust(top=titleAdj)

    plt.show()

    return fig



def aracnoPlot(groupDegree, edgeList, focus, percentaje, tags=[], xAxis=False, saveFig=False, 
               outPath=False, filterDgrRatio=False, pdf=False):
    '''
    Function to plot neigbourh interactions surrounding a bin of interest
    :param groupDegree: Dictionary with the cumulated degree of all the members in each group
        first key for group ID, and inside 'degree' for degree and 'members' for group members
    :param edgeList: list with all pairwise interactions
    :param focus: focus bin or bins range (last unit wont be used)
    :param percentaje: percentyle value above which we count the interactions (for each 
        interactions from the viewpoint bin towards the others)
    :param [] tags: tags fro title and file. First one is chromosome number and next one whatever 
        you want
    :param False xAxis: range to change deafult xAxis limits
    :param False saveFig: True if you want to store the plot in a PDF
    :param False outPath: False if you have already open a matplotlib.backends.backend_pdf.PdfPages
        instance, or a real PATH if you want the function to open it and close for you
    :param False filterDgrRatio: If integer, it sets the minimum degree ratio 
        (number of bins in group/cumulated degree) to take into account a group
    
    '''
    if len(tags) > 1:
        chrom = tags[0]
        tag = tags[1]
    # Open PDf if asked
    if saveFig == True:
        if outPath != False:
            if filterDgrRatio == False:
                pdf = matplotlib.backends.backend_pdf.PdfPages(outPath + \
                                                               'AracnoPlot_%s_chr%s_%s-%s.pdf' %(tag,
                                                                                                chrom, 
                                                                                                focus[0], 
                                                                                                focus[1]))
            else:
                pdf = matplotlib.backends.backend_pdf.PdfPages(outPath + \
                                                               'AracnoPlot_%s_chr%s_%s-%s_Flt%s.pdf' %(tag,
                                                                                                chrom, 
                                                                                                focus[0], 
                                                                                                focus[1],
                                                                                                filterDgrRatio))

    ## Plot
    fg, ax = plt.subplots(1, 1, figsize=(20, 8))
    plotHigh = 100
    stop = False  # stop if after filtering there is no data

    # If we are not going to filter
    if filterDgrRatio == False:
        # First add grey ellipses
        for k in edgeList:
            x1 = k[0]
            x2 = k[1]
            #thick = highPeaks[e[0]][e[1]] / emaxi
            if not (x1 in focus):
                pac = mpatches.Ellipse([(x2+x1)/2.0, 0], x2-x1, plotHigh*2, angle=0, fill=False, color='black', alpha=0.8)
                ax.add_patch(pac)

        # Then red ones
        for k in edgeList:
            x1 = k[0]
            x2 = k[1]
            #thick = highPeaks[e[0]][e[1]] / emaxi
            if (x1 in focus):
                pac = mpatches.Ellipse([(x2+x1)/2.0, 0], x2-x1, plotHigh*2, angle=0, fill=False, color='red')
            ax.add_patch(pac)

    # If we asked for a filtering
    else:
        # Get nodes to filter
        removed = []
        for gr in sorted(groupDegree.keys()):
            if groupDegree[gr]['degree']/float(len(groupDegree[gr]['members'])) < filterDgrRatio:
                removed += groupDegree[gr]['members']
                del groupDegree[gr]

        # get edges to filter
        todel = []
        for ned, ed in enumerate(edgeList):
            if ed[0] in removed or ed[1] in removed:
                todel.append(ned)    

        # Filter edges
        for i in range(len(todel) - 1, -1, -1):
            del edgeList[todel[i]]
               
        # if everything was filtered
        if len(edgeList) == 0:
            stop = True
            
        else:
            for k in edgeList:
                x1 = k[0]
                x2 = k[1]


                # Apply the filtering
                #thick = highPeaks[e[0]][e[1]] / emaxi
                if not (x1 in focus):
                    pac = mpatches.Ellipse([(x2+x1)/2.0, 0], x2-x1, plotHigh*2, 
                                           angle=0, fill=False, color='black', alpha=0.8)
                    ax.add_patch(pac)

            # Then red ones
            for k in edgeList:
                x1 = k[0]
                x2 = k[1]
                # Apply the filtering
                #thick = highPeaks[e[0]][e[1]] / emaxi
                if (x1 in focus):
                    pac = mpatches.Ellipse([(x2+x1)/2.0, 0], x2-x1, plotHigh*2, 
                                           angle=0, fill=False, color='red')
                ax.add_patch(pac)

    if stop == False:
        # get value ranges
        if xAxis == False:
            mini = min(min(ed) for ed in edgeList)
            maxi = max(max(ed) for ed in edgeList)
        else:
            mini = xAxis[0]
            maxi = xAxis[1]
        # set plot limits
        plt.xlim(mini, maxi)
        plt.ylim(0, plotHigh + plotHigh*0.1)

        # remove all axis data but down
        plt.tick_params(axis='both', left='off', top='off', right='off', bottom='on', labelleft='off', 
                        labeltop='off', labelright='off', labelbottom='on')
        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)
        ax.spines['left'].set_visible(False)

        # Will prepare degree data to be added in the plot
        # Find the mean
        for gd in groupDegree:
            # dont add the filtered groups
            mean = np.mean(groupDegree[gd]['members'])
            # Add text
            ypos = plotHigh + plotHigh*0.05
            plt.text(mean, ypos, groupDegree[gd]['degree'], horizontalalignment='center',
                    size=20)
            # Add lines delimitating each group
            # get smallest member
            minMemb = min(groupDegree[gd]['members'])
            maxMemb = max(groupDegree[gd]['members'])

            # Select altitude of the bracket
            up = ypos - plotHigh*0.05
            h=plotHigh*0.05
            # plot it
            plt.plot([minMemb , minMemb, maxMemb, maxMemb], [up, up+h, up+h, up], 
                                     lw=1.5, color='grey')

        # Dont allow exponential notation
        ax.get_xaxis().get_major_formatter().set_useOffset(False)
    
    else:
        plt.text(0.5, 0.5, 'Everything was filtered', horizontalalignment='center',
                        size=20)
        plt.xlim(0,1)
        plt.ylim(0,1)
            
            
    # Title and axis data
    if filterDgrRatio == False:
        plt.title('Filtered interactions of bins interacting with our focus bin/bins \
(thress %s) in %s\nred=interaction between bin and focus not filtered, \
grey=interction filtered' %(percentaje, tag))
    else:
        plt.title('Filtered interactions of bins interacting with our focus bin/bins \
(thress %s, flt %s) in %s\nred=interaction between bin and focus not filtered, \
grey=interction filtered' %(percentaje, filterDgrRatio, tag))
    plt.xlabel('Bin number in chromosome %s' %chrom)
    plt.show()

    # Store file
    if saveFig == True:
        pdf.savefig(fg , bbox_inches='tight')
        if outPath != False:
            pdf.close()
            

def getMatrixOrder(histdict, method='ward', metric='euclidean'):
    marks = sorted(histdict.keys())
    # get data in matrix format
    matrix = [[histdict[mark][k] for mark in marks ]for k in sorted(histdict[marks[0]].keys())]
    # normalize all columns by zscore values
    matrix = zscore(matrix, axis=0)
    #matrix = zip(*[[(v/float(sum(m))) if sum(m) else 0 for v in m] for m in zip(*matrix)])
    #matrix = np.array(matrix)
    
    # transpose matrix to get column clusters
    matrix2 = np.matrix.transpose(np.array(matrix))


    # Cluster it
    Y = sch.linkage(matrix2, method=method, metric=metric)
    Z2 = sch.dendrogram(Y, no_plot=True)
    markOrder = Z2['leaves']
    
    # return clustering order
    return markOrder

def squarePlot(histdict, region, title = '', prevOrder='', saveFig = True, minMax = "", zscore = False):
    # Generate some variables we will need
    marks = sorted(histdict.keys())
    nmark = len(marks)
    
    # get title
    if title == '':
        title = 'Marker intensities'
        
    # get data in matrix format
    matrix = [[histdict[mark][k] for mark in marks ]for k in sorted(histdict[marks[0]].keys())]
    if zscore == True:
        # normalize all columns by zscore values
        matrix = zscore(matrix, axis=0)
    #matrix = zip(*[[(v/float(sum(m))) if sum(m) else 0 for v in m] for m in zip(*matrix)])
    else:
        matrix = np.array(matrix)
    
    if prevOrder == '':
        markOrder = range(nmark)
        
    else:
        markOrder = prevOrder
    matrix = matrix[:,markOrder]

    # set colorbar limits
    mini = 100
    maxi = 0
    for i in matrix:
        for ii in i:
            mini = min(mini, ii)
            maxi = max(maxi, ii)
    colorLim = max(abs(mini), abs(maxi))
    colorLim = [-colorLim, colorLim]
    # check if there are given ones and are bigger
    if minMax != '':
        if (minMax[0] > colorLim[0]) or (minMax[1] < colorLim[1]):
            print 'There are values smaller/bigger than the provided range'
            print 'Default range will be used'
        
        # Just using the given values limit for asiggning color
        else:
            colorLim = minMax

    

    fig = plt.figure(figsize=(20, 10))
    plt.imshow(matrix, 
               cmap = "RdBu_r", aspect='auto', interpolation='none', origin='lower', 
               vmin=colorLim[0], vmax=colorLim[1])
    plt.xticks(range(nmark), [marks[i] for i in markOrder])
    plt.colorbar()

    
    # Add y labels tick
    ylabels = [int(i) for i in sorted(histdict[marks[0]].keys())]
    plt.yticks(range(len(ylabels)), ylabels)

    # Add labels
    plt.ylabel("Distance from center (nm)")
    plt.xlabel("ChIP mark Zscore")
    plt.title('%s to perifery in %s\n' %(title, region), size=20)
    
    
    #plt.colorbar(boundaries=np.linspace(-colorLim,colorLim,20))
    if saveFig == True:
        pdf.savefig( fig, bbox_inches='tight')
    plt.show()
    

    
# Newer version
def radialPlot(histdict, region, valRange, markOrder='', timesEach=7, nylabels=10, 
               title='', saveFig=True, minMax = "", oneCheese=False, color='RdBu_r',
               unEqualRadi = False, colorTrend='bimodal', fixedMinMax=False):
    

    '''
    Plot radial matrix with given data from selected point of the model outwards it
    :param histdict: Dictionary with two leves. First the data to separate into plot portions,
        second the value asociated to each radius in the plot
    :param nmark: Number of keys we have in the first level of histdict
    :param bins: Number of distance chunks we have (number of keys in the second level of histdict)
    :param valRange: Value ranges analyzed from start point to end point
    :param fixedMinMax: Wether the function can look for a range that contains or values (False), or 
        the range given in minMax cannot be changed (True)
    '''

    # Generate some variables we will need
    if oneCheese == False:
        marks = sorted(histdict.keys())
        nmark = len(marks)
        # Need to add one value more to bin so we see all data
        bins = len(histdict[marks[0]].keys()) + 1
    else:
        marks = ['']
        nmark = 1
        # Need to add one value more to bin so we see all data
        bins = len(histdict.keys()) + 1



    # get title
    if title == '':
        title = 'Marker intensities'
    # Generate some data...
    # Note that all of these are _2D_ arrays, so that we can use meshgrid
    # You'll need to "grid" your data to use pcolormesh if it's un-ordered points
    # theta es el angulo y r es la distancia en 'radio' a la que esta cada sub circunferencia
    portions = nmark * timesEach + 1
    theta, r = np.mgrid[0:2*np.pi:portions+0j, 0:1:bins+0j]
    z = np.random.random(theta.size).reshape(theta.shape)
    # Por algun motivo la ultima lista que se mete en pcolormesh no aparece, asi que hay q darle indices +1
    #aunque no existan
    # Lo mismo pasa cuando mira cuantos radios tiene cada quesito, siempre tendra en cuenta como que hay uno mas
    zs = []

    # If we are dealing with radius with different sizes
    if unEqualRadi == True:
        r = []
        for i in range(portions):
            r.append([v / valRange[-1] for v in valRange])
        r = np.array(r)

    # get new marks order
    if oneCheese == False:
        if markOrder != '':
            newMarks = [marks[i] for i in markOrder]
        else:
            newMarks = [marks[i] for i in range(nmark)]
        for n in range(nmark):
            values = [histdict[newMarks[n]][i] for i in sorted(histdict[newMarks[n]].keys())]
            for i in range(timesEach):
                zs.append(values)
    else:
        newMarks = ''
        values = [histdict[i] for i in sorted(histdict.keys())]
        for i in range(timesEach):
            zs.append(values)
    





    fig, ax2 = plt.subplots(ncols=1, subplot_kw=dict(projection='polar'), figsize=(10,10))
    #fig = plt.figure(figsize=(10,10))
    #plt.polar()
    
    infinites = []
    # Take color limit 
    if fixedMinMax == False:
        for n in range(0, (nmark * timesEach), timesEach):
            if sum(zs[n]) != 0:
                ## set colorbar limits
                # keep track of infinites
                for iz, z in enumerate(zs[n]):
                    if z ==  float('Inf') or z ==  float('-Inf'):
                        infinites.append([n, iz])
                # Find extrem values for color ranges
                mini = min([z for z in zs[n] if (z !=  float('Inf') and z !=  float('-Inf'))])
                maxi = max([z for z in zs[n] if (z !=  float('Inf') and z !=  float('-Inf'))])
                colorLim = max(abs(mini), abs(maxi))
                colorLim = [-colorLim, colorLim]
                # check if there are given ones and are bigger
                if minMax != '':
                    if (minMax[0] > colorLim[0]) or (minMax[1] < colorLim[1]):
                        print 'There are values smaller/bigger than the provided range'
                        print 'Default range will be used'
                        minMax = colorLim
                        #print colorLim
                    # Just using the given values limit for asiggning color
                    else:
                        colorLim = minMax
            else:
                colorLim = minMax
    else:
        colorLim = minMax
    if colorTrend == 'unimodal':
        colorLim[0] = 0     
    for n in range(0, (nmark * timesEach), timesEach):
        #if sum(zs[n]) != 0:
        ## Plot cheese
        plt.pcolormesh(theta[n:n+timesEach+1], r[n:n+timesEach+1], zs[n:n+timesEach], 
                       cmap=color, vmin=colorLim[0], vmax=colorLim[1], edgecolors='face')
                       #cmap=colors[n/timesEach], )
        ## Add vertical line separating portions
        plt.axvline(x=theta[n][0], color='black', alpha=0.3, linestyle='--')

    # If there is an infinite value we add an asterisk
    #print infinites
    #for infi in infinites:
    #    plt.plot((theta[infi[0]][0] + theta[infi[0] + timesEach][0])/2, 
    #             (r[0][infi[1]] + r[0][infi[1] + 1])/2, 
    #             '*', c = 'white', markersize = 20)
    #plt.plot(0.5178449428994164, 0.25, '*', c='white', markersize = 20)
    ## Add labels
    # get position for x labels
    angles = [i[0] for i in theta]
    # we should put the label more or less in the middle portion
    labpos = timesEach / 2
    angles = [ angles[n + labpos] for n in range(0, len(angles) - 1, timesEach)]
    # Add x labels  
    plt.xticks(angles, newMarks)


    #ax2.set_ylim([0, 1])

    ## Add y tick values
    # we remove the starting point from valRAnge
    #valRange = valRange[1:]
    # we need to create a range from 0 to one with same numbers as our bins
    if unEqualRadi == True:
        binrange = r[0][1:]
        # The plot will for sure have too many divisions
        fibonacci_numbers = [0, 1]
        for i in range(2,len(binrange) + 3):
            fibonacci_numbers.append(fibonacci_numbers[i-1]+fibonacci_numbers[i-2])
        # sacrilegious change
        fibonacci_numbers[3] = 0
        fibonacci_numbers[4] = 2
        binrangePos = [i for i in fibonacci_numbers[3:] if i <= len(binrange)]
        valRangeS = [int(valRange[1:][i]) if i in binrangePos else '' for i in range(0, len(valRange[1:]))]
        plt.yticks(binrange, valRangeS)
        alpha=0.3
    else:
        nytics = bins
        binrange = np.linspace(0,1,nytics,endpoint=True)[1:]
        # If the plot has to many divisions we need to show just a few axis
        steps = (bins / nylabels) + 1
        if bins > nylabels :
            binrangePos = [i for i in range(0, len(binrange), steps)]
            #binrangeS = [binrange[i] if i in binrangePos else '' for i in range(0, len(binrange))]
            valRangeS = [int(valRange[1:][i]) if i in binrangePos else '' for i in range(0, len(valRange[1:]))]
            plt.yticks(binrange, valRangeS)
        # If there are just a few divisions we just use the range values to delimitate y axis
        else:
            plt.yticks(binrange, valRange[1:])
        # transparency for radius lines
        alpha=0.3



    # Add LINES if we dont have to many divisions
    if len(valRange) < 100:
        plt.grid(axis='y', linestyle='--', linewidth=1, alpha=alpha)

    # Add colorbar
    if colorTrend == 'bimodal':
        cbar = plt.colorbar(ticks=[colorLim[0], 0, colorLim[1]], orientation='vertical', fraction=0.040, pad=0.15)
    else:
        cbar = plt.colorbar(ticks=[colorLim[0], colorLim[1]], orientation='vertical', fraction=0.040, pad=0.15)
    #cbar.ax.set_yticklabels(['Low', 'Average', 'High'])  # horizontal colorbar

    
    # title etc
    plt.title('%s to perifery in %s\n\n' %(title, region), size=20)

    # move positions of y labels
    if oneCheese == True:
        ax2.set_rlabel_position(90)

    
    if saveFig == True:
        pdf.savefig( fig, bbox_inches='tight')
    plt.show()

    return colorLim, infinites
    
# Function to know if we are dealing with a NaN
# works because NaN isn't equal to anything, even itself
def isNaN(num):
    return num != num
    
    
# Function to load coverage files info and stats
def covLoading(covFiles, regiones, resol, discrete=False):
    '''
    param False discrete: False if you want actual values, or threshold if you want
        0 if smaller or 1 if greater or equal
    '''

    notPresent = set()
    # Load files coverage
    covDict = {}
    for regi in regiones:
        covDict[regi] = {}
        for cfi in covFiles:
            marker = cfi.split('/')[-1].split('_')[0]
            marker = marker.replace('.', '')
            covDict[regi][marker] = []
        
            with open(cfi, 'r') as f:
                for line in f:
                    line = line.rstrip().split('\t')
                    # in line[0] we would have region id
                    if line[0] in covDict[regi][marker]:
                        if discrete == False:
                            covDict[regi][marker].append(float(line[1]))

                        else:
                            if float(line[1]) < discrete:
                                covDict[regi][marker].append(0)
                            else:
                                covDict[regi][marker].append(1)
                    else:
                        # just in case we want too see regions we havent add to analysis
                        notPresent.add(line[0])
    marks = sorted(covDict[regi].keys())         
    nmark = len(marks)

    # Check if region length is ok
    for regi in covDict.keys():
        for m in marks:
            reg = regiones[regi]
            regionStart = reg[1]
            regionEnd = reg[2]
            longi = ((regionEnd - regionStart) / resol) + 1
            if len(covDict[regi][m]) != longi:
                difference = len(covDict[regi][m]) - longi
                print 'Region %s has %s more/less positions in file %s' \
                %(regi, difference, m)
                # If more in file, we remove from the end #### CHANGE WHEN CORRECT FILES ###
                exit()
    # get coverage total average and standar deviation
    statsCov = {}
    for regi in covDict.keys():
        statsCov[regi] = {}
        for m in marks:
            # get total average
            statsCov[regi][m] = {'tAvg': np.mean(covDict[regi][m]), 'tSD':np.std(covDict[regi][m])}

    return covDict, statsCov, marks, nmark


def center_of_mass(self):
    """
    Gives the center of mass of a model

    :returns: the center of mass of a given model
    """
    r_x = sum(self['x'])/len(self)
    r_y = sum(self['y'])/len(self)
    r_z = sum(self['z'])/len(self)
    return dict((('x', r_x), ('y', r_y), ('z', r_z)))

def getBinDistribution(cdistDict, mini, modelRadius, resThres, listOfBins = '', 
                       maxRadi=False, groupBy = 'binNumber'):

    '''
    groupBy str 'binNumber': Options are 'binNumber' and 'density'
    '''
    # dictionary to show bins distribution
    binsDistdic = {}

    if listOfBins == '':
        listOfBins = cdistDict.keys()
        binsRep = Counter(listOfBins)
    else:
        binsRep = Counter(listOfBins)
        
    # get mean and std
    mean = np.mean(binsRep.values())
    std = np.std(binsRep.values())
        
    # compute model radius
    maximumDistance = resThres * (int(modelRadius/resThres) + 1)
    # Need to create this one in tha Density case to avoid white areas at the end
    maximumDistanceDensi = round(modelRadius + 1)
    valRange = range(0, int(maximumDistance) + 1, int(resThres))
    # set the limit of our search diameter
    if maxRadi == False:
        maxRadi2 = maximumDistance
        maxRadi2Densi = maximumDistanceDensi
    else:
        maxRadi2 = maxRadi
        maxRadi2Densi = maxRadi

    # If groupBy != 'density' we set valRangeTrimmed, otherwise will ve overwriten
    valRangeTrimmed = [i for i in valRange if i <= maxRadi2]
    if groupBy == 'density':
        # First of all we modify maximumDistance so its real
        valRangeTrimmed = [0]
        # get list with all cut sites
        # Since distances are from center of bin, we need to ad a radius of a bin
        #modelDiameter += resThres/2 ### !!!!!!!!! NO LO VEO NECESARIO

        # Now we need to compute de volume of the first sphere (with r equal to 
        #bin radius)
        firstVolume = (4 * np.pi * (resThres)**3) / 3
        valRangeTrimmed.append(resThres)
        # With info of the first volume we can multiply it for each radius in the plot
        #to obatain the distance for the next radius
        n = 2
        nextRadi = round((((n * firstVolume) * 3) / (4 * np.pi))**(1./3))
        # ULTIMO CAMBIO DE < A <=
        while nextRadi <= maxRadi2Densi:
            valRangeTrimmed.append(nextRadi)
            # Find in next radius away
            n += 1
            nextRadi = round((((n * firstVolume) * 3) / (4 * np.pi))**(1./3))


        # create the valrange we use to store distance data
        valRange = list(sorted(set(valRangeTrimmed + [round(modelRadius) + 1])))
        #print valRange
        # we create the value ranges we have prepared
        for mrange in valRange[1:]:
            #mrange = np.mean([v, valRange[nv + 1]])
            binsDistdic[mrange] = []

            
        # for each bin in our model for all values
        for nbin in listOfBins:
            for nb in cdistDict[nbin]:
                # we check between which range positions should our value be located and add it
                pos = [n+1 for n, va in enumerate(valRange) if va <= nb < valRange[n+1]][0]
                mrange = valRange[pos]

                binsDistdic[mrange].append(binsRep[nbin])

                # Add values to the continuous list (not binned)
                # Is a defaultdict, so if value doesnt exist is the same, will add it
                #histdictContinuous[k][nb] += covDict[k][region][nbin]

           
    # We adjust valRange and maxRadi so the last one points to the range value locating the limit
    #maxRadi2 = maxRadi2 - ((valRange[1] - valRange[0]) / 2)
    #print valRange
    else:
        
        # we create the value ranges we have prepared
        for nv, v in enumerate(valRange[:-1]):
            #mrange = np.mean([v, valRange[nv + 1]])
            mrange = valRange[nv + 1]
            binsDistdic[mrange] = []
            
        # for each bin in our model for all values
        for nbin in listOfBins:
            for nb in cdistDict[nbin]:
                # we check between which range positions should our value be located and add it
                # if we substract the ranging starting point to our value and divide it by the range
                #length we get the position in valRange where our value is located
                pos = int((nb - int(mini)) / (valRange[1] - valRange[0]))
                mrange = valRange[pos + 1]
                #mrange = np.mean([valRange[pos], valRange[pos + 1]])

                binsDistdic[mrange].append(binsRep[nbin])


    # If there is a maximum distance set from the point of interest we romeve the ranges above
    if maxRadi != False:
        for di in binsDistdic.keys():
            if di > maxRadi2:
                del binsDistdic[di]
  



    # Normalize data 
    for piece in binsDistdic.keys():
        #prev = histdict2[k][piece] 
        binsDistdic[piece] = sum([i for i in binsDistdic[piece]])
        # !!!!!!! if the result is a nan, means that there wasnt enough bins in this distance, so we change 
        # it for a 0
        if isNaN(binsDistdic[piece]):
            #print k, piece, 'removed'
            binsDistdic[piece] = 0

    if maxRadi == False:
        return binsDistdic, valRange
    else:
        return binsDistdic, list(sorted(set(valRangeTrimmed)))

# Function to get distance from interest point 
def getDistancesFromPoint(mods, cluster, interest='center'):
    '''
    :param mods: Clustered ensemble of models in tadbit StructuralModels object
    :param cluster: Integer indicating the models from wich cluster will
        be measured.
    :param 'center' interest: Point from wich we measure distances. If
        not set model center of mass will be set as default. If an integer
        is probided, the distances will be measured from the bin this integer
        relates to
    :returns: cdistDict, a dictionary with model bins as keys (integers) and
                    the distance from this bin, in all the models from the 
                    ensemble belonging to the selected cluster, to the interest
                    bin
              cdistDictMean, a dictionary like cdistDict but just containing the 
                    mean distances
              mini, a float with the minimum distance detected
              maxi, a float with the maximum distance detected
              
    
    '''
    # Here we store distances from interest point 
    # get cluster ids
    models = [mods[str(m)]['index'] for m in mods.clusters[cluster]]
    models = [mods[mdl] for mdl in models]

    # variables to know values range
    mini = 100
    maxi = 0

    cdistDict = {}
    # create bin indexes
    for nbin in range(len(models[0]['x'])):
        cdistDict[nbin] = []
    for mo in models:
        if interest == 'center':
            # get center of mass
            mcen = center_of_mass(mo)
        else:
            # In this case interest should be a bin
            mcen = {'x': mo['x'][interest], 'y': mo['y'][interest], 'z': mo['z'][interest]}

        # go for each bin
        for nbin, x in enumerate(mo['x']):
            pos = {'x': x, 'y': mo['y'][nbin], 'z': mo['z'][nbin]}
            # get distance 
            dist = sqrt((mcen['x'] - pos['x'])**2 + (mcen['y'] - pos['y'])**2 + (mcen['z'] - pos['z'])**2)
            cdistDict[nbin].append(dist)
            # store range
            mini = min(mini, dist)
            maxi = max(maxi, dist)



    # get distance dict with mean
    cdistDictMean = {}
    for nbin in cdistDict.keys():
        cdistDictMean[nbin] = np.mean(cdistDict[nbin])

    return cdistDict, cdistDictMean, mini, maxi


# Funtion to get the distribution from interest point of all bins
# Funtion to get the distribution from interest point of all bins
def getMarkerDistribution(cdistDict, cdistDictMean, covDict, statsCov, marks, region, 
                          mini, modelRadius, resThres, maxRadi=False, groupBy = 'binNumber',
                         discrete=False,method = 'divVolume', pval=0.01):


    
    '''
    zscore bool True: wether to use zscores or raw data instead
    groupBy str 'binNumber': Options are 'binNumber' and 'density'. binNumber stands for
        a radius distribution of fixed distance whereas 'density' will do all spheres with 
        mantaining equal volumes
    para 'zscore' method: Which method to use at the time to normalize. Can choose between
        zscore, contingency, divVolume and percentage
    param False discrete: NOT False if you want to build a contingency table and obtain the 
        odds
    '''
    # Create list to tell when there is no data
    noData = set()
    # Get histogram input
    histdict = {}
    # dictionary for pseudo normalization
    histdict2 = {}
    histdictMean = {}
    #histdictContinuous = {}
    #histdictContinuousMean = {}

    # compute model radius
    maximumDistance = resThres * (int(modelRadius/resThres) + 1)
    # Need to create this one in tha Density case to avoid white areas at the end
    maximumDistanceDensi = round(modelRadius + 1)
    valRange = range(0, int(maximumDistance) + 1, int(resThres))
    # set the limit of our search diameter
    if maxRadi == False:
        maxRadi2 = maximumDistance
        maxRadi2Densi = maximumDistanceDensi
    else:
        maxRadi2 = maxRadi
        maxRadi2Densi = maxRadi

    # If groupBy != 'density' we set valRangeTrimmed, otherwise will ve overwriten
    valRangeTrimmed = [i for i in valRange if i <= maxRadi2]
    if groupBy == 'density':
        # First of all we modify maximumDistance so its real
        valRangeTrimmed = [0]
        # get list with all cut sites
        # Since distances are from center of bin, we need to ad a radius of a bin
        #modelDiameter += resThres/2 ### !!!!!!!!! NO LO VEO NECESARIO

        # Now we need to compute de volume of the first sphere (with r equal to 
        #bin radius)
        firstVolume = (4 * np.pi * (resThres)**3) / 3
        valRangeTrimmed.append(resThres)
        # With info of the first volume we can multiply it for each radius in the plot
        #to obatain the distance for the next radius
        n = 2
        nextRadi = round((((n * firstVolume) * 3) / (4 * np.pi))**(1./3))
        # ULTIMO CAMBIO DE < A <=
        while nextRadi <= maxRadi2Densi:
            valRangeTrimmed.append(nextRadi)
            # Find in next radius away
            n += 1
            nextRadi = round((((n * firstVolume) * 3) / (4 * np.pi))**(1./3))


        # create the valrange we use to store distance data
        # if we have not reached the end of the model because it lies in the middle
        #of two radii checked
        if (round(modelRadius) + 1) > max(valRangeTrimmed):
            valRange = list(sorted(set(valRangeTrimmed + [round(modelRadius) + 1])))
        else:
            valRange = valRangeTrimmed
        #print valRange
        # in each marker
        for k in covDict.keys():
            histdict[k] = {}
            histdict2[k] = {}
            histdictMean[k] = {}
            #histdictContinuous[k] = defaultdict(int)
            #histdictContinuousMean[k] = defaultdict(int)

            # we create the value ranges we have prepared
            for mrange in valRange[1:]:
                #mrange = np.mean([v, valRange[nv + 1]])
                histdict[k][mrange] = 0
                histdict2[k][mrange] = []
                histdictMean[k][mrange] = 0

            # for each bin in our model for all values
            for nbin in cdistDict.keys():
                for nb in cdistDict[nbin]:
                    # we check between which range positions should our value be located and add it
                    pos = [n+1 for n, va in enumerate(valRange) if va <= nb < valRange[n+1]][0]
                    mrange = valRange[pos]

                    histdict[k][mrange] += covDict[k][nbin]
                    histdict2[k][mrange].append(covDict[k][nbin])

                    # Add values to the continuous list (not binned)
                    # Is a defaultdict, so if value doesnt exist is the same, will add it
                    #histdictContinuous[k][nb] += covDict[k][nbin]

                # Same with the mean distance values for a Bin in all models
                pos = [n+1 for n, va in enumerate(valRange) if va <= cdistDictMean[nbin] < valRange[n+1]][0]
                mrange = valRange[pos]
                histdictMean[k][mrange] += covDict[k][nbin]

                # Add values to the continuous list (not binned)
                # Is a defaultdict, so if value doesnt exist is the same, will add it
                #pos = cdistDictMean[nbin]
                #histdictContinuousMean[k][pos] += covDict[k][nbin]

    # We adjust valRange and maxRadi so the last one points to the range value locating the limit
    #maxRadi2 = maxRadi2 - ((valRange[1] - valRange[0]) / 2)
    #print valRange
    else:
        # in each marker
        for k in covDict.keys():
            histdict[k] = {}
            histdict2[k] = {}
            histdictMean[k] = {}
            #histdictContinuous[k] = defaultdict(int)
            #histdictContinuousMean[k] = defaultdict(int)

            # we create the value ranges we have prepared
            for nv, v in enumerate(valRange[:-1]):
                #mrange = np.mean([v, valRange[nv + 1]])
                mrange = valRange[nv + 1]
                histdict[k][mrange] = 0
                histdict2[k][mrange] = []
                histdictMean[k][mrange] = 0

            # for each bin in our model for all values
            for nbin in cdistDict.keys():
                for nb in cdistDict[nbin]:
                    # we check between which range positions should our value be located and add it
                    # if we substract the ranging starting point to our value and divide it by the range
                    #length we get the position in valRange where our value is located
                    pos = int((nb - int(mini)) / (valRange[1] - valRange[0]))
                    mrange = valRange[pos + 1]
                    #mrange = np.mean([valRange[pos], valRange[pos + 1]])

                    histdict[k][mrange] += covDict[k][nbin]
                    histdict2[k][mrange].append(covDict[k][nbin])

                    # Add values to the continuous list (not binned)
                    # Is a defaultdict, so if value doesnt exist is the same, will add it
                    #histdictContinuous[k][nb] += covDict[k][nbin]

                # Same with the mean distance values for a Bin in all models
                pos = int((cdistDictMean[nbin] - int(mini)) / (valRange[1] - valRange[0]))
                #histdictMean[k][np.mean([valRange[pos], valRange[pos + 1]])] += covDict[k][nbin]
                mrange = valRange[pos + 1]
                histdictMean[k][mrange] += covDict[k][nbin]

                # Add values to the continuous list (not binned)
                # Is a defaultdict, so if value doesnt exist is the same, will add it
                #pos = cdistDictMean[nbin]
                #histdictContinuousMean[k][pos] += covDict[k][nbin]
        #return histdict2

    # Get values distribution 
    #listMark = []
    #for m in histdict2.keys():
    #    tempi = []
    #    for rang in histdict2[m].keys():
    #        # Always PI - TB
    #        tempi.append(histdict2[m][rang])  
    #    listMark += tempi

    # If there is a maximum distance set from the point of interest we remove the ranges above
    if maxRadi != False:
        for k in marks:
            for di in histdict2[k].keys():
                if di > maxRadi2:
                    del histdict2[k][di]
                    del histdict[k][di]
                    del histdictMean[k][di]

            #for di in histdictContinuous[k].keys():
            #    if di > maxRadi:
            #        del histdictContinuous[k][di]
            #        del histdictContinuousMean[k][di]

    
    ## Normalize data in histdict2
    # obtain numer of models been analysed
    nmodels = len(cdistDict[cdistDict.keys()[0]])
    # start iterating over particles
    for k in marks:
        # get number of particles with positive and negative mark
        # number of particles with positive signal
        regiPositive = np.nansum(covDict[k]) * nmodels
        # number of particles with negative signal
        regiNegative = (len(covDict[k]) - np.nansum(covDict[k])) * nmodels
        for piece in histdict2[k].keys():
            if method == 'zscore':
                #prev = histdict2[k][piece] 
                histdict2[k][piece] = ((np.mean([i for i in histdict2[k][piece]]) - statsCov[k]['tAvg'])\
                                       / statsCov[k]['tSD'])
            # if we want odds ratio from a contingency table
            if discrete != False:
                if method == 'contingency':
                    # Get positive presence and negative signal number in the spherical shell
                    positive = np.nansum(histdict2[k][piece])
                    negative = len(histdict2[k][piece]) - positive
                    # get positive and negatives out of the shell
                    restPos = regiPositive - positive
                    restNeg = regiNegative - negative
                    contingencyTab = [[positive,negative],[restPos,restNeg]]
                    # get odds
                    #print contingencyTab, k, piece
                    oddsratio, pvalue = stats.fisher_exact(contingencyTab)
                    # convert to logarithm if significant and not 0
                    #if k == 'NKX61':
                    #print contingencyTab
                    #print piece, oddsratio, pvalue
                    if pvalue <= pval and oddsratio != 0:
                        oddsratio = np.log(oddsratio)
                    else:
                        oddsratio=0
                    # assign log Odds ratio value
                    #if k == 'NKX61':
                    #    print oddsratio
                    histdict2[k][piece] = oddsratio
                elif method == 'percentage':
                    # obtain percentage of positives in our spherical shell
                    shellPositives = np.nansum([i for i in histdict2[k][piece]])
                    wholePositives = np.nansum(covDict[k]) * nmodels
                    if wholePositives != 0:
                        #print shellPositives, wholePositives
                        histdict2[k][piece] = round((shellPositives / float(wholePositives)) * 100)
                    else:
                        histdict2[k][piece] = 0
                        noData.add(k)
            elif method == 'divVolume':
                # divide by volume and multiply result by 1000 to get higher values
                histdict2[k][piece] = (np.nansum([i for i in histdict2[k][piece]]) / firstVolume) * 1000

            # if the result is a nan, means that there wasnt enough bins in this distance, so we change 
            # it for a 0
            if isNaN(histdict2[k][piece]) or histdict2[k][piece] == float('Inf'):
                #print k, piece, 'removed'
                histdict2[k][piece] = 0
    
    print 'YOY SHOULD CHANGE THE PART THAT COLLAPSES RADIUS SMALLER IN DIFFERENCE THAN 1NM'
    print 'AT LEAST SHOW A WARNING OR STOP IN THERE'
    if maxRadi == False:
        return histdict2, valRange, noData, (histdict, histdictMean) # ,histdictContinuous, histdictContinuousMean)
    else:
        return histdict2, list(sorted(set(valRangeTrimmed))), noData, (histdict, histdictMean)
    