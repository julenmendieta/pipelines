import numpy as np
import copy
from matplotlib import pyplot as plt


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
                 fromOri = 0, fromEnd = 0, savefig = True):
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

    if savefig == True:
        pdf.savefig( f )
    plt.show()


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
                 fromOri = 0, fromEnd = 0, savefig = True, ylim1 = (-30, 30), ylim2 = (-4, 4)):
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

    if savefig == True:
        pdf.savefig( f )

    plt.show()

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
                 fromOri = 0, fromEnd = 0, savefig = True):
    fontSize = 15
    # propor is the proportion of size for normal plots against grey bar info plots
    # marksAll is to add the marks in all the plots, and MarksBar is to add marks just in the
    #top bar (is a list of lists)
    xa,ya,eMa,ePa,da=GetXYA(indirs[0] + 'accessibility.cluster%s.dat' %clust,diff=diff)
    xc,yc=GetXY(indirs[0] + 'consistency.cluster%s.dat' %clust,diff=diff)
    xi,yi,eMi,ePi,di=GetXY_smooth4(indirs[0] + 'interactionscluster%s.dat' %clust,diff=diff)

    xa2,ya2,eMa2,ePa2,da2=GetXYA(indirs[1] + 'accessibility.cluster%s.dat' %clust,diff=diff)
    xc2,yc2=GetXY(indirs[1] + 'consistency.cluster%s.dat' %clust,diff=diff)
    xi2,yi2,eMi2,ePi2,di2=GetXY_smooth4(indirs[1] + 'interactionscluster%s.dat' %clust,diff=diff)

    nplot = len(marksBar) + 3 #we alway have 3 plots at least
    #f, (ax3, ax1, ax2, ax5) = plt.subplots(4, sharex=True, sharey=False,figsize=figsize)
    f, allAxis = plt.subplots(nplot, sharex=True, sharey=False,figsize=figsize,
                              gridspec_kw = {'height_ratios':[1] *  ((nplot) - 3) + propor * 3})
    # change tic size
    #f.tick_params(axis='both', which='major', labelsize=10)
    #f.tick_params(axis='both', which='minor', labelsize=8)
    #gs1 = gridspec.GridSpec(figsize[0], figsize[1])
    #gs1.update(wspace=0.025, hspace=0.05)
    color1 = convertRGB([84,39,136])
    color2 = convertRGB([179,88,6])
    allAxis[-3].plot(xc, yc,'-', color = color1, lw=1.75)
    allAxis[-3].plot(xc2, yc2,'-', color = color2, lw=1.75)
    allAxis[-3].set_title('Consistency per particle', fontsize=fontSize)
    #allAxis[-3].axhline(y=np.mean(yc),c="blue")#,linewidth=0.2)
    allAxis[-3].set_ylabel("Consistency %",rotation=rotation, fontsize=fontSize)
    allAxis[-3].set_ylim(0, 100)

    allAxis[-2].plot(xa, ya,'-', color = color1, lw=1.75)
    allAxis[-2].plot(xa2, ya2,'-', color = color2, lw=1.75)
    allAxis[-2].set_title('Accessibility per particle', fontsize=fontSize)
    #allAxis[-2].axhline(y=np.mean(np.nan_to_num(ya)),c="blue")#,linewidth=0.2)
    allAxis[-2].set_ylabel("Accessibility" ,rotation=rotation,multialignment='center',
                          fontsize=fontSize)
    #allAxis[-2].fill_between(xa,eMa,ePa,alpha=0.5, edgecolor='#bebebe', facecolor='#bebebe' )
    allAxis[-2].set_ylim(0, 100)

    allAxis[-1].plot(xi, yi,'-', color = color1, lw=1.75)
    allAxis[-1].plot(xi2, yi2,'-', color = color2, lw=1.75)
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
    start, end = allAxis[-3].get_xlim()
    binrange = range(int(start), int(end) + 1)
    posrange = range(regionStart, regionEnd + resol, 5000)
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

    if savefig == True:
        pdf.savefig( f )
    plt.show()


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


