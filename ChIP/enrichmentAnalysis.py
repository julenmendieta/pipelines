from selenium import webdriver
from selenium.webdriver.common.keys import Keys

from selenium.webdriver.support.ui import WebDriverWait
from selenium.webdriver.support import expected_conditions as EC
from selenium.webdriver.common.by import By
from selenium.common.exceptions import TimeoutException

from collections import defaultdict
import pandas as pd
import numpy as np
import time
import os

import seaborn as sns
from matplotlib import pyplot as plt


def scrollToElem(driver, element):
    desired_y = (element.size['height'] / 2) + element.location['y']
    window_h = driver.execute_script('return window.innerHeight')
    window_y = driver.execute_script('return window.pageYOffset')
    current_y = (window_h / 2) + window_y
    scroll_y_by = desired_y - current_y

    driver.execute_script("window.scrollBy(0, arguments[0]);", scroll_y_by)
    
def download_wait(directory, timeout, nfiles=None):
    """
    Wait for downloads to finish with a specified timeout.

    Args
    ----
    directory : str
        The path to the folder where the files will be downloaded.
    timeout : int
        How many seconds to wait until timing out.
    nfiles : int, defaults to None
        If provided, also wait for the expected number of files.

    """
    seconds = 0
    dl_wait = True
    while dl_wait and seconds < timeout:
        time.sleep(1)
        dl_wait = False
        files = os.listdir(directory)
        if nfiles and len(files) != nfiles:
            dl_wait = True

        for fname in files:
            if fname.endswith('.crdownload'):
                dl_wait = True

        seconds += 1
    return seconds

def getGREATout(GREATout, infile, minGene=3, minGeneGO=5, maxGeneGO=500,
               executable_path='chromedriver', interactive=False,
               backGround=False):
    '''
    Function to acces GREAT GO webpage and get GO associated to a specific
        coordinate file
        Works with Selenium version 3.141.0, but it need to be ran interactively
        
    :param GREATout: output directory for the distance to TSS plot and the GO
        terms
    :param infile: path to the bed file with the coordinates to look into
    :param 3 minGene: Minimum number of genes from our list that should be
        in any GO term
    :param 5 minGeneGO: minimum number of genes in GO to be check
    :param 500 maxGeneGO: maximum number of genes in GO to be check
    :param 'chromedriver' executable_path: path to chromedriver executable
        (only name if is in PATH)
    :param False interactive: If you want to see the browser as changes are done
    
    '''
    options = webdriver.ChromeOptions()
    # do not open the web browser
    if not interactive:
        options.add_argument("--headless")
    #prefs = {"download.default_directory" : GREATout}
    #options.add_experimental_option("prefs", prefs)

    # options to change download location and automatically download pdf instead of opening
    options.add_experimental_option('prefs', {
    "download.default_directory": '%s' %GREATout, #Change default directory for downloads
    "download.prompt_for_download": False, #To auto download the file
    "download.directory_upgrade": True,
    "plugins.always_open_pdf_externally": True #It will not show PDF directly in chrome
    })

    # Using Chrome to access web
    driver = webdriver.Chrome(options=options, executable_path=executable_path)
    # Open the website
    driver.get('http://great.stanford.edu/public/html/')

    # select species
    # Select the id box (mm10)
    species_tick = driver.find_element_by_xpath('/html/body/div[2]/div[4]/div/form/fieldset/div[1]/div/ul/li[3]/label/input')
    # click on it
    species_tick.click()

    # uppload test regions file
    driver.find_element_by_id("fgFile").send_keys(infile)
    
    # Uppload background file if present
    if backGround != False:
        # Select the bed file box in background regions
        backg_tick = driver.find_element_by_xpath("/html/body/div[2]/div[4]/div/form/fieldset/div[3]/div/ul/li[2]/label/input")
        # click on it
        backg_tick.click()

        # uppload file
        driver.find_element_by_id("bgFile").send_keys(backGround)
    

    # submit
    submit_buttom = driver.find_element_by_id("submit_button")
    submit_buttom.click()

    # wait for page to load
    time.sleep(10)

    ## job page
    #get current tab handle
    mainTab = driver.current_window_handle
    
    # Open distances section if we used background
    if backGround != False:
        distances_button = driver.find_element_by_xpath("/html/body/div[2]/div[12]/div/h3/a/img")
        distances_button.click()

    # histogram with distances to TSS (opens new tab)
    xpathSection=13
    try:
        distTSS_file = driver.find_element_by_xpath(f"/html/body/div[2]/div[{xpathSection}]/div/div/div[2]/div[2]/div[2]/a")
        distTSS_file.click()
    except:
        xpathSection=14
        distTSS_file = driver.find_element_by_xpath(f"/html/body/div[2]/div[{xpathSection}]/div/div/div[2]/div[2]/div[2]/a")
        distTSS_file.click()
        print("Warning: Your set hits a large fraction of the genes in the genome, \
which often does not work well with the GREAT Significant by Both view due to a \
saturation of the gene-based hypergeometric test.")

    # get all tab Ids and switch to new one
    #chwd = driver.window_handles
    #driver.switch_to.window([c for c in chwd if c != mainTab])

    ## defined parameters (we will do clicks, so has to be in order from top to down)
    # in the future i might think about scrolling in serenity
    # open section
    paramSection_button = driver.find_element_by_id('global_controls_container')
    driver.execute_script("arguments[0].style.display = 'block';", paramSection_button)

    ## Set some parameters
    # min number of our genes in each GO
    minGene_button = driver.find_element_by_id('minAnnotFgHitGenes')
    minGene_button.clear()
    minGene_button.send_keys(minGene)
    driver.find_element_by_id('minAnnotFgHitGenesSet').click()

    # min number of gens in a GO to be used
    minGeneGO_button = driver.find_element_by_id('allMinAC')
    minGeneGO_button.clear()
    minGeneGO_button.send_keys(minGeneGO)
    driver.find_element_by_id('allACSet').click()
    # max number of gens in a GO to be used
    minGeneGO_button = driver.find_element_by_id('allMaxAC')
    minGeneGO_button.clear()
    minGeneGO_button.send_keys(maxGeneGO)
    driver.find_element_by_id('allACSet').click()
    
    ## Define GOs to show
    # This initially was before setting some parameters
    # ontology terms to show
    hideall_button = driver.find_element_by_id('hide_all_ontos_btn')
    scrollToElem(driver, hideall_button)
    hideall_button.click()

    # GO molecular function
    GOMol_button = driver.find_element_by_id('onto_cb_GOMolecularFunction')
    GOMol_button.click()
    # GO bio process
    GOBio_button = driver.find_element_by_id('onto_cb_GOBiologicalProcess')
    GOBio_button.click()
    
    
    # download all data
    downloadAll_button = driver.find_element_by_xpath(f'/html/body/div[2]/div[{xpathSection+1}]/div/h3/select/option[2]')
    scrollToElem(driver, downloadAll_button)
    downloadAll_button.click()
    # in headless mode we need to set some time or wont have any download
    download_wait(GREATout, 10, 2)
    #time.sleep(10)
    
    # close the browser
    driver.close()




################### Plot related

def greatToDf(GREATout, minGreatP, BinomFdrQ, HyperFdrQ, Ontology, 
             Desc, goID, BinomP, ObsGenes, TotalGenes,
             backGround=False):
    GreatOut = defaultdict(list)
    with open(f"{GREATout}/greatExportAll.tsv", 'r') as f:
        header = f.readline()
        header = f.readline()
        header = f.readline()
        header = f.readline()[2:-1].split('\t')
        
        for line in f:
            if not line.startswith('#'):
                line = line.split('\t')
                # Store only significnat terms
                store = False
                if backGround == False:
                    if ((float(line[BinomFdrQ]) < minGreatP) and 
                        (float(line[HyperFdrQ]) < minGreatP)):
                        GreatOut['Ontology'] += [line[Ontology]]
                        GreatOut['GO description'] += [line[Desc]]
                        GreatOut['ID'] += [line[goID]]
                        GreatOut['-log10(Binomial p value)'] += [-np.log10(float(line[BinomP]))]
                        GreatOut['% of observed genes'] += [round((int(line[ObsGenes]) / int(line[TotalGenes])) * 100, 2)]
                else:
                    if float(line[HyperFdrQ]) < minGreatP:
                        GreatOut['Ontology'] += [line[Ontology]]
                        GreatOut['GO description'] += [line[Desc]]
                        GreatOut['ID'] += [line[goID]]
                        GreatOut['-log10(Hyper FDR Q-Val)'] += [-np.log10(float(line[HyperFdrQ]))]

    df_great = pd.DataFrame.from_dict(GreatOut)
    return df_great


def plotGreat(df_onto2, title='', plot1='-log10(Binomial p value)',
             plot2 = False, yLabel='GO description',
             color=None, palette=None, figsize=(20, 8)):

    
    if plot2 == False:
        ncol = 1
        if figsize == 'auto':
            figsize = (10, max((0.35 * len(df_onto2.index)) + 1, 3))
    else:
        ncol = 2
        if figsize == 'auto':
            figsize = (20, max((0.35 * len(df_onto2.index)) + 1, 3))
    
    fig = plt.figure(figsize=figsize)
    #ax = fig.add_subplot(2, 1, 1)
    #ax2 = fig.add_subplot(2, 2, 1)

    ax1 = plt.subplot2grid((len(df_onto2.index), ncol), (0, 0), colspan=1, rowspan=len(df_onto2.index))
    sns.barplot(y=yLabel, x=plot1, data=df_onto2, orient='h', 
                       order=df_onto2[yLabel], ax=ax1, color=color, palette=palette)

    ax1.grid(zorder=0)
    ax1.set_xlabel(plot1)

    # set axis limit
    maxi = max(d for d in df_onto2[plot1] if np.isinf(d) == False)
    ax1.set_xlim(np.nanmin(df_onto2[plot1]) * 0.9, maxi * 1.01)

    if plot2 != False:
        ax2 = plt.subplot2grid((len(df_onto2.index), ncol), (0, 1), colspan=1,  rowspan=len(df_onto2.index))
        sns.barplot(y=yLabel, x=plot2, data=df_onto2, orient='h', 
                           order=df_onto2[yLabel], ax=ax2, color=color, palette=palette)

        ax2.set_xlim(0, 100)
        ax2.set_yticks([])
        ax2.grid(zorder=0)
        ax2.set_xlabel(plot2)

    if title != '':
        fig.suptitle(title)
    plt.show()

    return fig
