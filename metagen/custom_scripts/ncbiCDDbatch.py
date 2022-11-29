"""
Takes an input fasta file, splits it, and iteratively searchs the NCBI conserved
domain BATCH search website for hits & returns sequences with specified hits.

- Requires selenium to be installed
- Requires a specified webdriver (Chrome by default)
- Waits for your browser for a maximum of 10 minutes if browser freezes
- Requires you to modify the Keys.COMMAND if you use windows
- Make sure your version of chrome driver is compatible with your chrome

"""

#! /usr/local/bin/python3.10

import os
from itertools import islice
from selenium import webdriver
from selenium.webdriver.support.ui import WebDriverWait
from selenium.webdriver.support import expected_conditions as EC
from selenium.webdriver.common.by import By
from selenium.common.exceptions import TimeoutException
from selenium.webdriver.common.keys import Keys
import pyperclip

inputdir = input("\nEnter the input file directory with all fasta files: ")
maxseqinput = input(
    "\nEnter the maximum number of sequences to load in each entry (max 950): "
)


# Read input file fasta file function
def fasta2dict(fastafile):
    """Reads a fasta file and returns a dictionary"""
    try:
        with open(fastafile, "r", errors="ignore", encoding="UTF-8") as inputfile:
            content = inputfile.read()
            contentsplit = content.split(">")

            cleanlist = []
            for seql in contentsplit:
                seqr = repr(seql)
                clean = seqr.replace("\\n", "*NEWLINE*", 1)
                clean2 = clean.replace("\\n", "")
                clean3 = clean2.replace("'", "")
                cleanlist.append(clean3)

            seq_dict = {}
            for seqc in cleanlist:
                seqto = seqc
                seqsplit = seqto.split("*NEWLINE*")
                if len(seqsplit) > 1:
                    seqkey = ">" + seqsplit[0]
                    seq_dict[seqkey] = seqsplit[1]

            return seq_dict

    except FileNotFoundError:
        print("Could not open/read file: ", fastafile)


# split fasta dictionary into 950 sequence chunks for CDD batch results
def chunks(inputdict, maxsize=950):
    """split fasta dictionary into 950 sequence chunks for CDD batch results"""
    it_value = iter(inputdict)
    for _ in range(0, len(inputdict), maxsize):
        yield {k: inputdict[k] for k in islice(it_value, maxsize)}


# make new output directory
os.chdir(inputdir)
outputdir = os.path.join(inputdir, "CDsearch_batch_res")
FILECNTNAME = 1
while True:
    try:
        os.mkdir(outputdir)
        break
    except FileExistsError:
        outputdir += "_" + str(FILECNTNAME)
        FILECNTNAME += 1


# check to see if file has already been procssessed
# and if not process the file
for filename in os.listdir(inputdir):
    os.chdir(inputdir)
    # extract different filename components
    root, ext = os.path.splitext(filename)

    # check to see if file is fasta or fa format
    if filename == "CDD_batch_res":
        continue

    if ext in (".fa", ".fasta", ".fna", ".fas"):

        # read fasta file
        sequencedict = fasta2dict(filename)

        for sliceddict in chunks(sequencedict, int(maxseqinput)):

            # convert dictionary into string to input into NCBI search text box
            SEQUENCE = "\n".join(
                "\n".join((key, val)) for (key, val) in sliceddict.items()
            )

            # make sure to download the chromedriver.exe and place it in a known location
            # Needs to be in your PATH /usr/local/bin/


            # set the default download location to output directory
            chromeopt = webdriver.ChromeOptions()  # object of ChromeOptions
            prefs = {"download.default_directory": outputdir}
            chromeopt.add_experimental_option("prefs", prefs)
            # chromeopt.add_argument("--headless")
            # WINDOW_SIZE = "1920,1080"
            # chromeopt.add_argument("--window-size=%s" % WINDOW_SIZE)

            driver = webdriver.Chrome(options=chromeopt)
            driver.get("https://www.ncbi.nlm.nih.gov/Structure/bwrpsb/bwrpsb.cgi")

            DELAY = 10 * 60  # the MAXIMUM time that the program will wait
            # this textarea correspond to the box where you past the desired sequence
            TXTBOX = '//*[@id="frm_New_Search"]/div/table/tbody/tr/td/table/tbody/tr[1]/td[1]/div[2]/textarea'
            SRCHBTT = '//*[@id="frm_New_Search"]/div/table/tbody/tr/td/table/tbody/tr[2]/td/table/tbody/tr/td[3]/input'
            DWNLDHBTT = '//*[@id="tbl_DLPanel"]/tbody/tr[3]/td[3]/input'

            textarea = (
                WebDriverWait(driver, DELAY)
                .until(EC.visibility_of_element_located((By.XPATH, TXTBOX)))
                .find_elements(By.XPATH, TXTBOX)
            )[0]

            # this command send the sequence you want to search
            # pyperclip is used to paste lots of text very fast
            pyperclip.copy(SEQUENCE)
            # you need to change this to Keys.CONTROL for windows
            textarea.send_keys(Keys.COMMAND, "v")

            # the button variable searchs for the submit button
            try:
                subbutton = (
                    WebDriverWait(driver, DELAY)
                    .until(EC.visibility_of_element_located((By.XPATH, SRCHBTT)))
                    .find_element(By.XPATH, SRCHBTT)
                )
                subbutton.click()
            except TimeoutException as ex:
                ISRUNNING = 0
                print("Exception has been thrown. " + str(ex))
                driver.close()

            # Wait until the results page is loaded and download results
            try:
                downbutton = (
                    WebDriverWait(driver, DELAY)
                    .until(EC.visibility_of_element_located((By.XPATH, DWNLDHBTT)))
                    .find_element(By.XPATH, DWNLDHBTT)
                )
                downbutton.click()
            except TimeoutException as ex:
                ISRUNNING = 0
                print("Exception has been thrown. " + str(ex))
                driver.close()

            # close current window
            # driver.close()
