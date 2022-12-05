# Disentangling the Functional Role of Fungi in Cold Seep Sediment
Hi and welcome! Erfan here. This README.md is a step-by-step guide to reproduce the figures and tables in our study titled "Disentangling the Functional Role of Fungi in Cold Seep Sediment". Some of the steps can take a long time and high computation power. To help save your precious time, we've gone ahead and commented out those sections of the code, and instead directly upload the preprocessed output (denoted as [OPTIONAL]). Feel free to uncomment these code blocks and run the analysis for yourself!


* Note that our instructions are written for a linux/mac workbench, but can be modified easily for Windows and other OS.

____________________________________________

#### 1) Clone this GitHub repository.

1. Open your terminal
2. Enter the directory you'd like to download the files in
3. Make sure you have git installed on your device (https://www.atlassian.com/git/tutorials/install-git)
4. Paste the command ``` git clone https://github.com/erfanshekarriz/SeepFungiNet2022 ```

This should download all the data we have onto your device.

____________________________________________

#### 2) Explore the repository.

Our repository has the following tree structure:
- **data**: where you can find all your data.
  - **Data**: Phyloseq files to explore the microbiome data for yourself (highly encouraged) (https://joey711.github.io/phyloseq/), as well as other files here that we will use for analysis. We purposefully haven't included separate ASV/metadata/taxonomy data so that everything is clean and packaged for you using Phyloseq.
  - **graphs**: The OUTPUT of all figures in our paper.
  - **networks**: Final SPIEC-Easi and igrpah networks that you can use and explore. An added interactive cold seep network graph in .html format I highly recommend to check out.
- **R**: where all the R code is deposited. The name of the R codes allow you to track which figures they output! I've numbered them in the order which I've coded, but you can go in any order if you want, as long as the data is downloaded already with the ```git clone``` command.
- **renv**: file that is important for reproducibility.
- **renv.lock**: file that includes all the packages, version, sources, and in general anything else you would need to set up your environment exactly like ours to stably reproduce our results.
- **coldseepfungi.Rproj**: where you will click to begin your exciting journey of reproducing our results!
- [OPTIONAL] **metagen**: bash scripts written for metagenomic processing
- [OPTIONAL] **metatrans**: bash scripts written for metatranscriptomic processing

____________________________________________

#### 3) Setup environment using renv


First you have to install all the packages necessary for doing the analysis with two easy steps!
1. Open the R project file called ```coldseepfungi.Rproj``` (by double-clicking the file)
2. In the Console, type ```renv::restore()``` to install all packages necessary (might take a while so grab some coffee)

An exact instruction of how to use renv can be found at https://www.youtube.com/watch?v=yc7ZB4F_dc0&ab_channel=RiffomonasProject (at 15 mins 30 seconds)
Credits to Patrick Patrick D. Schloss!

____________________________________________

#### 4) Reproduce the figures!
In the ```Files``` panel enter the ```R``` directory.

Some of the steps can take a long time and high computation power. To help save your precious time, we've gone ahead and commented out those sections of the code, and instead directly upload the preprocessed output (denoted as [OPTIONAL]). Feel free to uncomment these code blocks and run the analysis for yourself! All figures are saved to the ./data/graphs or  ./data/graphs/supplementary directory

- [OPTIONAL] code you can skip because they mostly involve trivial 'cleaning' tasks not relevant to the output of our figures.
- [ANALYSIS] code that doesn't generate a figure, but an important file object that we will use later for analysis.
- [FIGXX] code that specifically generates a figure in our study.

These are the corresponding files and some notes:

1. **make_phyloseq.R** [OPTIONAL]: generate phyloseq files.

2. **alphadiversity_rar.R** [OPTIONAL] [ANALYSIS] [FIG2E] : performs rarefaction using the iNEXT package which can take an extremely long time. The blocks of code that are used to do the computation have been commented out, so you can still used it to generate the figures without having to wait for the output.

3. **phyloseq_split.R** [OPTIONAL]: splits the phyloseq file into seep and non-seep subsets to infer cold seep networks.

4. **phyloseq4network.R** [OPTIONAL]: filters phyloseq to only include taxa prevelant in 20% of samples (N=26)

5. **spieceasi_SLR.R** [ANALYSIS]: a versatile and complex pipeline that can generate SPIEC-Easi networks with either glasso, mb, or slr methods. It automatically turns your SPIEC.Easi object into an igraph weighted and non-weighted file and saves summary statistics of that network.  For SLR, it automatically does β parameterization using extended Bayesian Information Criterion mentioned in our study and needs some setting up. This is the only file that is difficult to reproduce only because it needs all the phyloseq files to be in the same directory as input. If you only have the 16S phyloseq in your directory, then it will make a bacterial-archaeal network. If you only have the 18S fungal phyloseq in your directory, then it will make a fungi-only network. If you have both, then it will make a multi-domain network (bacteria-archaea-fungi). If you don't understand the code please feel free to contact me! I would recommend you to skip this portion since the output files have already been generated for you for the next steps of analysis.

6. **filternetwork.R** [ANLYSIS]: filters the network to only includes ASVs that are present in both ROV1 and ROV2 (the most active cold seep sites). This only eliminates around 5 ASVs.

7. **networkattack.R** [FIG4A]: performs iterative attack and robustness calculation for the networks we've previously generated.

8. **networkeffic.R** [FIG4B]: calculates and plots the node-based efficiency of the BAF, BF, BA, and B network.

9. **domainspecBAF.R** [FIG4D]: looks at the properties of specific ASVs in the BAF network and plots them.

10. **keystonespec.R** [FIG4C]: performs and plots keystone species analysis.

11. **plotnetwork.R** [FIG2]: plots the different cold seep networks.

12. **percentfungi.R** [FIG2D]: plots the barchart for the relative abundance of fungi in all sites compared to other microbial eukaryotes.

13. **envdataplot.R** [FIG1A] [FIG1C]: plots the environmental and abiotic distributions across sites and depths.

14. **LDAeffec.R** [ANALYSIS] [FIG2A]: performs and plots LefSe results.

15. **robustPCA.R** [OPTIONAL]: evalute the robust PCA results of the different networks

16. **venndiagram.R** [FIG2C]: count the number and percentage of 18S fungal ASVs shared by different sites.

17. **phylumtree.R** [FIG2A]: plot and calucate the relative and CLR abundance of 18S fungal phylums.

18. **selba.R** [ANALYSIS] [FIG5A] [FIG5B]: performs and plots SELECTION OF BALANCES (selbal) to find fungal ASVs that correlate with high methane concentrations.

19. **visNetwork.R** [FIG5C]: make the interactive cold seep network in HTML format.

20. **globalmap.R** [FIG6A]: plot the global map of meta-omics data.

21. **metatransres.R** [FIG6B]: plot the barchart of metatranscriptomic bar chart.


____________________________________________

#### 5) Reproduce our meta-omics analysis

Reproducing our metaomics isn't as straight forward and needs a high performance computer and lots of patience, but we can assure you is fully possible. We use the VEBA suite (https://github.com/jolespin/veba) for all our analysis, so you only need to set up the correct version (v1.0.3) on your system, download our raw files, and get to working.

The VEBA suite ensures that each step is logged so you can track our analysis if you're ever lost.
The VEBA outputs can be accessed via our Figshare link in our manusript.

____________________________________________

### CITATION
