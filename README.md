# Pasteur_Internship

## Introduction

The main scope of my internship at the Hellenic Pasteur Institute involved analysing the polyreactive behaviour of Natural Antibodies (NAbs), in order to ascertain the mechanisms of epitopic valency through a combination of Amino Acid (AA) and biophysical descriptors. Publicly available NAb **CDR3H loop** and corresponding antigenic **epitope** sequence data were collected to compute their AA sequence biophysical descriptors (i.e. Natural Antibody Database). Sequence and biophysical descriptor data were then used in protracted clustering analysis, automated by the platform. Following multiple tests and modifications, a final design was decided upon with platform able to:

- Compute Composition Transition and Distribution (CTD) statistics on each protein sequence received; mimicing LISA statistics but offering larger-detail based-depth.
- Automate BLOSUM45 pairwise comparison between sequences.
- Identify O and N glycosylation sites, through regular expression formulae and small optimised algorithms.
- Parse foreign `csv` formats, merge to internal storage system, update database with computed statistics / metadata on all NAbs.
- Merge pairwise CTD and BLOSUM45 comparisons between NAb sequences using Similarity Network Fusion. 
- Perform HDBSCAN and Spectral Clustering; hyperparameters (i.e. cluster size) are chosen by the user. 
- Prints interactive graphs following clustering analysis, print operation logs easily accesible by the user.

The purpose of this design is to enable a user to use known NAb - to antigen -  sequence / metadata records and observe how **incomplete NAb records** (no specific antigen / lacking antigen metadata) cluster alongside the known examples. Done so, through mathematically inferred affinity between NAb's, the aim is to draw conclusions on unknown NAb valencies and with additional statistical tests, figure out the most influential determinants in NAb valency. 

# Pipeline GUI

On the welcome page, you are greeted with a request to upload an example `.csv` (ex. protein sequences with some metadata of your project), this example is scanned by the system which then provides you a list of data columns that you can include in your analysis. Importantly, to save time and effort for the user, instead of formatting your `.csv` file to match the style used by the platform's backend, your selections on the dropdown list will send instructions to the system on how to parse your project's `.csv` files.

### Welcome Page
![welcome_page](./Images/welcome_page_early.png)

Navigating the platform from the sidebar on the left, you can select the "*Clustering Analysis*" to show the central analysis page. On this page you have dropdown blocks, on which you can:

- **Process Raw Data**: Offers the user to designate their folder / directory where their project `.csv` files are stored, for the system to access and process automatically once the user commences analysis by clicking "*Execute Pipeline*". 
    - Note that the system can detect reruns, if the user introduces new files to their directory, select the "Rerun" option.
- **Separate CTD/BLOSUM Descriptor Analysis**: Compare and contrast the capabilities of BLOSUM based pairwise comparison between your NAbs and the CTD strategy developed for this project. You can perform a PROCRUSTES comparison between CTD and BLOSUM networks, as well as a Mantel test - comparing the separate CTD and BLOSUM affinity matrices.

Leveraging techniques used in Support Vector Machines for understanding protein sequences, i.e. CTD statistics, were used distinguish proteins and gain a low level "understanding" of their morphology and reactivity - without need for 3D data from ex. X-Ray diffraction. By generating normalised pairwise distance matrices of CDR3H, CTD statistics and BLOSUM45 comparisons, I convert the distances to affinities, then perform similarity network fusion to conjoin CTD and BLOSUM pairwise affinities to generate a single fused matrix for clustering analysis. All of this, wrapped in a Python based - local host server - graphical interface, will allow a user to cluster (HDBSCAN / Spectral) their CDR3H sequences, with parsing, computation, storage and visualisation tasks fully automated by an Object Oriented backend. Therefore, instead of imputing CDR3H valencies, one can observe the clusters to which CDR3H sequences with no known antigen belong to and draw conclusions from there.

### Projection of distance matrix to space
### Clustering upon distance projection
![cluster_example](./Images/clustering_early.png)

# Pipeline CLI (Under development)

Designed to be used in a server environment, the CLI version is an argsparse based alternative to the GUI version which wraps the Object-Oriented backend much the same. Although some features are kept absent (i.e. data visualisation, CTD-BLOSUM statistic comparison) due to a lack of time to implement them / debugging.

Issue commands with the following logic (UNIX):

`python pipeline_cli [option] --{argument} \path`

Observe the figure below for a detailed explanation of the pipeline's capabilities and inner-workings:

![pipeline_img](./Images/pipeline_outline.png)

# Overcoming Challenges: Note to the reader

Initially, the CLI interface was developed to faithfully reproduce the project's experimentation process, by automating data parsing, performing biophysical descriptor computations, data storage, SQL data injection and linking to the ML models tested. 

The data generated by the CLI platform were used to train machine learning models - aimed to predict the relationship between an antibody's characteristics and those of a corresponding antigen. Initial modelling experiments involved a `Multi-Layer Perceptron` with two hidden layers of 100 neurons each, using a `Leaky ReLU` activation function, and a linear activation function for the output layer. Further testing was conducted with a stacked regression model, combinining a `HistGradientBoostingRegressor` and `RidgeCV` as base estimators, with a final `RidgeCV` model working as a meta-regressor, aimed at learning how to impute antigenic sequence/biophysical descriptors from just antibody data. Additional tests with non-parametric models were attempted, but none could achieve a similar training ability (accounting for overfitting), such as the stacked regressor explained prior.

As it later became clear, the scarcity of publicly available data unreservedly hindered ML model training; with efforts to identify an optimal model that can handle learning from a dataset size of 1.2k records, proving futile. As such, any further experimentation was stopped, prompting for the project's methodology to be changed such that: the data processing pipeline becomes a centre point of development, not only automating every step of data preprocessing but also tackling data analysis itself (i.e. replacing standalone ML approaches). 

**Therefore, the objective shifted from imputing epitopes with high specificty, to a more pragmatic approach which seeks to approximate antibody epitopes through clustering techniques, overcoming the issue of restricted data availability, but accompanied with its own limitations.**