# Pasteur_Internship

## Preamble
This project analyses the polyreactive behaviour of Natural Antibodies (NAbs) to understand their epitopic valency based on a combination of Amino Acid (AA) and biophysical descriptors. Publicly available antibody and corresponding antigen sequence data were used to compute biophysical descriptors (i.e. Natural Antibody Database). These descriptors were then prepared to train machine learning models to predict the relationship between an antibody's characteristics and those of a corresponding antigen Initial modelling experiments involved a `Multi-Layer Perceptron` (MLP) with two hidden layers of 100 neurons each, using a `Leaky ReLU` activation function, and a linear activation function for the output layer. Further testing was conducted with a stacked regression model, combinining a `HistGradientBoostingRegressor` and `RidgeCV` as base estimators, with a final `RidgeCV` model working as a meta-regressor, aimed at learning how to impute antigenic sequence/biophysical descriptors from just antibody data. However, the scarcity of publicly available data was evident during ML training performance validation and no amount of tweaking could fix that issue. As such, any further experimentation was stopped, prompting for the project's methodology to be changed such that: the data processing pipeline becomes a centre point of development, not only automating every step of data preprocessing but also tackling data analysis itself (i.e. replacing standalone ML approaches). 

**Therefore, the objective shifts from predicting epitopes with high specificty, to a more pragmatic approach which seeks to approximate antibody epitopes through clustering techniques, overcoming the issue of restricted data availability.**

# Pipeline GUI

In this version you are greeted with a welcome on which you can select a number of columns from your source data files; columns to include in the platform's standardised analysis format. The column to the left of your webpage will give you the option to shift between the different components of the platform: namely, between the welcome page-issuing instructions-and the analysis dashboard on which all of the analysis takes place. 

### Welcome Page
![welcome_page](./Images/welcome_page_early.png)
### Projection of distance matrix to space
![umap_projection](./Images/distance_comp_projection.png)
### Clustering upon distance projection
![cluster_example](./Images/clustering_early.png)

# Pipeline CLI

Designed to be used in a server environment if needed, the CLI version was the first version of this piepline, developed to faithfully reproduce the project's experimentation process. Nonetheless, the design patterns and object oriented approach enable this application to be expanded or modified to one's needs.

Issue commands with the following logic (UNIX):

`python pipeline_cli [option] --{argument} \path`

Observe the figure below for a detailed explanation of the pipeline's capabilities and inner-workings:

![pipeline_img](./Images/pipeline_outline.png)

