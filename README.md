## Pasteur_Internship

This project analyses the polyreactive behaviour of Natural Antibodies (NAbs) to understand their epitopic valency based on a combination of Amino Acid (AA) and biophysical descriptors. Publicly available antibody and corresponding antigen sequence data were used to compute biophysical descriptors (i.e. Natural Antibody Database). These descriptors were then prepared to train machine learning models to predict the relationship between an antibody's characteristics and those of a corresponding antigen. Initial modelling experiments involved a Multi-Layer Perceptron (MLP) with two hidden layers of 100 neurons each, using a `Leaky ReLU` activation function, and a linear activation function for the output layer. Further testing was conducted with a stacked regression model, combinining a `HistGradientBoostingRegressor` and `RidgeCV` as base estimators, with a final `Ridge` model working as a meta-regressor, aimed at learning how to impute antigenic sequence/biophysical descriptors from just antibody data. However, the scarcity of publicly available data was evident during ML training performance validation and no amount of tweaking could fix that issue. As such, any further experimentation was stopped, prompting for the project's methodology to be changed such that: the data processing pipeline becomes a centre point of development, not only automating every step of the experimentation until ML training, but also tackling data analysis itself (i.e. replacing standalone ML approaches). 

Clustering algorithms (ex. Agglomerative) will be implemented in length (currently somewhat useful), to identify clonotypes in antibody data by means of analysing the outcoming clusters. 

# Pipeline GUI

More information will be added soon...otherwise try running the pipelineGUI file after cloning the repository and installing requirements..

# Pipeline CLI

Explained in internship report.
