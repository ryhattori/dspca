Demixed Subspace Principal Component Analysis (dsPCA)
===========================================

dsPCA decomposes high-dimensional data into low-dimensional orthogonal subspaces where pre-specified signals of interest are demixed into distinct subspaces. For example, when you want to separate signals related to variables A, B, and C in the original data, tsPCA decomposes the data into demixed subspaces A, B, C, and a subspace that is free of the signals related to A, B, and C. When tsPCA is applied to a neural population activity, tsPCA decomposes neural population activity into orthogonal subspaces where targeted task-related signals are demixed. tsPCA can identify subspaces for any types of variables, including continuous, discrete, and categorical variables.

## Installation
will make it available from PyPl

## How to use
dsPCA_demo.ipynb describes an example of the basic implementation.

