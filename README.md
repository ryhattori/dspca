Demixed Subspace Principal Component Analysis (dsPCA)
===========================================

dsPCA is a dimensionality reduction algorithm for high-dimensional data such as neural population activity. dsPCA describes high-dimensional data by finding best demixed linear coding axes for targeted variables and finding linear dimensions that concicely describe the remaining variance in the data. The demixed coding axes for targeted variables are identified with a supervised approach, and then the axes that capture the remaining variance in the data are identified with an unsupervised approach. Importantly, dsPCA completely removes the linear coding of target signals from data before the unsupervised axis identifications. Therefore, dsPCA decomposes a high-dimensional data into a subspace that fully captures all target signals in the data and a subspace that is free of all target signals.


## Installation
will make it available from PyPl

## How to use
dsPCA_demo.ipynb describes an example of the basic implementation.

