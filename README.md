Demixed Subspace Principal Component Analysis (dsPCA)
===========================================

dsPCA is a dimensionality reduction algorithm for high-dimensional data such as neural population activity. dsPCA describes high-dimensional data by finding best demixed linear coding axes for targeted variables and finding linear dimensions that concicely describe the remaining variance in the data. The demixed coding axes for targeted variables are identified with a supervised approach, and then the axes that capture the remaining variance in the data are identified with an unsupervised approach. Importantly, dsPCA completely removes the linear coding of target signals from data before the unsupervised axis identifications. Therefore, dsPCA decomposes a high-dimensional data into a subspace that fully captures all target signals in the data and a subspace that is free of all target signals.
The details of the algorithm can be found in **[Hattori and Komiyama, _Neuron_, 2021]**

## Installation
Run the following command to install the package.
~~~~
pip install dspca
~~~~

## How to use
dsPCA code is written for Python 3.x    
Please check an example demo in **dsPCA_demo.ipynb** for the basic implementation.

After installation, import dsPCA
~~~~
from dspca.dspca import dsPCA
~~~~
Provide data ([# of observations] X [# of dimensions to be reduced]) and targets ([# of observations] X [# of targets]) 
~~~~
dsPCA(data=data, targets=targets)
~~~~
The outputs of the function are in the order of
~~~~
projection_target_subspace
projection_targetfree_subspace
ax_targets, ax_targetfree
target_subspace_signal
targetfree_subspace_signal
target_subspace_var
targetfree_subspace_var
total_var
dot_target_ax
~~~~
projection_target_subspace
- Projection of data onto demixed target signal axes.

projection_targetfree_subspace
- Projection of population activity to the axes of target-free subspace.

ax_targets
- Axis vectors for the target signal subspace.

ax_targetfree
- Axis vectors for the target-free signal subspace.

target_subspace_signal
- Pearson correlation coefficient between the projected data along each target axis and the targeted task-related variables.

targetfree_subspace_signal
- Pearson correlation coefficient between the projected data along each target-free axis and the targeted task-related variables.

target_subspace_var
- Data variance along each taraget axis.

targetfree_subspace_var
- Data variance along each taraget-free axis.

total_var
- Total data variance of the original input data.

dot_target_ax
- Matrix with dot products between pairs of target axis vectors.

## Citation
Please cite the following paper:  
- Hattori, R. and Komiyama, T. (2021). Context-dependent persistency as a coding mechanism for robust and widely distributed value coding. _Neuron_, doi: 10.1016/j.neuron.2021.11.001.