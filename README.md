Demixed Subspace Principal Component Analysis (dsPCA)
===========================================

dsPCA is a dimensionality reduction algorithm for high-dimensional data such as neural population activity. dsPCA describes high-dimensional data by finding best demixed linear coding axes for targeted variables and finding linear dimensions that concicely describe the remaining variance in the data. The demixed coding axes for targeted variables are identified with a supervised approach, and then the axes that capture the remaining variance in the data are identified with an unsupervised approach. Importantly, dsPCA completely removes the linear coding of target signals from data before the unsupervised axis identifications. Therefore, dsPCA decomposes a high-dimensional data into a subspace that fully captures all target signals in the data and a subspace that is free of all target signals.


## Installation
will make it available from PyPl

## How to use
dsPCA code is written for Python 3.x    
Please check an example demo in dsPCA_demo.ipynb for the basic implementation.

After installation, import dsPCA
~~~~
from dspca import dsPCA
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
- Projections of population activity to demixed target signal axes ([Trials] x [target dimensions]). The dimension type is in the order of dQ, Qch, sQ.

projection_targetfree_subspace
- Projections of population activity to the axes of target-free subspace ([Trials] x [target-free dimensions]). The dimensions are ordered according to the amount of explained activity variance.

ax_targets
- Axis vectors for the target signal subspace.

ax_targetfree
- Axis vectors for the target-free signal subspace.

target_subspace_signal
- Pearson correlation coefficient between the activity along each target axis and the targeted task-related variables. The 1st dimension indicates target axis type, the 2nd dimension indicates target variables, the 3rd dimension specifies the correlation coefficient (0) or the p-value (1).

targetfree_subspace_signal
- Pearson correlation coefficient between the activity along each target-free axis and the targeted task-related variables. The 1st dimension indicates target-free axes, the 2nd dimension indicates target variables, the 3rd dimension specifies the correlation coefficient (0) or the p-value (1).

target_subspace_var
- Activity variance along each taraget axis

targetfree_subspace_var
- Activity variance along each taraget-free axis

total_var
- Total activity variance of the original input data (activity_mean)

dot_target_ax
- Matrix with dot products between pairs of target axis vectors.
