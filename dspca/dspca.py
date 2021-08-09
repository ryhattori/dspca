""" Demixed subspace Principal Component Analysis (dsPCA) """

# Author: Ryoma Hattori <rhattori0204@gmail.com>
#
# License: MIT License

from sklearn import linear_model
from sklearn.decomposition import PCA
from scipy.stats import pearsonr
import numpy as np
import numpy.linalg as la

def dsPCA(data, targets):
    model = linear_model.LinearRegression(n_jobs=-1, fit_intercept=True)
    n_target = targets.shape[1]
    n_dim_trial = data.shape[0]
    n_dim_realign = data.shape[1]

    # Fit multiple linear regression models to data
    coef_targets = np.zeros((n_dim_realign, n_target))
    for cell_id in range(n_dim_realign):
        coef_targets[cell_id, :] = model.fit(targets, data[:, cell_id]).coef_.squeeze()

    # Coding axes for the target signals
    ax_targets = np.zeros_like(coef_targets)
    for target_id in range(n_target):
        ax_targets[:, target_id] = coef_targets[:, target_id] / np.linalg.norm(coef_targets[:, target_id], ord=None)

    # Dot products between paired target axes
    dot_target_ax = np.zeros((n_target, n_target))
    for target1_id in range(n_target):
        for target2_id in range(n_target):
            dot_target_ax[target1_id, target2_id] = np.dot(ax_targets[:, target1_id], ax_targets[:, target2_id])

    # Projection of data onto each target axis
    projection_target_subspace = np.zeros((n_dim_trial, n_target))
    for target_id in range(n_target):
        projection_target_subspace[:, target_id] = np.matmul(data, ax_targets[:, target_id])

    # Projection of data onto target-free subspace
    ax_targets_temp = ax_targets.copy()
    data_for_projection = data.copy()
    ax_targetfree = np.identity(n_dim_realign)
    for target_id in range(n_target):
        qr_q, qr_r = np.linalg.qr(ax_targets_temp[:, target_id, None], mode='complete')
        projected_ax = np.matmul(ax_targets_temp.T, qr_q)
        projected = np.matmul(data_for_projection, qr_q)
        ax_targetfree = np.matmul(ax_targetfree, qr_q[:, 1:])
        data_for_projection = projected[:, 1:]
        ax_targets_temp = projected_ax[:, 1:].T

    # Re-arrange target-free axes with PCA, get target-free axes
    projection_targetfree_subspace_prepca = data_for_projection
    if n_dim_trial < (n_dim_realign - n_target):
        print('Trial number '+str(n_dim_trial)+' is smaller than (Dimension - # of targets) of '+str(n_dim_realign - n_target))
        print('# of dimensions for target-free subspace is set to the trial number '+str(n_dim_trial))
        pca_targetfree = PCA(n_components=n_dim_trial)
    else:
        pca_targetfree = PCA(n_components=(n_dim_realign - n_target))
    pca_targetfree.fit(projection_targetfree_subspace_prepca)
    projection_targetfree_subspace = np.matmul(projection_targetfree_subspace_prepca, pca_targetfree.components_.T)
    ax_targetfree = np.matmul(ax_targetfree, pca_targetfree.components_.T)

    # Variance along each axis
    target_subspace_var = np.var(projection_target_subspace, axis=0)
    targetfree_subspace_var = np.var(projection_targetfree_subspace, axis=0)
    total_var = np.sum(np.var(data, axis=0))

    # Pearson correlation coefficients between the data along each axis and the target signals
    target_subspace_signal = np.zeros((n_target, n_target, 2))
    targetfree_subspace_signal = np.zeros((projection_targetfree_subspace.shape[1], n_target, 2))
    for target_var_id in range(n_target):
        for target_ax_id in range(n_target):
            target_subspace_signal[target_ax_id, target_var_id, 0], target_subspace_signal[target_ax_id, target_var_id, 1] = pearsonr(projection_target_subspace[:, target_ax_id], targets[:, target_var_id])
        for target_ax_id in range(projection_targetfree_subspace.shape[1]):
            targetfree_subspace_signal[target_ax_id, target_var_id, 0], targetfree_subspace_signal[target_ax_id, target_var_id, 1] = pearsonr(projection_targetfree_subspace[:, target_ax_id], targets[:, target_var_id])

    return projection_target_subspace, projection_targetfree_subspace, ax_targets, ax_targetfree, target_subspace_signal, targetfree_subspace_signal, target_subspace_var, targetfree_subspace_var, total_var, dot_target_ax

