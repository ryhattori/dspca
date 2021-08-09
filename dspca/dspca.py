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

    coef_targets = np.zeros((n_dim_realign, n_target))
    for cell_id in range(n_dim_realign):
        coef_targets[cell_id, :] = model.fit(targets, data[:, cell_id]).coef_.squeeze()

    ax_targets = np.zeros_like(coef_targets)
    for target_id in range(n_target):
        ax_targets[:, target_id] = coef_targets[:, target_id] / np.linalg.norm(coef_targets[:, target_id], ord=None)

    dot_target_ax = np.zeros((n_target, n_target))
    for target1_id in range(n_target):
        for target2_id in range(n_target):
            dot_target_ax[target1_id, target2_id] = np.dot(ax_targets[:, target1_id], ax_targets[:, target2_id])

    projection_target_subspace = np.zeros((n_dim_trial, n_target))
    for target_id in range(n_target):
        projection_target_subspace[:, target_id] = np.matmul(data, ax_targets[:, target_id])

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

    target_subspace_var = np.var(projection_target_subspace, axis=0)
    targetfree_subspace_var = np.var(projection_targetfree_subspace, axis=0)
    total_var = np.sum(np.var(data, axis=0))

    target_subspace_signal = np.zeros((n_target, n_target, 2))
    targetfree_subspace_signal = np.zeros((projection_targetfree_subspace.shape[1], n_target, 2))
    for target_var_id in range(n_target):
        for target_ax_id in range(n_target):
            target_subspace_signal[target_ax_id, target_var_id, 0], target_subspace_signal[target_ax_id, target_var_id, 1] = pearsonr(projection_target_subspace[:, target_ax_id], targets[:, target_var_id])
        for target_ax_id in range(projection_targetfree_subspace.shape[1]):
            targetfree_subspace_signal[target_ax_id, target_var_id, 0], targetfree_subspace_signal[target_ax_id, target_var_id, 1] = pearsonr(projection_targetfree_subspace[:, target_ax_id], targets[:, target_var_id])

    return projection_target_subspace, projection_targetfree_subspace, ax_targets, ax_targetfree, target_subspace_signal, targetfree_subspace_signal, target_subspace_var, targetfree_subspace_var, total_var, dot_target_ax


# def dsPCA(data, targets, time_range):
#     model = linear_model.LinearRegression(n_jobs=-1, fit_intercept=True)
#     n_target = targets.shape[1]
#     n_dim_trial = data.shape[0]
#     n_dim_time = data.shape[1]
#     n_dim_realign = data.shape[2]
#
#     # Average
#     if len(time_range) == 1:
#         data_mean = data[:, time_range[0], :]
#     else:
#         data_mean = np.mean(data[:, time_range, :], axis=1)
#
#     coef_targets = np.zeros((n_dim_realign, n_target))
#     for cell_id in range(n_dim_realign):
#         coef_targets[cell_id, :] = model.fit(targets, data_mean[:, cell_id]).coef_.squeeze()
#
#     ax_targets = np.zeros_like(coef_targets)
#     for target_id in range(n_target):
#         ax_targets[:, target_id] = coef_targets[:, target_id] / np.linalg.norm(coef_targets[:, target_id], ord=None)
#
#     dot_target_ax = np.zeros((n_target, n_target))
#     for target1_id in range(n_target):
#         for target2_id in range(n_target):
#             dot_target_ax[target1_id, target2_id] = np.dot(ax_targets[:, target1_id], ax_targets[:, target2_id])
#
#     projection_targeted_subspace = np.zeros((n_dim_trial, n_dim_time, n_target))
#     for target_id in range(n_target):
#         projection_targeted_subspace[:, :, target_id] = np.matmul(data, ax_targets[:, target_id])
#
#     ax_targets_temp = ax_targets.copy()
#     data_for_projection = data.copy()
#     ax_targetfree = np.identity(n_dim_realign)
#     for target_id in range(n_target):
#         qr_q, qr_r = np.linalg.qr(ax_targets_temp[:, target_id, None], mode='complete')
#         projected_ax = np.matmul(ax_targets_temp.T, qr_q)
#         projected = np.matmul(data_for_projection, qr_q)
#         ax_targetfree = np.matmul(ax_targetfree, qr_q[:, 1:])
#         if len(time_range) == 1:
#             data_mean = projected[:, time_range[0], 1:]
#         else:
#             data_mean = np.mean(projected[:, time_range, 1:], axis=1)
#         data_for_projection = projected[:, :, 1:]
#         ax_targets_temp = projected_ax[:, 1:].T
#
#     projection_targetfree_subspace_prepca = data_for_projection
#     pca_targetfree = PCA(n_components=np.min([n_dim_trial, n_dim_realign - n_target]))
#     pca_targetfree.fit(data_mean)
#     projection_targetfree_subspace = np.matmul(projection_targetfree_subspace_prepca, pca_targetfree.components_.T)
#     ax_targetfree = np.matmul(ax_targetfree, pca_targetfree.components_.T)
#
#     if len(time_range) == 1:
#         target_subspace_var = np.var(projection_targeted_subspace[:, time_range[0], :], axis=0)
#         targetfree_subspace_var = np.var(projection_targetfree_subspace[:, time_range[0], :], axis=0)
#         total_var = np.sum(np.var(data[:, time_range[0], :], axis=0))
#     else:
#         target_subspace_var = np.var(np.mean(projection_targeted_subspace[:, time_range, :], axis=1), axis=0)
#         targetfree_subspace_var = np.var(np.mean(projection_targetfree_subspace[:, time_range, :], axis=1), axis=0)
#         total_var = np.sum(np.var(np.mean(data[:, time_range, :], axis=1), axis=0))
#
#     signal_targeted_ax = np.zeros((n_target, n_target, 2))
#     signal_targetfree_ax = np.zeros((projection_targetfree_subspace.shape[2], n_target, 2))
#     for target_var_id in range(n_target):
#         for target_ax_id in range(n_target):
#             if len(time_range) == 1:
#                 signal_targeted_ax[target_ax_id, target_var_id, 0], signal_targeted_ax[target_ax_id, target_var_id, 1] = pearsonr(projection_targeted_subspace[:, time_range[0], target_ax_id], targets[:, target_var_id])
#             else:
#                 signal_targeted_ax[target_ax_id, target_var_id, 0], signal_targeted_ax[target_ax_id, target_var_id, 1] = pearsonr(np.mean(projection_targeted_subspace[:, time_range, target_ax_id], axis=1), targets[:, target_var_id])
#         for target_ax_id in range(projection_targetfree_subspace.shape[2]):
#             if len(time_range) == 1:
#                 signal_targetfree_ax[target_ax_id, target_var_id, 0], signal_targetfree_ax[target_ax_id, target_var_id, 1] = pearsonr(projection_targetfree_subspace[:, time_range[0], target_ax_id], targets[:, target_var_id])
#             else:
#                 signal_targetfree_ax[target_ax_id, target_var_id, 0], signal_targetfree_ax[target_ax_id, target_var_id, 1] = pearsonr(np.mean(projection_targetfree_subspace[:, time_range, target_ax_id], axis=1), targets[:, target_var_id])
#
#     return projection_targeted_subspace, projection_targetfree_subspace, signal_targeted_ax, signal_targetfree_ax, ax_targets, ax_targetfree, target_subspace_var, targetfree_subspace_var, total_var, dot_target_ax
