import numpy as np
from matplotlib import colormaps as cm
from matplotlib.colors import ListedColormap
import os
from read_data import read_data_object
from matplotlib.patches import Circle, Rectangle
from defaults import *


def get_color(plot_var):
    redblue = cm['RdBu_r']
    only_red = ListedColormap(redblue(np.linspace(0.5, 1, 128)))
    color = only_red if plot_var == 'uv' else redblue
    return color


def get_title(plot_var):
    if plot_var == 'u':
        return r'$u$'
    elif plot_var == 'v':
        return r'$v$'
    elif plot_var == 'uv':
        return r'$\sqrt{u^2 + v^2}$'
    elif plot_var == 'T':
        return 'Temperature'
    elif plot_var == 'S':
        return 'Salinity'
    elif plot_var == 'p':
        return 'Pressure'


def set_nx_ny(data_blocks):
    nx = data_blocks.shape[2]
    ny = data_blocks.shape[1]
    return nx, ny


def set_axis(dx, dy, LX, LY, nx, ny):
    X = np.linspace(dx / 2 * (1 + cols_removed_beginning),
                    LX - dx / 2 * (1 + cols_removed_end), nx)
    Y = np.linspace(dy / 2, LY - dy / 2, ny)
    X, Y = np.meshgrid(X, Y, indexing='xy')
    return X, Y


def get_plot_args(dx, LX, LY, Z_MIN, Z_MAX, color):
    plot_args = {
        'cmap': color,
        'extent': [
            dx * cols_removed_beginning,
            LX - dx * cols_removed_end,
            0,
            LY],
        'vmin': Z_MIN,
        'vmax': Z_MAX,
        'interpolation': 'spline16'}
    return plot_args


def get_plot_args_quiver():
    plot_args = {
        'color': 'black',
        # 'scale': 0.01}
        'width': 0.002}
    return plot_args


def set_datablocks(plot_var, data_blocks_u, data_blocks_v,
                   data_blocks_T, data_blocks_S, data_blocks_p):
    if plot_var == 'u':
        data_blocks = data_blocks_u
    if plot_var == 'v':
        data_blocks = data_blocks_v
    if plot_var == 'uv':
        data_blocks = np.sqrt(data_blocks_u**2 + data_blocks_v**2)
    elif plot_var == 'T':
        data_blocks = data_blocks_T
    elif plot_var == 'S':
        data_blocks = data_blocks_S
    elif plot_var == 'p':
        data_blocks = data_blocks_p
    return data_blocks


def set_Z_max_min(data_blocks):
    Z_MAX = np.max(data_blocks)
    Z_MIN = np.min(data_blocks)
    # if w_on:
    #     Z_MAX = max(Z_MAX, -Z_MIN) / 4
    #     Z_MIN = -Z_MAX
    return Z_MAX, Z_MIN
