import numpy as np
import h5py
import os
from defaults import *

# def read_data_file(file_path):
#     """Read data from file. The file is of the form:
#     header_0
#     x_00 x_01 ... x_0m
#     x_10 x_11 ... x_1m
#     ...
#     x_n0 x_n1 ... x_nm

#     header_1
#     x_00 x_01 ... x_0m
#     x_10 x_11 ... x_1m
#     ...
#     x_n0 x_n1 ... x_nm

#     ...

#     header_k
#     x_00 x_01 ... x_0m
#     x_10 x_11 ... x_1m
#     ...
#     x_n0 x_n1 ... x_nm
#     """
#     headers = []
#     data_blocks = []

#     with open(file_path, 'r') as file:
#         lines = file.readlines()

#     i = 0
#     while i < len(lines):
#         line = lines[i].strip()
#         if line:
#             header = float(line)
#             headers.append(header)

#             i += 1
#             block_lines = []
#             while i < len(lines) and lines[i].strip():
#                 aux = list(map(float, lines[i].split()))
#                 block_lines.append(aux)
#                 i += 1
#             data_block = np.array(block_lines)
#             data_blocks.append(data_block)
#         i += 1

#     headers_array = np.array(headers)
#     data_blocks_array = np.array(data_blocks)

#     return headers_array, data_blocks_array


def get_num_frames(folder_path):
    num_frames = len([f for f in os.listdir(folder_path)
                      if os.path.isfile(os.path.join(folder_path, f))])
    return num_frames


def read_data_object(file_path):
    """
    Read data from file. The file is of the form:

    A1 x0
    A2 x1
    ...
    An xn
    """
    with open(file_path, 'r') as file:
        lines = file.readlines()

    data = []
    for line in lines:
        header, value = line.split()
        data.append(float(value))

    data_array = np.array(data)

    return data_array


def readSetupFromHDF5(hdf5_file):
    with h5py.File(hdf5_file, 'r') as file:
        # Accessing the first element and converting to scalar
        NX = int(file['NX'][()].item())
        NY = int(file['NY'][()].item())
        LX = file['L'][()].item()
        LY = file['H'][()].item()
        Pr = file['Pr'][()].item()
        Le = file['Le'][()].item()
        Ra_T = file['Ra_T'][()].item()
        R_rho = file['R_rho'][()].item()
        dx = file['dx'][()].item()
        dy = file['dy'][()].item()
        dt = file['dt'][()].item()
        Tfinal = file['Tfinal'][()].item()
        animation_on = bool(file['animation_on'][()].item())
        plot_var = file['plot_var'][()].item().decode('utf-8')
    return NX, NY, LX, LY, Pr, Le, Ra_T, R_rho, dx, dy, dt, Tfinal, animation_on, plot_var


def readSolutionFromHDF5(hdf5_file):
    with h5py.File(hdf5_file, 'r') as file:
        # NX = int(file['NX'][()].item())
        # NY = int(file['NY'][()].item())
        u = file['u'][:]
        v = file['v'][:]
        T = file['T'][:]
        S = file['S'][:]
        p = file['p'][:]
        t = file['t'][()].item()
    return t, u, v, T, S, p


def set_data(folder_path):
    num_frames = get_num_frames(folder_path)

    times = []
    data_blocks_u = []
    data_blocks_v = []
    data_blocks_T = []
    data_blocks_S = []
    data_blocks_p = []
    for i in range(num_frames):
        file_path = folder_path + 'sol_' + str(i) + '.h5'
        t, u, v, T, S, p = readSolutionFromHDF5(file_path)
        times.append(t)
        # remove ghost cells
        data_blocks_u.append(u[1:-1, 1:-1])
        data_blocks_v.append(v[1:-1, 1:-1])
        data_blocks_T.append(T[1:-1, 1:-1])
        data_blocks_S.append(S[1:-1, 1:-1])
        data_blocks_p.append(p[1:-1, 1:-1])

    times = np.array(times)
    data_blocks_u = np.array(data_blocks_u)
    data_blocks_v = np.array(data_blocks_v)
    data_blocks_T = np.array(data_blocks_T)
    data_blocks_S = np.array(data_blocks_S)
    data_blocks_p = np.array(data_blocks_p)

    arrays_to_process = [
        data_blocks_u,
        data_blocks_v,
        data_blocks_T,
        data_blocks_S,
        data_blocks_p]

    # Iterate over each array and remove columns
    for i, data_array in enumerate(arrays_to_process):
        # arrays_to_process[i] = np.delete(data_array,
        #                                   np.s_[:cols_removed_beginning],
        #                                   axis=1)
        # arrays_to_process[i] = np.delete(arrays_to_process[i],
        #                                   np.s_[-cols_removed_end:],
        #                                   axis=1)

        tmp = arrays_to_process[i].copy()

        arrays_to_process[i] = np.zeros(
            (len(tmp), tmp.shape[2], tmp.shape[1]))

        # transpose the data and reverse
        for j, data in enumerate(tmp):
            arrays_to_process[i][j] = data.T[::-1]

    # Update the original arrays with the modified ones
    data_blocks_u, data_blocks_v, data_blocks_T, data_blocks_S, data_blocks_p = arrays_to_process

    return times, data_blocks_u, data_blocks_v, data_blocks_T, data_blocks_S, data_blocks_p
