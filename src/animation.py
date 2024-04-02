import matplotlib.pyplot as plt
from matplotlib import pyplot as plt
from matplotlib.animation import FuncAnimation
import numpy as np
import sys
import time
from read_data import *
from defaults import *
from plot_functions import *


def animate(folder_path: str, save: bool = False):
    # Create a figure
    fig, ax = plt.subplots()

    # Read data
    num_frames = get_num_frames(folder_path)

    # set data
    times, data_blocks_u, data_blocks_v, data_blocks_T, data_blocks_S, data_blocks_p = set_data(
        folder_path)

    FPS = int(num_frames / duration) + 1  # +1 to ensure positivity
    interval = 1000 / FPS  # interval between frames in milliseconds

    data_blocks = set_datablocks(
        plot_var,
        data_blocks_u,
        data_blocks_v,
        data_blocks_T,
        data_blocks_S,
        data_blocks_p)

    # Print some information
    print("Length of data: ", num_frames)
    print("Number of frames: ", num_frames)
    print("Interval: ", interval)
    print("Expected duration: ", duration)

    # Create data
    nx, ny = set_nx_ny(data_blocks)
    X, Y = set_axis(dx, dy, LX, LY, nx, ny)

    Z = data_blocks[0]
    Z_MAX, Z_MIN = set_Z_max_min(data_blocks)

    color = get_color(plot_var)

    plot_args = get_plot_args(dx, LX, LY, Z_MIN, Z_MAX, color)
    plot_args_quiver = get_plot_args_quiver()

    text_args = {'x': 0.5, 'y': y_pos_text, 's': '', 'transform': ax.transAxes,
                 'fontsize': FONTSIZE_TIME, 'horizontalalignment': 'center'}
    time_text = ax.text(**text_args)

    # set title (omega or sqrt(u^2 + v^2)
    ax.set_title(get_title(plot_var), y=y_pos_title)

    # Axes
    ax.set_xlabel(X_LABEL)
    ax.set_ylabel(Y_LABEL)

    # Create plot
    # add quiver plot with u and v with 1/numth of the points
    num = 16
    uu, vv, XX, YY = data_blocks_u[0, ::num, ::num], data_blocks_v[0, ::num,
                                                                   ::num], X[::num, ::num], Y[::num, ::num]
    # uuvv = np.sqrt(uu**2 + vv**2)
    # uu, vv = uu/uuvv, vv/uuvv
    YY = YY[::-1]
    plot = [ax.imshow(Z, **plot_args),
            ax.quiver(XX, YY, uu, vv, **plot_args_quiver)]

    # Add a color bar which maps values to colors.
    fig.colorbar(plot[0], **colorbar_args)

    # Animation update function
    def init():
        return plot[0], time_text

    def update(frame):
        # Clear the previous frame
        # ax.clear()
        plot[0].remove()
        plot[1].remove()

        real_frame = int(frame)

        # Update the arrays
        Z = data_blocks[real_frame]
        uu, vv = data_blocks_u[real_frame, ::num,
                               ::num], data_blocks_v[real_frame, ::num, ::num]

        # Update the time text
        time_text.set_text('t = %.3f' % times[real_frame])

        # Update the plot
        plot[0] = ax.imshow(Z, **plot_args)
        plot[1] = ax.quiver(XX, YY, uu, vv, **plot_args_quiver)

        return plot[0], plot[1], time_text

    # Create the animation
    ani = FuncAnimation(fig, update, frames=num_frames,
                        interval=interval, blit=False, init_func=init)

    end_time = time.time()
    print("Total time for animating: ", int(
        UNIT_TIME * (end_time - start_time)), LABEL_TIME)

    if save:
        # Save the animation
        filename = 'data/videos/' + plot_var + '_Pr=' + str(Pr) + '_Le=' + str(
            Le) + '_Ra_T=' + str(Ra_T) + '_R_rho=' + str(R_rho) + '.mp4'
        ani.save(filename, writer='ffmpeg', dpi=300, fps=FPS)
    else:
        # Show the animation
        plt.show()


# count time
UNIT_TIME = 1000  # in seconds
LABEL_TIME = "ms"
start_time = time.time()
folder_results = 'output/results/'
folder_object = 'config/'

# Read data
NX, NY, LX, LY, Pr, Le, Ra_T, R_rho, dx, dy, dt, Tfinal, animation_on, plot_dt, plot_var = readSetupFromHDF5(
    folder_results + '../setup.h5')

if animation_on == False:
    print("Animation is off. Exiting...")
    sys.exit()

# if there is an argument, then save the animation
if len(sys.argv) > 1:
    save = True
    animate(folder_results, save)
else:
    animate(folder_results)
