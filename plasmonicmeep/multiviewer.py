"""
Slice viewer implemented using matplotlib.
From https://www.datacamp.com/community/tutorials/matplotlib-3d-volumetric-data
"""

from functools import partial
import matplotlib.pyplot as plt


def settitle(ax, fct, patt=None):
    """Sets the title for the specific axis."""
    if not patt:
        patt = "{:.0f} nm"
    ax.set_title(patt.format(fct(ax.index)))


def remove_keymap_conflicts(new_keys_set):
    """
    Removes key stroke conflicts from the rcParams dict.
    """
    for prop in plt.rcParams:
        if prop.startswith("keymap."):
            keys = plt.rcParams[prop]
            remove_list = set(keys) & new_keys_set
            for key in remove_list:
                keys.remove(key)


def multi_slice_viewer(volume, index_function=lambda x: x, resolution=None):
    """
    Entry function
    """

    # global COLORBAR  # pylint: disable=global-statement
    remove_keymap_conflicts({"j", "k"})
    fig, ax = plt.subplots(nrows=1, ncols=1)

    # Add custom attributes to the axis instead of global vars
    ax.volume = volume
    ax.cbar = None
    ax.index = volume.shape[-1] // 2

    ax.set_xlabel("x / px")
    ax.set_ylabel("y / px")

    # Plot the first frame
    img = ax.imshow(volume[:, :, ax.index], origin="lower")
    ax.cbar = plt.colorbar(img, ax=ax)

    # Partially apply function arguments
    processor = partial(process_key, img=img, idx_fct=index_function)
    fig.canvas.mpl_connect("key_press_event", processor)


def process_key(event, img, idx_fct):
    """
    Matplotlib connected callback function.
    """

    fig = event.canvas.figure
    ax = fig.axes[0]

    # Remove the old colorbar
    if ax.cbar is not None:
        ax.cbar.remove()

    if event.key == "k":
        _previous_slice(ax, idxf=idx_fct)

    elif event.key == "j":
        _next_slice(ax, idxf=idx_fct)

    # Build a new colorbar
    ax.cbar = plt.colorbar(img, ax=ax)

    img.autoscale()
    fig.canvas.draw()


def _previous_slice(ax, idxf=lambda x: x):
    """
    Render the previous slice
    """
    volume = ax.volume
    ax.index = (ax.index - 1) % volume.shape[-1]  # wrap around using %
    ax.images[0].set_array(volume[:, :, ax.index])
    settitle(ax=ax, fct=idxf)


def _next_slice(ax, idxf=lambda x: x):
    """
    Render the next slice
    """
    volume = ax.volume
    ax.index = (ax.index + 1) % volume.shape[-1]
    ax.images[0].set_array(volume[:, :, ax.index])
    settitle(ax=ax, fct=idxf)
