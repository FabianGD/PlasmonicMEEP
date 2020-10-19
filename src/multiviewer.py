"""
Slice viewer implemented using matplotlib.
From https://www.datacamp.com/community/tutorials/matplotlib-3d-volumetric-data

TODO!

"""


import matplotlib.pyplot as plt

AXIS = -1
COLORBAR = None
IMG = None
IDX_FCT = None


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


def multi_slice_viewer(volume, index_function=lambda x: x):
    global IMG
    global IDX_FCT
    IDX_FCT = index_function
    remove_keymap_conflicts({"j", "k"})
    fig, ax = plt.subplots()
    ax.volume = volume
    ax.index = volume.shape[AXIS] // 2
    IMG = ax.imshow(volume[:, :, ax.index])
    fig.canvas.mpl_connect("key_press_event", process_key)


def process_key(event):
    global COLORBAR
    fig = event.canvas.figure
    ax = fig.axes[0]

    if COLORBAR is not None:
        COLORBAR.remove()

    if event.key == "j":
        previous_slice(ax)

    elif event.key == "k":
        next_slice(ax)

    COLORBAR = plt.colorbar(IMG)
    # data = ax.images[0].get_array()

    IMG.autoscale()
    fig.canvas.draw()


def previous_slice(ax):
    volume = ax.volume
    ax.index = (ax.index - 1) % volume.shape[AXIS]  # wrap around using %
    ax.images[0].set_array(volume[:, :, ax.index])
    ax.set_title(IDX_FCT(ax.index))


def next_slice(ax):
    volume = ax.volume
    ax.index = (ax.index + 1) % volume.shape[AXIS]
    ax.images[0].set_array(volume[:, :, ax.index])
    ax.set_title(IDX_FCT(ax.index))
