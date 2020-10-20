"""
Calculate FDTD using MEEP of a plasmonic nanostructure
"""

import argparse
import sys

import matplotlib.pyplot as plt
import meep as mp
import numpy as np
from meep import materials

import model


def argparsing():
    """
    Argument parser
    """
    parser = argparse.ArgumentParser(
        description=(
            "Main FDTD module. Computes spectra and densities of a given model."
            "May be run in parallel using MPI."
        )
    )
    parser.add_argument(
        "-x",
        "--sx",
        type=float,
        default=1.0,
        help="The size of the box in µm in the x direction.",
    )
    parser.add_argument(
        "-y",
        "--sy",
        type=float,
        default=2.0,
        help="The size of the box in µm in the y direction.",
    )
    parser.add_argument(
        "-r",
        "--resolution",
        type=int,
        default=150,
        help="The resolution of the box (nr. of pixels per µm?)",
    )
    parser.add_argument(
        "-g",
        "--show-geometry",
        action="store_true",
        help="Plot only geometry and exit.",
    )
    parser.add_argument(
        "-o",
        "--output",
        default="data/",
        help="Output folder."
    )
    parser.add_argument(
        "--sinus",
        action="store_true",
        help="Create the sinus grating geometry for benchmarking."
    )

    return parser.parse_args()


def vec3_to_nparray(vec):
    """
    Utility to convert a Vector3 to a (3)-shaped np.ndarray
    """
    return np.asarray([vec.x, vec.y, vec.z])


def main():
    """
    Main computation
    """

    args = argparsing()

    # Inner size
    sizex = args.sx
    sizey = args.sy

    # Computational grid resolution in pixels per distance unit
    resolution = args.resolution
    pml_th = 30 / resolution

    # Full system size including PMLs
    fullx = sizex + 2 * pml_th
    fully = sizey + 2 * pml_th

    # The system cell and the waveguide material
    cell = mp.Vector3(fullx, fully, 0)
    mat = materials.Au

    #####################
    # run modes
    # geom - show geometry instead of actual simulation
    # saveref - save the reference file, needed for field enhancement calculations
    ####################

    geom = args.show_geometry
    saveref = True

    ####################
    # geometry, load AFM profile
    ####################

    metal_vert = model.create_sinus_grating(
        ampl=0.1,
        periodicity=0.5,
        thickness=0.04,
        resolution=resolution,
        sizex=fully,
        y=True,
    )

    ####################
    # source paraneters
    ####################

    # Frequency in  1 / µm. Speed of light c == 1
    # Size dimension a --> Time dimension a / c

    # Frequency has dim [c/a]
    cfreq = 1.5
    # Frequency width of the Gaussian pulse
    fwidth = 1.5

    # Field component to monitor
    comp = mp.Hz

    sources = [
        mp.Source(
            mp.GaussianSource(frequency=cfreq, fwidth=fwidth),
            size=mp.Vector3(0, fully, 0),
            component=comp,
            center=mp.Vector3(-sizex / 2, 0),
        )
    ]

    # Absorber on grating side because of field divergence at metal/pml interface
    pml_layers = [mp.PML(pml_th, direction=mp.X), mp.Absorber(pml_th, direction=mp.Y)]

    # empty cell for reference run
    geometry = []
    sim = mp.Simulation(
        cell_size=cell,
        boundary_layers=pml_layers,
        geometry=geometry,
        sources=sources,
        resolution=resolution,
        filename_prefix=args.output + "data",
        split_chunks_evenly=False,
    )

    # define monitors for further spectra calculation
    mon_height = sizey
    nfreq = 200

    # Small skip from PML layers
    mon_skip = sizex / 20

    # Flux regions are measurement areas that record the FT at the position at each time step.
    reflectance_fr = mp.FluxRegion(
        center=mp.Vector3(-sizex / 2 + mon_skip, 0, 0),
        size=mp.Vector3(0, mon_height, 0),
    )
    transmittance_fr = mp.FluxRegion(
        center=mp.Vector3(sizex / 2 - mon_skip, 0, 0), size=mp.Vector3(0, mon_height, 0)
    )

    # Add the flux regions to the simulation
    refl = sim.add_flux(cfreq, fwidth, nfreq, reflectance_fr)
    tran = sim.add_flux(cfreq, fwidth, nfreq, transmittance_fr)

    # The point records the field strength over time and gives the cancellation criterion
    point = mp.Vector3(sizex / 3, 0, 0)

    print("Starting reference run")

    #!######################
    #! Reference calculation
    #!######################

    if geom:
        sim.init_sim()
    else:
        if saveref:
            sim.run(
                mp.to_appended(
                    "ref",
                    mp.at_every(
                        0.5 / (cfreq + fwidth * 0.5),
                        mp.output_efield_x,
                        mp.output_efield_y,
                        mp.output_efield_z,
                    ),
                ),
                until_after_sources=mp.stop_when_fields_decayed(10, comp, point, 1e-2),
            )
        else:
            sim.run(
                until_after_sources=mp.stop_when_fields_decayed(10, comp, point, 1e-2)
            )

    # for normalization run, save flux fields data for reflection plane
    straight_refl_data = sim.get_flux_data(refl)
    # save incident power for transmission plane
    straight_tran_flux = mp.get_fluxes(tran)

    sim.reset_meep()

    #!#####################
    #! Material calculation
    #!#####################

    if geom:
        # overwrite mat to be correctly displayed
        # plot2D function does not handle correctly
        # dispersive materials
        mat = mp.Medium(epsilon=5)

    if args.sinus:
        geometry = [
            mp.Prism(
                metal_vert,
                height=100,
                center=mp.Vector3(0, 0, 0),
                axis=mp.Vector3(0, 0, 1),
                material=mat,
            )
        ]
    else:
        geometry = [
            *model.create_np_on_mirror(
                0.05, 0.005, mat, fullx, fully, y=True
            )
        ]

    sim = mp.Simulation(
        cell_size=cell,
        boundary_layers=pml_layers,
        geometry=geometry,
        sources=sources,
        resolution=resolution,
        filename_prefix=args.output + "data",
        split_chunks_evenly=False,
    )

    # same monitors
    refl = sim.add_flux(cfreq, fwidth, nfreq, reflectance_fr)
    tran = sim.add_flux(cfreq, fwidth, nfreq, transmittance_fr)
    sim.load_minus_flux_data(refl, straight_refl_data)

    print("Starting main run")

    if geom:
        vertices = np.asarray([vec3_to_nparray(vert) for vert in metal_vert])
        x, y, _ = vertices.T

        sim.plot2D(labels=True, fields=mp.Hz)

        if mp.am_master():
            plt.plot(x, y)
            plt.show()

        sys.exit()

    sim.run(
        mp.to_appended(
            "norm",
            mp.at_every(
                0.5 / (cfreq + fwidth * 0.5),
                mp.output_efield_x,
                mp.output_efield_y,
                mp.output_efield_z,
            ),
        ),
        until_after_sources=mp.stop_when_fields_decayed(10, comp, point, 1e-2),
    )

    refl_flux = np.asarray(mp.get_fluxes(refl))
    tran_flux = np.asarray(mp.get_fluxes(tran))
    straight_tran_flux = np.asarray(straight_tran_flux)

    flux_freqs = np.asarray(mp.get_flux_freqs(refl))

    wavelengths = 1 / flux_freqs
    refl_spectrum = -refl_flux / straight_tran_flux
    tran_spectrum = tran_flux / straight_tran_flux
    loss_spectrum = 1 - refl_spectrum - tran_spectrum

    if mp.am_master():
        plt.figure()
        plt.plot(wavelengths, refl_spectrum, "b-", label="reflectance")
        plt.plot(wavelengths, tran_spectrum, "r-", label="transmittance")
        plt.plot(wavelengths, loss_spectrum, "go-", label="loss")
        # plt.axis([5.0, 10.0, 0, 1])
        plt.xlabel("wavelength (μm)")
        plt.legend(loc="upper right")
        plt.grid(True)

        plt.savefig("spectra.svg")
        plt.show()


if __name__ == "__main__":
    main()
