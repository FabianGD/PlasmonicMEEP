"""
Calculate FDTD using MEEP of a plasmonic nanostructure
"""

import argparse
import sys
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt

import meep as mp
import numpy as np
from meep import materials

from .utils import append_attrs
from .model import create_sinus_grating, two_nps


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
        "--sizex",
        type=float,
        default=1.0,
        help="The size of the box in µm in the x direction.",
    )
    parser.add_argument(
        "-y",
        "--sizey",
        type=float,
        default=1.0,
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
        "-f",
        "--frequency",
        type=float,
        default=1.5,
        help=(
            "Change the central frequency of the incident laser field."
            "Frequency is given in units of 1/µm."
        ),
    )
    parser.add_argument(
        "-w",
        "--freq-width",
        type=float,
        default=1.5,
        help=(
            "Change the frequency width of the incident laser field."
            "Pulse with is given in units of µm."
        ),
    )

    parser.add_argument(
        "-g",
        "--show-geometry",
        action="store_true",
        help="Plot only geometry and exit.",
    )
    parser.add_argument("-o", "--output", default="data/", help="Output folder.")
    parser.add_argument(
        "--sinus",
        action="store_true",
        help="Create the sinus grating geometry for benchmarking.",
    )
    parser.add_argument(
        "--spectrum",
        "-s",
        action="store_true",
        help=(
            "At the end of the simulation, plot and save transmission/reflectance/loss spectra."
            "This might fail for MPI runs or on clusters."
        ),
    )

    return parser.parse_args()


def main():
    """
    Main computation
    """

    args = argparsing()

    # Set up output path, if not already existent.
    output_path = Path(args.output).resolve()
    if not output_path.is_dir():
        output_path.mkdir(parents=True, exist_ok=True)
    output = str(output_path)

    # Inner size
    sizex = args.sizex
    sizey = args.sizey

    # Computational grid resolution in pixels per distance unit
    resolution = args.resolution
    pml_th = 30 / resolution

    # Full system size including PMLs
    fullx = sizex + 2 * pml_th
    fully = sizey + 2 * pml_th

    # The system cell and the waveguide material
    cell = mp.Vector3(fullx, fully, 0)
    mat = materials.Au

    # Methods
    geom = args.show_geometry
    saveref = True

    # Frequency in  1 / µm. Speed of light c == 1
    # Size dimension a --> Time dimension a / c

    # Frequency has dim [c/a]
    cfreq = args.frequency
    # Frequency width of the Gaussian pulse
    fwidth = args.freq_width  # in units of µm

    # Field component to monitor
    comp = mp.Hz

    # List of sources, needed for the simulation
    sources = [
        mp.Source(
            mp.GaussianSource(frequency=cfreq, fwidth=fwidth),
            size=mp.Vector3(0, fully, 0),
            component=comp,
            center=mp.Vector3(-sizex / 2, 0),
        )
    ]

    # Absorber on grating side because of field divergence at metal/pml interface
    # pml_layers = [mp.PML(pml_th, direction=mp.X), mp.Absorber(pml_th, direction=mp.Y)]
    pml_layers = [mp.PML(pml_th, direction=mp.ALL)]

    # empty cell for reference run
    geometry = []
    sim = mp.Simulation(
        cell_size=cell,
        boundary_layers=pml_layers,
        geometry=geometry,
        sources=sources,
        resolution=resolution,
        split_chunks_evenly=False,
        force_complex_fields=True,
    )

    sim.use_output_directory(output)
    prefix = sim.get_filename_prefix()
    # Define monitors for further spectra calculation
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

    #######################
    # Reference calculation
    #######################

    dset = "ref"

    if geom:
        sim.init_sim()

    else:
        if saveref:
            sim.run(
                mp.to_appended(
                    dset,
                    mp.at_every(
                        0.5 / (cfreq + fwidth * 0.5),
                        mp.output_efield_x,
                        mp.output_efield_y,
                        mp.output_efield_z,
                    ),
                ),
                mp.to_appended(
                    dset + "-center",
                    mp.in_point(
                        mp.Vector3(),
                        mp.output_efield_x,
                        mp.output_efield_y,
                        mp.output_efield_z,
                    ),
                ),
                until=20
            )
        else:
            sim.run(until=20)

    # for normalization run, save flux fields data for reflection plane
    straight_refl_data = sim.get_flux_data(refl)
    # save incident power for transmission plane
    straight_tran_flux = mp.get_fluxes(tran)

    if mp.am_master():
        # This appends the calc attributes to the HDF5 file.
        append_attrs(
            output=output_path, prefix=prefix, dset=dset, cfreq=cfreq, fwidth=fwidth
        )

    sim.reset_meep()

    ######################
    # Material calculation
    ######################

    dset = "norm"

    if geom:
        # overwrite mat to be correctly displayed
        # plot2D function does not handle correctly
        # dispersive materials
        mat = mp.Medium(epsilon=5)

    if args.sinus:
        metal_vert = create_sinus_grating(
            ampl=0.1,
            periodicity=0.5,
            thickness=0.04,
            resolution=resolution,
            sizex=fully,
            y=True,
        )
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
        geometry = two_nps(
            radius=0.05, separation=0.005, center=mp.Vector3(), material=mat, y=True
        )

    sim = mp.Simulation(
        cell_size=cell,
        boundary_layers=pml_layers,
        geometry=geometry,
        sources=sources,
        resolution=resolution,
        split_chunks_evenly=False,
        force_complex_fields=True,
    )
    sim.use_output_directory(output)

    # same monitors
    refl = sim.add_flux(cfreq, fwidth, nfreq, reflectance_fr)
    tran = sim.add_flux(cfreq, fwidth, nfreq, transmittance_fr)
    sim.load_minus_flux_data(refl, straight_refl_data)

    print("Starting main run")

    if geom:
        sim.plot2D(labels=True)

        if mp.am_master():
            plt.show()

        sys.exit()

    sim.run(
        mp.to_appended(
            dset,
            mp.at_every(
                0.5 / (cfreq + fwidth * 0.5),
                mp.output_efield_x,
                mp.output_efield_y,
                mp.output_efield_z,
            ),
        ),
        mp.to_appended(
            dset + "-center",
            mp.in_point(
                mp.Vector3(),
                mp.output_efield_x,
                mp.output_efield_y,
                mp.output_efield_z,
            ),
        ),
        until=20
    )

    if mp.am_master():
        # This appends the calc attributes to the HDF5 file.
        append_attrs(
            output=output_path, prefix=prefix, dset=dset, cfreq=cfreq, fwidth=fwidth
        )

    refl_flux = np.asarray(mp.get_fluxes(refl))
    tran_flux = np.asarray(mp.get_fluxes(tran))
    straight_tran_flux = np.asarray(straight_tran_flux)
    flux_freqs = np.asarray(mp.get_flux_freqs(refl))

    wavelengths = 1 / flux_freqs
    refl_spectrum = -refl_flux / straight_tran_flux
    tran_spectrum = tran_flux / straight_tran_flux
    loss_spectrum = 1 - refl_spectrum - tran_spectrum

    if args.spectrum and mp.am_master():
        mpl.use("agg")

        fig, ax = plt.subplots()

        ax.plot(wavelengths, refl_spectrum, color="blue", label="reflectance")
        ax.plot(wavelengths, tran_spectrum, color="red", label="transmittance")
        ax.plot(wavelengths, loss_spectrum, color="green", label="loss")

        ax.set_xlabel("wavelength / μm")
        ax.legend(loc="upper right")
        ax.grid(True)

        fig.tight_layout()
        fig.savefig(output_path / "ReflTransLoss.pdf", dpi=300)

        if args.show_spectra:
            plt.show()


if __name__ == "__main__":
    main()
