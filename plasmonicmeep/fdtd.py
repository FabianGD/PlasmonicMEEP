"""
Calculate FDTD using MEEP of a plasmonic nanostructure
"""

import sys
import logging
from datetime import datetime
from pathlib import Path
from typing import Any, List

import matplotlib as mpl
import matplotlib.pyplot as plt

import meep as mp
import numpy as np
from meep import materials

from .utils import append_attrs, append_timegrid
from .model import two_nps
from .args import fdtd_argparsing


def gen_stepfuncs(
    dset: str,
    write_entire_cell: bool = True,
    write_single_point: bool = True,
    cell_timestep: float = 0.2,
    point: mp.Vector3 = mp.Vector3(),
) -> List[Any]:
    """Generate step functions to be used for MEEP runs.

    Args:
        dset (str): Data set name
        write_entire_cell (bool, optional): Whether to generate a step function that prints
            the field in the entire cell. Defaults to True.
        write_single_point (bool, optional): Whether to generate a step function that prints
            the field only at a specific point. Defaults to True.
        cell_timestep (float, optional): The timestep used for printing the entire field,
            in MEEP units of time. Defaults to 0.2.
        point (mp.Vector3, optional): The point in configuration space to read the field
            from. Defaults to mp.Vector3().

    Raises:
        ValueError: In case that both 'write_entire_cell' and 'write_single_point' are False.

    Returns:
        List[Any]: List of MEEP step functions.
    """

    step_functions = []

    if not write_entire_cell and not write_single_point:
        raise ValueError("You need to output something.")

    if write_entire_cell:
        step_functions.append(
            mp.to_appended(
                dset,
                mp.at_every(
                    cell_timestep,
                    mp.output_efield_x,
                    mp.output_efield_y,
                    mp.output_efield_z,
                ),
            )
        )

    if write_single_point:
        step_functions.append(
            mp.to_appended(
                dset + "-center",
                mp.in_point(
                    point,
                    mp.output_efield_x,
                    mp.output_efield_y,
                    mp.output_efield_z,
                ),
            )
        )

    return step_functions


def main():
    """
    Main computation
    """

    args = fdtd_argparsing()

    # TODO Possibility to change log levels.
    logging.basicConfig(level="INFO")

    # This will be used as a subfolder in the output folder
    subfolder = "_".join([datetime.now().strftime("%Y-%m-%dt%H-%M"), args.calc_name])

    # Set up output path, if not already existent.
    output_path = Path(args.output).expanduser().resolve() / subfolder
    if not output_path.is_dir():
        output_path.mkdir(parents=True, exist_ok=True)
    output = str(output_path)

    # Inner size
    sizex = args.sizex
    sizey = args.sizey

    # Computational grid resolution in pixels per distance unit
    resolution = args.resolution
    pml_th = 30 / resolution

    timestep = 0.5 / resolution
    n_timesteps = int(args.run_until / timestep)

    logging.info(
        "Calculated time step based on Courant factor and res: {}. "
        "Assuming MEEP wants to perform {:.0f} timesteps.".format(
            timestep, n_timesteps
        )
    )

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

    # Calculate the cell timestep
    if args.cfield_time_step == "auto":
        cell_timestep = 0.5 / (cfreq + fwidth * 0.5)
    else:
        cell_timestep = float(args.cfield_time_step)

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
        force_complex_fields=args.enable_complex,
        filename_prefix=args.calc_name,
    )

    sim.use_output_directory(output)

    # TODO: Prefix bug fix
    # Fix upstream bug that resets the prefix. The solution here saves us from
    # misery, however.
    prefix = sim.get_filename_prefix()

    if args.spectrum:
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
            center=mp.Vector3(sizex / 2 - mon_skip, 0, 0),
            size=mp.Vector3(0, mon_height, 0),
        )

        # Add the flux regions to the simulation
        refl = sim.add_flux(cfreq, fwidth, nfreq, reflectance_fr)
        tran = sim.add_flux(cfreq, fwidth, nfreq, transmittance_fr)

    logging.info("Starting reference run")

    #######################
    # Reference calculation
    #######################

    dset = "ref"

    if geom:
        sim.init_sim()

    else:
        if saveref:
            sim.run(
                *gen_stepfuncs(
                    dset,
                    write_entire_cell=not args.disable_cell_field,
                    write_single_point=not args.disable_single_point,
                    cell_timestep=cell_timestep,
                ),
                until=20,
            )
        else:
            sim.run(until=20)

    if args.spectrum:
        # for normalization run, save flux fields data for reflection plane
        straight_refl_data = sim.get_flux_data(refl)
        # save incident power for transmission plane
        straight_tran_flux = mp.get_fluxes(tran)

    if not args.disable_cell_field and mp.am_master():
        # This appends the calc attributes to the HDF5 file.
        append_attrs(
            folder=output_path, prefix=prefix, dset=dset, cfreq=cfreq, fwidth=fwidth
        )

    if not args.disable_single_point and mp.am_master():
        append_timegrid(
            folder=output_path,
            prefix=prefix,
            dset="-".join([dset, "center"]),
            n_timesteps=n_timesteps,
            timestep=timestep,
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
        force_complex_fields=args.enable_complex,
    )
    sim.use_output_directory(output)

    if args.spectrum:
        # same monitors
        refl = sim.add_flux(cfreq, fwidth, nfreq, reflectance_fr)
        tran = sim.add_flux(cfreq, fwidth, nfreq, transmittance_fr)
        sim.load_minus_flux_data(refl, straight_refl_data)

    logging.info("Starting main run")

    if geom:
        sim.plot2D(labels=True)

        if mp.am_master():
            plt.show()

        sys.exit()

    sim.run(
        *gen_stepfuncs(
            dset,
            write_entire_cell=not args.disable_cell_field,
            write_single_point=not args.disable_single_point,
            cell_timestep=cell_timestep,
        ),
        until=20,
    )

    if not args.disable_cell_field and mp.am_master():
        # This appends the calc attributes to the HDF5 file.
        append_attrs(
            folder=output_path, prefix=prefix, dset=dset, cfreq=cfreq, fwidth=fwidth
        )

    if not args.disable_single_point and mp.am_master():
        append_timegrid(
            folder=output_path,
            prefix=prefix,
            dset="-".join([dset, "center"]),
            n_timesteps=n_timesteps,
            timestep=timestep,
        )

    if args.spectrum:
        refl_flux = np.asarray(mp.get_fluxes(refl))
        tran_flux = np.asarray(mp.get_fluxes(tran))
        straight_tran_flux = np.asarray(straight_tran_flux)
        flux_freqs = np.asarray(mp.get_flux_freqs(refl))

        wavelengths = 1 / flux_freqs
        refl_spectrum = -refl_flux / straight_tran_flux
        tran_spectrum = tran_flux / straight_tran_flux
        loss_spectrum = 1 - refl_spectrum - tran_spectrum

        if mp.am_master():
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


if __name__ == "__main__":
    main()
