# Research project ideas

## Goals

- Overarching: Model plasmonic nanostructures using MEEP
  - Field enhancements, phase shifts etc. for small (and tiny!) cavities

## Todos

1. Introduce yourself to the literature of surface plasmons, nanostructure and FDTD
2. Get to know electrodynamics FDTD (finite difference time domain) simulations
3. Find analytical expressions for different particle geometries
4. Numerically simulate a range of geometries and extract field responses
5. If applicable, compare the numerical results with the theory
6. (Couple MEEP to pyQD to do coupled ED-QD simulations)

## Introductory problems (Write-your-own-code)

### Analytical solutions to the Mie scattering problem on spheres

We would like to compare the numerical scattering profile (which we will calculate using MEEP)
to the Mie-Lorenz scattering theory. Therefore it would be great to write the analytical solutions
to the Mie scattering problem of a plane wave with frequency $\omega_P$ on simple shapes such as a sphere with a radius $r$ in code.

> Whats the formula? --> May be too complex to start with.

### Build custom geometries

Start building custom geometries using MEEP, both based on `mp.Sphere`, `mp.Block`, etc.
and using custom geometries based on vertices (represented as `mp.Vector3`'s) and the
`mp.Prism` class.

Ideas:

- Two "wires" besides each other. This should simulate a micro- or nanocavity.
- Two triangles |> <| or <| |> facing.
- Put a sphere on a "surface", such as a slab of material.

### Calculate transmission and reflection spectra

Get to know the concept of flux regions and monitors.
For the geometries above, try to insert monitors, define flux regions and, lastly,
calculate transmission spectra using this setup. You might need a reference and normal run.
The MEEP code located in `./src/fdtd.py` has a similar setup.

### Find out how to calculate field enhancements

Use different types of materials, (gold, silicon dioxide, semiconductors, etc.)

## Literature

1. Ribeiro, R. F.; Martínez-Martínez, L. A.; Du, M.; Campos-Gonzalez-Angulo, J.; Yuen-Zhou, J.
   Polariton Chemistry: Controlling Molecular Dynamics with Optical Cavities. Chem. Sci. 2018,
   9 (30), 6325–6339. [DOI](https://doi.org/10.1039/C8SC01043A).
2. Herrera, F.; Owrutsky, J. Molecular Polaritons for Controlling Chemistry with Quantum Optics.
   J. Chem. Phys. 2020, 152 (10), 100902. [DOI](https://doi.org/10.1063/1.5136320).
3. [MEEP documentation](https://meep.readthedocs.io/en/latest/)
   - Especially the basic tutorials found
     [here](https://meep.readthedocs.io/en/latest/Python_Tutorials/Basics/)
4. [GitLab repo](https://gitlab.com/theoretical-chemistry-jena/quantum-dynamics/plasmonic-meep)
   with Python code and this document
5. Video introductions
   - [Andrey Bogdanov: Introduction to Green’s functions & scattering theory. Mie theory. Part 1.](https://www.youtube.com/watch?v=bUfcVTIlJz0)
   - [Andrey Bogdanov: Mie theory. Part 2](https://youtu.be/CaFyJRN_iYI)
   - [Kristina Frizyuk: Kerker effect. Mie theory. Part 3.](https://youtu.be/JblWhmOexy4)
