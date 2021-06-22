{ pkgs ? import ./nix/pkgs.nix }:

let
  devPython = pkgs.python3.withPackages(ps: with ps; [
    pkgs.qchem.python3.pkgs.meep
    numpy
    scipy
    matplotlib
    h5py
    joblib
    black
    pylint
    ipykernel
    ipympl
    pyyaml
    pint
    pycairo
    pygobject3

  ]);

  plasmonicPython = (import ./nix/default.nix {}).plasmonic-meep;

in
  {
    develop =
      with pkgs; mkShell {
        buildInputs = [
          which
          git

          gobject-introspection
          gtk3-x11
          # SSH is required for running with MPI
          openssh

          # The custom-linked python
          devPython
        ];
    };
    production =
      with pkgs; mkShell {
        buildInputs = [
          which
          git

          # SSH is required for running with MPI
          ssh

          # The custom-linked python
          plasmonicPython
        ];
    };
  }

