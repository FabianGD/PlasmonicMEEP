{ pkgs ? import ./nix/pkgs.nix }:

let
  devPython = pkgs.python3.withPackages(ps: with ps; [
    pkgs.qchem.python3.pkgs.meep
    numpy
    scipy
    matplotlib
    h5py
    pyyaml
    attrs
    joblib
    black
    pylint
    ipykernel
    ipympl
  ]);

  plasmonicPython = (import ./nix/default.nix {}).plasmonic-meep;

in
  {
    develop =
      with pkgs; mkShell {
        buildInputs = [
          which
          git

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

