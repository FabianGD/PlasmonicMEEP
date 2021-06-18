{ pkgs ? import ./nix/pkgs.nix }:

let
  customPython = pkgs.python3.withPackages(ps: with ps; [
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
  ]);
in
    with pkgs; mkShell {
        buildInputs = [
          which
          git
          customPython
        ];
    }
