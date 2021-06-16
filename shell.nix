{ pkgs ? import ./nix/pkgs.nix }:

let
  customPython = pkgs.python3.withPackages(ps: with ps; [
    pkgs.qchem.meep
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
            git
            customPython
        ];
    }
