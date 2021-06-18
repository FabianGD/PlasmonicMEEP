{ pkgs ? import ./nix/pkgs.nix }:

let
  customPython = pkgs.python3.withPackages(ps: with ps; [
    pkgs.qchem.meep
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

    pycairo
    pygobject3

  ]);
in
    with pkgs; mkShell {
        buildInputs = [
            git

            cairo

            gobjectIntrospection
            gtk3-x11

            customPython
        ];
    }
