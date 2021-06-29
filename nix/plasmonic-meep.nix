{ lib, buildPythonApplication, nix-gitignore
# Python deps
, numpy, scipy, matplotlib, h5py-mpi, joblib, meep, pandas
# Optional dependencies for development
, additionalDevDeps ? [ ] }:

buildPythonApplication rec {
    pname = "plasmonic-meep";
    version = "0.4.0";
    src = nix-gitignore.gitignoreSource [ ] ../.;

    nativeBuildInputs = additionalDevDeps ++ [ ];

    propagatedBuildInputs = [
        numpy
        scipy
        matplotlib
        pandas

        h5py-mpi
        meep

        joblib
    ];

  doCheck = false;

    meta = {
        description = "Set of scripts for calculation of plasmon resonance/electric field enhancement on different structures";
        licence = lib.licences.gpl3Only;
        homepage = "https://gitlab.com/theoretical-chemistry-jena/quantum-dynamics/plasmonic-meep";
    };
}
