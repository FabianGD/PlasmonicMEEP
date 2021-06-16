{ lib, buildPythonApplication, numpy, scipy, matplotlib, h5py, joblib, meep, pandas }:

buildPythonApplication rec {
    pname = "plasmonic-meep";
    version = "0.1.0";
    src = lib.cleanSource ../.;

    propagatedBuildInputs = [
        numpy
        scipy
        matplotlib
        pandas

        h5py
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
