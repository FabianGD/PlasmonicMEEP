{ lib, buildPythonApplication, nix-gitignore, openssh
# Python deps
, numpy, scipy, matplotlib, h5py-mpi, joblib, meep, pandas, pyyaml, pytestCheckHook
# Optional dependencies for development
, additionalDevDeps ? [ ] }:

buildPythonApplication rec {
    pname = "plasmonic-meep";
    version = "0.5.0";
    src = nix-gitignore.gitignoreSource [ ] ../.;

    nativeBuildInputs = additionalDevDeps ++ [ ];

    propagatedBuildInputs = [
        joblib
        h5py-mpi
        matplotlib
        meep
        numpy
        pandas
        scipy
        pyyaml
    ];

    doCheck = true;

    checkInputs = [
        pytestCheckHook

        # OpenSSH is necessary because of MEEP's MPI
        openssh
    ];
    pytestFlagsArray = [ "tests/" ];

    meta = {
        description = "Set of scripts for calculation of plasmon resonance/electric field enhancement on different structures";
        licence = lib.licences.gpl3Only;
        homepage = "https://gitlab.com/theoretical-chemistry-jena/quantum-dynamics/plasmonic-meep";
    };
}
