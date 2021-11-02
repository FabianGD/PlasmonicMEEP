{ lib, python3Packages, buildPythonApplication, nix-gitignore
, joblib, h5py-mpi, matplotlib, meep, numpy, pandas, pyyaml
, pytestCheckHook, openssh, additionalDevDeps ? [ ] }:

buildPythonApplication rec {
    pname = "plasmonic-meep";
    version = "0.5.1";
    src = nix-gitignore.gitignoreSource [ ] ../.;

    nativeBuildInputs = additionalDevDeps;

    propagatedBuildInputs = [
        joblib
        h5py-mpi
        matplotlib
        meep
        numpy
        pandas
        pyyaml
    ];

    doCheck = true;

    checkInputs = [
        pytestCheckHook

        # OpenSSH is necessary because of MEEP's MPI
        openssh
    ];
    pytestFlagsArray = [ "tests/" ];

    meta = with lib; {
        description = "Set of scripts for calculation of plasmon resonance/electric field enhancement on different structures";
        licence = lib.licences.gpl3Only;
        maintainers = [ maintainers.fabiangd ];
        homepage = "https://gitlab.com/theoretical-chemistry-jena/quantum-dynamics/plasmonic-meep";
    };
}
