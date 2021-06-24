{ pkgs ? import ./nix/pkgs.nix }:

let

in with pkgs;
  {
    production = mkShell {
      buildInputs = [
        # The custom-linked python
        ((import ./nix/default.nix {}).plasmonic-meep)

        which
        git
        # SSH is required for running with MPI
        openssh
      ];
    };

    # This allows for a local editable(!) dev install
    dev = pkgs.python3Packages.buildPythonPackage rec {

      name = "plasmonicmeep";
      src = ./.;
      version = "dev";

      nativeBuildInputs = [
        which
        git
        openssh
        mpi

        python3Packages.ipykernel
        python3Packages.black
        python3Packages.pylint
      ];

      propagatedBuildInputs = with pkgs.python3Packages; [

        pkgs.qchem.python3.pkgs.meep
        numpy
        scipy
        matplotlib
        h5py-mpi
        joblib
        pandas

        ipympl
      ];
    };
  }


