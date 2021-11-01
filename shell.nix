{ pkgs ? import ./nix/pkgs.nix }:

let

in with pkgs;
  {
    production = mkShell {
      buildInputs = [
        # The custom-linked python
        (import ./default.nix {}).plasmonic-meep

        which
        git
        # SSH is required for running with MPI
        openssh
      ];
    };

    # This allows for a local editable(!) dev install
    dev = pkgs.python3Packages.callPackage ./. {
      qchemPkgs = pkgs.qchem;
      additionalDevDeps = with python3Packages; [
        mpi
        openssh

        ipykernel
        black
        pylint
        ipympl
        pytest
      ];
    };
  }


