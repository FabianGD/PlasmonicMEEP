{ pkgs ? import ./nix/pkgs.nix }:

let

in with pkgs;
  {
    production = mkShell {
      buildInputs = [
        # The custom-linked python
        (import ./default.nix {
          meep = qchem.python3.pkgs.meep;
        }).plasmonic-meep

        which
        git
        # SSH is required for running with MPI
        openssh
      ];
    };

    # This allows for a local editable(!) dev install
    dev = pkgs.python3Packages.callPackage ./. {
      meep = qchem.python3.pkgs.meep;
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


