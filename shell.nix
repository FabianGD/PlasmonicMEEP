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
    dev = pkgs.python3Packages.callPackage ./nix/plasmonic-meep.nix {
      meep = pkgs.qchem.python3.pkgs.meep;
      additionalDevDeps = with python3Packages; [
        mpi
        openssh

        ipykernel
        black
        pylint
        ipympl
      ];
    };
  }


