{ pkgs ? import ./nix/pkgs.nix
, matplotlibrc ? import (import ./nix/sources.nix).matplotlibrc { inherit pkgs; } }:

with pkgs; {
  production = mkShell {
    buildInputs = [
      # The custom-linked python
      (import ./default.nix { meep = qchem.python3.pkgs.meep; })

      which
      git
      # SSH is required for running with MPI
      openssh
      matplotlibrc
    ];

    additionalShellHook = ''
      export MATPLOTLIBRC=${matplotlibrc}
    '';
  };

  # This allows for a local editable(!) dev install
  dev = pkgs.python3Packages.callPackage ./nix/plasmonic-meep.nix {
    meep = qchem.python3.pkgs.meep;
    additionalDevDeps = with python3Packages; [
      mpi
      openssh

      ipykernel
      black
      pylint
      ipympl
      pytest

      matplotlibrc
    ];
    additionalShellHook = ''
      export MATPLOTLIBRC=${matplotlibrc}
    '';
  };
}

