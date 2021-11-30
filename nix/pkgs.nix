let
  sources = import ./sources.nix;

  qchem = import sources.nixos-qchem;

  nixpkgs = import sources.nixpkgs {
    overlays = [ qchem ];

    config = {
      allowUnfree = true;
      qchem-config = {
        allowEnv = false;
        optAVX = true;
        srcurl = "https://troja.ipc3.uni-jena.de/nix-src";
      };
    };
  };

in nixpkgs
