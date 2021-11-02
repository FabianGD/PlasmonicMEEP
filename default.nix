{ pkgs ? import ./nix/pkgs.nix, meep, additionalDevDeps ? [ ] }:

with pkgs;
python3.pkgs.callPackage ./nix/plasmonic-meep.nix {
  inherit meep additionalDevDeps;
}
