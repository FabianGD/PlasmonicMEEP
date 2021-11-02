{ pkgs ? import ./nix/pkgs.nix, meep, additionalDevDeps ? [] }:

with pkgs; {
    plasmonic-meep = python3.pkgs.callPackage ./nix/plasmonic-meep.nix {
        inherit meep additionalDevDeps;
    };
}