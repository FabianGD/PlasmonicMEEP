{ pkgs ? import ./pkgs.nix }:

with pkgs; {
    plasmonic-meep = python3.pkgs.callPackage ./plasmonic-meep.nix {
        meep = qchem.python3.pkgs.meep;
    };
}