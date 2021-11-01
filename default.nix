{ pkgs ? import ./nix/pkgs.nix, qchemPkgs ? pkgs.qchem, additionalDevDeps ? [] }:

with pkgs; {
    plasmonic-meep = python3.pkgs.callPackage ./nix/plasmonic-meep.nix {
        meep = qchemPkgs.python3.pkgs.meep;
        inherit additionalDevDeps;
    };
}