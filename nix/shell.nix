{ pkgs ? import ./pkgs.nix }:

let
    plasmonic-meep = (import ./default.nix {}).plasmonic-meep;

in
    with pkgs; mkShell {
        buildInputs = [
            git

            python3Packages.pygobject
            python3Packages.black
            python3Packages.pylint

            plasmonic-meep

        ];
    }