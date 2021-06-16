{ pkgs ? import ./nix/pkgs.nix }:

let
    plasmonic-meep = (import ./nix/default.nix {}).plasmonic-meep;

in
    with pkgs; mkShell {
        buildInputs = [
            git

            python3Packages.black
            python3Packages.pylint

            plasmonic-meep

        ];
    }
