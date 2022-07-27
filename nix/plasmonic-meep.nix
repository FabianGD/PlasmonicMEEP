{ lib, buildPythonPackage, nix-gitignore, joblib, h5py-mpi
, matplotlib, meep, numpy, pandas, pyyaml, pytestCheckHook, openssh, tqdm
, additionalDevDeps ? [ ], additionalShellHook ? "" }:

buildPythonPackage rec {
  pname = "plasmonic-meep";
  version = "0.5.1";
  src = nix-gitignore.gitignoreSource [ ] ../.;

  nativeBuildInputs = additionalDevDeps;

  propagatedBuildInputs =
    [ joblib h5py-mpi matplotlib meep numpy pandas pyyaml tqdm ];

  doCheck = true;

  checkInputs = [
    pytestCheckHook

    # OpenSSH is necessary because of MEEP's MPI
    openssh
  ];
  pytestFlagsArray = [ "tests/" ];

  postShellHook = additionalShellHook;

  meta = with lib; {
    description =
      "Set of scripts for calculation of plasmon resonance/electric field enhancement on different structures";
    licence = licences.gpl3Only;
    maintainers = [ maintainers.fabiangd ];
    homepage =
      "https://gitlab.com/theoretical-chemistry-jena/quantum-dynamics/plasmonic-meep";
  };
}
