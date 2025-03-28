with import <nixpkgs> {};

pkgs.mkShell {
  buildInputs = with pkgs; [
    xorg.libX11
    libGL
    gcc
    cmake
    gfortran
    gsl
    doxygen
    pkgs.texlive.combined.scheme-full
    pkgs.biber
    sfml
    eigen
    fmt
  ];

  packages = [
  (pkgs.python3.withPackages(python-pkgs: [
    python-pkgs.jupyter
    python-pkgs.ipykernel
    python-pkgs.numpy
    python-pkgs.matplotlib
    python-pkgs.scipy
    python-pkgs.pandas
    python-pkgs.sympy
  ]))
  ];

  shellHook = ''
    echo 'This shell.nix file has been made for use in the Hydrodynamics and Holomorphicity project.' 
    '';

}
