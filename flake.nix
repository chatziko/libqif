{
  description = "libqif development environment";

  inputs = {
    nixpkgs.url = "github:NixOS/nixpkgs/nixos-unstable";
  };

  outputs = { self, nixpkgs }:
    let
      systems = [ "x86_64-linux" "aarch64-linux" "x86_64-darwin" "aarch64-darwin" ];
      forAllSystems = f: nixpkgs.lib.genAttrs systems (system: f nixpkgs.legacyPackages.${system});
    in {
      devShells = forAllSystems (pkgs:
        let
          # sphinx must run on the same python as the qif module (docs import it via autosummary).
          # Use the default python3: or-tools propagates python3 packages into PYTHONPATH.
          python = pkgs.python3.withPackages (ps: with ps; [
            numpy
            sphinx
          ]);

          # ortoolsConfig.cmake calls find_dependency on these. Take them from or-tools' own inputs
          # so that the versions match the ones it was built with (eg it pins older abseil/protobuf).
          ortoolsDeps = pkgs.lib.filter
            (d: builtins.elem (d.pname or "") [ "zlib" "bzip2" "abseil-cpp" "protobuf" "re2" "cbc" "highs" "eigen" "scipopt-scip" ])
            (pkgs.or-tools.buildInputs ++ pkgs.or-tools.propagatedBuildInputs);

          # or-tools also needs GLPK, but the released packages are built without it (the official or-tools binaries
          # don't need it). So we pass it only via GLPK_ROOT, which or-tools' FindGLPK uses but libqif's find_library doesn't.
          glpk = pkgs.lib.findFirst (d: (d.pname or "") == "glpk") null pkgs.or-tools.buildInputs;
        in {
          default = pkgs.mkShell {
            # build tools (the C/C++ compiler comes from mkShell's stdenv)
            nativeBuildInputs = with pkgs; [
              cmake
              pkg-config            # ortoolsConfig.cmake finds Clp via pkg-config
              python
              act                   # run the github actions locally
            ];

            buildInputs = with pkgs; [
              gmp
              gsl
              openblasCompat        # BLAS/LAPACK for armadillo. Not openblas: it uses 64-bit ints, armadillo expects 32-bit
              or-tools
            ] ++ ortoolsDeps;

            GLPK_ROOT = glpk;
          };
        });
    };
}
