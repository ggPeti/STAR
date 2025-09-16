{
  description = "Nix flake for building the STAR aligner from local sources";

  inputs = {
    nixpkgs.url = "github:NixOS/nixpkgs/nixos-unstable";
  };

  outputs = { self, nixpkgs }:
    let
      systems = [
        "x86_64-linux"
        "aarch64-linux"
        "x86_64-darwin"
        "aarch64-darwin"
      ];
      forAllSystems = f: nixpkgs.lib.genAttrs systems (system: f system);
    in {
      packages = forAllSystems (system:
        let
          pkgs = import nixpkgs { inherit system; };
        in {
          # Expose the STAR package, building from the local repository
          star = pkgs.callPackage ./nix/package.nix { src = ./source; };

          # Default package: STAR
          default = self.packages.${system}.star;
        });

      apps = forAllSystems (system:
        let
          pkg = self.packages.${system}.star;
        in {
          # Expose an app so `nix run` works
          star = {
            type = "app";
            program = "${pkg}/bin/STAR";
          };

          default = self.apps.${system}.star;
        });

      devShells = forAllSystems (system:
        let
          pkgs = import nixpkgs { inherit system; };
        in {
          # Simple dev shell with the tools needed to build via Makefile
          default = pkgs.mkShell {
            packages = with pkgs; [
              gnumake
              xxd
              zlib
            ];
          };
        });
    };
}
