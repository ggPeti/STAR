{ lib
, stdenv
, xxd
, zlib
, llvmPackages
, versionCheckHook
, nix-update-script
, src
}:

stdenv.mkDerivation {
  pname = "star";
  version = "2.7.11b";

  inherit src;

  postPatch = ''
    substituteInPlace Makefile --replace-fail "-std=c++11" "-std=c++14"
  '';

  nativeBuildInputs = [ xxd ];

  buildInputs = [ zlib ] ++ lib.optionals stdenv.isDarwin [ llvmPackages.openmp ];

  enableParallelBuilding = true;

  makeFlags = lib.optionals stdenv.hostPlatform.isAarch64 [ "CXXFLAGS_SIMD=" ];

  preBuild = lib.optionalString stdenv.isDarwin ''
    export CXXFLAGS="$CXXFLAGS -DSHM_NORESERVE=0"
  '';

  buildFlags = [
    "STAR"
    "STARlong"
  ];

  installPhase = ''
    runHook preInstall
    install -D STAR STARlong -t $out/bin
    runHook postInstall
  '';

  nativeInstallCheckInputs = [ versionCheckHook ];
  versionCheckProgram = "${placeholder "out"}/bin/STAR";
  versionCheckProgramArg = "--version";
  doInstallCheck = true;

  passthru.updateScript = nix-update-script { };

  meta = with lib; {
    description = "Spliced Transcripts Alignment to a Reference";
    longDescription = ''
      STAR (Spliced Transcripts Alignment to a Reference) is a fast RNA-seq
      read mapper, with support for splice-junction and fusion read detection.
    '';
    mainProgram = "STAR";
    homepage = "https://github.com/alexdobin/STAR";
    license = licenses.gpl3Plus;
    platforms = platforms.unix;
    maintainers = [ maintainers.arcadio ];
  };
}
