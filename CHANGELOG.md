# Changelog

## [8.1] - November 2017
### Added
- Electroweak one-loop corrections for $Z$, $t\bar t$ and di-jet production [1].
- $Z\gamma$ process at NNLO including anomalous couplings [2].
- $Z\gamma$ decay for H+jet production [2].
- $H$+2jet process with a finite top-quark mass and $H$+jet with finite
  top-quark mass effects [1]. This adds a new parameter [mtex] to the
  the input file, which specifies to which order in 1/mt^k (k=0,2,4) the finite
  part of the virtual corrections are computed.
- Support for random seeds by setting the seed value (previously [ij]) to 0.
- Support for boosted (as opposed to hadronic) definition of jettiness, which
  is now also used as default.
- Support for $p_T$ and rapidity ranges for most cuts in the input file.
- Native implementation of two PDF sets containing photons: mrstqed and
  CT14qed.
- EXPERIMENTAL: Integration routine that adaptively selects the cross section
  contribution with the largest integration uncertainty. It also adapts
  integration calls in the warmup phase and then continuously increases them to
  not run into grid bias problems. This ingration can be enabled by setting
  newIntegration to .true. in Need/mcfmmain.f. After each iteration a snapshot is
  saved for resumption. When readin is set to .true. in the gridinfo_logic
  common block, the integration will be resumed from any stage. Currently this
  integration mode does not stop the integration, but the user can at any point
  abort manually when the desired precision has been reached. It has only been
  thoroughly tested with $Z\gamma$ and $H+$jet. Other processes might need
  tweaks or specific initializations.

### Changed
- Fixed FROOT implementation.
- Upgraded QCDLoop to QCDLoop 2.
- Using new random number generator (Mersenne Twister) from libstdc++ (C++11).
- Fixed bug concerning integration uncertainty of histogram bins.

### Removed
- Removed flags virtonly,realonly,vanillafiles in input.DAT, as well as broken
  parameters for resuming integrations.
- Removed support for PDFLIB.
- Removed noomp version. If just a single thread is wanted, please
  set OMP_NUM_THREADS to 1.

[1] Campbell, Wackeroth, Zhou https://arxiv.org/abs/1608.03356
[2] Neumann, Williams https://arxiv.org/abs/1609.00367
[3] Campbell, Neumann, Williams https://arxiv.org/abs/1708.02925
