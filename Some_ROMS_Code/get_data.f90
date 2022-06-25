      SUBROUTINE get_data (ng)
!
!svn $Id: get_data.F 1099 2022-01-06 21:01:01Z arango $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2022 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This routine reads in forcing, climatology and other data from      !
!  NetCDF files.  If there is more than one time-record,  data is      !
!  loaded into global  two-time  record arrays. The interpolation      !
!  is carried elsewhere.                                               !
!                                                                      !
!  Currently, this routine is only executed in serial mode by the      !
!  main thread.                                                        !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_boundary
      USE mod_clima
      USE mod_forces
      USE mod_grid
      USE mod_iounits
      USE mod_ncparam
      USE mod_scalars
      USE mod_sources
      USE mod_stepping
!
      USE strings_mod, ONLY : FoundError
!
      implicit none
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng
!
!  Local variable declarations.
!
      logical :: Lprocess
      logical, save :: recordless = .FALSE.
      logical, dimension(3) :: update = (/ .FALSE., .FALSE., .FALSE. /)
!
      integer :: ILB, IUB, JLB, JUB
      integer :: LBi, UBi, LBj, UBj
      integer :: i, ic, my_tile
!
      character (len=*), parameter :: MyFile =                          &
     &  "ROMS/Nonlinear/get_data.F"
!
!  Lower and upper bounds for nontiled (global values) boundary arrays.
!
      my_tile=-1                           ! for global values
      ILB=BOUNDS(ng)%LBi(my_tile)
      IUB=BOUNDS(ng)%UBi(my_tile)
      JLB=BOUNDS(ng)%LBj(my_tile)
      JUB=BOUNDS(ng)%UBj(my_tile)
!
!  Lower and upper bounds for tiled arrays.
!
      LBi=LBOUND(GRID(ng)%h,DIM=1)
      UBi=UBOUND(GRID(ng)%h,DIM=1)
      LBj=LBOUND(GRID(ng)%h,DIM=2)
      UBj=UBOUND(GRID(ng)%h,DIM=2)
!
!-----------------------------------------------------------------------
!  Turn on input data time wall clock.
!-----------------------------------------------------------------------
!
      CALL wclock_on (ng, iNLM, 3, 90, MyFile)
!
!=======================================================================
!  Point Sources/Sinks time dependent data.
!=======================================================================
!
!  Point Source/Sink vertically integrated mass transport.
!
      IF (LuvSrc(ng).or.LwSrc(ng)) THEN
        CALL get_ngfld (ng, iNLM, idRtra, SSF(ng)%ncid,                 &
     &                  1, SSF(ng), recordless, update(1),              &
     &                  1, Nsrc(ng), 1, 2, 1, Nsrc(ng), 1,              &
     &                  SOURCES(ng) % QbarG)
        IF (FoundError(exit_flag, NoError, 109, MyFile)) RETURN
      END IF
!
!  Tracer Sources/Sinks.
!
      DO i=1,NT(ng)
        IF (LtracerSrc(i,ng)) THEN
          CALL get_ngfld (ng, iNLM, idRtrc(i), SSF(ng)%ncid,            &
     &                    1, SSF(ng), recordless, update(1),            &
     &                    1, Nsrc(ng), N(ng), 2, 1, Nsrc(ng), N(ng),    &
     &                    SOURCES(ng) % TsrcG(:,:,:,i))
          IF (FoundError(exit_flag, NoError, 125, MyFile)) RETURN
        END IF
      END DO
!
!=======================================================================
!  Read in forcing data from FORCING NetCDF file(s).
!=======================================================================
!
!  Set switch to process surface atmospheric fields.
!
      Lprocess=.TRUE.
!
!  Surface wind components.
!
      IF (Lprocess) THEN
        CALL get_2dfld (ng, iNLM, idUair, FRCncid(idUair,ng),           &
     &                  nFfiles(ng), FRC(1,ng), update(1),              &
     &                  LBi, UBi, LBj, UBj, 2, 1,                       &
     &                  GRID(ng) % rmask,                               &
     &                  FORCES(ng) % UwindG)
        IF (FoundError(exit_flag, NoError, 169, MyFile)) RETURN
!
        CALL get_2dfld (ng , iNLM, idVair, FRCncid(idVair,ng),          &
     &                  nFfiles(ng), FRC(1,ng), update(1),              &
     &                  LBi, UBi, LBj, UBj, 2, 1,                       &
     &                  GRID(ng) % rmask,                               &
     &                  FORCES(ng) % VwindG)
        IF (FoundError(exit_flag, NoError, 181, MyFile)) RETURN
      END IF
!
!  Surface air pressure.
!
      IF (Lprocess) THEN
        CALL get_2dfld (ng, iNLM, idPair, FRCncid(idPair,ng),           &
     &                  nFfiles(ng), FRC(1,ng), update(1),              &
     &                  LBi, UBi, LBj, UBj, 2, 1,                       &
     &                  GRID(ng) % rmask,                               &
     &                  FORCES(ng) % PairG)
        IF (FoundError(exit_flag, NoError, 232, MyFile)) RETURN
      END IF
!
!  Surface solar shortwave radiation.
!
      IF (Lprocess) THEN
        CALL get_2dfld (ng, iNLM, idSrad, FRCncid(idSrad,ng),           &
     &                  nFfiles(ng), FRC(1,ng), update(1),              &
     &                  LBi, UBi, LBj, UBj, 2, 1,                       &
     &                  GRID(ng) % rmask,                               &
     &                  FORCES(ng) % srflxG)
        IF (FoundError(exit_flag, NoError, 394, MyFile)) RETURN
      END IF
!
!  Surface downwelling longwave radiation.
!
      IF (Lprocess) THEN
        CALL get_2dfld (ng, iNLM, idLdwn, FRCncid(idLdwn,ng),           &
     &                  nFfiles(ng), FRC(1,ng), update(1),              &
     &                  LBi, UBi, LBj, UBj, 2, 1,                       &
     &                  GRID(ng) % rmask,                               &
     &                  FORCES(ng) % lrflxG)
        IF (FoundError(exit_flag, NoError, 451, MyFile)) RETURN
      END IF
!
!  Surface air temperature.
!
      IF (Lprocess) THEN
        CALL get_2dfld (ng, iNLM, idTair, FRCncid(idTair,ng),           &
     &                  nFfiles(ng), FRC(1,ng), update(1),              &
     &                  LBi, UBi, LBj, UBj, 2, 1,                       &
     &                  GRID(ng) % rmask,                               &
     &                  FORCES(ng) % TairG)
        IF (FoundError(exit_flag, NoError, 472, MyFile)) RETURN
      END IF
!
!  Surface air humidity.
!
      IF (Lprocess) THEN
        CALL get_2dfld (ng, iNLM, idQair, FRCncid(idQair,ng),           &
     &                  nFfiles(ng), FRC(1,ng), update(1),              &
     &                  LBi, UBi, LBj, UBj, 2, 1,                       &
     &                  GRID(ng) % rmask,                               &
     &                  FORCES(ng) % HairG)
        IF (FoundError(exit_flag, NoError, 492, MyFile)) RETURN
      END IF
!
!  Rain fall rate.
!
      IF (Lprocess) THEN
        CALL get_2dfld (ng, iNLM, idrain, FRCncid(idrain,ng),           &
     &                  nFfiles(ng), FRC(1,ng), update(1),              &
     &                  LBi, UBi, LBj, UBj, 2, 1,                       &
     &                  GRID(ng) % rmask,                               &
     &                  FORCES(ng) % rainG)
        IF (FoundError(exit_flag, NoError, 511, MyFile)) RETURN
      END IF
!
!=======================================================================
!  Read in open boundary conditions from BOUNDARY NetCDF file.
!=======================================================================
!
!  Free-surface open boundary conditions.
!
      IF (LprocessOBC(ng)) THEN
        IF (LBC(iwest,isFsur,ng)%acquire) THEN
          CALL get_ngfld (ng, iNLM, idZbry(iwest),                      &
     &                    BRYncid(idZbry(iwest),ng),                    &
     &                    nBCfiles(ng), BRY(1,ng),                      &
     &                    recordless, update(1),                        &
     &                    JLB, JUB, 1, 2, 0, Mm(ng)+1, 1,               &
     &                    BOUNDARY(ng) % zetaG_west)
          IF (FoundError(exit_flag, NoError, 875, MyFile)) RETURN
        END IF
!
        IF (LBC(ieast,isFsur,ng)%acquire) THEN
          CALL get_ngfld (ng, iNLM, idZbry(ieast),                      &
     &                    BRYncid(idZbry(ieast),ng),                    &
     &                    nBCfiles(ng), BRY(1,ng),                      &
     &                    recordless, update(1),                        &
     &                    JLB, JUB, 1, 2, 0, Mm(ng)+1, 1,               &
     &                    BOUNDARY(ng) % zetaG_east)
          IF (FoundError(exit_flag, NoError, 888, MyFile)) RETURN
        END IF
!
        IF (LBC(isouth,isFsur,ng)%acquire) THEN
          CALL get_ngfld (ng, iNLM, idZbry(isouth),                     &
     &                    BRYncid(idZbry(isouth),ng),                   &
     &                    nBCfiles(ng), BRY(1,ng),                      &
     &                    recordless, update(1),                        &
     &                    ILB, IUB, 1, 2, 0, Lm(ng)+1, 1,               &
     &                    BOUNDARY(ng) % zetaG_south)
          IF (FoundError(exit_flag, NoError, 901, MyFile)) RETURN
        END IF
!
        IF (LBC(inorth,isFsur,ng)%acquire) THEN
          CALL get_ngfld (ng, iNLM, idZbry(inorth),                     &
     &                    BRYncid(idZbry(inorth),ng),                   &
     &                    nBCfiles(ng), BRY(1,ng),                      &
     &                    recordless, update(1),                        &
     &                    ILB, IUB, 1, 2, 0, Lm(ng)+1, 1,               &
     &                    BOUNDARY(ng) % zetaG_north)
          IF (FoundError(exit_flag, NoError, 914, MyFile)) RETURN
        END IF
      END IF
!
!  2D momentum components open boundary conditions.
!
      IF (LprocessOBC(ng)) THEN
        IF (LBC(iwest,isUbar,ng)%acquire) THEN
          CALL get_ngfld (ng, iNLM, idU2bc(iwest),                      &
     &                    BRYncid(idU2bc(iwest),ng),                    &
     &                    nBCfiles(ng), BRY(1,ng),                      &
     &                    recordless, update(1),                        &
     &                    JLB, JUB, 1, 2, 0, Mm(ng)+1, 1,               &
     &                    BOUNDARY(ng) % ubarG_west)
          IF (FoundError(exit_flag, NoError, 934, MyFile)) RETURN
        END IF
!
        IF (LBC(iwest,isVbar,ng)%acquire) THEN
          CALL get_ngfld (ng, iNLM, idV2bc(iwest),                      &
     &                    BRYncid(idV2bc(iwest),ng),                    &
     &                    nBCfiles(ng), BRY(1,ng),                      &
     &                    recordless, update(1),                        &
     &                    JLB, JUB, 1, 2, 1, Mm(ng)+1, 1,               &
     &                    BOUNDARY(ng) % vbarG_west)
          IF (FoundError(exit_flag, NoError, 947, MyFile)) RETURN
        END IF
!
        IF (LBC(ieast,isUbar,ng)%acquire) THEN
          CALL get_ngfld (ng, iNLM, idU2bc(ieast),                      &
     &                    BRYncid(idU2bc(ieast),ng),                    &
     &                    nBCfiles(ng), BRY(1,ng),                      &
     &                    recordless, update(1),                        &
     &                    JLB, JUB, 1, 2, 0, Mm(ng)+1, 1,               &
     &                    BOUNDARY(ng) % ubarG_east)
          IF (FoundError(exit_flag, NoError, 960, MyFile)) RETURN
        END IF
!
        IF (LBC(ieast,isVbar,ng)%acquire) THEN
          CALL get_ngfld (ng, iNLM, idV2bc(ieast),                      &
     &                    BRYncid(idV2bc(ieast),ng),                    &
     &                    nBCfiles(ng), BRY(1,ng),                      &
     &                    recordless, update(1),                        &
     &                    JLB, JUB, 1, 2, 1, Mm(ng)+1, 1,               &
     &                    BOUNDARY(ng) % vbarG_east)
          IF (FoundError(exit_flag, NoError, 973, MyFile)) RETURN
        END IF
!
        IF (LBC(isouth,isUbar,ng)%acquire) THEN
          CALL get_ngfld (ng, iNLM, idU2bc(isouth),                     &
     &                    BRYncid(idU2bc(isouth),ng),                   &
     &                    nBCfiles(ng), BRY(1,ng),                      &
     &                    recordless, update(1),                        &
     &                    ILB, IUB, 1, 2, 1, Lm(ng)+1, 1,               &
     &                    BOUNDARY(ng) % ubarG_south)
          IF (FoundError(exit_flag, NoError, 986, MyFile)) RETURN
        END IF
!
        IF (LBC(isouth,isVbar,ng)%acquire) THEN
          CALL get_ngfld (ng, iNLM, idV2bc(isouth),                     &
     &                    BRYncid(idV2bc(isouth),ng),                   &
     &                    nBCfiles(ng), BRY(1,ng),                      &
     &                    recordless, update(1),                        &
     &                    ILB, IUB, 1, 2, 0, Lm(ng)+1, 1,               &
     &                    BOUNDARY(ng) % vbarG_south)
          IF (FoundError(exit_flag, NoError, 999, MyFile)) RETURN
        END IF
!
        IF (LBC(inorth,isUbar,ng)%acquire) THEN
          CALL get_ngfld (ng, iNLM, idU2bc(inorth),                     &
     &                    BRYncid(idU2bc(inorth),ng),                   &
     &                    nBCfiles(ng), BRY(1,ng),                      &
     &                    recordless, update(1),                        &
     &                    ILB, IUB, 1, 2, 1, Lm(ng)+1, 1,               &
     &                    BOUNDARY(ng) % ubarG_north)
          IF (FoundError(exit_flag, NoError, 1012, MyFile)) RETURN
        END IF
!
        IF (LBC(inorth,isVbar,ng)%acquire) THEN
          CALL get_ngfld (ng, iNLM, idV2bc(inorth),                     &
     &                    BRYncid(idV2bc(inorth),ng),                   &
     &                    nBCfiles(ng), BRY(1,ng),                      &
     &                    recordless, update(1),                        &
     &                    ILB, IUB, 1, 2, 0, Lm(ng)+1, 1,               &
     &                    BOUNDARY(ng) % vbarG_north)
          IF (FoundError(exit_flag, NoError, 1025, MyFile)) RETURN
        END IF
      END IF
!
!  3D momentum components open boundary conditions.
!
      IF (LprocessOBC(ng)) THEN
        IF (LBC(iwest,isUvel,ng)%acquire) THEN
          CALL get_ngfld (ng, iNLM, idU3bc(iwest),                      &
     &                    BRYncid(idU3bc(iwest),ng),                    &
     &                    nBCfiles(ng), BRY(1,ng),                      &
     &                    recordless, update(1),                        &
     &                    JLB, JUB, N(ng), 2, 0, Mm(ng)+1, N(ng),       &
     &                    BOUNDARY(ng) % uG_west)
          IF (FoundError(exit_flag, NoError, 1046, MyFile)) RETURN
        END IF
!
        IF (LBC(iwest,isVvel,ng)%acquire) THEN
          CALL get_ngfld (ng, iNLM, idV3bc(iwest),                      &
     &                    BRYncid(idV3bc(iwest),ng),                    &
     &                    nBCfiles(ng), BRY(1,ng),                      &
     &                    recordless, update(1),                        &
     &                    JLB, JUB, N(ng), 2, 1, Mm(ng)+1, N(ng),       &
     &                    BOUNDARY(ng) % vG_west)
          IF (FoundError(exit_flag, NoError, 1059, MyFile)) RETURN
        END IF
!
        IF (LBC(ieast,isUvel,ng)%acquire) THEN
          CALL get_ngfld (ng, iNLM, idU3bc(ieast),                      &
     &                    BRYncid(idU3bc(ieast),ng),                    &
     &                    nBCfiles(ng), BRY(1,ng),                      &
     &                    recordless, update(1),                        &
     &                    JLB, JUB, N(ng), 2, 0, Mm(ng)+1, N(ng),       &
     &                    BOUNDARY(ng) % uG_east)
          IF (FoundError(exit_flag, NoError, 1072, MyFile)) RETURN
        END IF
!
        IF (LBC(ieast,isVvel,ng)%acquire) THEN
          CALL get_ngfld (ng, iNLM, idV3bc(ieast),                      &
     &                    BRYncid(idV3bc(ieast),ng),                    &
     &                    nBCfiles(ng), BRY(1,ng),                      &
     &                    recordless, update(1),                        &
     &                    JLB, JUB, N(ng), 2, 1, Mm(ng)+1, N(ng),       &
     &                    BOUNDARY(ng) % vG_east)
          IF (FoundError(exit_flag, NoError, 1085, MyFile)) RETURN
        END IF
!
        IF (LBC(isouth,isUvel,ng)%acquire) THEN
          CALL get_ngfld (ng, iNLM, idU3bc(isouth),                     &
     &                    BRYncid(idU3bc(isouth),ng),                   &
     &                    nBCfiles(ng), BRY(1,ng),                      &
     &                    recordless, update(1),                        &
     &                    ILB, IUB, N(ng), 2, 1, Lm(ng)+1, N(ng),       &
     &                    BOUNDARY(ng) % uG_south)
          IF (FoundError(exit_flag, NoError, 1098, MyFile)) RETURN
        END IF
!
        IF (LBC(isouth,isVvel,ng)%acquire) THEN
          CALL get_ngfld (ng, iNLM, idV3bc(isouth),                     &
     &                    BRYncid(idV3bc(isouth),ng),                   &
     &                    nBCfiles(ng), BRY(1,ng),                      &
     &                    recordless, update(1),                        &
     &                    ILB, IUB, N(ng), 2, 0, Lm(ng)+1, N(ng),       &
     &                    BOUNDARY(ng) % vG_south)
          IF (FoundError(exit_flag, NoError, 1111, MyFile)) RETURN
        END IF
!
        IF (LBC(inorth,isUvel,ng)%acquire) THEN
          CALL get_ngfld (ng, iNLM, idU3bc(inorth),                     &
     &                    BRYncid(idU3bc(inorth),ng),                   &
     &                    nBCfiles(ng), BRY(1,ng),                      &
     &                    recordless, update(1),                        &
     &                    ILB, IUB, N(ng), 2, 1, Lm(ng)+1, N(ng),       &
     &                    BOUNDARY(ng) % uG_north)
          IF (FoundError(exit_flag, NoError, 1124, MyFile)) RETURN
        END IF
!
        IF (LBC(inorth,isVvel,ng)%acquire) THEN
          CALL get_ngfld (ng, iNLM, idV3bc(inorth),                     &
     &                    BRYncid(idV3bc(inorth),ng),                   &
     &                    nBCfiles(ng), BRY(1,ng),                      &
     &                    recordless, update(1),                        &
     &                    ILB, IUB, N(ng), 2, 0, Lm(ng)+1, N(ng),       &
     &                    BOUNDARY(ng) % vG_north)
          IF (FoundError(exit_flag, NoError, 1137, MyFile)) RETURN
        END IF
      END IF
!
!  Tracer variables open boundary conditions.
!
      IF (LprocessOBC(ng)) THEN
        DO i=1,NT(ng)
          IF (LBC(iwest,isTvar(i),ng)%acquire) THEN
            CALL get_ngfld (ng, iNLM, idTbry(iwest,i),                  &
     &                      BRYncid(idTbry(iwest,i),ng),                &
     &                      nBCfiles(ng), BRY(1,ng),                    &
     &                      recordless, update(1),                      &
     &                      JLB, JUB, N(ng), 2, 0, Mm(ng)+1, N(ng),     &
     &                      BOUNDARY(ng) % tG_west(:,:,:,i))
            IF (FoundError(exit_flag, NoError, 1158, MyFile)) RETURN
          END IF
        END DO
!
        DO i=1,NT(ng)
          IF (LBC(ieast,isTvar(i),ng)%acquire) THEN
            CALL get_ngfld (ng, iNLM, idTbry(ieast,i),                  &
     &                      BRYncid(idTbry(ieast,i),ng),                &
     &                      nBCfiles(ng), BRY(1,ng),                    &
     &                      recordless, update(1),                      &
     &                      JLB, JUB, N(ng), 2, 0, Mm(ng)+1, N(ng),     &
     &                      BOUNDARY(ng) % tG_east(:,:,:,i))
            IF (FoundError(exit_flag, NoError, 1173, MyFile)) RETURN
          END IF
        END DO
!
        DO i=1,NT(ng)
          IF (LBC(isouth,isTvar(i),ng)%acquire) THEN
            CALL get_ngfld (ng, iNLM, idTbry(isouth,i),                 &
     &                      BRYncid(idTbry(isouth,i),ng),               &
     &                      nBCfiles(ng), BRY(1,ng),                    &
     &                      recordless, update(1),                      &
     &                      ILB, IUB, N(ng), 2, 0, Lm(ng)+1, N(ng),     &
     &                      BOUNDARY(ng) % tG_south(:,:,:,i))
            IF (FoundError(exit_flag, NoError, 1188, MyFile)) RETURN
          END IF
        END DO
!
        DO i=1,NT(ng)
          IF (LBC(inorth,isTvar(i),ng)%acquire) THEN
            CALL get_ngfld (ng, iNLM, idTbry(inorth,i),                 &
     &                      BRYncid(idTbry(inorth,i),ng),               &
     &                      nBCfiles(ng), BRY(1,ng),                    &
     &                      recordless, update(1),                      &
     &                      ILB, IUB, N(ng), 2, 0, Lm(ng)+1, N(ng),     &
     &                      BOUNDARY(ng) % tG_north(:,:,:,i))
            IF (FoundError(exit_flag, NoError, 1203, MyFile)) RETURN
          END IF
        END DO
      END IF
!
!=======================================================================
!  Read in data from Climatology NetCDF file.
!=======================================================================
!
!  Free-surface climatology.
!
      IF (LsshCLM(ng)) THEN
        CALL get_2dfld (ng, iNLM, idSSHc, CLMncid(idSSHc,ng),           &
     &                  nCLMfiles(ng), CLM(1,ng), update(1),            &
     &                  LBi, UBi, LBj, UBj, 2, 1,                       &
     &                  GRID(ng) % rmask,                               &
     &                  CLIMA(ng) % sshG)
        IF (FoundError(exit_flag, NoError, 1229, MyFile)) RETURN
      END IF
!
!  2D momentum components climatology.
!
      IF (Lm2CLM(ng)) THEN
        CALL get_2dfld (ng, iNLM, idUbcl, CLMncid(idUbcl,ng),           &
     &                  nCLMfiles(ng), CLM(1,ng), update(1),            &
     &                  LBi, UBi, LBj, UBj, 2, 1,                       &
     &                  GRID(ng) % umask,                               &
     &                  CLIMA(ng) % ubarclmG)
        IF (FoundError(exit_flag, NoError, 1247, MyFile)) RETURN
!
        CALL get_2dfld (ng, iNLM, idVbcl, CLMncid(idVbcl,ng),           &
     &                  nCLMfiles(ng), CLM(1,ng), update(1),            &
     &                  LBi, UBi, LBj, UBj, 2, 1,                       &
     &                  GRID(ng) % vmask,                               &
     &                  CLIMA(ng) % vbarclmG)
        IF (FoundError(exit_flag, NoError, 1259, MyFile)) RETURN
      END IF
!
!  3D momentum components climatology.
!
      IF (Lm3CLM(ng)) THEN
        CALL get_3dfld (ng, iNLM, idUclm, CLMncid(idUclm,ng),           &
     &                  nCLMfiles(ng), CLM(1,ng), update(1),            &
     &                  LBi, UBi, LBj, UBj, 1, N(ng), 2, 1,             &
     &                  GRID(ng) % umask,                               &
     &                  CLIMA(ng) % uclmG)
        IF (FoundError(exit_flag, NoError, 1278, MyFile)) RETURN
!
        CALL get_3dfld (ng, iNLM, idVclm, CLMncid(idVclm,ng),           &
     &                  nCLMfiles(ng), CLM(1,ng), update(1),            &
     &                  LBi, UBi, LBj, UBj, 1, N(ng), 2, 1,             &
     &                  GRID(ng) % vmask,                               &
     &                  CLIMA(ng) % vclmG)
        IF (FoundError(exit_flag, NoError, 1290, MyFile)) RETURN
      END IF
!
!  Tracers variables climatology.
!
      ic=0
      DO i=1,NT(ng)
        IF (LtracerCLM(i,ng)) THEN
          ic=ic+1
          CALL get_3dfld (ng, iNLM, idTclm(i),                          &
     &                    CLMncid(idTclm(i),ng),                        &
     &                    nCLMfiles(ng), CLM(1,ng), update(1),          &
     &                    LBi, UBi, LBj, UBj, 1, N(ng), 2, 1,           &
     &                    GRID(ng) % rmask,                             &
     &                    CLIMA(ng) % tclmG(:,:,:,:,ic))
          IF (FoundError(exit_flag, NoError, 1312, MyFile)) RETURN
        END IF
      END DO
!
!-----------------------------------------------------------------------
!  Turn off input data time wall clock.
!-----------------------------------------------------------------------
!
      CALL wclock_off (ng, iNLM, 3, 1523, MyFile)
!
      RETURN
      END SUBROUTINE get_data
