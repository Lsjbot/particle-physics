CDECK  ID>, BLANKDEK.   
*CMZ :  3.15/09 07/04/92  14.40.14  by  Federico Carminati
*-- Author :
CDECK  ID>, NUMAIN. 
      PROGRAM NUMAIN
C
      COMMON/NUFLAG/FGEANT,KGEN,KFLUXE,KFLUXS,FWEI,
     *       IRNALG,IRNSAV,IEVTRS,IEVACT,NEVACT,
     *       KCERAL,FOPTIN,FDPST,FRSCAT,IVTXAL,FTISME,
     *       ZINVBS
      LOGICAL FGEANT,FWEI,FOPTIN,FDPST,FRSCAT,FTISME,ZINVBS
      COMMON/NUZZIO/IZUNII,IZUNIO,IOZH0,IOEZH,IOTZH,IOCPA,IOCP1,
     *      IOPMH,IOSDIR,IOSHHI,IOSHLO,IOLAST
      COMMON/NUFWIO/IFWIO,IWGEAN,IWZH0,IWEZH,IWTZH,IWCPA,IWCP1,IWPMH,
     +              IWLAST
      COMMON/NUHBIO/IHBUNI
C
C-----------------------------------------------
C Memory requirements for HBOOK, Geant, and our
C private Zebra structure.
C-----------------------------------------------
C
      PARAMETER ( LHMEM  = 400000 )
      PARAMETER ( LGMEM  = 100000 )
      PARAMETER ( LZMEM  = 600000 )
      PARAMETER ( IFFSZE = 1000   )
      COMMON/PAWC/HHH(LHMEM)
      COMMON/GCBANK/Q(LGMEM)
      COMMON/NUBANK/FENCE(16),QSPACE(LZMEM)
      DIMENSION LQN(LZMEM),IQN(LZMEM-8)
      DIMENSION QN(LZMEM-8)
      EQUIVALENCE ( QSPACE(1), LQN(1) )
      EQUIVALENCE ( IQN(1), QN(1), LQN(9) )
      COMMON/CFREAD/FFSPAC(IFFSZE)
C-----------------------------------------------
      COMMON/GCNUM/NMATE ,NVOLUM,NROTM,NTMED,NTMULT,NTRACK,NPART
     +            ,NSTMAX,NVERTX,NHEAD,NBIT
      COMMON /GCNUMX/ NALIVE,NTMSTO
C
      COMMON/GCFLAG/IDEBUG,IDEMIN,IDEMAX,ITEST,IDRUN,IDEVT,IEORUN
     +        ,IEOTRI,IEVENT,ISWIT(10),IFINIT(20),NEVENT,NRNDM(2)
      COMMON/GCFLAX/BATCH, NOLOG
      LOGICAL BATCH, NOLOG
C
      COMMON/GCLIST/NHSTA,NGET ,NSAVE,NSETS,NPRIN,NGEOM,NVIEW,NPLOT
     +       ,NSTAT,LHSTA(20),LGET (20),LSAVE(20),LSETS(20),LPRIN(20)
     +             ,LGEOM(20),LVIEW(20),LPLOT(20),LSTAT(20)
C
      PARAMETER ( MAXPM  = 250 )
      PARAMETER ( MXPMHI = 100*MAXPM )
      COMMON/NUCPM/PMINT(MAXPM),TFIRST(MAXPM),TAVE(MAXPM),NPMINT(MAXPM),
     +             NHITOT, IPMHIT(MXPMHI),
     +             TIMHIT(MXPMHI), CPHINT(MXPMHI), NPEINT(MXPMHI),
     +             PHPPMH(MXPMHI), ANGPMH(MXPMHI), BETPMH(MXPMHI),
     +             IPAHIT(MXPMHI)
      COMMON/NUPMPO/NNPM,X0PM(MAXPM),Y0PM(MAXPM),Z0PM(MAXPM),
     +              IORIPM(MAXPM), LAYSTR(MAXPM)
C     IORIPM specifies the orientation of the PM tubes:
C     IORIPM(i) = -1: PM tube number i is looking downward
C     IORIPM(i) = +1: PM tube number i is looking upward
C     LAYSTR(i) = the number of the layer or string that PM tube
C                 number i belongs to
      COMMON/NUGRID/ZSIZE,HGPOS,DISTX,DISTY,DISTZ,OFFLAY,OFFLX,OFFLY,
     *              OFFX,OFFY,NBINS,NLAY,NSIDE,NSTR,NPMSTR,IDCONF,
     *              ICONTY,PMDIAM,PMRAD,PMR2,PMAREA,XYFUZZ,
     *              VXPLSI, VXPLDI
      COMMON/NUFILE/DATFIL,RNDFIN,RNDFOU,HISFIL
      CHARACTER*48  DATFIL,RNDFIN,RNDFOU,HISFIL
      PARAMETER ( MXZBIN = 50 )
      PARAMETER ( MXTBIN = 40 )
      PARAMETER ( MXSHOW = 5000 )
      COMMON/NUSHLC/SHPHOT(0:MXZBIN,0:MXTBIN)
     *             ,SHT(0:MXZBIN,0:MXTBIN),SHT2(0:MXZBIN,0:MXTBIN)
     *             ,SHLATE(0:MXZBIN,0:MXTBIN)
     *             ,NBZ,NBTH,ZBINN(0:MXZBIN),THBINN(0:MXTBIN)
     *             ,EESHOW(MXSHOW),IDXSHO(MXSHOW)
     *             ,USESHL,NWBIN,INSLIB,NSH,CGEN,NZBMAX,NTBMAX
     *             ,QQCER2,QQSLIB,NPCER2,NPSLIB
     *             ,SHOFIL,CHDIR,CHNPAH,CHCPAH,CHNPAZ,CHCPAZ
      CHARACTER*48 SHOFIL
      CHARACTER*5 CHDIR
      CHARACTER*20 CHNPAH,CHCPAH,CHNPAZ,CHCPAZ
      LOGICAL USESHL,INSLIB,CGEN
C
      INTEGER IRSEED
C
C   'Random' dummies to be loaded:
      EXTERNAL EURRAN
      EXTERNAL EURRIN
      EXTERNAL GRNDM
      EXTERNAL GRANOR
      EXTERNAL GRNDMQ
      EXTERNAL NORRAN
      EXTERNAL NORRIN
      EXTERNAL NORRUT
      EXTERNAL NRAN
      EXTERNAL NRANIN
      EXTERNAL NRANUT
      EXTERNAL RAND
      EXTERNAL RANNOR
      EXTERNAL RLU
      EXTERNAL RNDM
      EXTERNAL RDMIN
      EXTERNAL RDMOUT
      EXTERNAL RNDM2
      EXTERNAL RD2IN
      EXTERNAL RD2OUT
      EXTERNAL RN32
      EXTERNAL RN32IN
      EXTERNAL RN32UT
C
C    GEANT dummies to be loaded:
      EXTERNAL GUTREV
      EXTERNAL GUKINE
      EXTERNAL GUTRAK
      EXTERNAL GUSTEP
      EXTERNAL GUOUT
      EXTERNAL GUPHAD
      EXTERNAL GUHADR
C
C   'Replace' routines to be loaded:
      EXTERNAL FUNRAN
C
C    LEPTO block data to be loaded:
      EXTERNAL LEPTOD
C
C    LUND block datas to be loaded:
      EXTERNAL LUDATA
      EXTERNAL LUEDAT
      EXTERNAL LUHDAT
C
C--> Call dummy routine to ensure correct routines and BLOCK DATAs
C--> will be loaded from library:
C
C   'Random' dummies:
      CALL NULOAD(EURRAN,'EURRAN')
      CALL NULOAD(EURRIN,'EURRIN')
      CALL NULOAD(GRNDM, 'GRNDM' )
      CALL NULOAD(GRANOR,'GRANOR')
      CALL NULOAD(GRNDMQ,'GRNDMQ')
      CALL NULOAD(NORRAN,'NORRAN')
      CALL NULOAD(NORRIN,'NORRIN')
      CALL NULOAD(NORRUT,'NORRUT')
      CALL NULOAD(NRAN,  'NRAN')
      CALL NULOAD(NRANIN,'NRANIN')
      CALL NULOAD(NRANUT,'NRANUT')
      CALL NULOAD(RAND,  'RAND'  )
      CALL NULOAD(RANNOR,'RANNOR')
      CALL NULOAD(RLU ,  'RLU'   )
      CALL NULOAD(RNDM,  'RNDM'  )
      CALL NULOAD(RDMIN, 'RDMIN' )
      CALL NULOAD(RDMOUT,'RDMOUT')
      CALL NULOAD(RNDM2, 'RNDM2' )
      CALL NULOAD(RD2IN, 'RD2IN' )
      CALL NULOAD(RD2OUT,'RD2OUT')
      CALL NULOAD(RN32,  'RN32'  )
      CALL NULOAD(RN32IN,'RN32IN')
      CALL NULOAD(RN32UT,'RN32UT')
C
C    GEANT dummies:
      CALL NULOAD(GUTREV,'GUTREV')
      CALL NULOAD(GUKINE,'GUKINE')
      CALL NULOAD(GUTRAK,'GUTRAK')
      CALL NULOAD(GUSTEP,'GUSTEP')
      CALL NULOAD(GUOUT, 'GUOUT' )
      CALL NULOAD(GUPHAD,'GUPHAD')
      CALL NULOAD(GUHADR,'GUHADR')
C
C   'Replace' routines:
      CALL NULOAD(FUNRAN,'FUNRAN')
C
C    LEPTO block data:
      CALL NULOAD(LEPTOD,'LEPTOD')
C
C    LUND block datas:
      CALL NULOAD(LUDATA,'LUDATA')
      CALL NULOAD(LUEDAT,'LUEDAT')
      CALL NULOAD(LUHDAT,'LUHDAT')
C
C
C--> Set default parameters
C
      CALL NUDATA
C
C--> Initialise Geant (has to be done before NUREAD)
C
      CALL GZEBRA(LGMEM)
      CALL GINIT
C
C--> Read data cards, ours and Geant's
C
      CALL NUREAD
C
C--> Now, after having read in our datacards, and before we do anything
C    else, let's initialise the random-number generation
C
      CALL NURINI
C
C--> Initialise HBOOK, negative LHMEM to cohabit with Geant
C
      CALL HLIMIT(-LHMEM)
C
C--> Initialise Geant (cont'd)
C
      CALL UGINIT
C
C--> Possibly book standard Geant histograms
C
      IF (FGEANT .AND. (NHSTA .GT. 0)) THEN
         CALL GBHSTA
      ENDIF
C
C--> Initialise event generator
C
      IF (KGEN .EQ. 2) THEN
         CALL NULPIN
      ELSEIF (KGEN .EQ. 3) THEN
         CALL NUPYIN
      ELSEIF (KGEN .EQ. 4) THEN
         CALL NUFRIN
      ENDIF
C
C--> NUSIM initialisation, including histogram booking
C
      CALL NUINIT
C
C--> Event loop
C
      IF (FGEANT) THEN
         CALL GRUN
      ELSE
         CALL NULOOP
      ENDIF
C
C--> NOTE: the number of actually simulated events may differ from
C    IEVENT in case we did restart the simulation from event number
C    IEVTRS > 1. - Therefore, write out IEVACT rather than IEVENT.
C
      IF (IEVACT .LT. 1E5) THEN
         WRITE(6,1) IEVACT
      ELSEIF (IEVACT .LT. 1E7) THEN
         WRITE(6,2) IEVACT
      ELSE
         WRITE(6,3) IEVACT
      ENDIF
   1  FORMAT(' '/' ',37('*')/
     +' ','NUMAIN: END OF RUN AFTER ',I5,' EVENTS'/' ',37('*')/' ')
   2  FORMAT(' '/' ',39('*')/
     +' ','NUMAIN: END OF RUN AFTER ',I7,' EVENTS'/' ',39('*')/' ')
   3  FORMAT(' '/' ',46('*')/
     +' ','NUMAIN: END OF RUN AFTER ',I, ' EVENTS'/' ',46('*')/' ')
C
      IF (IEVACT .NE. IEVENT)
     +WRITE(6,*) 'LAST SIMULATED EVENT WAS NUMBER ',IEVENT
C
C
C--> Save random-number seeds on file after each run, if anything
C    of this kind has been requested:
C
      IF (IRNSAV .EQ. 0) GOTO 923
C
      IF (IRNSAV .LT. 0) THEN
         REWIND(34,ERR=921)
         GOTO 922
  921    CONTINUE
         WRITE(6,4)
     +   '*** NUMAIN: ERROR REWINDING RANDOM-NUMBER OUTPUT FILE ',
     +   '***         (UNIT 34) AT END OF RUN'
   4     FORMAT(' '/' ',A/' ',A/' ')
         CLOSE(34)
         GOTO 923
  922    CONTINUE
      ENDIF
      CALL RNWRIT(34,1,IEVENT+1)
      CLOSE(34)
C
  923 CONTINUE
C
C ---------------------------------------------------------------------
C
C--> Termination
C
C    User termination:
      CALL NUSEND
C
      CALL HPRINT(0)
C
      IF (IHBUNI .GT. 0) THEN
C       Write out all histograms and ntuples
C      WRITE(6,*)'BEFORE IF USESHL ',USESHL
      IF ( USESHL ) THEN
C        WRITE(6,*)'CHPAH = ',CHNPAH,CHCPAH
        CALL RZNDIR(CHNPAZ,'P')
        CALL RZCDIR(CHCPAZ,'P')
        CALL RZNDIR(CHNPAH,' ')
        CALL RZCDIR(CHCPAH,' ')
        CALL RZNDIR(CHNPAH,'P')
        CALL RZCDIR(CHCPAH,'P')
      ENDIF
         CALL HROUT (0, ICYCLE, ' ')
C       Close the RZ file for 'NUSIM'
         CALL HLDIR (' ', ' ')
         CALL HREND ('NUSIM')
      ENDIF
C
C    GEANT termination:
      IF (FGEANT) THEN
         CALL GLAST
      ENDIF
C
C    Close Zebra files:
      CALL NUZEND
C
      STOP
      END
CDECK  ID>, NUREAD. 
      SUBROUTINE NUREAD
************************************************************************
*                                                                      *
* NUREAD (vsn 1.0)                                                     *
*                                                                      *
* This routine reads the run dependent values for the title constants  *
* from file, using FFREAD.                                             *
*                                                                      *
* Shamelessly stolen from DELSIM, routine SXREAD,                      *
* written by:  Fons Rademakers (i91@nikhefh.uucp), 13-Jan-1987         *
* modified by: Sverker Johansson, 13-April-1992                        *
*                                                                      *
************************************************************************
C
C
C-----------------------------------------------
C Memory requirements for HBOOK, Geant, and our
C private Zebra structure.
C-----------------------------------------------
C
      PARAMETER ( LHMEM  = 400000 )
      PARAMETER ( LGMEM  = 100000 )
      PARAMETER ( LZMEM  = 600000 )
      PARAMETER ( IFFSZE = 1000   )
      COMMON/PAWC/HHH(LHMEM)
      COMMON/GCBANK/Q(LGMEM)
      COMMON/NUBANK/FENCE(16),QSPACE(LZMEM)
      DIMENSION LQN(LZMEM),IQN(LZMEM-8)
      DIMENSION QN(LZMEM-8)
      EQUIVALENCE ( QSPACE(1), LQN(1) )
      EQUIVALENCE ( IQN(1), QN(1), LQN(9) )
      COMMON/CFREAD/FFSPAC(IFFSZE)
C-----------------------------------------------
      COMMON/NUFLAG/FGEANT,KGEN,KFLUXE,KFLUXS,FWEI,
     *       IRNALG,IRNSAV,IEVTRS,IEVACT,NEVACT,
     *       KCERAL,FOPTIN,FDPST,FRSCAT,IVTXAL,FTISME,
     *       ZINVBS
      LOGICAL FGEANT,FWEI,FOPTIN,FDPST,FRSCAT,FTISME,ZINVBS
      COMMON/NUGENC/E0MIN,E0MAX,TH0MIN,TH0MAX,CTHMIN,CTHMAX,PHI0,
     +              NPPROJ,NPTARG
      PARAMETER ( MAXZH = 500 )
      COMMON/NUZH/NZH,JZH,DISTZH,ZH(MAXZH),EZH(MAXZH),THEZH(MAXZH),
     +            CPATH(MAXZH),CPTOT,HISTZH,IHZH,FTHEZH,ETHLIM
      LOGICAL HISTZH,FTHEZH
      PARAMETER ( NLMAX = 4 )
      PARAMETER ( NBB = 160 )
      PARAMETER ( NPM = 10 )
      COMMON/NUPATC/CERPAT(NLMAX,NBB,NBB),TPAT(NLMAX,NBB,NBB)
     *        ,T2PAT(NLMAX,NBB,NBB)
     *        ,TFPAT(NLMAX,NBB,NBB),TLPAT(NLMAX,NBB,NBB)
     *        ,NPHTOT(NLMAX)
      COMMON/NUPM/PM(NLMAX,-NPM:NPM,-NPM:NPM)
     *           ,BINWID,PMTBIN,KBINCX,KBINCY
     *           ,NPHPM(NLMAX,-NPM:NPM,-NPM:NPM)
     *           ,TT(NLMAX,-NPM:NPM,-NPM:NPM)
     *           ,TSIG(NLMAX,-NPM:NPM,-NPM:NPM)
     *           ,TFL(NLMAX,-NPM:NPM,-NPM:NPM)
      COMMON/NUGRID/ZSIZE,HGPOS,DISTX,DISTY,DISTZ,OFFLAY,OFFLX,OFFLY,
     *              OFFX,OFFY,NBINS,NLAY,NSIDE,NSTR,NPMSTR,IDCONF,
     *              ICONTY,PMDIAM,PMRAD,PMR2,PMAREA,XYFUZZ,
     *              VXPLSI, VXPLDI
      COMMON/NUZZIO/IZUNII,IZUNIO,IOZH0,IOEZH,IOTZH,IOCPA,IOCP1,
     *      IOPMH,IOSDIR,IOSHHI,IOSHLO,IOLAST
      COMMON/NUFWIO/IFWIO,IWGEAN,IWZH0,IWEZH,IWTZH,IWCPA,IWCP1,IWPMH,
     +              IWLAST
      COMMON/NUHBIO/IHBUNI
      COMMON/NUZLOG/LOGLMZ,LOGLFZ,LOGLDZ
      PARAMETER ( MAXPM  = 250 )
      PARAMETER ( MXPMHI = 100*MAXPM )
      COMMON/NUCPM/PMINT(MAXPM),TFIRST(MAXPM),TAVE(MAXPM),NPMINT(MAXPM),
     +             NHITOT, IPMHIT(MXPMHI),
     +             TIMHIT(MXPMHI), CPHINT(MXPMHI), NPEINT(MXPMHI),
     +             PHPPMH(MXPMHI), ANGPMH(MXPMHI), BETPMH(MXPMHI),
     +             IPAHIT(MXPMHI)
      COMMON/NUPMPO/NNPM,X0PM(MAXPM),Y0PM(MAXPM),Z0PM(MAXPM),
     +              IORIPM(MAXPM), LAYSTR(MAXPM)
C     IORIPM specifies the orientation of the PM tubes:
C     IORIPM(i) = -1: PM tube number i is looking downward
C     IORIPM(i) = +1: PM tube number i is looking upward
C     LAYSTR(i) = the number of the layer or string that PM tube
C                 number i belongs to
      COMMON/NUTRIC/LTRIG, NPEMIN
      LOGICAL LTRIG
C     LTRIG  = Trigger flag; to be set .TRUE. if the current event has
C              generated a trigger
C     NPEMIN = Minimum # photo-electrons defining a signal in a PM tube
C              (to be used in the trigger algorithm)
      COMMON/NUMISC/START,  IPTYPE, PIN, EMU0, THEMU0, PHIMU0,
     +              ZTRMAX, IORIVX, IORITR, VOLD(7),
     +              KTMED,  REFIND, ATTLEN, CWLGT1, CWLGT2,
     +              CUTKC2, IPRHIS(3), TDIAM, INUTRA,
     +              ZCLAPP, DCLAPP
      INTEGER START
      COMMON /NUTIMC/ TRGATE
      REAL TRGATE
C     TRGATE = Trigger gate: time interval (in seconds) during which
C              the acqusition system is open for recording data; i.e.,
C              sensitive to PM-tube signals.
      COMMON /NUNOIC/ RTNOIS
      REAL RTNOIS
C     RTNOIS = Rate of noise hits in a PM tube (# hits/second).
      COMMON/NUFILE/DATFIL,RNDFIN,RNDFOU,HISFIL
      CHARACTER*48  DATFIL,RNDFIN,RNDFOU,HISFIL
      PARAMETER ( MXZBIN = 50 )
      PARAMETER ( MXTBIN = 40 )
      PARAMETER ( MXSHOW = 5000 )
      COMMON/NUSHLC/SHPHOT(0:MXZBIN,0:MXTBIN)
     *             ,SHT(0:MXZBIN,0:MXTBIN),SHT2(0:MXZBIN,0:MXTBIN)
     *             ,SHLATE(0:MXZBIN,0:MXTBIN)
     *             ,NBZ,NBTH,ZBINN(0:MXZBIN),THBINN(0:MXTBIN)
     *             ,EESHOW(MXSHOW),IDXSHO(MXSHOW)
     *             ,USESHL,NWBIN,INSLIB,NSH,CGEN,NZBMAX,NTBMAX
     *             ,QQCER2,QQSLIB,NPCER2,NPSLIB
     *             ,SHOFIL,CHDIR,CHNPAH,CHCPAH,CHNPAZ,CHCPAZ
      CHARACTER*48 SHOFIL
      CHARACTER*5 CHDIR
      CHARACTER*20 CHNPAH,CHCPAH,CHNPAZ,CHCPAZ
      LOGICAL USESHL,INSLIB,CGEN
      INTEGER ISHFIL(12)
      DATA ISHFIL/12*0/
C
      CHARACTER*16 FILNAM
      LOGICAL FEXIST
      CHARACTER*4 NAMESP
      DATA IFFIN/18/
      DATA IPRNT/6/
C
      INTEGER IDTFIL(12)
      DATA IDTFIL/12*0/
C
      INTEGER IRNFII(12)
      DATA IRNFII/12*0/
C
      INTEGER IRNFIO(12)
      DATA IRNFIO/12*0/
C
      INTEGER IHIFIL(12)
      DATA IHIFIL/12*0/
C
      DATA FILNAM/'SIMRUN:TIT      '/
*
*-- Initialize FFREAD
*
      CALL FFINIT(IFFSZE)
*
*-- Set FFREAD options
*
      CALL FFSET('LINP', IFFIN)
      CALL FFSET('LOUT', IPRNT)
      CALL FFSET('SIZE', 6    )
*
*-- Define user keys. They are in the same order as the title constants
*-- in NUDATA. For a description of the constants see NUDATA.
*
      CALL FFKEY('FGEANT', FGEANT   , 1  ,'LOGICAL')
      CALL FFKEY('IRNALG', IRNALG   , 1  ,'INTEGER')
      CALL FFKEY('IRNSAV', IRNSAV   , 1  ,'INTEGER')
      CALL FFKEY('IEVTRS', IEVTRS   , 1  ,'INTEGER')
      CALL FFKEY('KGEN  ', KGEN     , 1  ,'INTEGER')
      CALL FFKEY('NPPROJ', NPPROJ   , 1  ,'INTEGER')
      CALL FFKEY('KFLUXE', KFLUXE   , 1  ,'INTEGER')
      CALL FFKEY('KFLUXS', KFLUXS   , 1  ,'INTEGER')
      CALL FFKEY('E0MIN ', E0MIN    , 1  ,'REAL   ')
      CALL FFKEY('E0MAX ', E0MAX    , 1  ,'REAL   ')
      CALL FFKEY('TH0MIN', TH0MIN   , 1  ,'REAL   ')
      CALL FFKEY('TH0MAX', TH0MAX   , 1  ,'REAL   ')
      CALL FFKEY('PHI0  ', PHI0     , 1  ,'REAL   ')
      CALL FFKEY('IVTXAL', IVTXAL   , 1  ,'INTEGER')
      CALL FFKEY('FWEI  ', FWEI     , 1  ,'LOGICAL')
      CALL FFKEY('ZSIZE ', ZSIZE    , 1  ,'REAL   ')
      CALL FFKEY('HGPOS ', HGPOS    , 1  ,'REAL   ')
      CALL FFKEY('VXPLSI', VXPLSI,    1  ,'REAL   ')
      CALL FFKEY('VXPLDI', VXPLDI,    1  ,'REAL   ')
      CALL FFKEY('NZH   ', NZH      , 1  ,'INTEGER')
      CALL FFKEY('DISTZH', DISTZH   , 1  ,'REAL   ')
      CALL FFKEY('IHZH  ', IHZH     , 1  ,'INTEGER')
      CALL FFKEY('HISTZH', HISTZH   , 1  ,'LOGICAL')
      CALL FFKEY('KCERAL', KCERAL   , 1  ,'INTEGER')
      CALL FFKEY('ZTRMAX', ZTRMAX   , 1  ,'REAL   ')
      CALL FFKEY('FOPTIN', FOPTIN   , 1  ,'LOGICAL')
      CALL FFKEY('FTISME', FTISME   , 1  ,'LOGICAL')
      CALL FFKEY('RTNOIS', RTNOIS,    1  ,'REAL   ')
      CALL FFKEY('TRGATE', TRGATE,    1  ,'REAL   ')
      CALL FFKEY('REFIND', REFIND   , 1  ,'REAL   ')
      CALL FFKEY('FDPST ', FDPST    , 1  ,'LOGICAL')
      CALL FFKEY('FTHEZH', FTHEZH   , 1  ,'LOGICAL')
      CALL FFKEY('KTMED ', KTMED    , 1  ,'INTEGER')
      CALL FFKEY('ETHLIM', ETHLIM   , 1  ,'REAL   ')
      CALL FFKEY('IZUNII', IZUNII   , 1  ,'INTEGER')
      CALL FFKEY('IZUNIO', IZUNIO   , 1  ,'INTEGER')
      CALL FFKEY('IHBUNI', IHBUNI   , 1  ,'INTEGER')
      CALL FFKEY('LOGLMZ', LOGLMZ   , 1  ,'INTEGER')
      CALL FFKEY('LOGLFZ', LOGLFZ   , 1  ,'INTEGER')
      CALL FFKEY('LOGLDZ', LOGLDZ   , 1  ,'INTEGER')
      CALL FFKEY('ZINVBS', ZINVBS   , 1  ,'LOGICAL')
      CALL FFKEY('FRSCAT', FRSCAT   , 1  ,'LOGICAL')
      CALL FFKEY('IDCONF', IDCONF   , 1  ,'INTEGER')
      CALL FFKEY('NLAY  ', NLAY     , 1  ,'INTEGER')
      CALL FFKEY('NSIDE ', NSIDE    , 1  ,'INTEGER')
      CALL FFKEY('NSTR  ', NSTR     , 1  ,'INTEGER')
      CALL FFKEY('NPMSTR', NPMSTR   , 1  ,'INTEGER')
      CALL FFKEY('DISTX ', DISTX    , 1  ,'REAL   ')
      CALL FFKEY('DISTY ', DISTY    , 1  ,'REAL   ')
      CALL FFKEY('DISTZ ', DISTZ    , 1  ,'REAL   ')
      CALL FFKEY('NNPM  ', NNPM     , 1  ,'INTEGER')
      CALL FFKEY('PMDIAM', PMDIAM   , 1  ,'REAL   ')
      CALL FFKEY('XYFUZZ', XYFUZZ   , 1  ,'REAL   ')
      CALL FFKEY('ATTLEN', ATTLEN   , 1  ,'REAL   ')
      CALL FFKEY('CWLGT1', CWLGT1   , 1  ,'REAL   ')
      CALL FFKEY('CWLGT2', CWLGT2   , 1  ,'REAL   ')
      CALL FFKEY('NPEMIN', NPEMIN   , 1  ,'INTEGER')
      CALL FFKEY('IFWIO ', IFWIO    , 1  ,'INTEGER')
      CALL FFKEY('IWGEAN', IWGEAN   , 1  ,'INTEGER')
      CALL FFKEY('IWZH0 ', IWZH0    , 1  ,'INTEGER')
      CALL FFKEY('IWEZH ', IWEZH    , 1  ,'INTEGER')
      CALL FFKEY('IWTZH ', IWTZH    , 1  ,'INTEGER')
      CALL FFKEY('IWCPA ', IWCPA    , 1  ,'INTEGER')
      CALL FFKEY('IWCP1 ', IWCP1    , 1  ,'INTEGER')
      CALL FFKEY('IWPMH ', IWPMH    , 1  ,'INTEGER')
      CALL FFKEY('IPRHIS', IPRHIS(1), 3  ,'INTEGER')
      CALL FFKEY('DATFIL', IDTFIL(1),12  ,'MIXED')
      CALL FFKEY('RNDFIN', IRNFII(1),12  ,'MIXED')
      CALL FFKEY('RNDFOU', IRNFIO(1),12  ,'MIXED')
      CALL FFKEY('HISFIL', IHIFIL(1),12  ,'MIXED')
      CALL FFKEY('SHOFIL',ISHFIL(1),12,'MIXED')
      CALL FFKEY('USESHL',USESHL,1,'LOGICAL')
*
*     Set flag for lab name
      LABORA = 0
*
*-- Get user defined steering cards
*
      CALL NUUFFR
*
*-- Now read data cards, using Geant's card-reading routine
*-- so that we catch also Geant cards.
*
      CALL GFFGO
*
      WRITE(6,*) ' '
      IF (IRNALG .LE. 0) THEN
         WRITE(6,*) '*** NUREAD ***'
         WRITE(6,*) '*** Invalid random-number starting algorithm: ',
     +              'IRNALG = ',IRNALG
         WRITE(6,*) '*** Changing this to IRNALG = 1 (clock start)'
         IRNALG = 1
      ELSEIF (IRNALG .EQ. 1) THEN
         WRITE(6,*)
     +   'Chosen random-number starting algorithm: clock start'
      ELSEIF (IRNALG .EQ. 2) THEN
         WRITE(6,*)
     +   'Chosen random-number starting algorithm: restart'
         WRITE(6,*)
     +   'Simulation is to restart at event number ',IEVTRS
         IF (IEVTRS .LE. 0) THEN
            WRITE(6,*) '*** INVALID EVENT NUMBER!           ***'
            WRITE(6,*) '*** CHANGING THIS TO EVENT NUMBER 1 ***'
            IEVTRS = 1
         ENDIF
      ELSEIF (IRNALG .GT. 921350143) THEN
         WRITE(6,*) '*** NUREAD ***'
         WRITE(6,*)
     +   '*** Invalid start seed(s) for random-number generation: '
         WRITE(6,*) '*** IRNALG (= ISEED) = ',IRNALG
         WRITE(6,*) '*** Changing this to IRNALG (= ISEED) = 17560127'
         IRNALG = 17560127
      ELSEIF (IRNALG .GE. 3) THEN
         WRITE(6,*)
     +   'Chosen random-number starting algorithm: fixed start'
         WRITE(6,*) 'Start seed(s) = ',IRNALG
      ENDIF
*
      IF (IZUNII .GT. 0) THEN
         IF (IZUNIO .GT. 0) THEN
            WRITE(6,*) ' '
            WRITE(6,*) '*** NUREAD: ',
     +      'Both ZEBRA input and output requested, not allowed!'
            WRITE(6,*) '*** Output disabled.'
            WRITE(6,*) ' '
            IZUNIO = 0
         ENDIF
      ENDIF
*
      IF (E0MIN .GT. E0MAX) THEN
         WRITE(6,*) ' '
         WRITE(6,*) '*** NUREAD:'
         WRITE(6,*) '*** E0MIN > E0MAX:'
         WRITE(6,*) '*** E0MIN = ',E0MIN,'   E0MAX = ',E0MAX
         WRITE(6,*) '*** E0MAX changed to ',E0MIN
         E0MAX = E0MIN
      ENDIF
*
      IF (TH0MIN .GT. TH0MAX) THEN
         WRITE(6,*) ' '
         WRITE(6,*) '*** NUREAD:'
         WRITE(6,*) '*** TH0MIN > TH0MAX:'
         WRITE(6,*) '*** TH0MIN = ',TH0MIN,'   TH0MAX = ',TH0MAX
         WRITE(6,*) '*** TH0MAX changed to ',TH0MIN
         TH0MAX = TH0MIN
      ENDIF
*
      IF (KFLUXE .EQ. 1) THEN
         IF (E0MIN .EQ. E0MAX) THEN
            WRITE(6,*) ' '
            WRITE(6,*) '*** NUREAD:'
            WRITE(6,*) '*** KFLUXE, E0MIN, E0MAX = ',
     +                      KFLUXE, E0MIN, E0MAX
            WRITE(6,*) '*** I.e., monoenergetic energy distr. requested'
            WRITE(6,*) '*** ==> KFLUXE changed to 0'
            KFLUXE = 0
         ENDIF
      ENDIF
*
      IF ((IVTXAL .LT. 1) .OR. (IVTXAL .GT. 3)) THEN
         WRITE(6,*) ' '
         WRITE(6,*) '*** NUREAD: IVTXAL =',IVTXAL
         WRITE(6,*) '*** No such alternative; setting IVTXAL = 1'
         IVTXAL = 1
      ENDIF
*
      IF (KCERAL .EQ. 1) THEN
         IF (NLAY .GT. NLMAX) THEN
            WRITE(6,*) ' '
            WRITE(6,*) '*** NUREAD:'
            WRITE(6,*) '*** Too many layers!'
            WRITE(6,*) '*** NLAY adjusted to ',NLMAX
            NLAY = NLMAX
         ENDIF
      ENDIF
*
      IF (KCERAL .NE. 0) THEN
         IF (CWLGT1 .GE. CWLGT2) THEN
            WRITE(6,*) ' '
            WRITE(6,*) '*** NUREAD:'
            WRITE(6,*) '*** Impossible Cherenkov-wavelength interval!'
            IF (CWLGT1 .GT. CWLGT2) THEN
               WRITE(6,*) '*** CWLGT1 > CWLGT2: ',CWLGT1, CWLGT2
            ELSE
               WRITE(6,*) '*** CWLGT1 = CWLGT2: ',CWLGT1, CWLGT2
            ENDIF
            CWLGT2 = CWLGT1 + 100.
            WRITE(6,*) '*** Wavelength interval changed to be:'
            WRITE(6,*) '*** ',CWLGT1,'   < lambda <   ',CWLGT2
         ENDIF
      ENDIF
*
      IF (IVTXAL .EQ. 3) THEN
         VXPLSM = VXPLSI/100.
         VXPLDM = VXPLDI/100.
         WRITE(6,*) ' '
         WRITE(6,*)
     +   'Length of the side of the ''vertex plane''       = ',
     +    VXPLSM,'   metres'
         WRITE(6,*)
     +   'Distance of the ''vertex plane'' from the origin = ',
     +    VXPLDM,'   metres'
      ENDIF
*
*-- Transfer the name of the input or output file from IDTFIL to DATFIL
*
      IF (IDTFIL(1) .NE. 0) THEN
         CALL UHTOC(IDTFIL(1),4,DATFIL,48)
      ENDIF
      IF (IZUNII .GT. 0) THEN
         WRITE(6,101) 'Name of input data file: ', DATFIL
      ELSEIF (IZUNIO .GT. 0) THEN
         WRITE(6,101) 'Name of output data file: ', DATFIL
      ENDIF
*
*-- Transfer the name of the random-number file from IRNFII to RNDFIN
*
      IF (IRNFII(1) .NE. 0) THEN
         CALL UHTOC(IRNFII(1),4,RNDFIN,48)
      ENDIF
      IF (IRNALG .EQ. 2) THEN
         WRITE(6,101)
     +   'Name of the input file for random-number seeds:',
     +   RNDFIN
      ENDIF
*
*-- Transfer the name of the random-number file from IRNFIO to RNDFOU
*
      IF (IRNFIO(1) .NE. 0) THEN
         CALL UHTOC(IRNFIO(1),4,RNDFOU,48)
      ENDIF
      IF (IRNSAV .NE. 0) THEN
         WRITE(6,101)
     +   'Name of the output file for random-number seeds: ',
     +   RNDFOU
         IF (ABS(IRNSAV) .EQ. 1) THEN
            WRITE(6,*) 'Seeds will be written out after each event.'
         ELSE
            WRITE(6,*) 'Seeds will be written out after every ',
     +                 ABS(IRNSAV),'th event.'
         ENDIF
         IF (IRNSAV .GT. 0) THEN
            WRITE(6,*)
     +      'All saved seeds will be kept in the output file.'
         ELSE
            WRITE(6,*)
     +      'Only the last saved seeds will be kept in the output file.'
         ENDIF
      ENDIF
*
*-- Transfer the name of the histogram file from IHIFIL to HISFIL
*
      IF (IHIFIL(1) .NE. 0) THEN
         CALL UHTOC(IHIFIL(1),4,HISFIL,48)
      ENDIF
      IF (IHBUNI .GT. 0) THEN
         WRITE(6,101) 'Name of the output histogram/ntuple file: ',
     +   HISFIL
      ENDIF
*
  101 FORMAT(' '/' ',A/' ',A48/' ')
*
      IF ( ISHFIL(1) .NE. 0 ) THEN
        CALL UHTOC(ISHFIL(1),4,SHOFIL,48)
      ENDIF
C      WRITE(6,101)'Name of the shower library file : ',SHOFIL
*-- Give values to some constants
*
      CTHMIN = COS(3.141592654*TH0MIN/180.)
      CTHMAX = COS(3.141592654*TH0MAX/180.)
*
      PMRAD  = 0.5*PMDIAM
      PMR2   = PMRAD*PMRAD
      PMAREA = 3.141592654*PMR2
*
*
      IF(LABORA.NE.0) CALL UHTOC( LABORA , 4 , LABO   ,4)
  999 CONTINUE
*
      RETURN
      END
CDECK  ID>, NUZINI. 
      SUBROUTINE NUZINI
C
***********************************************************
*                                                         *
* NUZINI (vsn 1.0)                                        *
*                                                         *
* Create Zebra structure                                  *
*                                                         *
***********************************************************
C
C
C-----------------------------------------------
C Memory requirements for HBOOK, Geant, and our
C private Zebra structure.
C-----------------------------------------------
C
      PARAMETER ( LHMEM  = 400000 )
      PARAMETER ( LGMEM  = 100000 )
      PARAMETER ( LZMEM  = 600000 )
      PARAMETER ( IFFSZE = 1000   )
      COMMON/PAWC/HHH(LHMEM)
      COMMON/GCBANK/Q(LGMEM)
      COMMON/NUBANK/FENCE(16),QSPACE(LZMEM)
      DIMENSION LQN(LZMEM),IQN(LZMEM-8)
      DIMENSION QN(LZMEM-8)
      EQUIVALENCE ( QSPACE(1), LQN(1) )
      EQUIVALENCE ( IQN(1), QN(1), LQN(9) )
      COMMON/CFREAD/FFSPAC(IFFSZE)
C-----------------------------------------------
      COMMON/NUZLIN/IXSTOR,IXDIV1,IXDIV2,NLSTRU,NLREF,LROOT
      COMMON/NULINK/JZZH0,JZEZH,JZTZH,JZCPA,JZCP1,JZPMH
     *      ,JZSDIR,JZSHLO,JZSHHI,JZLAST
      PARAMETER ( MAXZH = 500 )
      COMMON/NUZH/NZH,JZH,DISTZH,ZH(MAXZH),EZH(MAXZH),THEZH(MAXZH),
     +            CPATH(MAXZH),CPTOT,HISTZH,IHZH,FTHEZH,ETHLIM
      LOGICAL HISTZH,FTHEZH
      COMMON/NUZZIO/IZUNII,IZUNIO,IOZH0,IOEZH,IOTZH,IOCPA,IOCP1,
     *      IOPMH,IOSDIR,IOSHHI,IOSHLO,IOLAST
      COMMON/NUFWIO/IFWIO,IWGEAN,IWZH0,IWEZH,IWTZH,IWCPA,IWCP1,IWPMH,
     +              IWLAST
      COMMON/NUHBIO/IHBUNI
      COMMON/NUZLOG/LOGLMZ,LOGLFZ,LOGLDZ
      COMMON/NUFILE/DATFIL,RNDFIN,RNDFOU,HISFIL
      CHARACTER*48  DATFIL,RNDFIN,RNDFOU,HISFIL
      PARAMETER ( MXZBIN = 50 )
      PARAMETER ( MXTBIN = 40 )
      PARAMETER ( MXSHOW = 5000 )
      COMMON/NUSHLC/SHPHOT(0:MXZBIN,0:MXTBIN)
     *             ,SHT(0:MXZBIN,0:MXTBIN),SHT2(0:MXZBIN,0:MXTBIN)
     *             ,SHLATE(0:MXZBIN,0:MXTBIN)
     *             ,NBZ,NBTH,ZBINN(0:MXZBIN),THBINN(0:MXTBIN)
     *             ,EESHOW(MXSHOW),IDXSHO(MXSHOW)
     *             ,USESHL,NWBIN,INSLIB,NSH,CGEN,NZBMAX,NTBMAX
     *             ,QQCER2,QQSLIB,NPCER2,NPSLIB
     *             ,SHOFIL,CHDIR,CHNPAH,CHCPAH,CHNPAZ,CHCPAZ
      CHARACTER*48 SHOFIL
      CHARACTER*5 CHDIR
      CHARACTER*20 CHNPAH,CHCPAH,CHNPAZ,CHCPAZ
      LOGICAL USESHL,INSLIB,CGEN
      DIMENSION KEYDIR(2)
C
      DATA LEN/900/
C
      NLSTRU = 10
      NLREF  = 10
      LIM2 = LZMEM/2
      CALL MZSTOR(IXSTOR,'NUSIM',' ',FENCE,LQN(1),LQN(NLSTRU+1),
     *         LQN(NLSTRU+NLREF+1),LQN(LIM2),LQN(LZMEM))
      CALL MZLOGL(IXSTOR,LOGLMZ)
C
      CALL MZLINK(IXSTOR,'/NULINK/',JZZH0,JZLAST,JZZH0)
C
      CALL MZFORM('ZH0 ','1I -F',IOZH0)
      CALL MZFORM('EZH ','1I -F',IOEZH)
      CALL MZFORM('TZH ','1I -F',IOTZH)
      CALL MZFORM('CPA ','1I -F',IOCPA)
      CALL MZFORM('CP1 ','1I -F',IOCP1)
      CALL MZFORM('PMH ','1I -F',IOPMH)
      CALL MZFORM('EOEV','1I',   IOLAST)
      CALL MZFORM('SDIR','1I 92F / 1F 1I',IOSDIR)
      CALL MZFORM('SHLO','3I -F',IOSHLO)
      CALL MZFORM('SHHI','3I -F',IOSHHI)
C
      CALL MZBOOK(IXSTOR,JZZH0,JZZH0,1,'ZH0 ',1,1,NZH+2,IOZH0,0)
      CALL MZBOOK(IXSTOR,JZLAST,JZLAST,1,'EOEV',1,1,1,IOLAST,0)
C
      IF (IZUNII .GT. 0) THEN
C
C
        CALL GOPEN(IZUNII,'IX',LEN,IERR)
        IF (IERR .NE. 0) GOTO 911
        CALL FZLOGL(IZUNII,LOGLFZ)
C
      ELSEIF (IZUNIO .GT. 0) THEN
C
C
        CALL GOPEN(IZUNIO,'OX',LEN,IERR)
        IF (IERR .NE. 0) GOTO 911
        CALL FZLOGL(IZUNIO,LOGLFZ)
C
      ENDIF
      CALL RZNDIR(CHNPAH,'R')
      CALL RZCDIR(CHCPAH,'R')
      CALL RZNDIR(CHNPAH,'P')
      CALL RZCDIR(CHCPAH,'P')
      LRECL = 0
      CHDIR = 'LUN61'
      CALL RZOPEN(61,CHDIR,SHOFIL,'X',LRECL,ISTAT)
      CALL RZFILE(61,CHDIR,'X')
      CALL RZLOGL(61,2)
      CALL RZLDIR(' ',' ')
      CALL RZSTAT(CHDIR,1,' ')
C
      KEYDIR(1) = -1
      KEYDIR(2) = -1
      ICYCLE = 9999
      CALL RZIN(IXSTOR,JZSDIR,2,KEYDIR,ICYCLE,' ')
C
      NWBIN = MXZBIN+MXTBIN+2
      DUM = NBINZ(0.0)
      DUM = NBINTH(0.0)
      DO 26 I=0,MXZBIN
        ZZZ = QN(JZSDIR+2+I)
C                IF ( ZZZ .NE. ZBINN(I)) THEN
C                  WRITE(6,*)'ZBINN CHANGED: ',ZZZ,ZBINN(I)
C                ENDIF
        ZBINN(I) = ZZZ
        IF ( ZZZ .GT. 0 ) THEN
          NZBMAX = I
        ENDIF
   26 CONTINUE
      NZBMAX = NZBMAX-1
C      WRITE(6,*)'NZBMAX = ',NZBMAX
      DO 28 I=0,MXTBIN
        TTT = QN(JZSDIR+3+MXZBIN+I)
        THBINN(I) = TTT
        IF ( TTT .GT. 0 ) THEN
          NTBMAX = I
        ENDIF
   28 CONTINUE
      NTBMAX = NTBMAX-1
C      WRITE(6,*)'NTBMAX = ',NTBMAX
C
      NSH = IQN(JZSDIR+1)
      IF ( NSH .GT. MXSHOW ) THEN
        NSH = MXSHOW
      ENDIF
      DO 35 I=1,NSH
        EESHOW(I) = QN(JZSDIR+NWBIN+2*I)
   35 CONTINUE
C
      CALL SORTZV(EESHOW,IDXSHO,NSH,1,0,0)
      DO 37 I=1,NSH
        EESHOW(I) = QN(JZSDIR+NWBIN+2*(IDXSHO(I)))
C        WRITE(6,*)'EESHOW= ',EESHOW(I)
   37 CONTINUE
C

C
CCC      CALL DZSHOW('Shower directory',IXSTOR,JZSDIR,'B',1,1,40)
C
      CALL RZNDIR(CHNPAZ,'R')
      CALL RZCDIR(CHCPAZ,'R')
      CALL RZNDIR(CHNPAZ,'P')
      CALL RZCDIR(CHCPAZ,'P')
      CALL RZNDIR(CHNPAH,' ')
      CALL RZCDIR(CHCPAH,' ')
C
C
      RETURN
C
 911  WRITE(6,*) ' '
      IF (IZUNII .GT. 0) THEN
         WRITE(6,*) '*** OPEN error in NUZINI: Unit = ',IZUNII
      ELSEIF (IZUNIO .GT. 0) THEN
         WRITE(6,*) '*** OPEN error in NUZINI: Unit = ',IZUNIO
      ENDIF
C
      IZUNII = 0
      IZUNIO = 0
C
      RETURN
      END
CDECK  ID>, NUHBOK. 
      SUBROUTINE NUHBOK
C
      COMMON/NUMISC/START,  IPTYPE, PIN, EMU0, THEMU0, PHIMU0,
     +              ZTRMAX, IORIVX, IORITR, VOLD(7),
     +              KTMED,  REFIND, ATTLEN, CWLGT1, CWLGT2,
     +              CUTKC2, IPRHIS(3), TDIAM, INUTRA,
     +              ZCLAPP, DCLAPP
      INTEGER START
      COMMON/NUFLAG/FGEANT,KGEN,KFLUXE,KFLUXS,FWEI,
     *       IRNALG,IRNSAV,IEVTRS,IEVACT,NEVACT,
     *       KCERAL,FOPTIN,FDPST,FRSCAT,IVTXAL,FTISME,
     *       ZINVBS
      LOGICAL FGEANT,FWEI,FOPTIN,FDPST,FRSCAT,FTISME,ZINVBS
      PARAMETER ( MAXZH = 500 )
      COMMON/NUZH/NZH,JZH,DISTZH,ZH(MAXZH),EZH(MAXZH),THEZH(MAXZH),
     +            CPATH(MAXZH),CPTOT,HISTZH,IHZH,FTHEZH,ETHLIM
      LOGICAL HISTZH,FTHEZH
      COMMON/NUGENC/E0MIN,E0MAX,TH0MIN,TH0MAX,CTHMIN,CTHMAX,PHI0,
     +              NPPROJ,NPTARG
      COMMON/NUGRID/ZSIZE,HGPOS,DISTX,DISTY,DISTZ,OFFLAY,OFFLX,OFFLY,
     *              OFFX,OFFY,NBINS,NLAY,NSIDE,NSTR,NPMSTR,IDCONF,
     *              ICONTY,PMDIAM,PMRAD,PMR2,PMAREA,XYFUZZ,
     *              VXPLSI, VXPLDI
      PARAMETER ( NLMAX = 4 )
      PARAMETER ( NBB = 160 )
      PARAMETER ( NPM = 10 )
      COMMON/NUPATC/CERPAT(NLMAX,NBB,NBB),TPAT(NLMAX,NBB,NBB)
     *        ,T2PAT(NLMAX,NBB,NBB)
     *        ,TFPAT(NLMAX,NBB,NBB),TLPAT(NLMAX,NBB,NBB)
     *        ,NPHTOT(NLMAX)
      COMMON/NUPM/PM(NLMAX,-NPM:NPM,-NPM:NPM)
     *           ,BINWID,PMTBIN,KBINCX,KBINCY
     *           ,NPHPM(NLMAX,-NPM:NPM,-NPM:NPM)
     *           ,TT(NLMAX,-NPM:NPM,-NPM:NPM)
     *           ,TSIG(NLMAX,-NPM:NPM,-NPM:NPM)
     *           ,TFL(NLMAX,-NPM:NPM,-NPM:NPM)
      PARAMETER ( MAXPM  = 250 )
      PARAMETER ( MXPMHI = 100*MAXPM )
      COMMON/NUCPM/PMINT(MAXPM),TFIRST(MAXPM),TAVE(MAXPM),NPMINT(MAXPM),
     +             NHITOT, IPMHIT(MXPMHI),
     +             TIMHIT(MXPMHI), CPHINT(MXPMHI), NPEINT(MXPMHI),
     +             PHPPMH(MXPMHI), ANGPMH(MXPMHI), BETPMH(MXPMHI),
     +             IPAHIT(MXPMHI)
      COMMON/NUPMPO/NNPM,X0PM(MAXPM),Y0PM(MAXPM),Z0PM(MAXPM),
     +              IORIPM(MAXPM), LAYSTR(MAXPM)
C     IORIPM specifies the orientation of the PM tubes:
C     IORIPM(i) = -1: PM tube number i is looking downward
C     IORIPM(i) = +1: PM tube number i is looking upward
C     LAYSTR(i) = the number of the layer or string that PM tube
C                 number i belongs to
      COMMON/NUZZIO/IZUNII,IZUNIO,IOZH0,IOEZH,IOTZH,IOCPA,IOCP1,
     *      IOPMH,IOSDIR,IOSHHI,IOSHLO,IOLAST
      COMMON/NUFWIO/IFWIO,IWGEAN,IWZH0,IWEZH,IWTZH,IWCPA,IWCP1,IWPMH,
     +              IWLAST
      COMMON/NUHBIO/IHBUNI
      COMMON/NUFILE/DATFIL,RNDFIN,RNDFOU,HISFIL
      CHARACTER*48  DATFIL,RNDFIN,RNDFOU,HISFIL
      PARAMETER (MXLYST = 99)
      PARAMETER (NNUVA1 = 11)
      PARAMETER (MXNUVA = NNUVA1 + 3*MXLYST)
      COMMON /NUNTUP/ NLYSTR, NNUVAR, NUNTVE(MXNUVA)
      REAL NUNTVE
C
C     MXLYST = Maximum # layers or strings that we can handle in our
C              NTUPLE
C     NNUVA1 = Number of variables for our NTUPLE, other than those
C              representing # firing PM tubes per layer or string
C     MXNUVA = Maximum # variables that we can handle in our NTUPLE
C     NLYSTR = Actual # layers or strings in the simulated set-up
C     NNUVAR = Actual # variables for our NTUPLE
C     NUNTVE = REAL vector containing the information that we fill our
C              NTUPLE with.
*
      PARAMETER ( MXZBIN = 50 )
      PARAMETER ( MXTBIN = 40 )
      PARAMETER ( MXSHOW = 5000 )
      COMMON/NUSHLC/SHPHOT(0:MXZBIN,0:MXTBIN)
     *             ,SHT(0:MXZBIN,0:MXTBIN),SHT2(0:MXZBIN,0:MXTBIN)
     *             ,SHLATE(0:MXZBIN,0:MXTBIN)
     *             ,NBZ,NBTH,ZBINN(0:MXZBIN),THBINN(0:MXTBIN)
     *             ,EESHOW(MXSHOW),IDXSHO(MXSHOW)
     *             ,USESHL,NWBIN,INSLIB,NSH,CGEN,NZBMAX,NTBMAX
     *             ,QQCER2,QQSLIB,NPCER2,NPSLIB
     *             ,SHOFIL,CHDIR,CHNPAH,CHCPAH,CHNPAZ,CHCPAZ
      CHARACTER*48 SHOFIL
      CHARACTER*5 CHDIR
      CHARACTER*20 CHNPAH,CHCPAH,CHNPAZ,CHCPAZ
      LOGICAL USESHL,INSLIB,CGEN
      PARAMETER ( MAXQQ  = 5000.)
      PARAMETER ( MAXTIM = 100.)
      PARAMETER ( MAXPS  = 5000)
      COMMON /MYSHOC/ XPRTRA(MAXPS), YPRTRA(MAXPS), ZPRTRA(MAXPS),
     +                DXPTRA(MAXPS), DYPTRA(MAXPS), DZPTRA(MAXPS),
     +                PEMAS1,PTOT1(MAXPS),PTOT2,IPS, IPSTRA,
     +                ZSHMAX,ZDIST(MAXPS),ETOT1(MAXPS),TOTPHO,
     +                DEGTHE(MAXQQ),
     +                GARMS, ELES,NEUS,MUS,NPBS,OTHERS,PIO,PIS
C
      CHARACTER*32 CHTIT1
      CHARACTER*37 CHTIT2
      CHARACTER*40 CHTIT3
C
      CHARACTER*6  CNUVAR(MXNUVA), CLSTOT, CLSDOW, CLSUPP,
     +                             CLYSTT, CLYSTD, CLYSTU
C
      DATA CNUVAR/'PIN   ','THETA ','PHI   ','DCLAPP','ZCLAPP',
     +            '      ','      ','      ',
     +            'NPMTOT','NPMDOW','NPMUPP',
     +     MXLYST*'      ',
     +     MXLYST*'      ',
     +     MXLYST*'      '/
C
C (These histograms are for 'numu.f' -- change as needed.)
C
C ----------------------------------------------------------------------
C Avoid using histogram IDs 1, 2, and 3 - they are used by NUPMIN for
C the histograms displaying the PM-tube positions (although those histo-
C grams should already have been deleted when NUHBOK is called).
C ----------------------------------------------------------------------
C
C--> Detector geometry, etc.
C
      IF (FGEANT .AND. (IVTXAL .EQ. 1)) THEN
         HGPOS1 = HGPOS + 0.1
         ZSIZE1 = ZSIZE + 0.1
         CALL HBOOK2(4,
     +   'Vertex position: y vs. x (Track enters from below)',
     +                    40,-HGPOS,HGPOS1,   40,-HGPOS,HGPOS1,   0.)
         CALL HBOOK2(5,
     +   'Vertex position: y vs. x (Track enters from above)',
     +                    40,-HGPOS,HGPOS1,   40,-HGPOS,HGPOS1,   0.)
         CALL HBOOK2(6,
     +   'Vertex position: z vs. y (Track enters from the left)',
     +                    40,-HGPOS,HGPOS1,   40,-ZSIZE,ZSIZE1,   0.)
         CALL HBOOK2(7,
     +   'Vertex position: z vs. y (Track enters from the right)',
     +                    40,-HGPOS,HGPOS1,   40,-ZSIZE,ZSIZE1,   0.)
         CALL HBOOK2(8,
     +   'Vertex position: z vs. x (Track enters from the front)',
     +                    40,-HGPOS,HGPOS1,   40,-ZSIZE,ZSIZE1,   0.)
         CALL HBOOK2(9,
     +   'Vertex position: z vs. x (Track enters from the rear)',
     +                    40,-HGPOS,HGPOS1,   40,-ZSIZE,ZSIZE1,   0.)
      ENDIF
C
C ----------------------------------------------------------------------
C
C--> Basic histograms, useful in most MC runs
C
      IF (FGEANT) THEN
         HGPOS1 = HGPOS + 0.1
         ZSIZE1 = ZSIZE + 0.1
        CALL HBOOK2(10,
     +  'Horizontal ground-pos (z = ZTRMAX): y vs. x',
     +  40,-HGPOS,HGPOS,   40,-HGPOS,HGPOS,   0.)
        CALL HBPRO(10,0.)
        CALL HBOOK1(20,'Log(E) at ground level (z = ZTRMAX)',
     +  40,-2.,6.,0.)
        CALL HBOOK1(30,'Z coordinate of secondary interactions',
     +  40,-ZSIZE,ZSIZE,0.)
        CALL HBOOK2(40,
     +  'Horizontal ground-pos (z = ZTRMAX): y vs. x -- muons',
     +  40,-HGPOS,HGPOS,   40,-HGPOS,HGPOS,   0.)
        CALL HBPRO(40,0.)
        CALL HBOOK2(50,
     +  'Horizontal ground-pos (z = ZTRMAX): y vs. x -- e-m',
     +  40,-HGPOS,HGPOS,   40,-HGPOS,HGPOS,   0.)
        CALL HBPRO(50,0.)
        CALL HBOOK2(60,
     +  'Horizontal ground-pos (z = ZTRMAX): y vs. x --  had',
     +  40,-HGPOS,HGPOS,   40,-HGPOS,HGPOS,   0.)
        CALL HBPRO(60,0.)
        CALL HBOOK1(70,'Long. shower profile (z at each step)',
     +  40,-ZSIZE,ZSIZE1,0.)
        CALL HBOOK2(80,'Shower shape (z vs. x at each step)',
     +  40,-HGPOS,HGPOS1,   40,-ZSIZE,ZSIZE1,   0.)
        IF (FDPST) THEN
           CALL HBOOK1(90,'Interaction mechanisms',40,0.5,40.5,0.)
        ENDIF
        CALL HBOOK1(100,'Theta(track) at ground level (z = ZTRMAX):',
     +  50,0.,1.25,0)
        CALL HBOOK1(110,
     +  'Radius (x**2 + y**2) at ground level (z = ZTRMAX):',
     +  50,0.,1.5*HGPOS,0)
        IF (KGEN .EQ. 1) THEN
           CALL HBOOK1(120,'Theta(muon) generated',40,0.,0.1,0.)
        ELSEIF (KGEN .EQ. 2) THEN
           CALL HBOOK1(120,'Theta(lepton) generated',40,0.,0.1,0.)
        ENDIF
        CALL HBOOK1(130,'E(mu)/E(nu) generated',50,0.,1.,0.)
        IF (KGEN .EQ. 1) THEN
           CALL HBOOK1(140,
CCC     +    'E(muon; z = ZTRMAX)/E(muon; origin)',    50,0.,1.,0.)
     +    'E(track; z = ZTRMAX)/E(track; origin)',    50,0.,1.,0.)
        ELSEIF (KGEN .EQ. 2) THEN
           CALL HBOOK1(140,
CCC     +    'E(muon; z = ZTRMAX)/E(parent neutrino)', 50,0.,1.,0.)
     +    'E(track; z = ZTRMAX)/E(parent neutrino)', 50,0.,1.,0.)
        ENDIF
        IF (FDPST) THEN
           CALL HBOOK1(150,'Cos(scattering angle)',40,0.8,1.0,0.)
           CALL HBOOK1(160,'dE/E in tracking step',40,0.,1.,  0.)
           CALL HBOOK1(161,'dE/E in tracking step',40,0.,0.04,0.)
           CALL HBOOK1(170,'Interaction mechanisms',40,0.5,40.5,0)
           CALL HBPROF(180,'dE/dx(E), mean',  40,0.,E0MAX,0.,1.,' ')
           CALL HBPROF(181,'dE/dx(E), spread',40,0.,E0MAX,0.,1.,'S')
           CALL HBOOK1(190,'Step length',40,0.,100.,0.)
           CALL HIDOPT(190,'LOGY')
        ENDIF
        CALL HBOOK1(200,'dE/dx (GeV/cm)',50,0.,0.10,0.)
      ENDIF
C
C ----------------------------------------------------------------------
C
C--> Step-by-step process histograms
C
      IF (FDPST) THEN
        CALL HBOOK1(802,'dE/E, muls',50,0.,0.10,0.)
        CALL HBOOK1(803,'dE/E, loss',50,0.,0.10,0.)
        CALL HBOOK1(806,'dE/E, pair',50,0.,0.10,0.)
        CALL HBOOK1(809,'dE/E, brem',50,0.,0.10,0.)
        CALL HBOOK1(821,'dE/E, munu',50,0.,0.10,0.)
        CALL HBOOK1(822,'dE/E, tofm',50,0.,0.10,0.)
        CALL HBOOK1(902,'Cos(scattering angle), muls',40,0.8,1.00,0.)
        CALL HBOOK1(903,'Cos(scattering angle), loss',40,0.8,1.00,0.)
        CALL HBOOK1(906,'Cos(scattering angle), pair',40,0.8,1.00,0.)
        CALL HBOOK1(909,'Cos(scattering angle), brem',40,0.8,1.00,0.)
        CALL HBOOK1(921,'Cos(scattering angle), munu',40,0.8,1.00,0.)
        CALL HBOOK1(922,'Cos(scattering angle), tofm',40,0.8,1.00,0.)
      ENDIF
C
C ----------------------------------------------------------------------
C
C--> Cerenkov photon generation histograms
C
      IF (KCERAL .EQ. 1) THEN
        CALL HBOOK1(210,'Radiating path/muon path',        40,0.,20.,0.)
        CALL HBPROF(220,'C-path/mupath(E), mean',40,0.,E0MAX,0.,20.,' ')
        CALL HBPROF(230,'C-path/mupath(E), spread',
     +                                           40,0.,E0MAX,0.,20.,'S')
        CALL HBOOK1(240,'Total # Cerenkov photons (in each layer)',
     +                                                    40,0.,20e6,0.)
        CALL HBOOK2(250,'Cerenkov vector: y vs. x',
     +                                           40,-1.,1.,40,-1.,1.,0.)
        CALL HBOOK1(260,'Cerenkov hit: R',              100,0.,HGPOS,0.)
        CALL HBOOK2(270,'Cerenkov hit: y vs. x',
     +                         40,-HGPOS,HGPOS,   40,-HGPOS,HGPOS,   0.)
C
CCC        NBINS = 2*NINT(0.02*ZSIZE)
        NBINS  = 160
        HHSIZE = 2000.0
        BINWID = 2.0*HHSIZE/FLOAT(NBINS)
        WRITE(6,*) 'NBINS,HHSIZE,BINWID = ',NBINS,HHSIZE,BINWID
C
        CHTIT1 = 'Hit pattern: y vs. x -- Layer xx'
        DO 10 ILAY = 1,NLAY
           WRITE(CHTIT1(31:32),'(I2)') ILAY
           CALL HBOOK2(280+ILAY,CHTIT1,
     +               NBINS,-HHSIZE,HHSIZE,   NBINS,-HHSIZE,HHSIZE,   0.)
   10   CONTINUE
C
        CALL HBOOK1(290,'Cos(Angle photon/PM)',             40,0.,1.,0.)
        CALL HIDOPT(290,'LOGY')
        CALL HBOOK1(300,'Arrival time of photons',       40,0.,400.0,0.)
        CALL HIDOPT(300,'LOGY')
        CALL HBOOK1(301,'Arrival time of photons',        40,0.,40.0,0.)
        CALL HIDOPT(301,'LOGY')
        CALL HBPROF(310,'Tof vs. R',            40,0.,2000.,0.,1E-5,'S')
        CALL HIDOPT(310,'LOGY')
        CALL HIDOPT(310,'PERR')
C
        CHTIT2 = 'Hit time pattern: y vs. x -- Layer xx'
        DO 20 ILAY = 1,NLAY
           WRITE(CHTIT2(36:37),'(I2)') ILAY
           CALL HBOOK2(320+ILAY,CHTIT2,
     +               NBINS,-HHSIZE,HHSIZE,   NBINS,-HHSIZE,HHSIZE,   0.)
   20   CONTINUE
C
        CHTIT3 = 'Hit time**2 pattern: y vs. x -- Layer xx'
        DO 30 ILAY = 1,NLAY
           WRITE(CHTIT3(39:40),'(I2)') ILAY
           CALL HBOOK2(330+ILAY,CHTIT3,
     +               NBINS,-HHSIZE,HHSIZE,   NBINS,-HHSIZE,HHSIZE,   0.)
   30   CONTINUE
C
        CALL HBOOK1(340,'Time-sigma in bin',               50,0.,25.,0.)
        CALL HBOOK1(350,'Time-average in bin',             50,0.,50.,0.)
        CALL HBOOK2(360,'Time vs R',    40,0.,2000.,   40,0.,100.,   0.)
        CALL HBOOK2(370,'Time-average vs. R',   50,0.,50.,40,0.,100.,0.)
        CALL HBOOK2(380,'Log10(W) vs. R',  40,0.,4000.,  40,-5.,0.,  0.)
        CALL HBOOK2(390,'J vs. I, t < 2*R',  30,0.,60.,  30,0.,60.,  0.)
        CALL HBOOK2(400,'J vs. I, t > 2*R',  30,0.,60.,  30,0.,60.,  0.)
C
        CALL HBOOK1(700,'t, r < 0.2',                   50, 0., 10., 0.)
        CALL HBOOK1(701,'t, r = 1.0 -- 1.2',            50, 0., 10., 0.)
        CALL HBOOK1(702,'t, r = 2.0 -- 2.2',            50, 0., 20., 0.)
        CALL HBOOK1(703,'t, r = 5.0 -- 5.2',            50, 0., 50., 0.)
        CALL HBOOK1(704,'t, r = 10 -- 11',              50, 0.,100., 0.)
        CALL HBOOK1(705,'t, r = 20 -- 21',              50, 0.,100., 0.)
      ENDIF
C
C ----------------------------------------------------------------------
C
C-->  Backscattering histograms, tracking
C
      IF ((NZH .GT. 0) .AND. FTHEZH) THEN
        CALL HBOOK1(500,'Theta, all',                     40,0.,0.4,0.)
        CALL HBOOK1(510,'Theta, last step',               40,0.,0.4,0.)
        CALL HBPROF(520,'Theta vs. energy',   40,0.,E0MAX,   0.,1.,' ')
        CALL HBPROF(530,'Theta vs. energy',   40,0.,E0MAX,   0.,1.,'S')
        CALL HBOOK2(540,'Theta vs. energy',   40,0.,E0MAX,40,0.,0.4,0.)
        CALL HBOOK2(550,'Theta vs distance',  40,0.5,80.5,40,0.,0.4,0.)
        CALL HBOOK2(610,'ThetaGMOL vs Emu',40,0.,E0MAX,40,0.01,0.41,0.)
        CALL HBOOK1(620,'MAX(E*Theta)',            40, 0., 20000.0, 0.)
      ENDIF
C
C ----------------------------------------------------------------------
C
C-->  Backscattering histograms, second stage
C
      IF (FRSCAT) THEN
        CALL HBOOK1(1010,'Theta$',40,0.0,2.0,0)
        CALL HBOOK1(1020,'Z$',40,0.,400.,0)
        CALL HBOOK2(1030,'Theta vs Z$',60,0.,600.,50,-0.5,+0.5,0)
        CALL HBOOK2(1040,'Theta vs Z$',60,0.,600.,50,0.0,1.0,0)
        CALL HBSLIY(1040,10,0)
        CALL HBOOK1(1050,'ZMAX$',40,0.,200.,0)
        CALL HBOOK1(1060,'Upwards from ZMAX$',40,0.,80.,0)
        CALL HBOOK1(1065,'MAX upwards from ZMAX$',40,0.,80.,0)
        CALL HBOOK2(1070,'ZMAX vs up$',40,0.,400.,60,0.,60.,0)
        CALL HBOOK1(1080,'MAX(E*THETA)$',40,0.,20000.0,0)
        CALL HBOOK2(1090,'MAX(ETH) vs Up$',40,0.,20000.,40,0.,40.,0)
        CALL HBOOK2(1100,'E0 vs up$',40,0.,3000.,60,0.,60.,0)
        CALL HBOOK1(1110,'Fraction of scattering in largest single$'
     *               ,50,0.,1.,0)
        CALL HBOOK2(1120,'Singlefrac vs singlesize$'
     *               ,20,0.,1.,20,0.,1.,0)
        CALL HBOOK2(1125,'Singlefrac vs TOTALsize$'
     *               ,20,0.,1.,20,0.,1.,0)
        CALL HBOOK1(1130,'JZMAX$',40,0.,200.,0)
        CALL HBOOK2(1140,'ZMAX vs X(ZMAX)$',40,0.,4000.,40,0.,400.,0)
        CALL HBOOK2(1150,'THMAX vs E0$',20,0.,3000.,65,0.,1.3,0)
        CALL HBOOK2(1160,'Theta vs E$',60,0.,3000.,40,0.,0.8,0)
        CALL HBOOK2(1170,'Z vs E$',60,0.,3000.,40,0.,400.,0)
        CALL HBOOK1(1180,'Up through surface, E.gt.2$',45,0.,90.,0)
        CALL HBOOK1(1190,'Up through surface, E.gt.5$',45,0.,90.,0)
        CALL HBOOK1(1200,'f(E)$',300,0.,3000.,0)
        CALL HIDOPT(0,'INTE')
        CALL HIDOPT(0,'STAT')
      ENDIF
C
C ----------------------------------------------------------------------
C
C-->  User histograms
C
      CALL NUUHB
C
C ----------------------------------------------------------------------
C
      IF (IHBUNI .LE. 0) THEN
C       Obviously we shouldn't open the HBOOK output file if we have
C       chosen not to write any HBOOK stuff out. - Also: it's not all
C       that useful to book an NTUPLE if we're not going to write it
C       out, so we'll skip the NTUPLE booking in that case.
         GOTO 1001
      ENDIF
C
C--> Open the histogram/ntuple output file:
C
      CALL HROPEN(IHBUNI,'NUSIM',HISFIL,'N',1024,IERR)
      IF (IERR .NE. 0) THEN
         WRITE(6,*) ' '
         WRITE(6,*) '*** NUHBOK:'
         WRITE(6,*) '*** WARNING: Return code',IERR,' from HROPEN'
      ENDIF
      CALL RZQUOT(60000)
C    Possibly use different names for different applications in
C    the future:
      CALL HMDIR('NTUPLE','S')
      CALL RZQUOT(50000)
      CALL RZNDIR(CHNPAH,'R')
      CALL RZCDIR(CHCPAH,'R')
      CALL RZNDIR(CHNPAH,'P')
      CALL RZCDIR(CHCPAH,'P')
C
C ----------------------------------------------------------------------
C
C-->  NTUPLE booking
C
      IF ((KCERAL .LE. 0) .OR. (KGEN .NE. 1)) THEN
C       The NTUPLE stuff has only been implemented for the case
C       KCERAL >= 0 (Cherenkov light generated) and KGEN = 1 (single
C       track) - it doesn't seem to make sense for other cases.
         GOTO 1001
      ENDIF
C
      IF (ICONTY .EQ. 1) THEN
         NLYSTR = NLAY
         CLYSTT = 'NPMLAx'
         CLYSTD = 'NPDLAx'
         CLYSTU = 'NPULAx'
         CLSTOT = 'NLYTOT'
         CLSDOW = 'NLYDOW'
         CLSUPP = 'NLYUPP'
      ELSEIF (ICONTY .EQ. 2) THEN
         NLYSTR = NSTR
         CLYSTT = 'NPMSTx'
         CLYSTD = 'NPDSTx'
         CLYSTU = 'NPUSTx'
         CLSTOT = 'NSTTOT'
         CLSDOW = 'NSTDOW'
         CLSUPP = 'NSTUPP'
      ELSEIF (ICONTY .EQ. 3) THEN
         WRITE(6,*) ' '
         WRITE(6,*) '*** NUHBOK:'
         WRITE(6,*) '*** ICONTY = ',ICONTY
         WRITE(6,*)
     +   '*** Neither string configuration nor layer configuration -'
         WRITE(6,*) '*** Assuming ',NSTR,' STRINGS for the NTUPLE'
         NLYSTR = NSTR
         CLYSTT = 'NPMSTx'
         CLYSTD = 'NPDSTx'
         CLYSTU = 'NPUSTx'
         CLSTOT = 'NSTTOT'
         CLSDOW = 'NSTDOW'
         CLSUPP = 'NSTUPP'
      ENDIF
C
      CNUVAR(NNUVA1 - 5) = CLSTOT
      CNUVAR(NNUVA1 - 4) = CLSDOW
      CNUVAR(NNUVA1 - 3) = CLSUPP
C
      NNUVAR = NNUVA1 + 3*NLYSTR
C
      IF (NLYSTR .GT. MXLYST) THEN
         WRITE(6,*) ' '
         WRITE(6,*) '*** NUHBOK:'
         WRITE(6,*) '*** Too many layers or strings: NLYSTR = ',NLYSTR
         WRITE(6,*) '*** Largest allowed number is',MXLYST
         WRITE(6,*) '*** NTUPLE not booked'
         GOTO 1001
      ENDIF
C
      NOFFS1 = NNUVA1 + NLYSTR
      NOFFS2 = NNUVA1 + 2*NLYSTR
      DO 100 ILS = 1,NLYSTR
         IF (ILS .LE. 9) THEN
            WRITE(CLYSTT(6:6),'(I1)') ILS
            WRITE(CLYSTD(6:6),'(I1)') ILS
            WRITE(CLYSTU(6:6),'(I1)') ILS
         ELSE
            WRITE(CLYSTT(5:6),'(I2)') ILS
            WRITE(CLYSTD(5:6),'(I2)') ILS
            WRITE(CLYSTU(5:6),'(I2)') ILS
         ENDIF
         CNUVAR(NNUVA1 + ILS) = CLYSTT
         CNUVAR(NOFFS1 + ILS) = CLYSTD
         CNUVAR(NOFFS2 + ILS) = CLYSTU
  100 CONTINUE
C
      CALL HBOOKN(9001, 'HITPAT', NNUVAR, '//NUSIM/NTUPLE', 1000,
     +            CNUVAR)
C
C ----------------------------------------------------------------------
C
 1001 CONTINUE
C
      RETURN
      END
CDECK  ID>, NUZEND. 
      SUBROUTINE NUZEND
C
      COMMON/NUFLAG/FGEANT,KGEN,KFLUXE,KFLUXS,FWEI,
     *       IRNALG,IRNSAV,IEVTRS,IEVACT,NEVACT,
     *       KCERAL,FOPTIN,FDPST,FRSCAT,IVTXAL,FTISME,
     *       ZINVBS
      LOGICAL FGEANT,FWEI,FOPTIN,FDPST,FRSCAT,FTISME,ZINVBS
      COMMON/NUZZIO/IZUNII,IZUNIO,IOZH0,IOEZH,IOTZH,IOCPA,IOCP1,
     *      IOPMH,IOSDIR,IOSHHI,IOSHLO,IOLAST
      COMMON/NUFWIO/IFWIO,IWGEAN,IWZH0,IWEZH,IWTZH,IWCPA,IWCP1,IWPMH,
     +              IWLAST
      COMMON/NUHBIO/IHBUNI
      COMMON/NULINK/JZZH0,JZEZH,JZTZH,JZCPA,JZCP1,JZPMH
     *      ,JZSDIR,JZSHLO,JZSHHI,JZLAST
      PARAMETER ( MXZBIN = 50 )
      PARAMETER ( MXTBIN = 40 )
      PARAMETER ( MXSHOW = 5000 )
      COMMON/NUSHLC/SHPHOT(0:MXZBIN,0:MXTBIN)
     *             ,SHT(0:MXZBIN,0:MXTBIN),SHT2(0:MXZBIN,0:MXTBIN)
     *             ,SHLATE(0:MXZBIN,0:MXTBIN)
     *             ,NBZ,NBTH,ZBINN(0:MXZBIN),THBINN(0:MXTBIN)
     *             ,EESHOW(MXSHOW),IDXSHO(MXSHOW)
     *             ,USESHL,NWBIN,INSLIB,NSH,CGEN,NZBMAX,NTBMAX
     *             ,QQCER2,QQSLIB,NPCER2,NPSLIB
     *             ,SHOFIL,CHDIR,CHNPAH,CHCPAH,CHNPAZ,CHCPAZ
      CHARACTER*48 SHOFIL
      CHARACTER*5 CHDIR
      CHARACTER*20 CHNPAH,CHCPAH,CHNPAZ,CHCPAZ
      LOGICAL USESHL,INSLIB,CGEN
      DIMENSION KEYDIR(2)
      DATA KEYDIR/-1,-1/
C
      IERR = 0
C
      IF (IZUNII .GT. 0) THEN
         CALL GCLOSE(IZUNII,IERR)
      ELSEIF (IZUNIO .GT. 0) THEN
         CALL GCLOSE(IZUNIO,IERR)
      ENDIF
C
      IF (IERR .NE. 0) THEN
         WRITE(6,1) '*** NUZEND: Error code ',IERR,'from GCLOSE ***'
      ENDIF
    1 FORMAT(' '/' ',A,I5,A)
C
      CALL RZOUT(IXSTOR,JZSDIR,KEYDIR,ICYCLE,'W')
      CALL RZEND(CHDIR)
      CLOSE(61)
      RETURN
      END
CDECK  ID>, NU1MU.  
      SUBROUTINE NU1MU
*
* Give a muon to Geant
*
* Use GSVERT and GSKINE
* GSVERT stores vertex parameters,
* Parameters: array(x,y,z) of vertex, beam track number,
*             target track number, user array, new vertex number
* GSKINE stores long life track params.
* Parameters: comp. of momentum, type of particle, vertex number,
*             user array, track number
*
      COMMON/NUMULT/ NGAM,NELE,NCHA,NNEUT,NMUON,MGTYPE(100)
     *              ,MGAM,MELE,MCHA,MNEUT,MMUON,FCONV,FHADR
      LOGICAL FCONV,FHADR
      COMMON/NUGRID/ZSIZE,HGPOS,DISTX,DISTY,DISTZ,OFFLAY,OFFLX,OFFLY,
     *              OFFX,OFFY,NBINS,NLAY,NSIDE,NSTR,NPMSTR,IDCONF,
     *              ICONTY,PMDIAM,PMRAD,PMR2,PMAREA,XYFUZZ,
     *              VXPLSI, VXPLDI
      COMMON/NUMISC/START,  IPTYPE, PIN, EMU0, THEMU0, PHIMU0,
     +              ZTRMAX, IORIVX, IORITR, VOLD(7),
     +              KTMED,  REFIND, ATTLEN, CWLGT1, CWLGT2,
     +              CUTKC2, IPRHIS(3), TDIAM, INUTRA,
     +              ZCLAPP, DCLAPP
      INTEGER START
      COMMON/NUFLAG/FGEANT,KGEN,KFLUXE,KFLUXS,FWEI,
     *       IRNALG,IRNSAV,IEVTRS,IEVACT,NEVACT,
     *       KCERAL,FOPTIN,FDPST,FRSCAT,IVTXAL,FTISME,
     *       ZINVBS
      LOGICAL FGEANT,FWEI,FOPTIN,FDPST,FRSCAT,FTISME,ZINVBS
      COMMON/NUGENC/E0MIN,E0MAX,TH0MIN,TH0MAX,CTHMIN,CTHMAX,PHI0,
     +              NPPROJ,NPTARG
      PARAMETER ( MAXZH = 500 )
      COMMON/NUZH/NZH,JZH,DISTZH,ZH(MAXZH),EZH(MAXZH),THEZH(MAXZH),
     +            CPATH(MAXZH),CPTOT,HISTZH,IHZH,FTHEZH,ETHLIM
      LOGICAL HISTZH,FTHEZH
C
C Pick up GEANT constants - pi and c, for example:
      COMMON/GCONST/PI,TWOPI,PIBY2,DEGRAD,RADDEG,CLIGHT,BIG,EMASS
      COMMON/GCONSX/EMMU,PMASS,AVO
C
C
      DIMENSION VERT(3), PLAB(4),
     +          VTXPLA(3), VTXROT(3,3), VXRINV(3,3), I1(3)
C
      EQUIVALENCE (VTXROT(1,1), VXRINV(1,1))
C
      INTEGER NVXTRA
      DATA NVXTRA/0/
C
      SAVE NVXTRA, VTXROT, VXRINV
C
      IF (KFLUXE .EQ. 0) THEN
         EMU0 = E0MIN
      ELSEIF (KFLUXE .EQ. 1) THEN
         EMU0  = E0MIN + RNUNIF(JJ)*(E0MAX-E0MIN)
      ELSEIF (KFLUXE .EQ. 2) THEN
         WRITE(6,*) ' '
         WRITE(6,*) '*** NU1MU: KFLUXE =',KFLUXE
         WRITE(6,*) '*** ==> Power law for energy distr.: ',
     +              'NOT IMPLEMENTED'
         RETURN
      ELSEIF (KFLUXE .EQ. 3) THEN
C       The weight returned by FEFLUX is not yet used (remains to
C       be implemented).
         EWGT = FEFLUX(EMU0)
      ELSEIF ( KFLUXE .EQ. 4 ) THEN
C       (* Flat in log(E) *)
        EMU0 = E0MIN*EXP(RNUNIF(JJ)*(LOG(E0MAX)-LOG(E0MIN)))
      ENDIF
      PINSQ = (EMU0 + EMMU)*(EMU0 - EMMU)
        PINSQ = (EMU0 + EMASS)*(EMU0 - EMASS)
      PIN   = SQRT(PINSQ)
C
      IF (KFLUXS .EQ. 0) THEN
         THEMU = TH0MIN
         PHI   = PHI0
         THEMU = THEMU*DEGRAD
         PHI   = PHI*DEGRAD
      ELSEIF (KFLUXS .EQ. 1) THEN
         CTHEMU = CTHMIN + RNUNIF(PIN)*(CTHMAX-CTHMIN)
         THEMU  = ACOS(CTHEMU)
         PHI    = 360.*RNUNIF(THEMU)
         PHI    = PHI*DEGRAD
      ELSEIF (KFLUXS .EQ. 2) THEN
         WRITE(6,*) ' '
         WRITE(6,*) '*** NU1MU: KFLUXS =',KFLUXS
         WRITE(6,*) '*** ==> Atmospheric NEUTRINOS!'
         WRITE(6,*) '*** NU1MU should not be called in this case!'
CCC         THEMU = THEMU*DEGRAD
CCC         PHI   = PHI*DEGRAD
         RETURN
      ELSEIF (KFLUXS .EQ. 3) THEN
         WRITE(6,*) ' '
         WRITE(6,*) '*** NU1MU: KFLUXS =',KFLUXS
         WRITE(6,*) '*** ==> Atmospheric muons: NOT IMPLEMENTED'
CCC         THEMU = THEMU*DEGRAD
CCC         PHI   = PHI*DEGRAD
         RETURN
      ELSEIF (KFLUXS .EQ. 4) THEN
C       The weight returned by FSFLUX is not yet used (remains to
C       be implemented).
         ANGWGT = FSFLUX(THEMU,PHI)
         THEMU  = THEMU*DEGRAD
         PHI    = PHI*DEGRAD
      ENDIF
C
      STHEMU = SIN(THEMU)
      CTHEMU = COS(THEMU)
      SPHIMU = SIN(PHI)
      CPHIMU = COS(PHI)
C
      PLAB(1) = STHEMU*CPHIMU*PIN
      PLAB(2) = STHEMU*SPHIMU*PIN
      PLAB(3) = CTHEMU*PIN
C
CCC      PTMU = SQRT(PLAB(1)*PLAB(1) + PLAB(2)*PLAB(2))
C
CCC      CALL HFILL(120,THEMU,0.,1.)
C
      THEMU0 = THEMU
      PHIMU0 = PHI
C
C       ELECTRON:
        IPTYPE = 3
C if we want to generate muons's direction is along the vertics:
C+ADD,P=NUDECK,D=NU1MU,C=288.
C        VERT(1) = 0.
C        VERT(2) = 0.
C        VERT(3) = 0.
C
C    Generate the primary-vertex position:
C    -------------------------------------
C
      IF (IVTXAL .EQ. 1) THEN
C
C       First generate a random point inside the GEANT volume:
         X0TR = -HGPOS + RNUNIF(PHI)*2.*HGPOS
         Y0TR = -HGPOS + RNUNIF(X0TR)*2.*HGPOS
         Z0TR = -ZSIZE + RNUNIF(Y0TR)*2.*ZSIZE
C       Extrapolate the track backwards from this point to the sides of
C       the GEANT volume, and check trough which surface it enters:
C       -- z-plane:
         IF (PLAB(3) .GT. 0.) THEN
            ZWALL = -ZSIZE
            IVXFIL = 4
         ELSE
            ZWALL =  ZSIZE
            IVXFIL = 5
         ENDIF
         STRACK = (Z0TR - ZWALL)*PIN/PLAB(3)
         IF (STRACK .LT. 0.) WRITE(6,*) '*** NU1MU: STRACK = ',STRACK
         XWALL = X0TR - STRACK*PLAB(1)/PIN
         YWALL = Y0TR - STRACK*PLAB(2)/PIN
         IF (((-HGPOS .LE. XWALL) .AND. (XWALL .LE. HGPOS))  .AND.
     +       ((-HGPOS .LE. YWALL) .AND. (YWALL .LE. HGPOS))) THEN
            CALL HFILL(IVXFIL,XWALL,YWALL,1.)
            GOTO 10
         ENDIF
C       -- x-plane:
         IF (PLAB(1) .GT. 0.) THEN
            XWALL = -HGPOS
            IVXFIL = 6
         ELSE
            XWALL =  HGPOS
            IVXFIL = 7
         ENDIF
         STRACK = (X0TR - XWALL)*PIN/PLAB(1)
         IF (STRACK .LT. 0.) WRITE(6,*) '*** NU1MU: STRACK = ',STRACK
         YWALL = Y0TR - STRACK*PLAB(2)/PIN
         ZWALL = Z0TR - STRACK*PLAB(3)/PIN
         IF (((-HGPOS .LE. YWALL) .AND. (YWALL .LE. HGPOS))  .AND.
     +       ((-ZSIZE .LE. ZWALL) .AND. (ZWALL .LE. ZSIZE))) THEN
            CALL HFILL(IVXFIL,YWALL,ZWALL,1.)
            GOTO 10
         ENDIF
C       -- y-plane:
         IF (PLAB(2) .GT. 0.) THEN
            YWALL = -HGPOS
            IVXFIL = 8
         ELSE
            YWALL =  HGPOS
            IVXFIL = 9
         ENDIF
         STRACK = (Y0TR - YWALL)*PIN/PLAB(2)
         IF (STRACK .LT. 0.) WRITE(6,*) '*** NU1MU: STRACK = ',STRACK
         XWALL = X0TR - STRACK*PLAB(1)/PIN
         ZWALL = Z0TR - STRACK*PLAB(3)/PIN
         IF (((-HGPOS .LE. XWALL) .AND. (XWALL .LE. HGPOS))  .AND.
     +       ((-ZSIZE .LE. ZWALL) .AND. (ZWALL .LE. ZSIZE))) THEN
            CALL HFILL(IVXFIL,XWALL,ZWALL,1.)
            GOTO 10
         ENDIF
C
C       If we have come here, then something must be wrong:
         WRITE(6,*) ' '
         WRITE(6,*) '*** NU1MU: Error in positioning of vertex! ***'
         XWALL = 0.
         YWALL = 0.
         ZWALL = 0.
C
   10    CONTINUE
         VERT(1) = XWALL
         VERT(2) = YWALL
         VERT(3) = ZWALL
C
      ELSEIF (IVTXAL .EQ. 2) THEN
C
         RRR = ZSIZE - 1.0
 192     CONTINUE
         VERT(1) = -RRR*(PLAB(1)/PIN)
         VERT(2) = -RRR*(PLAB(2)/PIN)
         VERT(3) = -RRR*(PLAB(3)/PIN)
         IF  (ABS(VERT(1)) .GT. HGPOS) THEN
            RRR = ABS(RRR*HGPOS/VERT(1))
            GOTO 192
         ENDIF
         IF (ABS(VERT(2)) .GT. HGPOS) THEN
            RRR = ABS(RRR*HGPOS/VERT(2))
            GOTO 192
         ENDIF
         IF (ABS(VERT(3)) .GT. ZSIZE) THEN
            RRR = ABS(RRR*ZSIZE/VERT(3))
            GOTO 192
         ENDIF
C
      ELSEIF (IVTXAL .EQ. 3) THEN
C
C     - The American AMANDA vertex algorithm -
C
C       First generate randomly a point inside a 400*400 m**2 square,
C       in a plane perpendicular to the track and 200 metres away from
C       the centre of the detector. We do this in a 'local' coordinate
C       system, with the Z axis along the trajectory of the particle
C       and centred at the origin of the 'master' coordinate system.
         VXPLS2 = VXPLSI/2.
         VTXPLA(1) = -VXPLS2 + VXPLSI*RNUNIF(PHI)
         VTXPLA(2) = -VXPLS2 + VXPLSI*RNUNIF(THEMU)
         VTXPLA(3) = -VXPLDI
C
         IF (THEMU .LE. 1.0E-5) THEN
C          Special treatment of the case of very small angles: then
C          we don't need to do any rotation.
            VERT(1) = VTXPLA(1)
            VERT(2) = VTXPLA(2)
            VERT(3) = VTXPLA(3)
            GOTO 193
         ENDIF
C
C       Then translate this point into the 'master' coordinate system.
C       We use the same angles and rotation matrix as those defined by
C       Qin Sun in NUCER2 - only here it's the inverse of that rotation
C       matrix that we need. Note that if we are using fixed angles,
C       then we need only do the matrix inversion once. Note also that
C       Qin's ALPHA corresponds to our PHI.
         IF (NVXTRA .EQ. 0) THEN
C          Note the trick with NVXTRA: Its value is changed from zero
C          ONLY if we're using FIXED angles.
            IF (KFLUXS .EQ. 0) NVXTRA = 1
C          With these angles, the inverse transformation goes like:
C           X1 =  CTHEMU*CPHIMU*X + CTHEMU*SPHIMU*Y - STHEMU*Z
C           Y1 = -SPHIMU*X        + CPHIMU*Y
C           Z1 =  STHEMU*CPHIMU*X + STHEMU*SPHIMU*Y + CTHEMU*Z
C          I.e., the inverse of our rotation matrix is:
            VXRINV(1,1) =  CTHEMU*CPHIMU
            VXRINV(1,2) =  CTHEMU*SPHIMU
            VXRINV(1,3) = -STHEMU
            VXRINV(2,1) = -SPHIMU
            VXRINV(2,2) =  CPHIMU
            VXRINV(2,3) =  0.
            VXRINV(3,1) =  STHEMU*CPHIMU
            VXRINV(3,2) =  STHEMU*SPHIMU
            VXRINV(3,3) =  CTHEMU
C          Now invert this inverted matrix (MATIN1 is in KERNLIB) -
C          note that VTXROT and VXRINV are EQUIVALENCEd:
            CALL MATIN1(VTXROT,3,3,3,0,I1,NINDEX,DETERM)
            IF (NINDEX .GT. 0) THEN
               WRITE(6,*) ' '
               WRITE(6,*) '*** NU1MU:'
               WRITE(6,*) '*** Problem with matrix inversion (MATIN1)'
               WRITE(6,*) '*** NINDEX = ',NINDEX
            ENDIF
         ENDIF
C       Now, finally, transform our vertex coordinates to the
C       'master' reference system:
         VERT(1) =  VTXROT(1,1)*VTXPLA(1) +
     +              VTXROT(1,2)*VTXPLA(2) +
     +              VTXROT(1,3)*VTXPLA(3)
         VERT(2) =  VTXROT(2,1)*VTXPLA(1) +
     +              VTXROT(2,2)*VTXPLA(2) +
     +              VTXROT(2,3)*VTXPLA(3)
         VERT(3) =  VTXROT(3,1)*VTXPLA(1) +
     +              VTXROT(3,2)*VTXPLA(2) +
     +              VTXROT(3,3)*VTXPLA(3)
C
  193    CONTINUE
C
      ELSE
C
         WRITE(6,*) ' '
         WRITE(6,*) '*** NU1MU:'
         WRITE(6,*) '*** Undefined vertex algorithm:'
         WRITE(6,*) '*** IVTXAL = ',IVTXAL
         VERT(1) = 0.
         VERT(2) = 0.
         VERT(3) = 0.
C
      ENDIF
C
      IF ((ABS(VERT(1)) .GT. HGPOS) .OR. (ABS(VERT(2)) .GT. HGPOS) .OR.
     +    (ABS(VERT(3)) .GT. ZSIZE)) THEN
         WRITE(6,*) ' '
         WRITE(6,*) '*** NU1MU:'
         WRITE(6,*)
     +   '*** Vertex position is outside of the detector volume!'
         WRITE(6,*) '*** X(vert) = ',VERT(1)
         WRITE(6,*) '*** Y(vert) = ',VERT(2)
         WRITE(6,*) '*** z(vert) = ',VERT(3)
      ENDIF
C
      CALL GSVERT(VERT,0,0,UBUF,0,NVTX)
      IF (NVTX .EQ. 0) THEN
         WRITE(6,*) '*** NU1MU: Failure to store vertex ***'
         RETURN
      ENDIF
      IORIVX = NVTX
        CALL HFILL(2001, VERT(1), VERT(2), 1.)
CC      CALL HFILL(2002, VERT(3), 0.,      1.)
C
      CALL GSKINE(PLAB, IPTYPE, NVTX, UBUF, 0, NT)
      IF (NT .EQ. 0) THEN
         WRITE(6,*) '*** NU1MU: Failure to store track ***'
         RETURN
      ENDIF
      IORITR = NT
C
C    INUTRA is our own track counter (GEANT's ITRA doesn't properly
C    take into account the tracks in the temporary stack):
      INUTRA = 0
C
C    Find the point of closest approach between the track and the z-axis:
      IF (THEMU .GT. 1.0E-5) THEN
         CTANN = CTHEMU/STHEMU
         ZCLAPP = VERT(3) - CTANN*(VERT(1)*CPHIMU + VERT(2)*SPHIMU)
         DCLAPP = ABS(VERT(1)*SPHIMU - VERT(2)*CPHIMU)
      ELSE
         ZCLAPP = 9999999999.
         DCLAPP = SQRT(VERT(1)*VERT(1) + VERT(2)*VERT(2))
      ENDIF
C
      RETURN
      END
CDECK  ID>, NUCER2. 
      SUBROUTINE NUCER2(PVEC,PARMAS,X00,Y00,Z00,TLEN,T0,OK)
C----------------------------------------------------------------------
C Called from GUSTEP to generate Cerenkov light, with parametrization
C algorithm.
C
C Transforms the coordinates of each PMT in the Lab to the system of the
C incoming muon, with the new z-axis along the direction of motion of the
C muon.
C Calculates the time when the photons arrive at each PMT.
C Calculates the efficiency of photons of each PMT.
C
C   Authors:  Qin Sun & Sverker Johansson                    June, 1992
C_______________________________________________________________________
C
C--- (X0PM,Y0PM,Z0PM) in /NUPMPO/ (in +KEEP,NUCPM) are the coordinates of
C--- each PMT in the Lab system.
C--- (X, Y, Z) are the origin of the new coordinates which are parallell
C--- moved to the A-th random point at the muon ray.
C--- (X1,Y1,Z1) are the coordinates of each PMT in the rotated system of
C--- the incoming muon.
C--- (X00,Y00,Z00) is a random point on the muon ray.
C--- PVEC(3) is the momentum vector of the muon; the unit is GeV/c.
C--- PARMAS  is the particle mass (GeV/c**2).
C--- ATTLEN  in /NUMISC/ is the attenuation length.
C--- TLEN is the lenght of the track segment being considered (cm).
C--- T0   is the starting time (s).
C
      INTEGER OK
      REAL PVEC(3)
C
      COMMON/NUFLAG/FGEANT,KGEN,KFLUXE,KFLUXS,FWEI,
     *       IRNALG,IRNSAV,IEVTRS,IEVACT,NEVACT,
     *       KCERAL,FOPTIN,FDPST,FRSCAT,IVTXAL,FTISME,
     *       ZINVBS
      LOGICAL FGEANT,FWEI,FOPTIN,FDPST,FRSCAT,FTISME,ZINVBS
      PARAMETER ( MAXPM  = 250 )
      PARAMETER ( MXPMHI = 100*MAXPM )
      COMMON/NUCPM/PMINT(MAXPM),TFIRST(MAXPM),TAVE(MAXPM),NPMINT(MAXPM),
     +             NHITOT, IPMHIT(MXPMHI),
     +             TIMHIT(MXPMHI), CPHINT(MXPMHI), NPEINT(MXPMHI),
     +             PHPPMH(MXPMHI), ANGPMH(MXPMHI), BETPMH(MXPMHI),
     +             IPAHIT(MXPMHI)
      COMMON/NUPMPO/NNPM,X0PM(MAXPM),Y0PM(MAXPM),Z0PM(MAXPM),
     +              IORIPM(MAXPM), LAYSTR(MAXPM)
C     IORIPM specifies the orientation of the PM tubes:
C     IORIPM(i) = -1: PM tube number i is looking downward
C     IORIPM(i) = +1: PM tube number i is looking upward
C     LAYSTR(i) = the number of the layer or string that PM tube
C                 number i belongs to
      COMMON/NUGRID/ZSIZE,HGPOS,DISTX,DISTY,DISTZ,OFFLAY,OFFLX,OFFLY,
     *              OFFX,OFFY,NBINS,NLAY,NSIDE,NSTR,NPMSTR,IDCONF,
     *              ICONTY,PMDIAM,PMRAD,PMR2,PMAREA,XYFUZZ,
     *              VXPLSI, VXPLDI
      COMMON/NUMISC/START,  IPTYPE, PIN, EMU0, THEMU0, PHIMU0,
     +              ZTRMAX, IORIVX, IORITR, VOLD(7),
     +              KTMED,  REFIND, ATTLEN, CWLGT1, CWLGT2,
     +              CUTKC2, IPRHIS(3), TDIAM, INUTRA,
     +              ZCLAPP, DCLAPP
      INTEGER START
C
      PARAMETER ( MXZBIN = 50 )
      PARAMETER ( MXTBIN = 40 )
      PARAMETER ( MXSHOW = 5000 )
      COMMON/NUSHLC/SHPHOT(0:MXZBIN,0:MXTBIN)
     *             ,SHT(0:MXZBIN,0:MXTBIN),SHT2(0:MXZBIN,0:MXTBIN)
     *             ,SHLATE(0:MXZBIN,0:MXTBIN)
     *             ,NBZ,NBTH,ZBINN(0:MXZBIN),THBINN(0:MXTBIN)
     *             ,EESHOW(MXSHOW),IDXSHO(MXSHOW)
     *             ,USESHL,NWBIN,INSLIB,NSH,CGEN,NZBMAX,NTBMAX
     *             ,QQCER2,QQSLIB,NPCER2,NPSLIB
     *             ,SHOFIL,CHDIR,CHNPAH,CHCPAH,CHNPAZ,CHCPAZ
      CHARACTER*48 SHOFIL
      CHARACTER*5 CHDIR
      CHARACTER*20 CHNPAH,CHCPAH,CHNPAZ,CHCPAZ
      LOGICAL USESHL,INSLIB,CGEN
C Pick up GEANT constants - pi and c, for example:
      COMMON/GCONST/PI,TWOPI,PIBY2,DEGRAD,RADDEG,CLIGHT,BIG,EMASS
      COMMON/GCONSX/EMMU,PMASS,AVO
C
C
C Pick up some other GEANT variables:
      COMMON/GCKINE/IKINE,PKINE(10),ITRA,ISTAK,IVERT,IPART,ITRTYP
     +      ,NAPART(5),AMASS,CHARGE,TLIFE,VERT(3),PVERT(4),IPAOLD
C
C
C Fine structure constant:
      PARAMETER (ALPHA = 1./137.036)
C
      DIMENSION X(MAXPM),   Y(MAXPM),   Z(MAXPM),
     +          X1(MAXPM),  Y1(MAXPM),  Z1(MAXPM),  D(MAXPM),
     +          X1SQ(MAXPM),Y1SQ(MAXPM),Z1SQ(MAXPM),DSQ(MAXPM)
      DIMENSION TPM(MAXPM), PMI(MAXPM)
C
      DATA ICALL/0/
C
      SAVE
C
      IF (ICALL .EQ. 0) THEN
         BETAC  = 1./REFIND
         CINMED = CLIGHT/REFIND
         TWOPIA = TWOPI*ALPHA
         CWTERM = 1./CWLGT1 - 1./CWLGT2
C       (ICALL  = 1 done further down.)
      ENDIF
C
CCC      WRITE(6,*) 'PVEC = ',PVEC
C
      CALL VZERO(TPM, NNPM)
C
C---- Transform the coordinate.
C
C---- First translation:
C
      DO 60 I = 1,NNPM
         X(I) = X0PM(I) - X00
         Y(I) = Y0PM(I) - Y00
         Z(I) = Z0PM(I) - Z00
   60 CONTINUE
C
C---- Then rotation:
C
C----   The momentum of the muon is P, and its total energy is ETOT
C----   The Cherenkov angle is ANG = ACOS(1/BETA*REFIND)
C
      PSQUAR = PVEC(1)*PVEC(1) + PVEC(2)*PVEC(2) + PVEC(3)*PVEC(3)
      P      = SQRT(PSQUAR)
      IF (P .LT. CUTKC2) GOTO 999
C
      ETOT  = SQRT(PSQUAR + PARMAS*PARMAS)
      BETA  = P/ETOT
      IF (BETA .LE. BETAC) GOTO 999
C
      CXBETA = CLIGHT*BETA
C
      COSANG = 1.0/(BETA*REFIND)
      ANG    = ACOS(COSANG)
      SINANG = SIN(ANG)
      TANANG = TAN(ANG)
C
      PT = SQRT(PVEC(1)*PVEC(1) + PVEC(2)*PVEC(2))
      THETA = ATAN2(PT,PVEC(3))
      IF (PT .GT. 1.0E-3) THEN
         ALFA = ATAN2(PVEC(2),PVEC(1))
      ELSE
         ALFA = 0.
      ENDIF
C
      COSTHE = COS(THETA)
      SINTHE = SIN(THETA)
      COSALF = COS(ALFA)
      SINALF = SIN(ALFA)
C
      DO 70 I = 1,NNPM
        X1(I) =  COSTHE*COSALF*X(I) +
     +           COSTHE*SINALF*Y(I) -
     -           SINTHE*Z(I)
        Y1(I) = -SINALF*X(I) + COSALF*Y(I)
        Z1(I) =  SINTHE*COSALF*X(I) +
     +           SINTHE*SINALF*Y(I) +
     +           COSTHE*Z(I)
        X1SQ(I) = X1(I)*X1(I)
        Y1SQ(I) = Y1(I)*Y1(I)
        Z1SQ(I) = Z1(I)*Z1(I)
C----   D(I) is the distance from each PMT to the muon ray
CCC        D(I) = SQRT(X1(I)*X1(I) + Y1(I)*Y1(I))
        DSQ(I) = X1SQ(I) + Y1SQ(I)
        D(I)   = SQRT(DSQ(I))
CCC        WRITE(6,100) X1(I),Y1(I),Z1(I),X0PM(I),Y0PM(I),Z0PM(I)
CCC  100   FORMAT(1X,'X1,Y1,Z1 = ',3F9.1,'  X0,Y0,Z0 = ',3F9.1)
CCC        WRITE(6,*) 'D(I)=       ',D(I)
   70 CONTINUE
C
C----   CALCULATE THE NUMBER OF PHOTONS RECEIVED BY EACH PMT
C
C    Number of Cherenkov photons, with wavelength between CWLGT1
C    and CWLGT2, emitted per centimetre along the (muon) track
C    (the factor 1.E7 is a conversion from nm to cm):
      REFBET = REFIND*BETA
      BNTERM = (1. - 1./REFBET)*(1. + 1./REFBET)
      QQ = TWOPIA*BNTERM*CWTERM
      QQ = QQ*1.E7
      IF (QQ .LT. 0.) THEN
         WRITE(6,*) ' '
         WRITE(6,*) '*** NUCER2: Negative # photos per cm: QQ = ',QQ
         WRITE(6,*) '*** (Setting QQ = 0.)'
         QQ = 0.
      ENDIF
C
C    How much of this track segment is actually 'seen' by a PM tube?
C    (Assume, for simplicity - and for the time being - spherical tubes.)
      TSEG = PMDIAM/SINANG
      IF (TSEG .GT. TLEN) TSEG = TLEN
C
C    Number of Cherenkov photons along this track segment:
      QQSEGM = QQ*TSEG
      IF (ICALL .EQ. 0) THEN
         WRITE(6,*) ' '
         WRITE(6,*) 'NUCER2: Number of Cherenkov photons emitted ',
     +              '(as a first example):'
         WRITE(6,*) '* Per centimetre: ',QQ
         WRITE(6,*) '* Along this track segment: ',QQSEGM
         WRITE(6,*) 'Length of this track segment is ',TLEN,' cm.;'
         WRITE(6,*) 'length of track segment seen by a PM tube is ',
     +              TSEG,' cm.;'
         WRITE(6,*) 'wavelength interval is ',
     +               CWLGT1,'  -  ',CWLGT2,' nm.'
         ICALL = 1
      ENDIF
C
      TANASQ = TANANG*TANANG
      TSEGB2 = TSEG/2.
C
      DO 110 I = 1,NNPM
         ZAIMPM = Z1(I) - D(I)/TANANG
         ZAIM1  = ZAIMPM - TSEGB2
         ZAIM2  = ZAIMPM + TSEGB2
         IF ((ZAIM2 .GE. 0.) .AND. (ZAIM1 .LE. TLEN)) THEN
C           The above condition means:
C           Is this PM tube within acceptance (in Z1) for Cherenkov
C           photons emitted from the track segment being considered?
C           (If yes, we're here.)
            IF (ZAIM2 .GT. TLEN) THEN
               DTSEG = TLEN - ZAIM1
               FRASEG = DTSEG/TSEG
            ELSEIF (ZAIM1 .LT. 0.) THEN
               DTSEG = ZAIM2
               FRASEG = DTSEG/TSEG
            ELSE
               FRASEG = 1.
            ENDIF
            IF (FRASEG .GT. 1.) THEN
               WRITE(6,*) ' '
               WRITE(6,*) '*** NUCER2: FRASEG = ',FRASEG
            ENDIF
C          (FRASEG is the fraction of the track segment actually 'seen'
C           by  t h i s  particular PM tube.)
            QQPM = QQSEGM*FRASEG
            IF (D(I) .GT. PMRAD) THEN
               RAD1 = D(I) - PMRAD
               RAD2 = D(I) + PMRAD
               CAREA1 = PI*RAD1*RAD1
               CAREA2 = PI*RAD2*RAD2
               RIAREA = CAREA2 - CAREA1
               PMI(I) = QQPM*PMAREA/RIAREA
            ELSE
               PMI(I) = QQPM
            ENDIF
C          Calculate the length of the photon path, from the radiating
C          track to the PM tube - not that the attenuation over this
C          distance will be taken into account later (in the subroutine
C          NUOMEF):
            PHPATH = D(I)/SINANG
C          Add to the total 'intensity' for this PM tube ('intensity'
C          here means # Cherenkov photons emitted from the track in
C          the direction towards the tube):
            PMINT(I) = PMINT(I) + PMI(I)
         ELSE
            PMI(I) = 0.
         ENDIF
C
         IF (PMI(I) .LE. 0.) GOTO 110
C
      QQCER2 = QQCER2+PMI(I)
C       Now find a reasonable INTEGER number for the # Cherenkov
C       photons reaching the PM tube:
         CALL POISSN(PMI(I),NPHTMP,IERR)
C
         IF (NPHTMP .LE. 0) GOTO 110
C
C----   CALCULATE THE TIME OF ARRIVAL OF EACH PHOTON AT EACH PMT
C
         TPHOT  = PHPATH/CINMED
         TMU    = ZAIMPM/CXBETA
         TPM(I) = T0 + TPHOT + TMU
CCC         WRITE(6,*) 'I, TPM(I) = ',I,TPM(I)
C       Note: The smearing of the time response is (potentially) done
C       further down in this routine, when the number of photo-electrons
C       has been calculated.
C
C----   CALCULATE THE EFFICIENCY OF THE PHOTON IN EACH PMT
C
C       Calculate the angle between the photon path and the Z axis in lab.
         IF (D(I) .GT. PMRAD) THEN
            COSLAM = (-X1(I)*SINTHE + COSTHE*D(I)/TANANG) /
     /               SQRT(X1SQ(I) + Y1SQ(I) + DSQ(I)/TANASQ)
            IF (COSLAM .LT. -1.) THEN
               WRITE(6,830) COSLAM, -1.
               COSLAM = -1.
            ELSEIF (COSLAM .GT. 1.) THEN
               WRITE(6,830) COSLAM, 1.
               COSLAM = 1.
            ENDIF
            RLAMBD = ACOS(COSLAM)
         ELSE
            COSLAM = 1.
            RLAMBD = 0.
         ENDIF
         IF (IORIPM(I) .GT. 0) RLAMBD = PI - RLAMBD
C
C       Now for the calculation of the number of photo-electrons...
C
         NPELEC = 0
C
C       Take into account the efficiency of the PM tube, including the
C       wavelength-dependent attenuation of the Cherenkov light, the
C       wavelength-dependent quantum efficiency, and the geometrical
C       acceptance of the tube:
         CALL NUOMEF(I,RLAMBD,PHPATH,BETA,NPHTMP,NPELEC)
C
C      (The output from the above subroutine call is NPELEC.)
C
C     - Note: The jump goes to statement label 111, not 110.
         IF (NPELEC .LE. 0) GOTO 111
      NPCER2 = NPCER2+NPELEC
C
         IF (FTISME) THEN
C          Find the smearing of the time response of the PM tube to
C          the signal giving rise to these NPELEC photo-electrons:
            CALL NUTSME(I,NPELEC,TSMEAR)
            TPM(I) = TPM(I) + TSMEAR
         ENDIF
C
C       Calculate the 'average' arrival time of the Cherenkov photons:
         IF (TAVE(I) .EQ. 0.) THEN
            TAVE(I) = TPM(I)
         ELSE
            TAVE(I) = (TAVE(I)*FLOAT(NPMINT(I)) + TPM(I)*FLOAT(NPELEC))
     /               /FLOAT(NPMINT(I) + NPELEC)
         ENDIF
C
         IF (TPM(I) .LT. TFIRST(I)) THEN
            TFIRST(I) = TPM(I)
         ENDIF
C
         NPMINT(I) = NPMINT(I) + NPELEC
C
C     - Continue here, even if NPELEC = 0.
  111    CONTINUE
C
C     - Use this condition, instead of simply NPELEC = 0, for skipping
C       the 'hit by hit' stuff. The reason for doing this is that we
C       want to be able to play around with the efficiecies 'off-line'.
C       (Note: For 'hits' with 0 photo-electrons there will be no smea-
C       ring of the time response in the 'hit-by-hit' vector TIMHIT,
C       even if FTISME is set = TRUE.)
         IF ((PMI(I) .LE. 0.1) .AND. (NPELEC .LE. 0)) GOTO 110
C
C       Fill the 'hit by hit' vectors (note that zero skipping is
C       effective, because of the condition on PMI(I) above):
         IF (NHITOT .LT. MXPMHI) THEN
            NHITOT = NHITOT + 1
            IPMHIT(NHITOT) = I
            TIMHIT(NHITOT) = TPM(I)
            CPHINT(NHITOT) = PMI(I)
            NPEINT(NHITOT) = NPELEC
            PHPPMH(NHITOT) = PHPATH
            ANGPMH(NHITOT) = RLAMBD
            BETPMH(NHITOT) = BETA
            IPAHIT(NHITOT) = IPART
         ELSE
            WRITE(6,*) ' '
            WRITE(6,*) '*** NUCER2: More than ',MXPMHI,' PM-tube hits!'
         ENDIF
C
  110 CONTINUE
C
  830 FORMAT(' '/' ','*** ERROR in NUCER2:'/' ','*** COSLAM = ',F7.3/
     +           ' ','*** (Adjusted to ',F4.1,')')
C
  999 CONTINUE
C
      RETURN
      END
CDECK  ID>, GUKINE. 
      SUBROUTINE GUKINE
* Throw in some particles.
*
* Called in GRUN
* Use GSVERT and GSKINE
* GSVERT stores vertex parameters,
* Parameters: array(x,y,z) of vertex, beam track number,
*             target track number, user array, new vertex number
* GSKINE stores long life track params.
* Parameters: comp. of momentum, type of particle, vertex number,
*             user array, track number
*     8 = pi+, 9 = pi-, 14 = proton
C
      COMMON/LUJETS/N,K(4000,5),P(4000,5),V(4000,5)
      INTEGER N,K
      REAL P,V
      SAVE /LUJETS/
*
      COMMON/NUMISC/START,  IPTYPE, PIN, EMU0, THEMU0, PHIMU0,
     +              ZTRMAX, IORIVX, IORITR, VOLD(7),
     +              KTMED,  REFIND, ATTLEN, CWLGT1, CWLGT2,
     +              CUTKC2, IPRHIS(3), TDIAM, INUTRA,
     +              ZCLAPP, DCLAPP
      INTEGER START
      COMMON/NUMULT/ NGAM,NELE,NCHA,NNEUT,NMUON,MGTYPE(100)
     *              ,MGAM,MELE,MCHA,MNEUT,MMUON,FCONV,FHADR
      LOGICAL FCONV,FHADR
      COMMON/NUGENC/E0MIN,E0MAX,TH0MIN,TH0MAX,CTHMIN,CTHMAX,PHI0,
     +              NPPROJ,NPTARG
      COMMON/NUFLAG/FGEANT,KGEN,KFLUXE,KFLUXS,FWEI,
     *       IRNALG,IRNSAV,IEVTRS,IEVACT,NEVACT,
     *       KCERAL,FOPTIN,FDPST,FRSCAT,IVTXAL,FTISME,
     *       ZINVBS
      LOGICAL FGEANT,FWEI,FOPTIN,FDPST,FRSCAT,FTISME,ZINVBS
      COMMON/GCFLAG/IDEBUG,IDEMIN,IDEMAX,ITEST,IDRUN,IDEVT,IEORUN
     +        ,IEOTRI,IEVENT,ISWIT(10),IFINIT(20),NEVENT,NRNDM(2)
      COMMON/GCFLAX/BATCH, NOLOG
      LOGICAL BATCH, NOLOG
C
C
      DIMENSION PBEAM(4),VERT(3),PVERT(3),UBUF(10)
      DIMENSION WMIN(10)
      DIMENSION PLAB(4)
CCC      SAVE JJEV         <--- Use IEVENT in GCFLAG instead!
CCC      DATA JJEV / 0 /   <--- Use IEVENT in GCFLAG instead!
*
CCC      JJEV = JJEV + 1   <--- Use IEVENT in GCFLAG instead!
C
      INTEGER ICALL
      SAVE    ICALL
      DATA    ICALL/0/
C
      INTEGER IEVONE, IEVLAS
C
      IF (ICALL .EQ. 0) THEN
         ICALL = 1
         IEVACT = 0
         NEVACT = NEVENT
         IF ((IRNALG .EQ. 2) .AND. (IEVTRS .GT. 1)) THEN
            IEVONE = IEVTRS
            IEVLAS = NEVENT + (IEVTRS - 1)
            WRITE(6,*) ' '
            WRITE(6,*) ' '
            WRITE(6,*) 'GUKINE (first call):'
            WRITE(6,*) 'Restart of simulation requested.'
            WRITE(6,*) 'First event to be (re-)simulated is number ',
     +                  IEVTRS
            WRITE(6,*) 'Number of events requested is ',NEVENT
            WRITE(6,*) '==> Last event should be number ',IEVLAS
            WRITE(6,*) '==> Changing IEVENT from ',IEVENT,' to ',IEVTRS
            WRITE(6,*) '    and      NEVENT from ',NEVENT,' to ',IEVLAS
            WRITE(6,*) ' '
            IEVENT = IEVTRS
            NEVENT = IEVLAS
         ENDIF
      ENDIF
C
      IEVACT = IEVACT + 1
C
      IORIVX = 0
      IORITR = 0
C
C    Included here 'for safety' (INUTRA is our own track counter):
      INUTRA = 0
C
* *> Generate an event
      IF (KGEN .EQ. 1) THEN
         CALL NU1MU
      ELSEIF (KGEN .EQ. 2) THEN
         PIN = E0MIN + RNUNIF(IEVENT)*(E0MAX-E0MIN)
         CALL LEPTO(10,0.,0.,PIN,0.,0.,0.,2)
*       Interface to GEANT
         CALL NULEMU
      ENDIF
C
C    Note: for the tests on # simulated events below we should use
C    the number of actually simulated events (IEVACT), rather than
C    the event number (IEVENT), since these might differ in case we
C    did restart the simulation from event # IEVTRS > 1.
C
      IF (IEVACT .LE. 100) THEN
         WRITE(6,*) ' '
         WRITE(6,*) 'Event number ',IEVENT
         IF ((IEVACT .EQ. 100) .AND. (IEVENT .NE. NEVENT)) THEN
            WRITE(6,*) '(Last printing of event number)'
         ENDIF
         IF (IEVACT .LE. 10) THEN
            IF ((KGEN .EQ. 2) .OR. (KGEN .EQ. 3)) THEN
               CALL LULIST(1)
            ENDIF
            CALL GPVERT(0)
            CALL GPKINE(0)
         ENDIF
      ENDIF
C
* *> Occasionally print current histograms during the run
*    (Don't print them if we're about to simulate the last event)
C
      IF (IEVENT .EQ. NEVENT) GOTO 100
C
      IF (IEVACT .LE. 1000) THEN
         IF (MOD(IEVACT,IPRHIS(1)) .EQ. 0) THEN
            WRITE(6,*) ' '
            IF (IEVACT .GT. 100) WRITE(6,*) 'Event number ',IEVENT
            WRITE(6,*) 'GUKINE: Printing histograms so far'
            CALL HPRINT(0)
         ENDIF
      ELSEIF (IEVACT .LE. 5000) THEN
         IF (MOD(IEVACT,IPRHIS(2)) .EQ. 0) THEN
            WRITE(6,*) ' '
            IF (IEVACT .GT. 100) WRITE(6,*) 'Event number ',IEVENT
            WRITE(6,*) 'GUKINE: Printing histograms so far'
            CALL HPRINT(0)
         ENDIF
      ELSE
         IF (MOD(IEVACT,IPRHIS(3)) .EQ. 0) THEN
            WRITE(6,*) ' '
            IF (IEVACT .GT. 100) WRITE(6,*) 'Event number ',IEVENT
            WRITE(6,*) 'GUKINE: Printing histograms so far'
            CALL HPRINT(0)
         ENDIF
      ENDIF
C
  100 CONTINUE
C
      NGAM  = 0
      NELE  = 0
      MELE  = 0
      NMUON = 0
      MMUON = 0
      NCHA  = 0
      NNEUT = 0
      FCONV = .FALSE.
      FHADR = .FALSE.
C
      RETURN
      END
CDECK  ID>, GUTREV. 
      SUBROUTINE GUTREV
C
C +CDE,NUMULT.    ---   Not needed since the ZEBRA stuff was removed.
      PARAMETER ( MAXZH = 500 )
      COMMON/NUZH/NZH,JZH,DISTZH,ZH(MAXZH),EZH(MAXZH),THEZH(MAXZH),
     +            CPATH(MAXZH),CPTOT,HISTZH,IHZH,FTHEZH,ETHLIM
      LOGICAL HISTZH,FTHEZH
C +CDE,NUMISC.    ---   Not needed since the ZEBRA stuff was removed.
      COMMON/NUFLAG/FGEANT,KGEN,KFLUXE,KFLUXS,FWEI,
     *       IRNALG,IRNSAV,IEVTRS,IEVACT,NEVACT,
     *       KCERAL,FOPTIN,FDPST,FRSCAT,IVTXAL,FTISME,
     *       ZINVBS
      LOGICAL FGEANT,FWEI,FOPTIN,FDPST,FRSCAT,FTISME,ZINVBS
      PARAMETER ( NLMAX = 4 )
      PARAMETER ( NBB = 160 )
      PARAMETER ( NPM = 10 )
      COMMON/NUPATC/CERPAT(NLMAX,NBB,NBB),TPAT(NLMAX,NBB,NBB)
     *        ,T2PAT(NLMAX,NBB,NBB)
     *        ,TFPAT(NLMAX,NBB,NBB),TLPAT(NLMAX,NBB,NBB)
     *        ,NPHTOT(NLMAX)
      COMMON/NUPM/PM(NLMAX,-NPM:NPM,-NPM:NPM)
     *           ,BINWID,PMTBIN,KBINCX,KBINCY
     *           ,NPHPM(NLMAX,-NPM:NPM,-NPM:NPM)
     *           ,TT(NLMAX,-NPM:NPM,-NPM:NPM)
     *           ,TSIG(NLMAX,-NPM:NPM,-NPM:NPM)
     *           ,TFL(NLMAX,-NPM:NPM,-NPM:NPM)
C +CDE,GCFLAG.    ---   Not needed since the ZEBRA stuff was removed.
C +CDE,NULINK.    ---   Not needed since the ZEBRA stuff was removed.
C +CDE,NUZLIN.    ---   Not needed since the ZEBRA stuff was removed.
C +CDE,NUZZIO.    ---   Not needed since the ZEBRA stuff was removed.
C +CDE,NUBANK.    ---   Not needed since the ZEBRA stuff was removed.
      PARAMETER ( MAXPM  = 250 )
      PARAMETER ( MXPMHI = 100*MAXPM )
      COMMON/NUCPM/PMINT(MAXPM),TFIRST(MAXPM),TAVE(MAXPM),NPMINT(MAXPM),
     +             NHITOT, IPMHIT(MXPMHI),
     +             TIMHIT(MXPMHI), CPHINT(MXPMHI), NPEINT(MXPMHI),
     +             PHPPMH(MXPMHI), ANGPMH(MXPMHI), BETPMH(MXPMHI),
     +             IPAHIT(MXPMHI)
      COMMON/NUPMPO/NNPM,X0PM(MAXPM),Y0PM(MAXPM),Z0PM(MAXPM),
     +              IORIPM(MAXPM), LAYSTR(MAXPM)
C     IORIPM specifies the orientation of the PM tubes:
C     IORIPM(i) = -1: PM tube number i is looking downward
C     IORIPM(i) = +1: PM tube number i is looking upward
C     LAYSTR(i) = the number of the layer or string that PM tube
C                 number i belongs to
C +CDE,NUGRID.    ---   Not needed since the ZEBRA stuff was removed.
      COMMON /NUTIMC/ TRGATE
      REAL TRGATE
C     TRGATE = Trigger gate: time interval (in seconds) during which
C              the acqusition system is open for recording data; i.e.,
C              sensitive to PM-tube signals.
      COMMON /NUNOIC/ RTNOIS
      REAL RTNOIS
C     RTNOIS = Rate of noise hits in a PM tube (# hits/second).
      PARAMETER ( MAXQQ  = 5000.)
      PARAMETER ( MAXTIM = 100.)
      PARAMETER ( MAXPS  = 5000)
      COMMON /MYSHOC/ XPRTRA(MAXPS), YPRTRA(MAXPS), ZPRTRA(MAXPS),
     +                DXPTRA(MAXPS), DYPTRA(MAXPS), DZPTRA(MAXPS),
     +                PEMAS1,PTOT1(MAXPS),PTOT2,IPS, IPSTRA,
     +                ZSHMAX,ZDIST(MAXPS),ETOT1(MAXPS),TOTPHO,
     +                DEGTHE(MAXQQ),
     +                GARMS, ELES,NEUS,MUS,NPBS,OTHERS,PIO,PIS
      PARAMETER ( MXZBIN = 50 )
      PARAMETER ( MXTBIN = 40 )
      PARAMETER ( MXSHOW = 5000 )
      COMMON/NUSHLC/SHPHOT(0:MXZBIN,0:MXTBIN)
     *             ,SHT(0:MXZBIN,0:MXTBIN),SHT2(0:MXZBIN,0:MXTBIN)
     *             ,SHLATE(0:MXZBIN,0:MXTBIN)
     *             ,NBZ,NBTH,ZBINN(0:MXZBIN),THBINN(0:MXTBIN)
     *             ,EESHOW(MXSHOW),IDXSHO(MXSHOW)
     *             ,USESHL,NWBIN,INSLIB,NSH,CGEN,NZBMAX,NTBMAX
     *             ,QQCER2,QQSLIB,NPCER2,NPSLIB
     *             ,SHOFIL,CHDIR,CHNPAH,CHCPAH,CHNPAZ,CHCPAZ
      CHARACTER*48 SHOFIL
      CHARACTER*5 CHDIR
      CHARACTER*20 CHNPAH,CHCPAH,CHNPAZ,CHCPAZ
      LOGICAL USESHL,INSLIB,CGEN
      COMMON/NUMISC/START,  IPTYPE, PIN, EMU0, THEMU0, PHIMU0,
     +              ZTRMAX, IORIVX, IORITR, VOLD(7),
     +              KTMED,  REFIND, ATTLEN, CWLGT1, CWLGT2,
     +              CUTKC2, IPRHIS(3), TDIAM, INUTRA,
     +              ZCLAPP, DCLAPP
      INTEGER START
C
      REAL RTNOGA, RTNOEV, T0GATE, TFHIT, TNOIPE
      SAVE RTNOGA, RTNOEV
C
      INTEGER NNOIPE
C
      LOGICAL YESHIT
C
      INTEGER ICALL
      SAVE ICALL
      DATA ICALL/0/
C
      IF (ICALL .EQ. 0) THEN
         ICALL = 1
C       Rate of noise in a PM tube during the duration of the trigger
C       gate:
         RTNOGA = RTNOIS*TRGATE
C       Rate of noise in ANY PM tube during the duration of the trigger
C       gate:
         RTNOEV = RTNOGA*FLOAT(NNPM)
         IF ((RTNOIS .GT. 0.) .AND. (TRGATE .GT. 0.)) THEN
            WRITE(6,*) ' '
            WRITE(6,*) 'GUTREV: Initialising noise simulation.'
            WRITE(6,*)
     +      'Noise rate per PM tube (single photo-electron) = ',
     +      RTNOIS,' Hz'
         ELSE
C          (One could think of a case where the user has set both
C          RTNOEV and TRGATE negative in order to turn of the noise
C          simulation; that's why I do it this way.)
            RTNOEV = 0.
            WRITE(6,*) ' '
            WRITE(6,*) 'GUTREV: Noise simulation NOT activated.'
         ENDIF
      ENDIF
C
      JZH = 1
      DO 10 I = 1,NZH
         EZH(I)   = -1.0
         CPATH(I) =  0.0
         THEZH(I) =  0.0
 10   CONTINUE
      CPTOT = 0.0
C
      IF (KCERAL .EQ. 1) THEN
         DO 15 I = 1,NBB
            DO 17 J = 1,NBB
               DO 19 ILAY = 1,4
                  TFPAT(ILAY,I,J) =  999999.9
                  TLPAT(ILAY,I,J) = -999999.9
 19            CONTINUE
 17         CONTINUE
 15      CONTINUE
      ENDIF
C
      IF (KCERAL .NE. 0) THEN
         CALL VZERO( PMINT,  NNPM )
         CALL VZERO( TAVE,   NNPM )
         CALL VZERO( NPMINT, NNPM )
         CALL VFILL( TFIRST, NNPM, 999999.9 )
      ENDIF
C
      IF (KCERAL .EQ. 2) THEN
         IF (NHITOT .GT. 0) THEN
            CALL VZERO( IPMHIT, NHITOT )
            CALL VZERO( TIMHIT, NHITOT )
            CALL VZERO( CPHINT, NHITOT )
            CALL VZERO( NPEINT, NHITOT )
            CALL VZERO( PHPPMH, NHITOT )
            CALL VZERO( ANGPMH, NHITOT )
            CALL VZERO( BETPMH, NHITOT )
            CALL VZERO( IPAHIT, NHITOT )
            NHITOT = 0
         ENDIF
      ENDIF
            ZSHMAX = -999999.
            IPS = 0
            CALL VZERO(XPRTRA,6*MAXPS)
            TOTPHO = 0.
            QQCER2 = 0.
            QQSLIB = 0.
            NPCER2 = 0
            NPSLIB = 0
            GARMS  = 0.
            ELES   = 0.
            NEUS   = 0.
            MUS    = 0.
            PIO    = 0.
            PIS    = 0.
            NPBS   = 0.
            OTHERS = 0.
C
      CALL GTREVE
C
C        CALL HFILL (4445,STIME,0.,1.)
C      write(6,*)'NPCER2,NPSLIB = ', NPCER2,NPSLIB
            CALL HFILL(4466,QQCER2,0.,1.)
            CALL HFILL(4467,QQSLIB,0.,1.)
            CALL HFILL(4468,QQCER2+QQSLIB,0.,1.)
            CALL HFILL(4470,FLOAT(NPCER2),0.,1.)
            CALL HFILL(4471,FLOAT(NPSLIB),0.,1.)
            CALL HFILL(4472,FLOAT(NPSLIB+NPCER2),0.,1.)
C ---------------------------------------------------------------------
C
C    Now for a little while we'll be simulating noise in the PM tubes...
C
      IF (KCERAL .LE. 0) THEN
C       If we're not simulating any Cherenkov light, then we probably
C       don't want to simulate any PM-tube noise, either:
         GOTO 999
      ENDIF
C
      IF (RTNOEV .LE. 0.) THEN
C      (Here we're giving the user a possibility to switch off the noise
C       simulation.)
         GOTO 999
      ENDIF
C
C    Total number of 'single noise photo-electrons' in ALL PM tubes
C    in the current event:
      NNOIPE = 0
      CALL POISSN(RTNOEV,NNOIPE,IERR)
      IF (IERR .NE. 0) THEN
         WRITE(6,*) ' '
         WRITE(6,*) '*** GUTREV: Error code ',IERR,' from POISSN ',
     +              'in noise simulation'
         GOTO 999
      ENDIF
C
      IF (NNOIPE .LE. 0) GOTO 999
C
C    For the timing of the noise photo-electrons: find which PM tube
C    fired first (if any); later we'll assume that the hit in that
C    tube triggered the acqusition system, i.e., opened the trigger
C    gate:
C NOTE: IN REALITY WE WILL REQUIRE A COINCIDENCE BETWEEN TWO OR MORE
C PM TUBES FOR OPENING THE TRIGGER GATE; THIS SHOULD BE TAKEN INTO
C ACCOUNT HERE. NOTE ALSO THAT WE  M I G H T  THEN GET A CONICIDENCE
C OF FALSE HITS, OR OF A COMBINATION OF FALSE AND TRUE HITS.
      T0GATE = 0.
      YESHIT = .FALSE.
      TFHIT  = 99999999.
      DO 101 I = 1,NNPM
         IF (NPMINT(I) .LE. 0) GOTO 101
         YESHIT = .TRUE.
         IF (TFIRST(I) .LT. TFHIT) TFHIT = TFIRST(I)
  101 CONTINUE
      IF (YESHIT) T0GATE = TFHIT
C NOTE: THIS IS FAR FROM BEING PERFECT - T0GATE = 0 IS REALLY QUITE
C ARBITRARY IN THE CASE WHEN THERE ARE NO TRUE HITS!
C
C    Distribute the noise photo-electrons among the PM tubes:
C
      DO 102 I = 1,NNOIPE
C
C       Which PM tube is noisy?
         XNPM = RNUNIF(I)*FLOAT(NNPM)
         INPM = INT(XNPM) + 1
         IF ((INPM .LT. 1) .OR. (INPM .GT. NNPM)) THEN
            WRITE(6,*) ' '
            WRITE(6,*) '*** GUTREV: Noise in non-existing PM tube: ',
     +                 'number ',INPM
            GOTO 102
         ENDIF
C
C       When did this false photo-electron appear?
         TNOIPE = T0GATE + RNUNIF(INPM)*TRGATE
C
C       Now add this false photo-electron to the 'CPA' bank:
         IF (TAVE(INPM) .EQ. 0.) THEN
            TAVE(INPM) = TNOIPE
         ELSE
            TAVE(INPM) = (TAVE(INPM)*FLOAT(NPMINT(INPM)) + TNOIPE)
     /               /FLOAT(NPMINT(INPM) + 1)
         ENDIF
         IF (TNOIPE .LT. TFIRST(INPM)) TFIRST(INPM) = TNOIPE
         NPMINT(INPM) = NPMINT(INPM) + 1
C
C       Next, add it to the 'PMH' bank (the -1:s indicate a noise hit):
         IF (KCERAL .NE. 2) GOTO 102
         IF (NHITOT .LT. MXPMHI) THEN
            NHITOT = NHITOT + 1
            IPMHIT(NHITOT) = INPM
            TIMHIT(NHITOT) = TNOIPE
            CPHINT(NHITOT) = -1.
            NPEINT(NHITOT) =  1
            PHPPMH(NHITOT) = -1.
            ANGPMH(NHITOT) = -1.
            BETPMH(NHITOT) = -1.
            IPAHIT(NHITOT) = -1
         ELSE
            WRITE(6,*) ' '
            WRITE(6,*) '*** GUTREV: Error adding noise hits -'
            WRITE(6,*) '*** more than ',MXPMHI,' PM-tube hits!'
         ENDIF
C
  102 CONTINUE
C
C    (End of noise simulation)
C
C ---------------------------------------------------------------------
C
C   --- ZEBRA stuff moved to GUOUT! ---
C
  999 CONTINUE
C
      RETURN
      END
CDECK  ID>, GUOUT.  
      SUBROUTINE GUOUT
C
      COMMON/NUMULT/ NGAM,NELE,NCHA,NNEUT,NMUON,MGTYPE(100)
     *              ,MGAM,MELE,MCHA,MNEUT,MMUON,FCONV,FHADR
      LOGICAL FCONV,FHADR
      PARAMETER ( MAXZH = 500 )
      COMMON/NUZH/NZH,JZH,DISTZH,ZH(MAXZH),EZH(MAXZH),THEZH(MAXZH),
     +            CPATH(MAXZH),CPTOT,HISTZH,IHZH,FTHEZH,ETHLIM
      LOGICAL HISTZH,FTHEZH
      COMMON/NUMISC/START,  IPTYPE, PIN, EMU0, THEMU0, PHIMU0,
     +              ZTRMAX, IORIVX, IORITR, VOLD(7),
     +              KTMED,  REFIND, ATTLEN, CWLGT1, CWLGT2,
     +              CUTKC2, IPRHIS(3), TDIAM, INUTRA,
     +              ZCLAPP, DCLAPP
      INTEGER START
      COMMON/NUFLAG/FGEANT,KGEN,KFLUXE,KFLUXS,FWEI,
     *       IRNALG,IRNSAV,IEVTRS,IEVACT,NEVACT,
     *       KCERAL,FOPTIN,FDPST,FRSCAT,IVTXAL,FTISME,
     *       ZINVBS
      LOGICAL FGEANT,FWEI,FOPTIN,FDPST,FRSCAT,FTISME,ZINVBS
      PARAMETER ( NLMAX = 4 )
      PARAMETER ( NBB = 160 )
      PARAMETER ( NPM = 10 )
      COMMON/NUPATC/CERPAT(NLMAX,NBB,NBB),TPAT(NLMAX,NBB,NBB)
     *        ,T2PAT(NLMAX,NBB,NBB)
     *        ,TFPAT(NLMAX,NBB,NBB),TLPAT(NLMAX,NBB,NBB)
     *        ,NPHTOT(NLMAX)
      COMMON/NUPM/PM(NLMAX,-NPM:NPM,-NPM:NPM)
     *           ,BINWID,PMTBIN,KBINCX,KBINCY
     *           ,NPHPM(NLMAX,-NPM:NPM,-NPM:NPM)
     *           ,TT(NLMAX,-NPM:NPM,-NPM:NPM)
     *           ,TSIG(NLMAX,-NPM:NPM,-NPM:NPM)
     *           ,TFL(NLMAX,-NPM:NPM,-NPM:NPM)
      COMMON/GCFLAG/IDEBUG,IDEMIN,IDEMAX,ITEST,IDRUN,IDEVT,IEORUN
     +        ,IEOTRI,IEVENT,ISWIT(10),IFINIT(20),NEVENT,NRNDM(2)
      COMMON/GCFLAX/BATCH, NOLOG
      LOGICAL BATCH, NOLOG
C
      COMMON/GCLIST/NHSTA,NGET ,NSAVE,NSETS,NPRIN,NGEOM,NVIEW,NPLOT
     +       ,NSTAT,LHSTA(20),LGET (20),LSAVE(20),LSETS(20),LPRIN(20)
     +             ,LGEOM(20),LVIEW(20),LPLOT(20),LSTAT(20)
C
      COMMON/NULINK/JZZH0,JZEZH,JZTZH,JZCPA,JZCP1,JZPMH
     *      ,JZSDIR,JZSHLO,JZSHHI,JZLAST
      COMMON/NUZLIN/IXSTOR,IXDIV1,IXDIV2,NLSTRU,NLREF,LROOT
      COMMON/NUZZIO/IZUNII,IZUNIO,IOZH0,IOEZH,IOTZH,IOCPA,IOCP1,
     *      IOPMH,IOSDIR,IOSHHI,IOSHLO,IOLAST
      COMMON/NUFWIO/IFWIO,IWGEAN,IWZH0,IWEZH,IWTZH,IWCPA,IWCP1,IWPMH,
     +              IWLAST
      COMMON/NUHBIO/IHBUNI
C
C-----------------------------------------------
C Memory requirements for HBOOK, Geant, and our
C private Zebra structure.
C-----------------------------------------------
C
      PARAMETER ( LHMEM  = 400000 )
      PARAMETER ( LGMEM  = 100000 )
      PARAMETER ( LZMEM  = 600000 )
      PARAMETER ( IFFSZE = 1000   )
      COMMON/PAWC/HHH(LHMEM)
      COMMON/GCBANK/Q(LGMEM)
      COMMON/NUBANK/FENCE(16),QSPACE(LZMEM)
      DIMENSION LQN(LZMEM),IQN(LZMEM-8)
      DIMENSION QN(LZMEM-8)
      EQUIVALENCE ( QSPACE(1), LQN(1) )
      EQUIVALENCE ( IQN(1), QN(1), LQN(9) )
      COMMON/CFREAD/FFSPAC(IFFSZE)
C-----------------------------------------------
      COMMON/NUZLOG/LOGLMZ,LOGLFZ,LOGLDZ
      PARAMETER ( MAXPM  = 250 )
      PARAMETER ( MXPMHI = 100*MAXPM )
      COMMON/NUCPM/PMINT(MAXPM),TFIRST(MAXPM),TAVE(MAXPM),NPMINT(MAXPM),
     +             NHITOT, IPMHIT(MXPMHI),
     +             TIMHIT(MXPMHI), CPHINT(MXPMHI), NPEINT(MXPMHI),
     +             PHPPMH(MXPMHI), ANGPMH(MXPMHI), BETPMH(MXPMHI),
     +             IPAHIT(MXPMHI)
      COMMON/NUPMPO/NNPM,X0PM(MAXPM),Y0PM(MAXPM),Z0PM(MAXPM),
     +              IORIPM(MAXPM), LAYSTR(MAXPM)
C     IORIPM specifies the orientation of the PM tubes:
C     IORIPM(i) = -1: PM tube number i is looking downward
C     IORIPM(i) = +1: PM tube number i is looking upward
C     LAYSTR(i) = the number of the layer or string that PM tube
C                 number i belongs to
      COMMON/NUGRID/ZSIZE,HGPOS,DISTX,DISTY,DISTZ,OFFLAY,OFFLX,OFFLY,
     *              OFFX,OFFY,NBINS,NLAY,NSIDE,NSTR,NPMSTR,IDCONF,
     *              ICONTY,PMDIAM,PMRAD,PMR2,PMAREA,XYFUZZ,
     *              VXPLSI, VXPLDI
      COMMON/NUFILE/DATFIL,RNDFIN,RNDFOU,HISFIL
      CHARACTER*48  DATFIL,RNDFIN,RNDFOU,HISFIL
      COMMON/NUTRIC/LTRIG, NPEMIN
      LOGICAL LTRIG
C     LTRIG  = Trigger flag; to be set .TRUE. if the current event has
C              generated a trigger
C     NPEMIN = Minimum # photo-electrons defining a signal in a PM tube
C              (to be used in the trigger algorithm)
      PARAMETER (MXLYST = 99)
      PARAMETER (NNUVA1 = 11)
      PARAMETER (MXNUVA = NNUVA1 + 3*MXLYST)
      COMMON /NUNTUP/ NLYSTR, NNUVAR, NUNTVE(MXNUVA)
      REAL NUNTVE
C
C     MXLYST = Maximum # layers or strings that we can handle in our
C              NTUPLE
C     NNUVA1 = Number of variables for our NTUPLE, other than those
C              representing # firing PM tubes per layer or string
C     MXNUVA = Maximum # variables that we can handle in our NTUPLE
C     NLYSTR = Actual # layers or strings in the simulated set-up
C     NNUVAR = Actual # variables for our NTUPLE
C     NUNTVE = REAL vector containing the information that we fill our
C              NTUPLE with.
*
      COMMON /MYCOUN/ NTRGEV(5), NQ(100)
      PARAMETER ( MAXQQ  = 5000.)
      PARAMETER ( MAXTIM = 100.)
      PARAMETER ( MAXPS  = 5000)
      COMMON /MYSHOC/ XPRTRA(MAXPS), YPRTRA(MAXPS), ZPRTRA(MAXPS),
     +                DXPTRA(MAXPS), DYPTRA(MAXPS), DZPTRA(MAXPS),
     +                PEMAS1,PTOT1(MAXPS),PTOT2,IPS, IPSTRA,
     +                ZSHMAX,ZDIST(MAXPS),ETOT1(MAXPS),TOTPHO,
     +                DEGTHE(MAXQQ),
     +                GARMS, ELES,NEUS,MUS,NPBS,OTHERS,PIO,PIS
      PARAMETER (MAXMEC=30)
      COMMON/GCTRAK/VECT(7),GETOT,GEKIN,VOUT(7),NMEC,LMEC(MAXMEC)
     + ,NAMEC(MAXMEC),NSTEP ,MAXNST,DESTEP,DESTEL,SAFETY,SLENG
     + ,STEP  ,SNEXT ,SFIELD,TOFG  ,GEKRAT,UPWGHT,IGNEXT,INWVOL
     + ,ISTOP ,IGAUTO,IEKBIN, ILOSL, IMULL,INGOTO,NLDOWN,NLEVIN
     + ,NLVSAV,ISTORY
C
      COMMON/GCONST/PI,TWOPI,PIBY2,DEGRAD,RADDEG,CLIGHT,BIG,EMASS
      COMMON/GCONSX/EMMU,PMASS,AVO
C
C+ADD,P=NUGU,D=GUOUT,C=284.
C        CALL HFILL(4442,ZSHMAX,0.,1.)
CC        CALL HFILL(4443,ZSTMAX,0.,1.)
C         WRITE(6,*)'ZSHMAX,ZSTMAX=',ZSHMAX,ZSTMAX
C          WRITE(6,*)'beta1,clight    =',beta1,clight
C          WRITE(6,*)'STIME    =',STIME
C
      DIMENSION DUMPAT(NBB,NBB)
C
      INTEGER IRSEED, ICALL
C
      SAVE ICALL
      DATA ICALL/0/
C
C    For the NTUPLE stuff:
      INTEGER IHITLS(MXLYST), IDOWLS(MXLYST), IUPPLS(MXLYST)
C
C ---------------------------------------------------------------------
C
C Skip most of the stuff in this routine if we're not running GEANT:
C
      IF (.NOT. FGEANT) GOTO 1000
C
C ---------------------------------------------------------------------
C
C    Possibly fill some standard GEANT histograms:
      IF (NHSTA .GT. 0) THEN
         CALL GFHSTA
      ENDIF
C
      IF (NZH .GT. 0) THEN
         JZMAX = 0
         DO 25 I = 1,NZH
            IF (EZH(I) .GT. 0) THEN
               JZMAX = I
            ENDIF
 25      CONTINUE
C
         WRITE(6,*) 'JZEZH = ',JZEZH
         CALL MZBOOK(IXSTOR,JZEZH,JZEZH,1,'EZH ',1,1,JZMAX+3,IOEZH,0)
         WRITE(6,*) 'JZEZH = ',JZEZH
         IQN(JZEZH+1) = JZMAX
          QN(JZEZH+2) = PIN
          QN(JZEZH+3) = EMU0
         DO 27 I = 1,JZMAX
            QN(JZEZH+3+I) = EZH(I)
 27      CONTINUE
         IF (LOGLDZ .GE. 2) THEN
            CALL DZSHOW('EZH AFTER STORE:',IXSTOR,JZEZH,' ',0,0,1,20)
         ENDIF
C
         IF (FTHEZH) THEN
            ETH = 0.0
            DO 30 I = 1,JZMAX
               IF (EZH(I)*THEZH(I) .GT. ETH) THEN
                  ETH = EZH(I)*THEZH(I)
               ENDIF
               CALL HFILL(500, THEZH(I), 0., 1.)
               CALL HFILL(520, EZH(I),   THEZH(I), 1.)
               CALL HFILL(530, EZH(I),   THEZH(I), 1.)
               CALL HFILL(540, EZH(I),   THEZH(I), 1.)
               CALL HFILL(550, FLOAT(I), THEZH(I), 1.)
 30         CONTINUE
            IF (JZMAX .GT. 0.) THEN
               CALL HFILL(510,THEZH(JZMAX),0.,1.)
            ENDIF
            CALL HFILL(620,ETH,0.,1.)
C
            IF (ETHLIM .LT. 0.) THEN
C-->           (If ETHLIM negative, then the condition below
C-->           becomes IF ( ABS(ETH) .LT. ABS(ETHLIM) )... )
               ETH = -ETH
            ENDIF
            IF (ETH .GT. ETHLIM) THEN
               WRITE(6,100) PIN,EMU0,JZMAX
               WRITE(6,200) (EZH(I),I=1,JZMAX)
               WRITE(6,200) (1000.0*THEZH(I),I=1,JZMAX)
            ENDIF
C
            CALL MZBOOK(IXSTOR,JZTZH,JZTZH,1,'TZH ',1,1,JZMAX+2,IOTZH,0)
            IQN(JZTZH+1) = JZMAX
             QN(JZTZH+2) = ETHLIM
            DO 37 I = 1,JZMAX
               QN(JZTZH+2+I) = THEZH(I)
 37         CONTINUE
C
         ELSE
            WRITE(6,*) 'NZH,EZH(1) = ',NZH,EZH(1)
            WRITE(6,100) PIN,EMU0,JZMAX
            WRITE(6,200) (EZH(I),I=1,JZMAX)
         ENDIF
      ENDIF
C
 100  FORMAT(1X,2F10.2,I3)
 200  FORMAT(8F10.2)
C
      IF (KCERAL .EQ. 1) THEN
         CPPTOT = 0.0
         ZMAX   = 0.0
         DZ = ZH(2) - ZH(1)
CCC         WRITE(6,*) 'DZ = ',DZ
         DO 20 I = 3,NZH
            CPPTOT = CPPTOT + CPATH(I)
CCC            WRITE(6,*) 'I,CPATH,EZH = ',I,CPATH(I),EZH(I)
            IF (EZH(I) .GT. 0.) THEN
               CALL HFILL(220,EZH(I),CPATH(I)/DZ,1.)
               CALL HFILL(230,EZH(I),CPATH(I)/DZ,1.)
               ZMAX = ZH(I)
            ENDIF
 20      CONTINUE
         WRITE(6,*) 'EMU0,CPTOT,CPPTOT,NPHTOT(4),NBINS,ATTLEN = ',
     +               EMU0,CPTOT,CPPTOT,NPHTOT,NBINS,ATTLEN
C
         IF (NBINS .NE. NBB) THEN
            WRITE(6,*) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
            WRITE(6,*) '  NBB != NBINS : ',NBB,NBINS
            WRITE(6,*) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
            STOP
         ENDIF
         DO 45 ILAY = 1,NLAY
            WRITE(6,*) 'Unpacking layer ',ILAY
C
            CALL HUNPAK(280+ILAY,DUMPAT,'HIST',0)
            DO 47 I = 1,NBB
               DO 48 J = 1,NBB
                  CERPAT(ILAY,I,J) = DUMPAT(I,J)
 48            CONTINUE
 47         CONTINUE
            CALL HUNPAK(320+ILAY,DUMPAT,'HIST',0)
            DO 57 I = 1,NBB
               DO 58 J = 1,NBB
                  TPAT(ILAY,I,J) = DUMPAT(I,J)
 58            CONTINUE
 57         CONTINUE
            CALL HUNPAK(330+ILAY,DUMPAT,'HIST',0)
            DO 67 I = 1,NBB
               DO 68 J = 1,NBB
                  T2PAT(ILAY,I,J) = DUMPAT(I,J)
 68            CONTINUE
 67         CONTINUE
C
            CALL NUHIT(PMTOT,NPMTOT,ILAY)
C
            DO 70 I = 1,NBB
               DO 75 J = 1,NBB
                  IF ( CERPAT(ILAY,I,J) .GT. 0.1 ) THEN
                     TTT = TPAT(ILAY,I,J)/CERPAT(ILAY,I,J)
                     V = (T2PAT(ILAY,I,J)-(TPAT(ILAY,I,J)**2)
     +                   /CERPAT(ILAY,I,J))/CERPAT(ILAY,I,J)
                     IF (V .LT. 0.) THEN
                        WRITE(6,*) 'V,TPAT,T2PAT,CERPAT = ',
     +                              V,TPAT(ILAY,I,J),T2PAT(ILAY,I,J),
     +                              CERPAT(ILAY,I,J)
                        V = 1.0
                     ENDIF
                     CALL HFILL(340,SQRT(V),0.,CERPAT(ILAY,I,J))
                     CALL HFILL(350,TTT,    0.,CERPAT(ILAY,I,J))
                     R = SQRT((I-0.5*NBB)**2 + (J-0.5*NBB)**2)
                     CALL HFILL(370,R,TTT,CERPAT(ILAY,I,J))
                     IF (TTT .LT. 2.*R) THEN
                        CALL HFILL(390,FLOAT(I),FLOAT(J),1.)
                     ELSE
                        CALL HFILL(400,FLOAT(I),FLOAT(J),1.)
                     ENDIF
                  ELSE
                     TTT = 0.0
                     V   = 1.0
                  ENDIF
                  TPAT(ILAY,I,J)  = TTT
                  T2PAT(ILAY,I,J) = V
 75            CONTINUE
 70         CONTINUE
CCC            CALL HROUT(280,ICYCLE,' ')
CCC            WRITE(6,*) 'GUOUT: AFTER HROUT; ICYCLE = ',ICYCLE
            CALL HRESET(280+ILAY,' ')
            CALL HRESET(320+ILAY,' ')
            CALL HRESET(330+ILAY,' ')
            IF (ZMAX .GT. 0.) THEN
               CALL HFILL(210,CPPTOT/ZMAX,0.,1.)
            ENDIF
            CALL HFILL(240,FLOAT(NPHTOT(ILAY)),0.,1.)
 45      CONTINUE
         WRITE(6,*) 'GUOUT: Output starting'
         NBBSQ = NBB*NBB
         NW = 5 + 5*NLAY*NBBSQ
         WRITE(6,*) 'NW for CP1 = ',NW
         CALL MZBOOK(IXSTOR,JZCP1,JZCP1,1,'CP1 ',1,1,NW,IOCP1,0)
CCC
CCC      NOTE: PPP, THETA, AND PHI HAVE NOT BEEN INITIALISED! THEY'RE
CCC      PART OF SVERKER J'S CODE, AND I DON'T KNOW WHAT THEY SHOULD
CCC      BE SET TO. AS WE'RE CURRENTLY NOT USING SVERKER'S CODE, I
CCC      HAVE CHOSEN THE FOLLOWING SIMPLE-MINDED WAY OF GETTING RID
CCC      OF SOME COMPILATION WARNINGS. - IF WE START USING SVERKER'S
CCC      CODE AGAIN, WE SHOULD OF COURSE FIND A MORE APPROPRIATE SOLU-
CCC      TION TO THE PROBLEM.
CCC                                                  B. E. 930506.
         PPP   = 0.
         THETA = 0.
         PHI   = 0.
         IF (ICALL .EQ. 0) THEN
            ICALL = 1
            WRITE(6,*)' '
            WRITE(6,*)'************************************************'
            WRITE(6,*)'***              GUOUT: WARNING!             ***'
            WRITE(6,*)'*** UNINITIALISED VARIABLES: PPP, THETA, PHI ***'
            WRITE(6,*)'*** SETTING THEM ALL TO ZERO                 ***'
            WRITE(6,*)'*** (THIS MESSAGE IS ONLY PRINTED ONCE.)     ***'
            WRITE(6,*)'************************************************'
         ENDIF
CCC
         QN(JZCP1+1)  = PPP
         QN(JZCP1+2)  = THETA
         QN(JZCP1+3)  = PHI
         QN(JZCP1+4)  = CPTOT
         IQN(JZCP1+5) = NPHTOT(ILAY)
         NJK = NLAY*NBBSQ
         DO 52 ILAY = 1,NLAY
            ILAYM1 = ILAY - 1
            DO 53 I = 1,NBB
               DO 54 J = 1,NBB
                  IJK = ILAYM1*NBBSQ + (I-1)*NBB + J
                  QN(JZCP1 + 5 + IJK        ) = CERPAT(ILAY,I,J)
                  QN(JZCP1 + 5 + IJK +   NJK) =   TPAT(ILAY,I,J)
                  QN(JZCP1 + 5 + IJK + 2*NJK) =  T2PAT(ILAY,I,J)
                  QN(JZCP1 + 5 + IJK + 3*NJK) =  TFPAT(ILAY,I,J)
                  QN(JZCP1 + 5 + IJK + 4*NJK) =  TLPAT(ILAY,I,J)
 54            CONTINUE
 53         CONTINUE
 52      CONTINUE
      ENDIF
C    Remark: The calculation of # photo-electrons for KCERAL = 2
C    (using POISSN) has been moved to NUCER2.
C
C***
C
      IF (KCERAL .NE. 0) THEN
CCC         WRITE(6,*) 'NNPM  = ',NNPM
CCC         WRITE(6,*) 'PMINT = ',PMINT
         CALL MZBOOK(IXSTOR,JZCPA,JZCPA,1,'CPA ',1,1,4*NNPM+1,IOCPA,0)
         IQN(JZCPA+1) = NNPM
         DO 317 I = 1,NNPM
            QN(JZCPA+1+I)         = PMINT(I)
            QN(JZCPA+1+NNPM+I)    = TAVE(I)
            QN(JZCPA+1+2*NNPM+I)  = TFIRST(I)
            QN(JZCPA+1+3*NNPM+I)  = FLOAT(NPMINT(I))
 317     CONTINUE
         IF (LOGLDZ .GE. 2) THEN
CCC            CALL DZSHOW('CPA AFTER BOOK:',IXSTOR,JZCPA,' ',0,0,1,82)
            CALL DZSHOW('CPA AFTER BOOK:',IXSTOR,JZCPA,' ',0,0,0,0)
         ENDIF
      ENDIF
C
      IF (KCERAL .EQ. 2) THEN
         CALL MZBOOK(IXSTOR,JZPMH,JZPMH,1,'PMH ',1,1,8*NHITOT+1,IOPMH,0)
         IQN(JZPMH+1) = NHITOT
         IJ = 1
         DO 318 I = 1,NHITOT
            QN(JZPMH+IJ+1) = FLOAT(IPMHIT(I))
            QN(JZPMH+IJ+2) = TIMHIT(I)
            QN(JZPMH+IJ+3) = CPHINT(I)
            QN(JZPMH+IJ+4) = NPEINT(I)
            QN(JZPMH+IJ+5) = PHPPMH(I)
            QN(JZPMH+IJ+6) = ANGPMH(I)
            QN(JZPMH+IJ+7) = BETPMH(I)
            QN(JZPMH+IJ+8) = FLOAT(IPAHIT(I))
            IJ = IJ + 8
 318     CONTINUE
         IF (LOGLDZ .GE. 2) THEN
            CALL DZSHOW('PMH AFTER BOOK:',IXSTOR,JZPMH,' ',0,0,0,0)
         ENDIF
      ENDIF
C
C ---------------------------------------------------------------------
C
C Enter here if we're not running GEANT:
C
 1000 CONTINUE
C
C ---------------------------------------------------------------------
C
C    Check if this event generated a trigger:
      IF (KCERAL .NE. 0) THEN
C       Special trigger routine for the effective-area studies:
            CALL MYTRIG
      ENDIF
C
C    Then fill the NTUPLE:
            IF( (KCERAL .GT. 0).AND.(KGEN .EQ. 1))
     +THEN
         NNUV11 = NNUVA1 - 1
         NNUV12 = NNUVA1 - 2
         NNUV13 = NNUVA1 - 3
         NNUV14 = NNUVA1 - 4
         NNUV15 = NNUVA1 - 5
         IOFFS1 = NNUVA1 + NLYSTR
         IOFFS2 = NNUVA1 + 2*NLYSTR
         NUNTVE(1) = PIN
         NUNTVE(2) = THEMU0
         NUNTVE(3) = PHIMU0
         NUNTVE(4) = DCLAPP
         NUNTVE(5) = ZCLAPP
         CALL VZERO(NUNTVE(NNUV15), (6 + 3*NLYSTR))
         CALL VZERO(IHITLS(1), NLYSTR)
         CALL VZERO(IDOWLS(1), NLYSTR)
         CALL VZERO(IUPPLS(1), NLYSTR)
            SPMALL = 0.
            SPMNQU = 0.
            NNQUI1 = 0
            NPMALL = 0
            NPMNQU = 0
            NNQUI2 = 0
         DO 500 I = 1,NNPM
            IF (NPMINT(I) .GE. NPEMIN) THEN
C             Total number of hit PM tubes in this event:
               NUNTVE(NNUV12) = NUNTVE(NNUV12) + 1.
               IF (IORIPM(I) .LT. 0) THEN
C                Total # hit downward-looking PM:s in this event:
                  NUNTVE(NNUV11) = NUNTVE(NNUV11) + 1.
               ELSEIF (IORIPM(I) .GT. 0) THEN
C                Total # hit upward-looking PM:s in this event:
                  NUNTVE(NNUVA1) = NUNTVE(NNUVA1) + 1.
               ENDIF
C             Number of layers/strings in this event with at least one
C             PM-tube hit:
               IF (IHITLS(LAYSTR(I)) .LE. 0) THEN
                  NUNTVE(NNUV15) = NUNTVE(NNUV15) + 1.
                  IHITLS(LAYSTR(I)) = 1
               ENDIF
               IF (IORIPM(I) .LT. 0) THEN
C                Number of layers/strings in this event with at least
C                one downward-looking PM-tube hit:
                  IF (IDOWLS(LAYSTR(I)) .LE. 0) THEN
                     NUNTVE(NNUV14) = NUNTVE(NNUV14) + 1.
                     IDOWLS(LAYSTR(I)) = 1
                  ENDIF
               ELSEIF (IORIPM(I) .GT. 0) THEN
C                Number of layers/strings in this event with at least
C                one upward-looking PM-tube hit:
                  IF (IUPPLS(LAYSTR(I)) .LE. 0) THEN
                     NUNTVE(NNUV13) = NUNTVE(NNUV13) + 1.
                     IUPPLS(LAYSTR(I)) = 1
                  ENDIF
               ENDIF
C             Total number of hit PM tubes in this layer/string in
C             this event:
               IND = NNUVA1 + LAYSTR(I)
               NUNTVE(IND) = NUNTVE(IND) + 1.
               IF (IORIPM(I) .LT. 0) THEN
C                Total number of hit downward-looking PM tubes in this
C                layer/string in this event:
                  INDDOW = IOFFS1 + LAYSTR(I)
                  NUNTVE(INDDOW) = NUNTVE(INDDOW) + 1.
               ELSEIF (IORIPM(I) .GT. 0) THEN
C                Total number of hit upward-looking PM tubes in this
C                layer/string in this event:
                  INDUPP = IOFFS2 + LAYSTR(I)
                  NUNTVE(INDUPP) = NUNTVE(INDUPP) + 1.
               ENDIF
            ENDIF
  501       CONTINUE
            CALL HFILL(4001,PMINT(I),0.,1.)
            SPMALL = SPMALL + PMINT(I)
            RPMINT = FLOAT(NPMINT(I))
            CALL HFILL(4005,RPMINT,0.,1.)
            NPMALL = NPMALL + NPMINT(I)
            IF (NPMINT(I) .GT. 0) THEN
              NNQUI2 = NNQUI2 + 1
              NPMNQU = NPMNQU + NPMINT(I)
              CALL HFILL(4006,RPMINT,0.,1.)
            ENDIF
            IF (PMINT(I) .GT. 0.) THEN
              NNQUI1 = NNQUI1 + 1
              SPMNQU = SPMNQU + PMINT(I)
              RPHPMI = RPMINT/PMINT(I)
              CALL HFILL(4009,RPHPMI,0.,1.)
              IF (NPMINT(I) .GT. 0) CALL HFILL(4010,RPHPMI,0.,1.)
            ENDIF
  500    CONTINUE
         CALL HFN(9001, NUNTVE)
            AVPIAL = SPMALL/FLOAT(NNPM)
C       write(6,*)'spmall,float(nnpm)', spmall, float(nnpm)
            AVPEAL = FLOAT(NPMALL)/FLOAT(NNPM)
            CALL HFILL(4003,AVPIAL,0.,1.)
            CALL HFILL(4007,AVPEAL,0.,1.)
            IF (NNQUI1 .GT. 0) THEN
              AVPINQ = SPMNQU/FLOAT(NNQUI1)
              CALL HFILL(4004,AVPINQ,0.,1.)
            ENDIF
            IF (NNQUI2 .GT. 0) THEN
              AVPENQ = FLOAT(NPMNQU)/FLOAT(NNQUI2)
              CALL HFILL(4008,AVPENQ,0.,1.)
            ELSE
              CALL HFILL(4008,0.,  0.,1.)
            ENDIF
C+ADD,P=NUGU,D=GUOUT,C=380.
C        WRITE(6,*)' Particles generated:'
C        write(6,*) ' garmmas          = ', garms
C        write(6,*) ' e+ and e-        = ', eles
C        write(6,*) ' neutrinos        = ', neus
C        write(6,*) ' mu+ and mu-      = ', mus
C        write(6,*) ' pi0              = ', pio
C        write(6,*) ' pi+ and pi-      = ', pis
C        write(6,*) ' n and p and pbar = ', npbs
C        write(6,*) ' others           = ', others
CCCC      CALL NUSHOU
C        WRITE(6,*)'NPMALL = ',NPMALL
C        WRITE(6,*)'--------------------NPMINT: -------------'
C        WRITE(6,733)NPMINT
C 733    FORMAT(1X,10I7)
C        WRITE(6,*)'-----------------------------------------'
      ENDIF
C
C--> ZEBRA output:
C
      CALL NUZOUT
C
C
C--> Save random-number seeds on file after this event, if requested:
C
      IF (IRNSAV .EQ. 0) GOTO 999
C
      IF (ABS(IRNSAV) .GE. 2) THEN
C       Note: the following test should be done on the number of events
C       we have actually simulated (IEVACT); this might differ from the
C       event number (IEVENT) in case we did restart the simulation from
C       event # IEVTRS > 1.
         IF (MOD(IEVACT,ABS(IRNSAV)) .NE. 0) GOTO 999
      ENDIF
C
      IF (IRNSAV .LT. 0) THEN
         REWIND(34,ERR=901)
         GOTO 902
  901    CONTINUE
         WRITE(6,*) ' '
         WRITE(6,*) '*** GUOUT:'
         WRITE(6,*)
     +   '*** ERROR REWINDING RANDOM-NUMBER OUTPUT FILE (UNIT 34)'
         WRITE(6,*)
     +   '*** AT THE END OF EVENT NUMBER ',IEVENT
         GOTO 999
  902    CONTINUE
      ENDIF
      CALL RNWRIT(34,1,IEVENT+1)
C
  999 CONTINUE
C
      RETURN
      END
CDECK  ID>, GUTRAK. 
      SUBROUTINE GUTRAK
      COMMON/NUMISC/START,  IPTYPE, PIN, EMU0, THEMU0, PHIMU0,
     +              ZTRMAX, IORIVX, IORITR, VOLD(7),
     +              KTMED,  REFIND, ATTLEN, CWLGT1, CWLGT2,
     +              CUTKC2, IPRHIS(3), TDIAM, INUTRA,
     +              ZCLAPP, DCLAPP
      INTEGER START
      COMMON/NUTRAK/PVEC(3),AAMASS,XM,YM,ZM,TLENG,T0,ISTEP,BBETA,PTOT
      COMMON/NUFLAG/FGEANT,KGEN,KFLUXE,KFLUXS,FWEI,
     *       IRNALG,IRNSAV,IEVTRS,IEVACT,NEVACT,
     *       KCERAL,FOPTIN,FDPST,FRSCAT,IVTXAL,FTISME,
     *       ZINVBS
      LOGICAL FGEANT,FWEI,FOPTIN,FDPST,FRSCAT,FTISME,ZINVBS
      COMMON/GCKINE/IKINE,PKINE(10),ITRA,ISTAK,IVERT,IPART,ITRTYP
     +      ,NAPART(5),AMASS,CHARGE,TLIFE,VERT(3),PVERT(4),IPAOLD
C
      PARAMETER ( MAXPM  = 250 )
      PARAMETER ( MXPMHI = 100*MAXPM )
      COMMON/NUCPM/PMINT(MAXPM),TFIRST(MAXPM),TAVE(MAXPM),NPMINT(MAXPM),
     +             NHITOT, IPMHIT(MXPMHI),
     +             TIMHIT(MXPMHI), CPHINT(MXPMHI), NPEINT(MXPMHI),
     +             PHPPMH(MXPMHI), ANGPMH(MXPMHI), BETPMH(MXPMHI),
     +             IPAHIT(MXPMHI)
      COMMON/NUPMPO/NNPM,X0PM(MAXPM),Y0PM(MAXPM),Z0PM(MAXPM),
     +              IORIPM(MAXPM), LAYSTR(MAXPM)
C     IORIPM specifies the orientation of the PM tubes:
C     IORIPM(i) = -1: PM tube number i is looking downward
C     IORIPM(i) = +1: PM tube number i is looking upward
C     LAYSTR(i) = the number of the layer or string that PM tube
C                 number i belongs to
      PARAMETER ( MAXQQ  = 5000.)
      PARAMETER ( MAXTIM = 100.)
      PARAMETER ( MAXPS  = 5000)
      COMMON /MYSHOC/ XPRTRA(MAXPS), YPRTRA(MAXPS), ZPRTRA(MAXPS),
     +                DXPTRA(MAXPS), DYPTRA(MAXPS), DZPTRA(MAXPS),
     +                PEMAS1,PTOT1(MAXPS),PTOT2,IPS, IPSTRA,
     +                ZSHMAX,ZDIST(MAXPS),ETOT1(MAXPS),TOTPHO,
     +                DEGTHE(MAXQQ),
     +                GARMS, ELES,NEUS,MUS,NPBS,OTHERS,PIO,PIS
      COMMON/GCONST/PI,TWOPI,PIBY2,DEGRAD,RADDEG,CLIGHT,BIG,EMASS
      COMMON/GCONSX/EMMU,PMASS,AVO
C
      PARAMETER (MAXMEC=30)
      COMMON/GCTRAK/VECT(7),GETOT,GEKIN,VOUT(7),NMEC,LMEC(MAXMEC)
     + ,NAMEC(MAXMEC),NSTEP ,MAXNST,DESTEP,DESTEL,SAFETY,SLENG
     + ,STEP  ,SNEXT ,SFIELD,TOFG  ,GEKRAT,UPWGHT,IGNEXT,INWVOL
     + ,ISTOP ,IGAUTO,IEKBIN, ILOSL, IMULL,INGOTO,NLDOWN,NLEVIN
     + ,NLVSAV,ISTORY
C
      PARAMETER ( MXZBIN = 50 )
      PARAMETER ( MXTBIN = 40 )
      PARAMETER ( MXSHOW = 5000 )
      COMMON/NUSHLC/SHPHOT(0:MXZBIN,0:MXTBIN)
     *             ,SHT(0:MXZBIN,0:MXTBIN),SHT2(0:MXZBIN,0:MXTBIN)
     *             ,SHLATE(0:MXZBIN,0:MXTBIN)
     *             ,NBZ,NBTH,ZBINN(0:MXZBIN),THBINN(0:MXTBIN)
     *             ,EESHOW(MXSHOW),IDXSHO(MXSHOW)
     *             ,USESHL,NWBIN,INSLIB,NSH,CGEN,NZBMAX,NTBMAX
     *             ,QQCER2,QQSLIB,NPCER2,NPSLIB
     *             ,SHOFIL,CHDIR,CHNPAH,CHCPAH,CHNPAZ,CHCPAZ
      CHARACTER*48 SHOFIL
      CHARACTER*5 CHDIR
      CHARACTER*20 CHNPAH,CHCPAH,CHNPAZ,CHCPAZ
      LOGICAL USESHL,INSLIB,CGEN
      COMMON/GCTATI/NSEC,KIND(30),EN(30),PL(30),PT(30),THETA(30),PHI(30)
      INTEGER NSEC,KIND
      REAL    EN,PL,PT,THETA,PHI
C
C
      IF (KCERAL .EQ. 2) THEN
         ISTEP = 0
         TLENG = 0.0
      ENDIF
C
      START = 2
C
C    INUTRA is our own track counter (GEANT's ITRA doesn't properly
C    take into account the tracks in the temporary stack):
      INUTRA = INUTRA + 1
C
CC        WRITE(6,*)'New track!'
CC        CALL GPCXYZ
            INSLIB = .FALSE.
            IF (INUTRA .GT. 1) THEN
C          IPSTRA = -1
              DO 111 I = 1,IPS
                DISTSQ = (VECT(1)-XPRTRA(I))*(VECT(1)-XPRTRA(I)) +
     +               (VECT(2)-YPRTRA(I))*(VECT(2)-YPRTRA(I)) +
     +               (VECT(3)-ZPRTRA(I))*(VECT(3)-ZPRTRA(I))
                IF (DISTSQ .LT. 1.) THEN
                  IPSTRA = I
                  GOTO 112
                ENDIF
  111         CONTINUE
  112         CONTINUE
C          IF (IPSTRA .EQ. -1) THEN
C            WRITE(6,*) '*** GUTRAK: IPSTRA = -1'
C            IPSTRA = 1
C         ELSE
C            WRITE(6,*) 'IPSTRA,DXPTRA = ',IPSTRA,DXPTRA(IPSTRA)
C          ENDIF
            ENDIF
C          CALL HFILL(2003,FLOAT(IPART),0.,1.)
            PTRACK = PVERT(1)*PVERT(1) + PVERT(2)*PVERT(2) + PVERT(3)
     &        *PVERT(3)
            PTRACK = SQRT(PTRACK)
            CALL HFILL(2004,FLOAT(IPART),PTRACK,1.)
C          CALL HFILL(2005,FLOAT(IPART),PTRACK,1.)
            IF (IPART .EQ. 1) THEN
C       Gammas
              IDFILL = 2006
              IFILL = 3006
              GARMS = GARMS + 1
            ELSEIF ((IPART .EQ. 2) .OR. (IPART .EQ. 3)) THEN
C       e+, e-
              IDFILL = 2007
              IFILL = 3007
              ELES   = ELES + 1
            ELSEIF (IPART .EQ. 4) THEN
C       Neutrinos
              IDFILL = 2008
              IFILL = 3008
              NEUS   = NEUS + 1
            ELSEIF ((IPART .EQ. 5) .OR. (IPART .EQ. 6)) THEN
C       mu+, mu-
              IDFILL = 2009
              IFILL = 3009
              MUS    = MUS + 1
            ELSEIF (IPART .EQ. 7) THEN
C       Pi0
              IDFILL = 2010
              IFILL = 3010
              PIO    = PIO + 1
            ELSEIF ((IPART .EQ. 8) .OR. (IPART .EQ. 9)) THEN
C       pi+, pi-
              IDFILL = 2011
              IFILL = 3011
              PIS    = PIS + 1
            ELSEIF ((IPART .EQ. 13) .OR. (IPART .EQ. 14) .OR. (IPART .
     &        EQ.
     &        15)) THEN
C       n, p, pbar
              IDFILL = 2012
              IFILL = 3012
              NPBS   = NPBS + 1
            ELSE
C       Others
              IDFILL = 2013
              IFILL = 3013
              OTHERS = OTHERS + 1
            ENDIF
C        CALL HFILL(IFILL,SLENG,0.,1.)
            CALL HFILL(IDFILL,PTRACK,0.,1.)
C ******************* ANGLE'S DISTRIBUTION , FIRST AND SECONDARY ************
C
            IF ( INUTRA .EQ. 1) THEN
C            WRITE (6,*) ' INUTRA .EQ.1 !'
            ELSE IF ((IPART .EQ. 2) .OR. (IPART .EQ. 3)) THEN
C Three components of momentum of the primary track:
              PX  = DXPTRA(IPSTRA)*PTOT1(IPSTRA)
              PY  = DYPTRA(IPSTRA)*PTOT1(IPSTRA)
              PZ  = DZPTRA(IPSTRA)*PTOT1(IPSTRA)
C Three components of momentum of the current track:
              CX  = PVERT(1)
              CY  = PVERT(2)
              CZ  = PVERT(3)

              CPX = PX*CX
              CPY = PY*CY
              CPZ = PZ*CZ

              SPP = SQRT(PX*PX+PY*PY+PZ*PZ)
              SCC = SQRT(CX*CX+CY*CY+CZ*CZ)
C
C            WRITE(6,*) 'IPSTRA,PX,PY,PZ,P=',IPSTRA,PX,PY,PZ,SPP
C            WRITE(6,*) 'PVERT,CURRENT:',PVERT(1),PVERT(2),
C     +                 PVERT(3),SCC
C
              COSANG = (CPX+CPY+CPZ)/(SPP*SCC)
              IF( ABS(COSANG ).GE. 1.) THEN
                COSANG = SIGN(1.,COSANG)
              ENDIF
              SHWANG = ACOS(COSANG)
C
C Transforms radius to degree:
C
              DEGANG = SHWANG*RADDEG
C
              CALL HFILL (4441, DEGANG,0.,1.)
C
            ENDIF
C *****************************************************************************
C
C
      CALL GTRACK
C
      IF (KCERAL .EQ. 2) THEN
         IF (TLENG .GT. 0.1) THEN
            PVEC(1) = PVEC(1)*PTOT
            PVEC(2) = PVEC(2)*PTOT
            PVEC(3) = PVEC(3)*PTOT
            IF ( USESHL .AND. INSLIB ) THEN
              CONTINUE
CCC          CGEN = .TRUE.
C          IF (INUTRA .EQ. 1) THEN
C            CALL NUCER2(PVEC,AAMASS,XM,YM,ZM,TLENG,T0,IOK)
C          ELSE
C            IF ( INSLIB ) THEN
C              CALL GPCXYZ
C              WRITE(6,*)'ERROR! Twice in NUSLIB.'
C            ENDIF
C            CALL NUSLIB(PVEC,AAMASS,XM,YM,ZM,TLENG,T0,IOK)
CCC          ENDIF
C940921            ELSEIF (INUTRA .EQ. 1) THEN
            ELSE
              CALL NUCER2(PVEC,AAMASS,XM,YM,ZM,TLENG,T0,IOK)
            ENDIF
         ENDIF
      ENDIF
C
      RETURN
      END
CDECK  ID>, GUSTEP. 
      SUBROUTINE GUSTEP
C
      PARAMETER (MAXMEC=30)
      COMMON/GCTRAK/VECT(7),GETOT,GEKIN,VOUT(7),NMEC,LMEC(MAXMEC)
     + ,NAMEC(MAXMEC),NSTEP ,MAXNST,DESTEP,DESTEL,SAFETY,SLENG
     + ,STEP  ,SNEXT ,SFIELD,TOFG  ,GEKRAT,UPWGHT,IGNEXT,INWVOL
     + ,ISTOP ,IGAUTO,IEKBIN, ILOSL, IMULL,INGOTO,NLDOWN,NLEVIN
     + ,NLVSAV,ISTORY
C
      COMMON/GCKINE/IKINE,PKINE(10),ITRA,ISTAK,IVERT,IPART,ITRTYP
     +      ,NAPART(5),AMASS,CHARGE,TLIFE,VERT(3),PVERT(4),IPAOLD
C
      INTEGER MXGKIN
      PARAMETER (MXGKIN=100)
      COMMON/GCKING/KCASE,NGKINE,GKIN(5,MXGKIN),
     +                           TOFD(MXGKIN),IFLGK(MXGKIN)
      INTEGER       KCASE,NGKINE ,IFLGK
      REAL          GKIN,TOFD
C
      COMMON/NUMULT/ NGAM,NELE,NCHA,NNEUT,NMUON,MGTYPE(100)
     *              ,MGAM,MELE,MCHA,MNEUT,MMUON,FCONV,FHADR
      LOGICAL FCONV,FHADR
      COMMON/NUMISC/START,  IPTYPE, PIN, EMU0, THEMU0, PHIMU0,
     +              ZTRMAX, IORIVX, IORITR, VOLD(7),
     +              KTMED,  REFIND, ATTLEN, CWLGT1, CWLGT2,
     +              CUTKC2, IPRHIS(3), TDIAM, INUTRA,
     +              ZCLAPP, DCLAPP
      INTEGER START
      PARAMETER ( MAXZH = 500 )
      COMMON/NUZH/NZH,JZH,DISTZH,ZH(MAXZH),EZH(MAXZH),THEZH(MAXZH),
     +            CPATH(MAXZH),CPTOT,HISTZH,IHZH,FTHEZH,ETHLIM
      LOGICAL HISTZH,FTHEZH
      COMMON/NUFLAG/FGEANT,KGEN,KFLUXE,KFLUXS,FWEI,
     *       IRNALG,IRNSAV,IEVTRS,IEVACT,NEVACT,
     *       KCERAL,FOPTIN,FDPST,FRSCAT,IVTXAL,FTISME,
     *       ZINVBS
      LOGICAL FGEANT,FWEI,FOPTIN,FDPST,FRSCAT,FTISME,ZINVBS
      COMMON/NUTRAK/PVEC(3),AAMASS,XM,YM,ZM,TLENG,T0,ISTEP,BBETA,PTOT
      COMMON/GCONST/PI,TWOPI,PIBY2,DEGRAD,RADDEG,CLIGHT,BIG,EMASS
      COMMON/GCONSX/EMMU,PMASS,AVO
C
      COMMON /MYCOUN/ NTRGEV(5), NQ(100)
      PARAMETER ( MAXQQ  = 5000.)
      PARAMETER ( MAXTIM = 100.)
      PARAMETER ( MAXPS  = 5000)
      COMMON /MYSHOC/ XPRTRA(MAXPS), YPRTRA(MAXPS), ZPRTRA(MAXPS),
     +                DXPTRA(MAXPS), DYPTRA(MAXPS), DZPTRA(MAXPS),
     +                PEMAS1,PTOT1(MAXPS),PTOT2,IPS, IPSTRA,
     +                ZSHMAX,ZDIST(MAXPS),ETOT1(MAXPS),TOTPHO,
     +                DEGTHE(MAXQQ),
     +                GARMS, ELES,NEUS,MUS,NPBS,OTHERS,PIO,PIS
      PARAMETER ( MXZBIN = 50 )
      PARAMETER ( MXTBIN = 40 )
      PARAMETER ( MXSHOW = 5000 )
      COMMON/NUSHLC/SHPHOT(0:MXZBIN,0:MXTBIN)
     *             ,SHT(0:MXZBIN,0:MXTBIN),SHT2(0:MXZBIN,0:MXTBIN)
     *             ,SHLATE(0:MXZBIN,0:MXTBIN)
     *             ,NBZ,NBTH,ZBINN(0:MXZBIN),THBINN(0:MXTBIN)
     *             ,EESHOW(MXSHOW),IDXSHO(MXSHOW)
     *             ,USESHL,NWBIN,INSLIB,NSH,CGEN,NZBMAX,NTBMAX
     *             ,QQCER2,QQSLIB,NPCER2,NPSLIB
     *             ,SHOFIL,CHDIR,CHNPAH,CHCPAH,CHNPAZ,CHCPAZ
      CHARACTER*48 SHOFIL
      CHARACTER*5 CHDIR
      CHARACTER*20 CHNPAH,CHCPAH,CHNPAZ,CHCPAZ
      LOGICAL USESHL,INSLIB,CGEN
      COMMON/GCTATI/NSEC,KIND(30),EN(30),PL(30),PT(30),THETA(30),PHI(30)
      INTEGER NSEC,KIND
      REAL    EN,PL,PT,THETA,PHI
C
C
C
        LOGICAL STOSEC
C
        REAL BETACR
        INTEGER ICALL
        SAVE BETACR, ICALL
        DATA ICALL/0/
C
        IF (ICALL .EQ. 0) THEN
          BETACR = 1.0/REFIND
          ICALL = 1
        ENDIF
C
C     In this simulation, we are not interested in particles with
C     energy below the Cherenkov threshold (only consider charged
C     particles at this point):
        IF (CHARGE .NE. 0.) THEN
          BETA = VECT(7)/GETOT
          IF (BETA .LE. BETACR) THEN
            ISTOP = 1
            RETURN
          ENDIF
        ENDIF
C
C   Here we calculate the total length and the time that the SHOWER takes:

        IF (CHARGE .NE. 0.) THEN
          IF ( VECT(3)- ZPRTRA(1) .GT. ZSHMAX ) THEN
            ZSHMAX = VECT(3) - ZPRTRA(1)
          ENDIF
          IF ( VECT(3) .GT. ZSTMAX ) THEN
            ZSTMAX = VECT(3)
          ENDIF
CC        CALL GPCXYZ
        ENDIF
      IF ((CHARGE .NE. 0.) .AND. (KCERAL .GT. 0)) THEN
         BETA = VECT(7)/GETOT
         IF (BETA .GT. 1.0/REFIND) THEN
            IF (KCERAL .EQ. 1) THEN
               CALL NUCER1
            ELSEIF (KCERAL .EQ. 2) THEN
               IF (ISTEP .EQ. 0) THEN
                  PVEC(1) = VECT(4)
                  PVEC(2) = VECT(5)
                  PVEC(3) = VECT(6)
                  PTOT    = VECT(7)
                  XM = VECT(1)
                  YM = VECT(2)
                  ZM = VECT(3)
                  T0 = TOFG
                  AAMASS = AMASS
                  BBETA  = BETA
                  TLENG  = STEP
                  ISTEP  = 1
        IF ( USESHL ) THEN
          CGEN = .FALSE.
          PVEC(1) = PVEC(1)*PTOT
          PVEC(2) = PVEC(2)*PTOT
          PVEC(3) = PVEC(3)*PTOT
          IF ( INSLIB ) THEN
            CALL GPCXYZ
            WRITE(6,*)'ERROR! Twice in NUSLIB.'
          ENDIF
          CALL NUSLIB(PVEC,AAMASS,XM,YM,ZM,TLENG,T0,IOK)
        ENDIF
               ELSE
                  TLENG = TLENG + STEP
                  ISTEP = ISTEP + 1
                  CGEN  = .FALSE.
                  IF (BBETA-BETA .GT. 0.05) THEN
                     CGEN = .TRUE.
                  ELSE
                     CTH = PVEC(1)*VECT(4) + PVEC(2)*VECT(5) +
     +                     PVEC(3)*VECT(6)
                     IF (CTH .LT. 0.99) THEN
                        CGEN = .TRUE.
                     ENDIF
                  ENDIF
        IF (TLENG .LT. 5.0) THEN
                     CGEN = .FALSE.
                  ENDIF
          IF ( CGEN .OR. USESHL ) THEN
                     PVEC(1) = PVEC(1)*PTOT
                     PVEC(2) = PVEC(2)*PTOT
                     PVEC(3) = PVEC(3)*PTOT
            IF ( USESHL ) THEN
              IF (INUTRA .EQ. 1) THEN
                IF ( CGEN ) THEN
                  CALL NUCER2(PVEC,AAMASS,XM,YM,ZM,TLENG,T0,IOK)
                ENDIF
              ELSE
                IF ( INSLIB ) THEN
                  CALL GPCXYZ
                  WRITE(6,*)'ERROR! Twice in NUSLIB.'
                ENDIF
                CALL NUSLIB(PVEC,AAMASS,XM,YM,ZM,TLENG,T0,IOK)
                IF ( INSLIB ) THEN
                  write(6,*)'Not supposed to come here.'
                ENDIF
              ENDIF
              PVEC(1) = PVEC(1)/PTOT
              PVEC(2) = PVEC(2)/PTOT
              PVEC(3) = PVEC(3)/PTOT
            ELSEIF ( CGEN ) THEN
              CALL NUCER2(PVEC,AAMASS,XM,YM,ZM,TLENG,T0,IOK)
            ENDIF
            IF ( CGEN ) THEN
                     TLENG = 0.0
                     PVEC(1) = VECT(4)
                     PVEC(2) = VECT(5)
                     PVEC(3) = VECT(6)
                     PTOT    = VECT(7)
                     XM = VECT(1)
                     YM = VECT(2)
                     ZM = VECT(3)
                     T0 = TOFG
                     AAMASS = AMASS
                     BBETA  = BETA
                  ENDIF
            ENDIF
               ENDIF
            ENDIF
         ENDIF
      ENDIF
C
      IF (NGKINE .GT. 0) THEN
         CALL HFILL(30,VECT(3),0.,1.)
            STOSEC = .FALSE.
C         IF ( INUTRA .NE. 1 ) GOTO 10
         DO 10 I=1,NGKINE
            IP = GKIN(5,I)
            IF (IP .EQ. 1) THEN
              CALL GSKING(I)
              STOSEC = .TRUE.
            ELSE
              PSECON = GKIN(1,I)*GKIN(1,I) + GKIN(2,I)*GKIN(2,I) +
     +               GKIN(3,I)*GKIN(3,I)
              PSECON = SQRT(PSECON)
              BETASE = PSECON/GKIN(4,I)
              IF (BETASE .GT. BETACR) THEN
                CALL GSKING(I)
                STOSEC = .TRUE.
              ENDIF
            ENDIF
 10      CONTINUE
            IF ((INUTRA .EQ. 1) .AND. STOSEC) THEN
C          WRITE(6,*) 'FIRST IPS=',IPS
              IF ( IPS .LT. MAXPS ) THEN
                IPS = IPS+1
                XPRTRA(IPS) = VECT(1)
                YPRTRA(IPS) = VECT(2)
                ZPRTRA(IPS) = VECT(3)
                DXPTRA(IPS) = VECT(4)
                DYPTRA(IPS) = VECT(5)
                DZPTRA(IPS) = VECT(6)
                PTOT1(IPS)  = VECT(7)
                ETOT1(1)  = SQRT(PEMAS1*PEMAS1+
     +                    PTOT1(1)*PTOT1(1))
                BETA1 = PTOT1(1)/ETOT1(1)
C           WRITE(6,*) 'IPS=',IPS
C           write(6,*)'dxptra,y,z',dxptra(IPS),
C     +               dyptra(IPS),dzptra(IPS)
              ELSE
C            WRITE(6,*)'IPS = ',IPS,' !!!'
              ENDIF
              PEMAS1 = AMASS
CCC            PTOT1  = PTOT
C            write(6,*)'ptot1=',ptot1(ips)
C            if ( ips .gt. 1 ) then
C              IF ( PTOT1(IPS-1)-PTOT1(IPS) .GT. 0.2*PTOT1(IPS-1)) THEN
C                CALL GPCXYZ
C              ENDIF
C            ENDIF
C       WRITE(6,*) 'store primary parameters:'
C       WRITE(6,*) 'VECT: 4,5,6,7 =', VECT(4),VECT(5),VECT(6),VECT(7)
CCC   TEMPORARY SOLUTION: STOP THE PRIMARY TRACK IF IT HAS GENERATED
CCC   ANY SECONDARIES:
C            ISTOP = 1
            ENDIF
      ENDIF
C
      RETURN
      END
CDECK  ID>, NUSDAT. 
            SUBROUTINE NUSDAT
C
      COMMON/NUFLAG/FGEANT,KGEN,KFLUXE,KFLUXS,FWEI,
     *       IRNALG,IRNSAV,IEVTRS,IEVACT,NEVACT,
     *       KCERAL,FOPTIN,FDPST,FRSCAT,IVTXAL,FTISME,
     *       ZINVBS
      LOGICAL FGEANT,FWEI,FOPTIN,FDPST,FRSCAT,FTISME,ZINVBS
      COMMON/NUZZIO/IZUNII,IZUNIO,IOZH0,IOEZH,IOTZH,IOCPA,IOCP1,
     *      IOPMH,IOSDIR,IOSHHI,IOSHLO,IOLAST
      COMMON/NUFWIO/IFWIO,IWGEAN,IWZH0,IWEZH,IWTZH,IWCPA,IWCP1,IWPMH,
     +              IWLAST
      COMMON/NUHBIO/IHBUNI
C
C-----------------------------------------------
C Memory requirements for HBOOK, Geant, and our
C private Zebra structure.
C-----------------------------------------------
C
      PARAMETER ( LHMEM  = 400000 )
      PARAMETER ( LGMEM  = 100000 )
      PARAMETER ( LZMEM  = 600000 )
      PARAMETER ( IFFSZE = 1000   )
      COMMON/PAWC/HHH(LHMEM)
      COMMON/GCBANK/Q(LGMEM)
      COMMON/NUBANK/FENCE(16),QSPACE(LZMEM)
      DIMENSION LQN(LZMEM),IQN(LZMEM-8)
      DIMENSION QN(LZMEM-8)
      EQUIVALENCE ( QSPACE(1), LQN(1) )
      EQUIVALENCE ( IQN(1), QN(1), LQN(9) )
      COMMON/CFREAD/FFSPAC(IFFSZE)
C-----------------------------------------------
      COMMON/GCNUM/NMATE ,NVOLUM,NROTM,NTMED,NTMULT,NTRACK,NPART
     +            ,NSTMAX,NVERTX,NHEAD,NBIT
      COMMON /GCNUMX/ NALIVE,NTMSTO
C
      COMMON/GCFLAG/IDEBUG,IDEMIN,IDEMAX,ITEST,IDRUN,IDEVT,IEORUN
     +        ,IEOTRI,IEVENT,ISWIT(10),IFINIT(20),NEVENT,NRNDM(2)
      COMMON/GCFLAX/BATCH, NOLOG
      LOGICAL BATCH, NOLOG
C
      COMMON/GCLIST/NHSTA,NGET ,NSAVE,NSETS,NPRIN,NGEOM,NVIEW,NPLOT
     +       ,NSTAT,LHSTA(20),LGET (20),LSAVE(20),LSETS(20),LPRIN(20)
     +             ,LGEOM(20),LVIEW(20),LPLOT(20),LSTAT(20)
C
      PARAMETER ( MAXPM  = 250 )
      PARAMETER ( MXPMHI = 100*MAXPM )
      COMMON/NUCPM/PMINT(MAXPM),TFIRST(MAXPM),TAVE(MAXPM),NPMINT(MAXPM),
     +             NHITOT, IPMHIT(MXPMHI),
     +             TIMHIT(MXPMHI), CPHINT(MXPMHI), NPEINT(MXPMHI),
     +             PHPPMH(MXPMHI), ANGPMH(MXPMHI), BETPMH(MXPMHI),
     +             IPAHIT(MXPMHI)
      COMMON/NUPMPO/NNPM,X0PM(MAXPM),Y0PM(MAXPM),Z0PM(MAXPM),
     +              IORIPM(MAXPM), LAYSTR(MAXPM)
C     IORIPM specifies the orientation of the PM tubes:
C     IORIPM(i) = -1: PM tube number i is looking downward
C     IORIPM(i) = +1: PM tube number i is looking upward
C     LAYSTR(i) = the number of the layer or string that PM tube
C                 number i belongs to
      COMMON/NUGRID/ZSIZE,HGPOS,DISTX,DISTY,DISTZ,OFFLAY,OFFLX,OFFLY,
     *              OFFX,OFFY,NBINS,NLAY,NSIDE,NSTR,NPMSTR,IDCONF,
     *              ICONTY,PMDIAM,PMRAD,PMR2,PMAREA,XYFUZZ,
     *              VXPLSI, VXPLDI
      COMMON/NUFILE/DATFIL,RNDFIN,RNDFOU,HISFIL
      CHARACTER*48  DATFIL,RNDFIN,RNDFOU,HISFIL
C
C My own COMMONs, defined above:
      COMMON /MYCOUN/ NTRGEV(5), NQ(100)
C
            RETURN
            END
CDECK  ID>, NUSINI. 
            SUBROUTINE NUSINI
C
      COMMON/NUFLAG/FGEANT,KGEN,KFLUXE,KFLUXS,FWEI,
     *       IRNALG,IRNSAV,IEVTRS,IEVACT,NEVACT,
     *       KCERAL,FOPTIN,FDPST,FRSCAT,IVTXAL,FTISME,
     *       ZINVBS
      LOGICAL FGEANT,FWEI,FOPTIN,FDPST,FRSCAT,FTISME,ZINVBS
      COMMON/NUZZIO/IZUNII,IZUNIO,IOZH0,IOEZH,IOTZH,IOCPA,IOCP1,
     *      IOPMH,IOSDIR,IOSHHI,IOSHLO,IOLAST
      COMMON/NUFWIO/IFWIO,IWGEAN,IWZH0,IWEZH,IWTZH,IWCPA,IWCP1,IWPMH,
     +              IWLAST
      COMMON/NUHBIO/IHBUNI
C
C-----------------------------------------------
C Memory requirements for HBOOK, Geant, and our
C private Zebra structure.
C-----------------------------------------------
C
      PARAMETER ( LHMEM  = 400000 )
      PARAMETER ( LGMEM  = 100000 )
      PARAMETER ( LZMEM  = 600000 )
      PARAMETER ( IFFSZE = 1000   )
      COMMON/PAWC/HHH(LHMEM)
      COMMON/GCBANK/Q(LGMEM)
      COMMON/NUBANK/FENCE(16),QSPACE(LZMEM)
      DIMENSION LQN(LZMEM),IQN(LZMEM-8)
      DIMENSION QN(LZMEM-8)
      EQUIVALENCE ( QSPACE(1), LQN(1) )
      EQUIVALENCE ( IQN(1), QN(1), LQN(9) )
      COMMON/CFREAD/FFSPAC(IFFSZE)
C-----------------------------------------------
      COMMON/GCNUM/NMATE ,NVOLUM,NROTM,NTMED,NTMULT,NTRACK,NPART
     +            ,NSTMAX,NVERTX,NHEAD,NBIT
      COMMON /GCNUMX/ NALIVE,NTMSTO
C
      COMMON/GCFLAG/IDEBUG,IDEMIN,IDEMAX,ITEST,IDRUN,IDEVT,IEORUN
     +        ,IEOTRI,IEVENT,ISWIT(10),IFINIT(20),NEVENT,NRNDM(2)
      COMMON/GCFLAX/BATCH, NOLOG
      LOGICAL BATCH, NOLOG
C
      COMMON/GCLIST/NHSTA,NGET ,NSAVE,NSETS,NPRIN,NGEOM,NVIEW,NPLOT
     +       ,NSTAT,LHSTA(20),LGET (20),LSAVE(20),LSETS(20),LPRIN(20)
     +             ,LGEOM(20),LVIEW(20),LPLOT(20),LSTAT(20)
C
      PARAMETER ( MAXPM  = 250 )
      PARAMETER ( MXPMHI = 100*MAXPM )
      COMMON/NUCPM/PMINT(MAXPM),TFIRST(MAXPM),TAVE(MAXPM),NPMINT(MAXPM),
     +             NHITOT, IPMHIT(MXPMHI),
     +             TIMHIT(MXPMHI), CPHINT(MXPMHI), NPEINT(MXPMHI),
     +             PHPPMH(MXPMHI), ANGPMH(MXPMHI), BETPMH(MXPMHI),
     +             IPAHIT(MXPMHI)
      COMMON/NUPMPO/NNPM,X0PM(MAXPM),Y0PM(MAXPM),Z0PM(MAXPM),
     +              IORIPM(MAXPM), LAYSTR(MAXPM)
C     IORIPM specifies the orientation of the PM tubes:
C     IORIPM(i) = -1: PM tube number i is looking downward
C     IORIPM(i) = +1: PM tube number i is looking upward
C     LAYSTR(i) = the number of the layer or string that PM tube
C                 number i belongs to
      COMMON/NUGRID/ZSIZE,HGPOS,DISTX,DISTY,DISTZ,OFFLAY,OFFLX,OFFLY,
     *              OFFX,OFFY,NBINS,NLAY,NSIDE,NSTR,NPMSTR,IDCONF,
     *              ICONTY,PMDIAM,PMRAD,PMR2,PMAREA,XYFUZZ,
     *              VXPLSI, VXPLDI
      COMMON/NUFILE/DATFIL,RNDFIN,RNDFOU,HISFIL
      CHARACTER*48  DATFIL,RNDFIN,RNDFOU,HISFIL
C
C My own COMMONs, defined above:
      COMMON /MYCOUN/ NTRGEV(5), NQ(100)
C
C--> My counters of # triggered events
C
            CALL VZERO(NTRGEV,5)
C
            RETURN
            END
CDECK  ID>, NUSEND. 
            SUBROUTINE NUSEND
C
      COMMON/NUFLAG/FGEANT,KGEN,KFLUXE,KFLUXS,FWEI,
     *       IRNALG,IRNSAV,IEVTRS,IEVACT,NEVACT,
     *       KCERAL,FOPTIN,FDPST,FRSCAT,IVTXAL,FTISME,
     *       ZINVBS
      LOGICAL FGEANT,FWEI,FOPTIN,FDPST,FRSCAT,FTISME,ZINVBS
      COMMON/NUZZIO/IZUNII,IZUNIO,IOZH0,IOEZH,IOTZH,IOCPA,IOCP1,
     *      IOPMH,IOSDIR,IOSHHI,IOSHLO,IOLAST
      COMMON/NUFWIO/IFWIO,IWGEAN,IWZH0,IWEZH,IWTZH,IWCPA,IWCP1,IWPMH,
     +              IWLAST
      COMMON/NUHBIO/IHBUNI
C
C-----------------------------------------------
C Memory requirements for HBOOK, Geant, and our
C private Zebra structure.
C-----------------------------------------------
C
      PARAMETER ( LHMEM  = 400000 )
      PARAMETER ( LGMEM  = 100000 )
      PARAMETER ( LZMEM  = 600000 )
      PARAMETER ( IFFSZE = 1000   )
      COMMON/PAWC/HHH(LHMEM)
      COMMON/GCBANK/Q(LGMEM)
      COMMON/NUBANK/FENCE(16),QSPACE(LZMEM)
      DIMENSION LQN(LZMEM),IQN(LZMEM-8)
      DIMENSION QN(LZMEM-8)
      EQUIVALENCE ( QSPACE(1), LQN(1) )
      EQUIVALENCE ( IQN(1), QN(1), LQN(9) )
      COMMON/CFREAD/FFSPAC(IFFSZE)
C-----------------------------------------------
      COMMON/GCNUM/NMATE ,NVOLUM,NROTM,NTMED,NTMULT,NTRACK,NPART
     +            ,NSTMAX,NVERTX,NHEAD,NBIT
      COMMON /GCNUMX/ NALIVE,NTMSTO
C
      COMMON/GCFLAG/IDEBUG,IDEMIN,IDEMAX,ITEST,IDRUN,IDEVT,IEORUN
     +        ,IEOTRI,IEVENT,ISWIT(10),IFINIT(20),NEVENT,NRNDM(2)
      COMMON/GCFLAX/BATCH, NOLOG
      LOGICAL BATCH, NOLOG
C
      COMMON/GCLIST/NHSTA,NGET ,NSAVE,NSETS,NPRIN,NGEOM,NVIEW,NPLOT
     +       ,NSTAT,LHSTA(20),LGET (20),LSAVE(20),LSETS(20),LPRIN(20)
     +             ,LGEOM(20),LVIEW(20),LPLOT(20),LSTAT(20)
C
      PARAMETER ( MAXPM  = 250 )
      PARAMETER ( MXPMHI = 100*MAXPM )
      COMMON/NUCPM/PMINT(MAXPM),TFIRST(MAXPM),TAVE(MAXPM),NPMINT(MAXPM),
     +             NHITOT, IPMHIT(MXPMHI),
     +             TIMHIT(MXPMHI), CPHINT(MXPMHI), NPEINT(MXPMHI),
     +             PHPPMH(MXPMHI), ANGPMH(MXPMHI), BETPMH(MXPMHI),
     +             IPAHIT(MXPMHI)
      COMMON/NUPMPO/NNPM,X0PM(MAXPM),Y0PM(MAXPM),Z0PM(MAXPM),
     +              IORIPM(MAXPM), LAYSTR(MAXPM)
C     IORIPM specifies the orientation of the PM tubes:
C     IORIPM(i) = -1: PM tube number i is looking downward
C     IORIPM(i) = +1: PM tube number i is looking upward
C     LAYSTR(i) = the number of the layer or string that PM tube
C                 number i belongs to
      COMMON/NUGRID/ZSIZE,HGPOS,DISTX,DISTY,DISTZ,OFFLAY,OFFLX,OFFLY,
     *              OFFX,OFFY,NBINS,NLAY,NSIDE,NSTR,NPMSTR,IDCONF,
     *              ICONTY,PMDIAM,PMRAD,PMR2,PMAREA,XYFUZZ,
     *              VXPLSI, VXPLDI
      COMMON/NUFILE/DATFIL,RNDFIN,RNDFOU,HISFIL
      CHARACTER*48  DATFIL,RNDFIN,RNDFOU,HISFIL
C
C My own COMMONs, defined above:
      COMMON /MYCOUN/ NTRGEV(5), NQ(100)
C
C For the effective-area studies:
            REAL EFFARE(5), EFAERR(5)
C
C--> Note: the number of actually simulated events (IEVACT) may differ
C    from the last event number (IEVENT) in case we did restart the
C    simulation from event number IEVTRS > 1.
C
C--> My counter of # triggered events
C
            WRITE(6,*) ' '
            WRITE(6,*) 'Number of triggered events = ',NTRGEV(1)
            WRITE(6,*) ' '
C
C ---- Number of triggered events, and (possibly) effective areas: ----
C
            IF (IVTXAL .NE. 3) GOTO 200
C
            VXPLSM = VXPLSI/100.
            VXPLAR = VXPLSM*VXPLSM
            DO 100 I = 1, 5
              REVENT    = FLOAT(IEVACT)
              EFFARE(I) = FLOAT(NTRGEV(I))/REVENT
              EFFARE(I) = EFFARE(I)*VXPLAR
              SQROOT    = SQRT(FLOAT(NTRGEV(I)))
              EFAERR(I) = SQROOT/REVENT
              EFAERR(I) = EFAERR(I)*VXPLAR
  100       CONTINUE
C
  200       CONTINUE
C
            IF (IVTXAL .EQ. 3) THEN
              WRITE(6,101)
            ELSE
              WRITE(6,102)
            ENDIF
  101       FORMAT(
     (/' ',4X,81('=')/
     / ' ',4X,'Number of triggers and effective area (in m**2) ',
     + 'for different trigger conditions:'/
     + ' ',4X,81('=')//
     + ' ',7X,'Trigger:',5X,'6-3',12X,'7-3',12X,'8-3',12X,'9-3',
     + 11X,'10-3'/
     / ' ',7X,80('-'))
  102       FORMAT(
     (/' ',4X,52('=')/
     / ' ',4X,'Number of triggers for different trigger conditions:'/
     + ' ',4X,52('=')//
     + ' ',7X,'Trigger:',5X,'6-3',12X,'7-3',12X,'8-3',12X,'9-3',
     + 11X,'10-3'/
     / ' ',7X,80('-'))
C
            IF (IVTXAL .EQ. 3) THEN
              WRITE(6,103) (NTRGEV(I), I=1,5), (EFFARE(I),EFAERR(I),
     &          I=1,5)
            ELSE
              WRITE(6,104) (NTRGEV(I), I=1,5)
            ENDIF
  103       FORMAT(' ',1X,' # trg.:',5(9X,I6)/
     /        ' ',2X,85('-')/
     /        ' ',1X,' eff. area:',5(F8.0,'+-',F5.0)/
     /        ' ',2X,85('-')//)
  104       FORMAT(' ',1X,' # trg.:',5(9X,I6)/' ',2X,85('-')//)
C
            RETURN
            END
CDECK  ID>, NUUHB.  
            SUBROUTINE NUUHB
      COMMON/NUMISC/START,  IPTYPE, PIN, EMU0, THEMU0, PHIMU0,
     +              ZTRMAX, IORIVX, IORITR, VOLD(7),
     +              KTMED,  REFIND, ATTLEN, CWLGT1, CWLGT2,
     +              CUTKC2, IPRHIS(3), TDIAM, INUTRA,
     +              ZCLAPP, DCLAPP
      INTEGER START
      COMMON/NUFLAG/FGEANT,KGEN,KFLUXE,KFLUXS,FWEI,
     *       IRNALG,IRNSAV,IEVTRS,IEVACT,NEVACT,
     *       KCERAL,FOPTIN,FDPST,FRSCAT,IVTXAL,FTISME,
     *       ZINVBS
      LOGICAL FGEANT,FWEI,FOPTIN,FDPST,FRSCAT,FTISME,ZINVBS
      PARAMETER ( MAXZH = 500 )
      COMMON/NUZH/NZH,JZH,DISTZH,ZH(MAXZH),EZH(MAXZH),THEZH(MAXZH),
     +            CPATH(MAXZH),CPTOT,HISTZH,IHZH,FTHEZH,ETHLIM
      LOGICAL HISTZH,FTHEZH
      COMMON/NUGENC/E0MIN,E0MAX,TH0MIN,TH0MAX,CTHMIN,CTHMAX,PHI0,
     +              NPPROJ,NPTARG
      COMMON/NUGRID/ZSIZE,HGPOS,DISTX,DISTY,DISTZ,OFFLAY,OFFLX,OFFLY,
     *              OFFX,OFFY,NBINS,NLAY,NSIDE,NSTR,NPMSTR,IDCONF,
     *              ICONTY,PMDIAM,PMRAD,PMR2,PMAREA,XYFUZZ,
     *              VXPLSI, VXPLDI
      PARAMETER ( NLMAX = 4 )
      PARAMETER ( NBB = 160 )
      PARAMETER ( NPM = 10 )
      COMMON/NUPATC/CERPAT(NLMAX,NBB,NBB),TPAT(NLMAX,NBB,NBB)
     *        ,T2PAT(NLMAX,NBB,NBB)
     *        ,TFPAT(NLMAX,NBB,NBB),TLPAT(NLMAX,NBB,NBB)
     *        ,NPHTOT(NLMAX)
      COMMON/NUPM/PM(NLMAX,-NPM:NPM,-NPM:NPM)
     *           ,BINWID,PMTBIN,KBINCX,KBINCY
     *           ,NPHPM(NLMAX,-NPM:NPM,-NPM:NPM)
     *           ,TT(NLMAX,-NPM:NPM,-NPM:NPM)
     *           ,TSIG(NLMAX,-NPM:NPM,-NPM:NPM)
     *           ,TFL(NLMAX,-NPM:NPM,-NPM:NPM)
      PARAMETER ( MAXPM  = 250 )
      PARAMETER ( MXPMHI = 100*MAXPM )
      COMMON/NUCPM/PMINT(MAXPM),TFIRST(MAXPM),TAVE(MAXPM),NPMINT(MAXPM),
     +             NHITOT, IPMHIT(MXPMHI),
     +             TIMHIT(MXPMHI), CPHINT(MXPMHI), NPEINT(MXPMHI),
     +             PHPPMH(MXPMHI), ANGPMH(MXPMHI), BETPMH(MXPMHI),
     +             IPAHIT(MXPMHI)
      COMMON/NUPMPO/NNPM,X0PM(MAXPM),Y0PM(MAXPM),Z0PM(MAXPM),
     +              IORIPM(MAXPM), LAYSTR(MAXPM)
C     IORIPM specifies the orientation of the PM tubes:
C     IORIPM(i) = -1: PM tube number i is looking downward
C     IORIPM(i) = +1: PM tube number i is looking upward
C     LAYSTR(i) = the number of the layer or string that PM tube
C                 number i belongs to
      COMMON/GCFLAG/IDEBUG,IDEMIN,IDEMAX,ITEST,IDRUN,IDEVT,IEORUN
     +        ,IEOTRI,IEVENT,ISWIT(10),IFINIT(20),NEVENT,NRNDM(2)
      COMMON/GCFLAX/BATCH, NOLOG
      LOGICAL BATCH, NOLOG
C
      COMMON/GCTATI/NSEC,KIND(30),EN(30),PL(30),PT(30),THETA(30),PHI(30)
      INTEGER NSEC,KIND
      REAL    EN,PL,PT,THETA,PHI
C
C
            WRITE(6,*) ' '
            WRITE(6,*) 'NUUHB: Booking user histograms.'
            WRITE(6,*) ' '
C
            IF (FGEANT) THEN
              CALL HBOOK2(2001,'Horizontal vertex position: y vs. x',
     +                      40,-HGPOS,HGPOS,   40,-HGPOS,HGPOS,   0.)
C            CALL HBOOK1(2002,'Vertical vertex position (z)',
C     +                            50, -1.25*ZSIZE, -0.75*ZSIZE,   0.)
C            CALL HBOOK1(2003,'Particle type in GUTRAK', 50, 0., 50.,
C     &          0.)
              CALL HBOOK2(2004,'GUTRAK: Momentum vs. Particle type',
     +                              50, 0., 30.,   50, 0., E0MAX, 0.)
C            CALL HBOOK2(2005,
C     +
C     + 'GUTRAK: Also Momentum vs. Particle type (different scale)',
C     +          50, 0., 30.,   50, 0., 300., 0.)
              CALL HBOOK1(2006,'GUTRAK: Momentum (gammas)',
     +                                          50, 0., E0MAX/5., 0.)
              CALL HBOOK1(2007,'GUTRAK: Momentum (e+, e-)',
     +                                          50, 0., E0MAX/5., 0.)
              CALL HBOOK1(2008,'GUTRAK: Momentum (neutrinos)',
     +                                          50, 0., E0MAX/5., 0.)
              CALL HBOOK1(2009,'GUTRAK: Momentum (mu+, mu-)',
     +                                          50, 0., E0MAX/10.,0.)
              CALL HBOOK1(2010,'GUTRAK: Momentum (pi0)',
     +                                          50, 0., E0MAX/10.,0.)
              CALL HBOOK1(2011,'GUTRAK: Momentum (pi+, pi-)',
     +                                          50, 0., E0MAX/10.,0.)
              CALL HBOOK1(2012,'GUTRAK: Momentum (n, p, pbar)',
     +                                          50, 0., E0MAX/10.,0.)
              CALL HBOOK1(2013,'GUTRAK: Momentum (others)',
     +                                          50, 0., E0MAX/10.,0.)
C**********************************************************************
              CALL HBOOK1(4461,'Time in NUSLIB',40,0.,1.0,0)
              CALL HBOOK1(4462,'Energy in NUSLIB',40,0.,100.,0)
              CALL HIDOPT(4462,'LOGY')
              CALL HBOOK1(4463,'Shower use vs LOG(E)',60,-1.,2.,0)
              CALL HBOOK1(4464,'Shower use vs index',40,0.,1000.,0)
              CALL HBOOK1(4465,'Frac in NUSLIB',40,0.,1.,0)
              CALL HBOOK1(4466,'QQCER2',45,0.,1E5,0)
              CALL HBOOK1(4467,'QQSLIB',50,0.,9000.,0)
              CALL HBOOK1(4468,'QQSLIB+QQCER2',45,0.,1E5,0)
              CALL HBOOK1(4470,'NPCER2',50,0.,500.,0)
              CALL HBOOK1(4471,'NPSLIB',50,0.,500.,0)
              CALL HBOOK1(4472,'NPCER2+NPSLIB',50,0.,200.,0)
            ENDIF
C
            IF (KCERAL .GT. 0) THEN
C            CALL HBOOK1(5551,' CHERENKOV PHOTONS',0.,5E6,0.)
              CALL HBOOK1(6671,'Case in NUCERL$',20,0.5,20.5,0)
              CALL HBOOK1(6672,'QQQ in NUCERL$',50,0.,5E5,0)
              CALL HBOOK1(6673,'PHPATH in NUCERL$',40,0.,4000.,0)
              CALL HBOOK1(6674,'RLAMBD in NUCERL$',45,0.,90.,0)
              CALL HBPROF(6675,'QQQ vs JCASE$',20,0.5,20.5,0.,5E7,' ')
              CALL HBOOK1(6676,'Photons/cm in shower$',40,0.,4000.,0)
              CALL HBPROF(6677,'Dist vs JCASE$',20,0.5,20.5,0.,500000.,
     &          ' ')
              CALL HBOOK1(6678,'NTH*-NTL*',10,-0.5,9.5,0)
              CALL HBOOK1(6679,'FRAC when upper chosen',40,0.,1.,0)
              CALL HBOOK1(6680,'FRAC when lower chosen',40,0.,1.,0)
C
              CALL HBOOK1(4001,
     +          '# Cherenkov photons in the PM tubes: Event by event',
     +          50, 0., 200.,   0.)
              CALL HIDOPT(4001,'LOGY')
              IF (KCERAL .EQ. 2) THEN
C              CALL HBOOK1(4002,
C     +            '# Cherenkov photons in the PM tubes: Track by track',
C     +            50, 0., 200.,   0.)
C              CALL HIDOPT(4002,'LOGY')
              ENDIF
              CALL HBOOK1(4003,
     +
     +'Average # Cher. photons in the PMs: Event by event; All PMs',
     +      50, 0., 500.,   0.)
              CALL HIDOPT(4003,'LOGY')
              CALL HBOOK1(4004,
     +
     +'Av. # Cher. ph. in the PMs: Event by event; All struck PMs',
     +      50, 0., 500.,   0.)
              CALL HIDOPT(4004,'LOGY')
              CALL HBOOK1(4005,
     +
     +'# photo-electrons in the PM tubes: Event by event; All tubes',
     +      70, 0., 350.,   0.)
              CALL HIDOPT(4005,'LOGY')
              CALL HBOOK1(4006,
     +
     +'# photo-electrons in the PMs: Ev. by ev.; Firing PMs only',
     +          90, 0., 450.,   0.)
              CALL HIDOPT(4006,'LOGY')
              CALL HBOOK1(4007,
     +
     +'Average # photo-electrons: Event by event; All PM tubes',
     +          60, 0., 300.,   0.)
              CALL HIDOPT(4007,'LOGY')
              CALL HBOOK1(4008,
     +
     +'Average # photo-electrons: Event by event; Firing PMs only',
     +      60, 0., 300.,   0.)
              CALL HIDOPT(4008,'LOGY')
              CALL HBOOK1(4009,
     +
     +'(# photo-electrons) / (# Cherenkov photons) - All struck PMs',
     +      70, -0.5, 3.,   0.)
              CALL HBOOK1(4010,
     +
     +'(# photo-electrons) / (# Cher. photons) - Firing PMs only',
     +          70, -0.5, 3.,   0.)
            ENDIF
C
C          CALL HBOOK1(4711,'Attenuation length (wavelength-dep.)',
C     +25, 0., 2500.,   0.)
C
C          CALL HBOOK1(4712,'Quantum efficiency (wavelength-dep.)',
C     +50, 0., 0.5,   0.)
C
            CALL HBOOK1(4718,'Total # hits in MYTRIG',40,0.,800.,0)
            CALL HBOOK1(4719,'# hit PM tubes in MYTRIG',40,-0.5,79.5,0)
            CALL HBOOK1(4720,'# hit strings in MYTRIG',11,-0.5,10.5,0)
C
            RETURN
            END
CDECK  ID>, NUUFFR. 
            SUBROUTINE NUUFFR
      PARAMETER ( MAXPM  = 250 )
      PARAMETER ( MXPMHI = 100*MAXPM )
      COMMON/NUCPM/PMINT(MAXPM),TFIRST(MAXPM),TAVE(MAXPM),NPMINT(MAXPM),
     +             NHITOT, IPMHIT(MXPMHI),
     +             TIMHIT(MXPMHI), CPHINT(MXPMHI), NPEINT(MXPMHI),
     +             PHPPMH(MXPMHI), ANGPMH(MXPMHI), BETPMH(MXPMHI),
     +             IPAHIT(MXPMHI)
      COMMON/NUPMPO/NNPM,X0PM(MAXPM),Y0PM(MAXPM),Z0PM(MAXPM),
     +              IORIPM(MAXPM), LAYSTR(MAXPM)
C     IORIPM specifies the orientation of the PM tubes:
C     IORIPM(i) = -1: PM tube number i is looking downward
C     IORIPM(i) = +1: PM tube number i is looking upward
C     LAYSTR(i) = the number of the layer or string that PM tube
C                 number i belongs to
C
            WRITE(6,*) ' '
            WRITE(6,*) 'NUUFFR: Defining user data cards.'
C
            RETURN
            END
C+DECK,NUPMUS. *********************************************************
            SUBROUTINE NUPMUS
C------------------------------------------------------
C  My own routine to initialise the PM-tube positions.
C
C  This routine will set up the PM tubes in a simple
C  string geometry, with nine strings each carrying
C  three strings.
C  Hence the values of NNPM, NSTR and NPMSTR as speci-
C  fied in the datacards will be disregarded.
C
C  The vertical distance between the PM tubes on each
C  string will be taken as DISTZ (as specified in the
C  datacards), and the whole detector will be placed
C  symmetrically around Z = 0.
C  The horizontal coordinates of the strings are speci-
C  fied below in this routine.
C
C  The numbering of the PM tubes will basically go from
C  below and upwards, to make the numbering similar to
C  that for layer detectors (ICONTY = 1). Thus the num-
C  bering does not go string by string.
C------------------------------------------------------
C
      PARAMETER ( MAXPM  = 250 )
      PARAMETER ( MXPMHI = 100*MAXPM )
      COMMON/NUCPM/PMINT(MAXPM),TFIRST(MAXPM),TAVE(MAXPM),NPMINT(MAXPM),
     +             NHITOT, IPMHIT(MXPMHI),
     +             TIMHIT(MXPMHI), CPHINT(MXPMHI), NPEINT(MXPMHI),
     +             PHPPMH(MXPMHI), ANGPMH(MXPMHI), BETPMH(MXPMHI),
     +             IPAHIT(MXPMHI)
      COMMON/NUPMPO/NNPM,X0PM(MAXPM),Y0PM(MAXPM),Z0PM(MAXPM),
     +              IORIPM(MAXPM), LAYSTR(MAXPM)
C     IORIPM specifies the orientation of the PM tubes:
C     IORIPM(i) = -1: PM tube number i is looking downward
C     IORIPM(i) = +1: PM tube number i is looking upward
C     LAYSTR(i) = the number of the layer or string that PM tube
C                 number i belongs to
      COMMON/NUGRID/ZSIZE,HGPOS,DISTX,DISTY,DISTZ,OFFLAY,OFFLX,OFFLY,
     *              OFFX,OFFY,NBINS,NLAY,NSIDE,NSTR,NPMSTR,IDCONF,
     *              ICONTY,PMDIAM,PMRAD,PMR2,PMAREA,XYFUZZ,
     *              VXPLSI, VXPLDI
C
            DATA NNSTR, NNPMST, NNNPM /9, 3, 27/
C
            REAL XSTRNG(9), YSTRNG(9)
            DATA XSTRNG/-500., -500., -500., 0., 0., 0., 500., 500.,
     &        500./
            DATA YSTRNG/-500., 0., 500., -500., 0., 500., -500., 0.,
     &        500./
C
            ICONTY = 2
C
C     Check if we override any datacards:
            IF (NNSTR .NE. NSTR) THEN
              WRITE(6,*) '*** NUPMUS: ***'
              WRITE(6,*) 'Number of strings set to',NNSTR
              WRITE(6,*) '***************'
            ENDIF
            IF (NNPMST .NE. NPMSTR) THEN
              WRITE(6,*) '*** NUPMUS: ***'
              WRITE(6,*) 'Number of PM tubes per string set to',NNPMST
              WRITE(6,*) '***************'
            ENDIF
            IF (NNNPM .NE. NNPM) THEN
              WRITE(6,*) '*** NUPMUS: ***'
              WRITE(6,*) 'Number of PM tubes set to',NNNPM
              WRITE(6,*) '***************'
            ENDIF
C
            NSTR   = NNSTR
            NPMSTR = NNPMST
            NNPM   = NNNPM
C
C     Set # 'layers' equal to # PM tubes per string:
            NLAY = NPMSTR
C
C      WRITE(6,*) 'NNPM,NSTR,NPMSTR = ',NNPM,NSTR,NPMSTR
C
            DO 10 IZPOS = 1,NPMSTR
              DO 20 ISTR = 1,NSTR
                IJK = (IZPOS - 1)*NSTR + ISTR
                X0PM(IJK) = XSTRNG(ISTR)
                Y0PM(IJK) = YSTRNG(ISTR)
                Z0PM(IJK) = (IZPOS - 0.5*NPMSTR - 0.5)*DISTZ
                IORIPM(IJK) = -1
                LAYSTR(IJK) = ISTR
C            WRITE(6,*) 'IJK, X0PM,Y0PM,Z0PM = ',
C     +                  IJK,X0PM(IJK),Y0PM(IJK),Z0PM(IJK)
   20         CONTINUE
   10       CONTINUE
C-->  Consistency checks:
C
            RETURN
            END
CDECK  ID>, FSFLUX. 
            FUNCTION FSFLUX(THETA,PHI)
      COMMON/NUGENC/E0MIN,E0MAX,TH0MIN,TH0MAX,CTHMIN,CTHMAX,PHI0,
     +              NPPROJ,NPTARG
      COMMON/GCONST/PI,TWOPI,PIBY2,DEGRAD,RADDEG,CLIGHT,BIG,EMASS
      COMMON/GCONSX/EMMU,PMASS,AVO
C
C
            EXTERNAL FCOSSQ
            REAL FCOSSQ
C
            REAL FSPACE(100)
            SAVE FSPACE, ICALL
            DATA ICALL/0/
C
            IF (ICALL .EQ. 0) THEN
              ANGMIN = TH0MIN*DEGRAD
              ANGMAX = TH0MAX*DEGRAD
              CALL FUNPRE(FCOSSQ,FSPACE,ANGMIN,ANGMAX)
              ICALL = 1
              CALL HBOOK1(4711,'FSFLUX: Theta distribution',
     +               60, TH0MIN, TH0MAX, 0.)
              CALL HBOOK1(4712,'FSFLUX: Phi distribution',
     +               36, 0., 360., 0.)
C            CALL HBOOK2(4713,'FSFLUX: Theta versus Phi',
C     +               36, 0., 360., 60, TH0MIN, TH0MAX, 0.)
              CALL HMINIM(4711,0.)
              CALL HMINIM(4712,0.)
            ENDIF
C
            CALL FUNRAN(FSPACE,XRAN)
            THETA = XRAN
C
C        PHI = TWOPI*RNDM(THETA)
            PHI = TWOPI*RNUNIF(THETA)
C
            THETA = THETA*RADDEG
            PHI   = PHI*RADDEG
            CALL HFILL(4711,THETA,0., 1.)
            CALL HFILL(4712,PHI,0.,   1.)
C        CALL HFILL(4713,PHI,THETA,1.)
C
            FSFLUX = 1.
C
            RETURN
            END
            FUNCTION FCOSSQ(ANGLE)
            REAL ANGLE
            COSANG = COS(ANGLE)
            COS2 = COSANG*COSANG
            FCOSSQ = COS2
            RETURN
            END
CDECK  ID>, MYTRIG. 
            SUBROUTINE MYTRIG
C------------------------------------------------------------------------
C
C  This is an special trigger routine for the effective-area studies.
C  It based on NUTRG2, rather than on NUTRIG.
C  The special thing with this routine is that it applies five diffe-
C  rent trigger conditions, instead of only one.
C  It still defines a 'standard' trigger, though: it is the usual
C  '8-3' one.
C
C  An in NUTRIG, a 'hit' in a PM tube is defined by # photo-electrons in
C  that tube exceeding the value of NPEMIN in /NUTRIC/.
C
C  So far this routine is in a very preliminary stage, with only the
C  ANANDA trigger implemented. For other detector configurations, the
C  trigger condition is considered as being always fulfilled.
C
C                                       Bengt Erlandsson 15/1/93
C
C------------------------------------------------------------------------
C
      COMMON/NUGRID/ZSIZE,HGPOS,DISTX,DISTY,DISTZ,OFFLAY,OFFLX,OFFLY,
     *              OFFX,OFFY,NBINS,NLAY,NSIDE,NSTR,NPMSTR,IDCONF,
     *              ICONTY,PMDIAM,PMRAD,PMR2,PMAREA,XYFUZZ,
     *              VXPLSI, VXPLDI
      PARAMETER ( MAXPM  = 250 )
      PARAMETER ( MXPMHI = 100*MAXPM )
      COMMON/NUCPM/PMINT(MAXPM),TFIRST(MAXPM),TAVE(MAXPM),NPMINT(MAXPM),
     +             NHITOT, IPMHIT(MXPMHI),
     +             TIMHIT(MXPMHI), CPHINT(MXPMHI), NPEINT(MXPMHI),
     +             PHPPMH(MXPMHI), ANGPMH(MXPMHI), BETPMH(MXPMHI),
     +             IPAHIT(MXPMHI)
      COMMON/NUPMPO/NNPM,X0PM(MAXPM),Y0PM(MAXPM),Z0PM(MAXPM),
     +              IORIPM(MAXPM), LAYSTR(MAXPM)
C     IORIPM specifies the orientation of the PM tubes:
C     IORIPM(i) = -1: PM tube number i is looking downward
C     IORIPM(i) = +1: PM tube number i is looking upward
C     LAYSTR(i) = the number of the layer or string that PM tube
C                 number i belongs to
      COMMON/NUTRIC/LTRIG, NPEMIN
      LOGICAL LTRIG
C     LTRIG  = Trigger flag; to be set .TRUE. if the current event has
C              generated a trigger
C     NPEMIN = Minimum # photo-electrons defining a signal in a PM tube
C              (to be used in the trigger algorithm)
      COMMON /MYCOUN/ NTRGEV(5), NQ(100)
C
            PARAMETER (MAXSTR = 100)
            INTEGER NHITPM, IHISTR(MAXSTR), NHISTR
C
C    Minimum # PM tubes required in the AMANDA trigger:
            DATA MINPM/8/
C    Minimum # strings required in the AMANDA trigger:
            DATA MINSTR/3/
C
            DATA ICALL/0/
C
            SAVE MINPM, MINSTR, ICALL
C
            IF (IDCONF .NE. 4) THEN
              LTRIG = .TRUE.
              GOTO 999
            ENDIF
C
C    Consider the AMANDA case in the rest of this routine!
C
            IF (ICALL .EQ. 0) THEN
              ICALL = 1
              IF (MINPM  .GT. NNPM) MINPM  = NNPM
              IF (MINSTR .GT. NSTR) MINSTR = NSTR
              IF (NPEMIN .GT. 1) THEN
                WRITE(6,1) MINPM, NPEMIN, MINSTR
              ELSEIF (NPEMIN .EQ. 1) THEN
                WRITE(6,2) MINPM, NPEMIN, MINSTR
              ELSE
                WRITE(6,3) NPEMIN
              ENDIF
            ENDIF
C
C    NPEMIN = 0 means that the trigger condition has been deactivated:
            IF (NPEMIN .LE. 0) THEN
              LTRIG = .TRUE.
              GOTO 999
            ENDIF
C
            IF (NSTR .GT. MAXSTR) THEN
              WRITE(6,4) NSTR, MAXSTR
              LTRIG = .TRUE.
              GOTO 999
            ENDIF
C
            LTRIG = .FALSE.
C
            CALL VZERO(IHISTR,NSTR)
            NHITPM = 0
            NHISTR = 0
C
C    Loop over all PM tubes, checking how many of them have been hit
C    (with # photo-electrons above threshold), and counting the number
C    of hit strings:
C
            NHTOT = 0
            DO 10 IPM = 1,NNPM
              IF (IORIPM(IPM) .GT. 0) GOTO 10
              NHTOT = NHTOT + NPMINT(IPM)
              IF (NPMINT(IPM) .GE. NPEMIN) THEN
                NHITPM = NHITPM + 1
                ISTR = LAYSTR(IPM)
                IF ((ISTR .LE. 0) .OR. (ISTR .GT. NSTR)) THEN
                  WRITE(6,5) ISTR
                  GOTO 10
                ENDIF
                IF (IHISTR(ISTR) .LE. 0) THEN
                  IHISTR(ISTR) = 1
                  NHISTR = NHISTR + 1
                ENDIF
              ENDIF
   10       CONTINUE
C
            CALL HFILL(4718,FLOAT(NHTOT),0.,1.)
            CALL HFILL(4719,FLOAT(NHITPM),0.,1.)
            CALL HFILL(4720,FLOAT(NHISTR),0.,1.)
C
C    Check the various trigger conditions:
C    1)  '6-3':
            IF ((NHITPM .GE.  6) .AND. (NHISTR .GE. 3))
     +        NTRGEV(1) = NTRGEV(1) + 1
C    2)  '7-3':
            IF ((NHITPM .GE.  7) .AND. (NHISTR .GE. 3))
     +        NTRGEV(2) = NTRGEV(2) + 1
C    3)  '8-3':
            IF ((NHITPM .GE.  8) .AND. (NHISTR .GE. 3))
     +        NTRGEV(3) = NTRGEV(3) + 1
C    4)  '9-3':
            IF ((NHITPM .GE.  9) .AND. (NHISTR .GE. 3))
     +        NTRGEV(4) = NTRGEV(4) + 1
C    5) '10-3':
            IF ((NHITPM .GE. 10) .AND. (NHISTR .GE. 3))
     +        NTRGEV(5) = NTRGEV(5) + 1
C
C    Finally, check the 'standard' trigger condition:
CCC      LTRIG = ((NHITPM .GE. MINPM) .AND. (NHISTR .GE. MINSTR))
            LTRIG = ((NHITPM .GE. 8) .AND. (NHISTR .GE. 3))
C
    1       FORMAT(' '/' ',46('*')/' ','*',44X,'*'/
     +        ' ','*   AMANDA trigger:',26X,'*'/
     +        ' ','*   ',15('-'),26X,'*'/' ','*',44X,'*'/
     +        ' ','*   Require at least',I5,' PM tubes',11X,'*'/
     +        ' ','*   with at least',5X,I3,' photo-electrons,   *'/
     +        ' ','*   distributed over at least',I4,' strings',4X,'*'/
     +        ' ','*',44X,'*'/' ',46('*')/' ')
    2       FORMAT(' '/' ',45('*')/' ','*',43X,'*'/
     +        ' ','*   AMANDA trigger:',25X,'*'/
     +        ' ','*   ',15('-'),25X,'*'/' ','*',43X,'*'/
     +        ' ','*   Require at least',I5,' PM tubes',10X,'*'/
     +        ' ','*   with at least',5X,I3,' photo-electron,   *'/
     +        ' ','*   distributed over at least',I4,' strings',3X,'*'/
     +        ' ','*',43X,'*'/' ',45('*')/' ')
    3       FORMAT(' '/' ',41('*')/' ','*',39X,'*'/
     +        ' ','*   AMANDA trigger:',21X,'*'/
     +        ' ','*   ',15('-'),21X,'*'/' ','*',39X,'*'/
     +        ' ','*   DEACTIVATED, because NPEMIN =',I4,3X,'*'/
     +        ' ','*',39X,'*'/' ',41('*')/' ')
    4       FORMAT(' '/' ','*** MYTRIG: ***'/
     +        ' ','*** Too many strings:',I5/
     +        ' ','*** Max # strings allowed by this routine is',I5/
     +        ' ','*** ==> Might be necessary to change in this code'/
     +        ' ','*** (Setting triger flag = .TRUE. for this event)')
    5       FORMAT(' '/' ','*** MYTRIG: ***'/
     +        ' ','*** Hit in non-existent string: number',I5)
C
  999       CONTINUE
C
            RETURN
            END
CDECK  ID>, NUPMAM. 
            SUBROUTINE NUPMAM
C----------------------------------------------------------------------
C
C  This routine sets up the PM-tube positions according to the AMANDA
C  geometry.
C
C  --- Special version for the `DOWN 20 PMT' set-up. ---
C
C  This is a string geometry, with the number of strings and the number
C  of PM tubes per string fixed by the AMANDA specifications. Thus the
C  values of NSTR, NNPM, and NPMSTR, as specified in the datacards, will
C  (to some extent) be disregarded. (Where needed, they will be adjusted
C  to the appropriate AMANDA values.)
C
C  The full AMANDA set-up is anticipated to have nine strings. Thus here
C  we allow for 1 - 10 strings (specified by datacard NSTR; the default
C  is taken to be 10). String number 1 is always taken to be that in the
C  centre.
C
C  The AMANDA detector will possibly have both downward-looking and up-
C  ward-looking PM tubes. In this implementation, the downward-looking
C  tubes and the upward-looking tubes are installed 'on an equal foot-
C  ing', although the former are more numerous.
C
C  The vertical distance between the PM tubes on each string is also
C  fixed by the AMANDA geometry. The value of DISTZ will be adjusted to
C  the distance between the tubes.
C
C  The central string will be placed at (x,y) = (0,0). As said above,
C  this is taken to be string number 1.
C
C  In reality, the AMANDA detector array is going to be placed approxi-
C  mately one kilometre below ground level. This fact is not taken into
C  account here, where the detector will instead be centred at z = 0.
C
C  All strings are assumed to be perfectly vertical.
C
C  The numbering of the PM tubes will go string by string, from below
C  and upwards.
C
C                                       Bengt Erlandsson 3/11/92
C
C----------------------------------------------------------------------
      PARAMETER ( MAXPM  = 250 )
      PARAMETER ( MXPMHI = 100*MAXPM )
      COMMON/NUCPM/PMINT(MAXPM),TFIRST(MAXPM),TAVE(MAXPM),NPMINT(MAXPM),
     +             NHITOT, IPMHIT(MXPMHI),
     +             TIMHIT(MXPMHI), CPHINT(MXPMHI), NPEINT(MXPMHI),
     +             PHPPMH(MXPMHI), ANGPMH(MXPMHI), BETPMH(MXPMHI),
     +             IPAHIT(MXPMHI)
      COMMON/NUPMPO/NNPM,X0PM(MAXPM),Y0PM(MAXPM),Z0PM(MAXPM),
     +              IORIPM(MAXPM), LAYSTR(MAXPM)
C     IORIPM specifies the orientation of the PM tubes:
C     IORIPM(i) = -1: PM tube number i is looking downward
C     IORIPM(i) = +1: PM tube number i is looking upward
C     LAYSTR(i) = the number of the layer or string that PM tube
C                 number i belongs to
      COMMON/NUGRID/ZSIZE,HGPOS,DISTX,DISTY,DISTZ,OFFLAY,OFFLX,OFFLY,
     *              OFFX,OFFY,NBINS,NLAY,NSIDE,NSTR,NPMSTR,IDCONF,
     *              ICONTY,PMDIAM,PMRAD,PMR2,PMAREA,XYFUZZ,
     *              VXPLSI, VXPLDI
      COMMON/NUMISC/START,  IPTYPE, PIN, EMU0, THEMU0, PHIMU0,
     +              ZTRMAX, IORIVX, IORITR, VOLD(7),
     +              KTMED,  REFIND, ATTLEN, CWLGT1, CWLGT2,
     +              CUTKC2, IPRHIS(3), TDIAM, INUTRA,
     +              ZCLAPP, DCLAPP
      INTEGER START
C
C # strings, # PM tubes per string, and total # PM tubes:
            PARAMETER (NNSTR  = 6)
            PARAMETER (NNPMST = 13)
            PARAMETER (NNNPM  = NNSTR*NNPMST)
C
C Vertical spacing between the PM tubes (in cm.):
            PARAMETER (ZSPACE = 1500.)
C
C Distance between the central string and those at the periphery (cm):
C +SELF,IF=-60MRAD.
            PARAMETER (RADSTR = 5000.)
C
C Diameter of the AMANDA PM tubes (in cm.):
            PARAMETER (AMPMDI = 19.)
C


C
C Attenuation length in the tracking medium (ice) of AMANDA, in cm.
C (included here for a consistency check):
CCC      PARAMETER (ATTAMA = 1200.)
C
CCC      INTEGER IORIAM(NNPMST)
CCC      DATA IORIAM/ -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
CCC     -             -1,-1,-1,-1,-1,-1,-1,-1,-1,-1/
C
            REAL XSTRNG(NNSTR), YSTRNG(NNSTR)
C
            PI = 4.*ATAN(1.)
C
            WRITE(6,101)
     &        'NUPMAM: initialising AMANDA PM-tube positions.'
C
C     Check if we override any datacards:
            IF (NSTR .NE. NNSTR) THEN
              WRITE(6,102)
     +          'Wrong number of strings requested:',NSTR,
     +          'Number of strings set to',  NNSTR
              NSTR = NNSTR
            ENDIF


            IF ((NPMSTR .NE. 13) .AND. (NPMSTR .NE. 20)) THEN
              WRITE(6,104)
     +          'Number of PM tubes per string set to',NNPMST
              NPMSTR = NNPMST
            ENDIF
            IF (NNPM .NE. (NSTR*NPMSTR)) THEN
              WRITE(6,105)
     +          'Number of PM tubes set to',NSTR*NPMSTR
              NNPM = NSTR*NPMSTR
            ENDIF
            IF ((DISTZ .NE. 1500.) .AND. (DISTZ .NE. 2000.)) THEN
              WRITE(6,106) 'DISTZ changed to',ZSPACE,' cm.'
              DISTZ  = ZSPACE
            ENDIF
            IF (AMPMDI .NE. PMDIAM) THEN
              WRITE(6,107) 'PMDIAM changed to',AMPMDI,' cm.'
            ENDIF
C
            WRITE(6,109)
C
            IF ((NPMSTR .NE. 13) .AND. (NPMSTR .NE. 20)) THEN
              WRITE(6,104)
     +          'Number of PM tubes per string set to',NNPMST
              NPMSTR = NNPMST
            ENDIF
            IF (NNPM .NE. (NSTR*NPMSTR)) THEN
              WRITE(6,105)
     +          'Number of PM tubes set to',NSTR*NPMSTR
              NNPM = NSTR*NPMSTR
            ENDIF
            IF ((DISTZ .NE. 1500.) .AND. (DISTZ .NE. 2000.)) THEN
              WRITE(6,106) 'DISTZ changed to',ZSPACE,' cm.'
              DISTZ  = ZSPACE
            ENDIF
            IF (AMPMDI .NE. PMDIAM) THEN
              WRITE(6,107) 'PMDIAM changed to',AMPMDI,' cm.'
            ENDIF
C
            WRITE(6,109)
C
  101       FORMAT(//' ',54('*')/' ','*',52X,'*'/' ','*   ',A,'   *'/
     +         ' ','*   ',46('='),'   *')
  102       FORMAT(' ','*',52X,'*'/' ','*   *** ',A,I4,' ***',11X,'*'/
     +        ' ','*   *** ',A,2X,I4,' ***',11X,'*')
  103       FORMAT(' ','*',52X,'*'/' ','*   ',A,21X,'*'/
     +        ' ','*   ',A,I4,25X,'*')
  104       FORMAT(' ','*',52X,'*'/' ','*   *** ',A,I4,' *** *')
  105       FORMAT(' ','*',52X,'*'/' ','*   *** ',A,I4,' ***',12X,'*')
  106       FORMAT(' ','*',52X,'*'/' ','*   *** ',A,F9.2,A,' ***',12X,
     &        '*')
  107       FORMAT(' ','*',52X,'*'/' ','*   *** ',A,F6.2,A,' ***',14X,
     &        '*')
  108       FORMAT(' ','*',52X,'*'/' ','*',3X,45('*'),4X,'*'/
     +        ' ','*   ',A,F10.3,A,4X,'*'/' ','*   ',A,4X,F10.3,A,4X,
     +        '*'/
     +        ' ','*   ',A,8X,'*'/' ','*',3X,45('*'),4X,'*')
  109       FORMAT(' ','*',52X,'*'/' ',54('*')/)
C
            PMDIAM = AMPMDI
            PMRAD  = 0.5*PMDIAM
            PMR2   = PMRAD*PMRAD
            PMAREA = PI*PMR2
C
C    Set # 'layers' equal to # PM tubes per string:
            NLAY = NPMSTR
C
C
            WRITE(6,*) ' '
            WRITE(6,*) 'NUPMAM: NNPM, NSTR = ',NNPM,NSTR
            WRITE(6,*) '        NPMSTR     = ',NPMSTR
C
C    Angular spacing between the strings along the periphery:
            PHISTR = 2.*PI/FLOAT(NNSTR-1)
C
C    Spacing between the strings along the periphery (in cm.):
            DPERST = 2.*RADSTR*SIN(PHISTR/2.)
C
C    Configuration type:
C    -------------------

            ICONTY = 2
C
C    Now set up the string positions:
C    --------------------------------
C
            XSTRNG(1) = 0.
            YSTRNG(1) = 0.
            PHI = -PHISTR/2.
            DO 10 ISTR = 2,NNSTR
              PHI = PHI + PHISTR
              XSTRNG(ISTR) = RADSTR*COS(PHI)
              YSTRNG(ISTR) = RADSTR*SIN(PHI)
   10       CONTINUE

C    Now set up the PM-tube positions:
C    ---------------------------------
C
C    Vertical extension of the detector array:
            ZDETEX = FLOAT(NPMSTR - 1)*DISTZ
            IF ((ZDETEX/2.) .GT. ZSIZE) THEN
              WRITE(6,*) ' '
              WRITE(6,*) '*** NUPMAM: WARNING! ***'
              WRITE(6,*) '*** GEANT volume is not large enough! ***'
            ENDIF
C
            ZPMUND =  -ZDETEX/2.
C
            IJK = 0
            DO 20 ISTR = 1,NSTR
              WRITE(6,201) 'String number',ISTR,':'

              ZPOS = ZPMUND - DISTZ
              DO 30 IZ = 1,NPMSTR
                IJK  = IJK  + 1
                ZPOS = ZPOS + DISTZ
                X0PM(IJK)   = XSTRNG(ISTR)
                Y0PM(IJK)   = YSTRNG(ISTR)
                Z0PM(IJK)   = ZPOS
CCC            IORIPM(IJK) = IORIAM(IZ)
                IORIPM(IJK) = -1
                LAYSTR(IJK) = ISTR
                IF (IORIPM(IJK) .LT. 0) THEN
                  WRITE(6,202) 'IJK, X0PM,Y0PM,Z0PM = ',
     +                       IJK,X0PM(IJK),Y0PM(IJK),Z0PM(IJK)
                ELSE
                  WRITE(6,203) 'IJK, X0PM,Y0PM,Z0PM = ',
     +                       IJK,X0PM(IJK),Y0PM(IJK),Z0PM(IJK),' *'
                ENDIF
   30         CONTINUE
   20       CONTINUE
C
  201       FORMAT(' '/' ',A,I3,A)
  202       FORMAT(' ',A,I4,3(2X,F10.3))
  203       FORMAT(' ',A,I4,3(2X,F10.3),A)
C
            RETURN
            END


CDECK  ID>, NUOMEF. 
            SUBROUTINE NUOMEF(IPM,ANGLE,PHPATH,BETA,NCPHOT,NPELEC)
C------------------------------------------------------------------------
C
C  *** THIS ROUTINE IS SO FAR MAINLY ADAPTED TO WORK WELL FOR AMANDA ***
C  *** SIMULATIONS.                                                  ***
C
C  This is a more sophisticated routine than NUPMEF. 'OM' in the subrou-
C  tine name stands for 'Optical Module' (instead of Photo-Multiplier).
C
C  Note that NUOMEF is a SUBROUTINE, whereas NUPMEF is a REAL FUNCTION.
C
C  Just like NUPMEF, the routine is expected to change rather substan-
C  tially in the future, developing gradually.
C
C  ANGLE is the angle between the photon path and the PM-tube axis,
C  defined such that 0 < ANGLE < pi/2 for a photon hitting the PM tube
C  from 'above', i.e., from the side where it is the most likely to be
C  'seen' by the PM tube.
C
C  Not taken into account in this version is, for example, the point of
C  impact on the optical module.
C
C  In the future presumably each PM tube in the detector should be as-
C  signed an integer number indicating its type, and QUANEF below should
C  depend on the type of tube.
C
C  Of course the user may provide his/her/its own version of NUOMEF,
C  where he/she/it includes whatever effects he/she/it likes.
C
C                                               Bengt Erlandsson 9/11/92
C
C  Angular efficiency updated (for AMANDA) 930118 / B. E.
C  Updated again for AMANDA                930129 / B. E.
C
C------------------------------------------------------------------------
C
C Input:  IPM    = PM-tube number (for future use)
C         ANGLE  = Angle between photon path and PM-tube axis
C         PHPATH = Length of photon path from radiating track to PM tube
C         BETA   = Beta of the radiating particle - NOT USED IN THE CUR-
C                  RENT IMPLEMENTATION, WHERE INSTEAD BETA = 1 IS ASSUMED.
C         NCPHOT = Average # (INTEGER) Cherenkov photons (in the wave-
C                 length interval considered) reaching the PM tube
C
C Output: NPELEC = Number of photo-electrons generated
C
            INTEGER IPM,   NCPHOT, NPELEC
            REAL    ANGLE, PHPATH, BETA
C
      COMMON/NUMISC/START,  IPTYPE, PIN, EMU0, THEMU0, PHIMU0,
     +              ZTRMAX, IORIVX, IORITR, VOLD(7),
     +              KTMED,  REFIND, ATTLEN, CWLGT1, CWLGT2,
     +              CUTKC2, IPRHIS(3), TDIAM, INUTRA,
     +              ZCLAPP, DCLAPP
      INTEGER START
      COMMON/NUGRID/ZSIZE,HGPOS,DISTX,DISTY,DISTZ,OFFLAY,OFFLX,OFFLY,
     *              OFFX,OFFY,NBINS,NLAY,NSIDE,NSTR,NPMSTR,IDCONF,
     *              ICONTY,PMDIAM,PMRAD,PMR2,PMAREA,XYFUZZ,
     *              VXPLSI, VXPLDI
      COMMON/GCONST/PI,TWOPI,PIBY2,DEGRAD,RADDEG,CLIGHT,BIG,EMASS
      COMMON/GCONSX/EMMU,PMASS,AVO
C
      PARAMETER ( MAXPM  = 250 )
      PARAMETER ( MXPMHI = 100*MAXPM )
      COMMON/NUCPM/PMINT(MAXPM),TFIRST(MAXPM),TAVE(MAXPM),NPMINT(MAXPM),
     +             NHITOT, IPMHIT(MXPMHI),
     +             TIMHIT(MXPMHI), CPHINT(MXPMHI), NPEINT(MXPMHI),
     +             PHPPMH(MXPMHI), ANGPMH(MXPMHI), BETPMH(MXPMHI),
     +             IPAHIT(MXPMHI)
      COMMON/NUPMPO/NNPM,X0PM(MAXPM),Y0PM(MAXPM),Z0PM(MAXPM),
     +              IORIPM(MAXPM), LAYSTR(MAXPM)
C     IORIPM specifies the orientation of the PM tubes:
C     IORIPM(i) = -1: PM tube number i is looking downward
C     IORIPM(i) = +1: PM tube number i is looking upward
C     LAYSTR(i) = the number of the layer or string that PM tube
C                 number i belongs to
C                  (Or maybe already taken into account in the calling
C                  routine?)
C
C     ===========================================================
C     Note that in this routine the unit for wavelength is always
C     nanometres.
C     ===========================================================
C
C    Fine structure constant:
            PARAMETER (ALPHA = 1./137.036)
C
C    For the calculation of angular efficiency. The numbers used are
C    based on values for the effective area of the Thorn EMI 8" PM
C    tubes as a function of angle (in 5-degree bins from 0 to 175
C    degrees), as reported by Tim Miller at Berkeley (see e-mails from
C    him of Jan. 7 and Jan. 18, 1993). These numbers take into account
C    both the correction for the Fresnel transmission in the ice-glass
C    interface and the correction for the shielding by hose clamps, etc.
C    Here we use these numbers to obtain the angular efficiency by de-
C    fining
C           ANGEFF = [Eff. area (angle)] / [Total PM-tube area]
C    which is probably the appropriate thing to do, since in NUCER2 we
C    assume spherical PM tubes of a certain size. (Note that we should
C    not divide by PMAREA, but by the total gemoetric photo-cathode area
C    of the particular PM tube that the values of AREANW were obtained
C    for; this area is the parameter AREMI8.)
C    The effective areas were given for certain incident angles (0, 5,
C    10 degrees, etc.), not for certain angular bins. I define the bins
C    by taking the mean of nearby values - the lower limit of -2.5 deg-
C    rees is of course not physical, but is only used a tool for simpli-
C    fying the code.
            PARAMETER (NBIAEF = 36)
            PARAMETER (AREMI8 = 270.)
            REAL AREANW(NBIAEF), ANGEFV(NBIAEF), ANGEBI, ANGEMI, ANGEMX
            DATA ANGEMI, ANGEMX /-2.5, 177.5/
            DATA AREANW / 299.045435, 297.498211, 293.954656, 287.
     &        357333,
     +              278.639762, 268.922298, 258.396415, 247.743225,
     +              234.673816, 221.741833, 208.187864, 194.998804,
     +              181.245966, 166.955534, 152.432254, 138.275631,
     +              125.436807, 111.934269,  98.596700,  85.448248,
     +              73.012225,  60.795377 , 49.870054,  41.480492,
     +              32.831791,  24.854928,  18.050132,  12.299419,
     +              7.741681,   3.841417,   1.553673,   0.319095,
     +              0.013892,   0.0,        0.0,        0.0/
C
C    For the estimation of quantum efficiencty vs. wavelength:
            PARAMETER (NWBINQ = 31)
            PARAMETER (WLMINQ = 295.)
            PARAMETER (WLMAXQ = 605.)
            REAL WLBINQ
            REAL QE(NWBINQ)
            DATA QE/0.010, 0.220, 0.265, 0.290, 0.295, 0.300, 0.300, 0.
     &        300,
     +        0.295, 0.290, 0.285, 0.275, 0.270, 0.260, 0.250, 0.240,
     +        0.225, 0.215, 0.200, 0.180, 0.160, 0.140, 0.120, 0.095,
     +        0.080, 0.065, 0.050, 0.035, 0.030, 0.020, 0.015/
C
C
C    Standard attenuation lengths for south pole ice, as measured by
C    Warren. Note that Warren lists values for every 0.1e-5 cm. from
C    4.0e-5 cm. to 6.0e-5 cm., but only for every 0.5e-5 cm. for shorter
C    wavelengths. The values for 3.1e-5, 3.2e-5, 3.3e-5, 3.4e-5, 3.6e-5,
C    3.7e-5, 3.8e-5, and 3.9e-5 are interpolated.
C    (Note: We use wavelengths specified in nanometres, not cm.)
            PARAMETER (NWBINA = 31)
            PARAMETER (WLMINA = 295.)
            PARAMETER (WLMAXA = 605.)
            REAL WLBINA
C    (Note: The values of the attenuation length are given in centi-
C    metres!)
            REAL ATTVEC(NWBINA)
            DATA ATTVEC/  433.743489,  479.010026,  530.516480,  590.
     &        125074,
     +              658.305122,  739.764016,  806.982815,  881.546846,
     +              966.116274, 1062.849798, 1174.575233, 1299.871057,
     +              1478.873374, 1645.111201, 1833.198309, 2325.315742,
     +              2392.525304, 2412.994314, 2329.096743, 2190.615804,
     +              2083.179897, 1896.472465, 1830.986083, 1660.474811,
     +              1466.615526, 1407.318638, 1354.510163, 1288.612474,
     +              1143.448856,  962.104683,  833.271959/

            DATA ICALL/0/
C
            SAVE
C
            IF (ICALL .EQ. 0) THEN
C
              TWOPIA = TWOPI*ALPHA
              CWTERM = 1./CWLGT1 - 1./CWLGT2
C
C       Initialise the values for the angular efficiencies. Note that
C       the effective areas AREANW were obtained for 8" tubes, so it
C       is this area, and not PMAREA, that should be used when conver-
C       ting AREANW to angular-efficiency values. This area is AREMI8,
C       which is not exactly equal to pi*(0.5*8")**2.
C       Note that the angular bins are converted to radians.
              AREINV = 1./AREMI8
              CALL VSCALE(AREANW,AREINV,ANGEFV,NBIAEF)
              ANGEMI = ANGEMI*DEGRAD
              ANGEMX = ANGEMX*DEGRAD
              ANGEBI = (ANGEMX - ANGEMI)/FLOAT(NBIAEF)
C
C
C       Bin width for the calculation of quantum efficiency:
              WLBINQ = (WLMAXQ - WLMINQ)/FLOAT(NWBINQ)
C
              IF (IDCONF .EQ. 4) THEN
C
C          OK, let's restrict the use of the wavelength-dependent
C          attenuation length to the AMANDA case:
                WRITE(6,*) ' '
                WRITE(6,*) '*** NUOMEF:'
                WRITE(6,*)
     +            '*** The value of',ATTLEN,
     +            ' cm. for the attenuation length,'
                WRITE(6,*)
     +
     +'*** as specified in datacards (or as default value), ',
     +
     +'*** will be taken as the  p e a k  attenuation length.'
                WRITE(6,*) ' '
                ATTFAC = ATTLEN/2412.994314
                DO 1 I = 1,NWBINA
                  ATTVEC(I) = ATTVEC(I)*ATTFAC
    1           CONTINUE
C
                WLBINA = (WLMAXA - WLMINA)/FLOAT(NWBINA)
C
              ENDIF
C
              ICALL  = 1
C
            ENDIF
C
C    Don't waste time on a trivial case:
            IF (NCPHOT .LE. 0) THEN
              NPELEC = 0
              RETURN
            ENDIF
C
C    Next calculate the angular efficiency. This should be independent
C    of wavelength.
C
C    Preliminary angular efficiency of AMANDA tubes:
            IF (ANGLE .LT. -0.001) THEN
C       (I didn't test on exactly zero, because I want to allow for
C       some presicion error. At zero degrees we have the HIGHEST
C       angular efficiency, not the lowest.)
              ANGDEG = ANGLE*RADDEG
              WRITE(6,*) ' '
              WRITE(6,*) '*** NUOMEF:'
              WRITE(6,*) '*** Angle less than zero: ',ANGDEG,
     &          '   degrees!'
              ANGEFF = ANGEFV(1)
            ELSEIF (ANGLE .GT. PI) THEN
              ANGDEG = ANGLE*RADDEG
              WRITE(6,*) ' '
              WRITE(6,*) '*** NUOMEF:'
              WRITE(6,*) '*** Angle greater than 180 degrees: ',ANGDEG,
     +              '   degrees!'
              ANGEFF = ANGEFV(NBIAEF)
            ELSEIF (ANGLE .GT. ANGEMX) THEN
              ANGEFF = ANGEFV(NBIAEF)
            ELSE
              IAEBIN = INT((ANGLE-ANGEMI)/ANGEBI) + 1
              ANGEFF = ANGEFV(IAEBIN)
            ENDIF
C
            IF (IDCONF .EQ. 4) THEN
C        Reduced angular efficiency due to shielding by hose clamps,
C        etc. (AMANDA):
              ANGEFF = ANGEFF*0.89


            ENDIF
C
C    Note: We might well have ANGEFF > 1, because the values of AREANW
C    include a correction factor for the QUANTUM EFFICIENCY (for this
C    Fresnel thing with the ice-glass interface). That is, we should
C    NOT do anything of the kind "IF (ANGEFF .GT. 1.) ANGEFF = 1."
C    here.
C
C    For the generation of wavelengths of the Cherenkov photons:
            REFBET = REFIND*BETA
            BNTERM = (1. - 1./REFBET)*(1. + 1./REFBET)
            CWLCON = TWOPIA*BNTERM*1.E7
            CWLINT = CWLCON*CWTERM
C
C    Loop over the NCPHOT Cherenkov photons, and check for each of them
C    if it actually managed to reach the PM tube and produce a photo-
C    electron:
C
            NPELEC = 0
C
            DO 10 I = 1,NCPHOT
C
C       Wavelength of this photon
C       (the formula is basically 1/lambda**2):
              WLGDEN = CWLCON - RNUNIF(I)*CWLINT*CWLGT1
              WLGT   = CWLGT1*CWLCON/WLGDEN
C
C       Attenuation length:
C       Non-AMANDA case:
              IF (IDCONF .NE. 4) THEN
C          Constant attenuation length for the non-AMANDA cases:
                ATTLGT = ATTLEN
                GOTO 11
              ENDIF
C       AMANDA case: Attenuation length of the Cherenkov light at this
C       wavelength.
              IF (WLGT .LT. WLMINA) THEN
C          Poor extrapolation outside measured wavelength range - ought
C          to be improved!
                ATTLGT = ATTVEC(1)
              ELSEIF (WLGT .GT. WLMAXA) THEN
C          Poor extrapolation outside measured wavelength range - ought
C          to be improved!
                ATTLGT = ATTVEC(NWBINA)
              ELSE
                IABIN = INT((WLGT-WLMINA)/WLBINA) + 1
                ATTLGT = ATTVEC(IABIN)
              ENDIF
C
   11         CONTINUE
C       (Back to the general case.)
C
C       Attenuation of the Cherenkov light: does the photon reach the
C       PM tube or not?
              THISLG = -ATTLGT*LOG(RNUNIF(WLGT))
              IF (THISLG .LT. PHPATH) THEN
C          No, it doesn't.
                GOTO 10
              ENDIF
C
C       Quantum efficiency of the PM tube at this wavelength:
              IF ((WLGT .LT. WLMINQ) .OR. (WLGT .GT. WLMAXQ)) THEN
                QUANEF = 0.
              ELSE
                IQBIN = INT((WLGT-WLMINQ)/WLBINQ) + 1
                QUANEF = QE(IQBIN)
              ENDIF
C
              IF (QUANEF .LE. 0.) THEN
                GOTO 10
              ENDIF
C
C       Now take into account the angular efficiency:
              EFFTOT = QUANEF*ANGEFF
C
C      (I don't think we need to do anything of the kind
C       "IF (EFFTOT .GT. 1.) EFFTOT = 1." here - it should
C       have no effect, as far as I can see.)
C
C       Is the PM tube efficient in this particular case?
              THISEF = RNUNIF(THISLG)
              IF (THISEF .GT. EFFTOT) THEN
C          No, it isn't.
                GOTO 10
              ENDIF
C
C       If we're here, then this photon does actually generate a
C       photo-electron:
              NPELEC = NPELEC + 1
C
   10       CONTINUE
C
C    That should be it!
C
            RETURN
            END
CDECK  ID>, NUSLIB. 
            SUBROUTINE NUSLIB(PVEC,PARMAS,X,Y,Z,TLENG,T0,OK)
C----------------------------------------------------------------------
C Spacial routine for reading the shower library.
C Shower library is generated in 'NUMAIN_ZLIBF.CRA' , subroutine:
C NUSHOU.
C_______________________________________________________________________
C
C
C-----------------------------------------------
C Memory requirements for HBOOK, Geant, and our
C private Zebra structure.
C-----------------------------------------------
C
      PARAMETER ( LHMEM  = 400000 )
      PARAMETER ( LGMEM  = 100000 )
      PARAMETER ( LZMEM  = 600000 )
      PARAMETER ( IFFSZE = 1000   )
      COMMON/PAWC/HHH(LHMEM)
      COMMON/GCBANK/Q(LGMEM)
      COMMON/NUBANK/FENCE(16),QSPACE(LZMEM)
      DIMENSION LQN(LZMEM),IQN(LZMEM-8)
      DIMENSION QN(LZMEM-8)
      EQUIVALENCE ( QSPACE(1), LQN(1) )
      EQUIVALENCE ( IQN(1), QN(1), LQN(9) )
      COMMON/CFREAD/FFSPAC(IFFSZE)
C-----------------------------------------------
      COMMON/NUZLIN/IXSTOR,IXDIV1,IXDIV2,NLSTRU,NLREF,LROOT
      COMMON/NULINK/JZZH0,JZEZH,JZTZH,JZCPA,JZCP1,JZPMH
     *      ,JZSDIR,JZSHLO,JZSHHI,JZLAST
      PARAMETER ( MXZBIN = 50 )
      PARAMETER ( MXTBIN = 40 )
      PARAMETER ( MXSHOW = 5000 )
      COMMON/NUSHLC/SHPHOT(0:MXZBIN,0:MXTBIN)
     *             ,SHT(0:MXZBIN,0:MXTBIN),SHT2(0:MXZBIN,0:MXTBIN)
     *             ,SHLATE(0:MXZBIN,0:MXTBIN)
     *             ,NBZ,NBTH,ZBINN(0:MXZBIN),THBINN(0:MXTBIN)
     *             ,EESHOW(MXSHOW),IDXSHO(MXSHOW)
     *             ,USESHL,NWBIN,INSLIB,NSH,CGEN,NZBMAX,NTBMAX
     *             ,QQCER2,QQSLIB,NPCER2,NPSLIB
     *             ,SHOFIL,CHDIR,CHNPAH,CHCPAH,CHNPAZ,CHCPAZ
      CHARACTER*48 SHOFIL
      CHARACTER*5 CHDIR
      CHARACTER*20 CHNPAH,CHCPAH,CHNPAZ,CHCPAZ
      LOGICAL USESHL,INSLIB,CGEN
      PARAMETER ( MAXQQ  = 5000.)
      PARAMETER ( MAXTIM = 100.)
      PARAMETER ( MAXPS  = 5000)
      COMMON /MYSHOC/ XPRTRA(MAXPS), YPRTRA(MAXPS), ZPRTRA(MAXPS),
     +                DXPTRA(MAXPS), DYPTRA(MAXPS), DZPTRA(MAXPS),
     +                PEMAS1,PTOT1(MAXPS),PTOT2,IPS, IPSTRA,
     +                ZSHMAX,ZDIST(MAXPS),ETOT1(MAXPS),TOTPHO,
     +                DEGTHE(MAXQQ),
     +                GARMS, ELES,NEUS,MUS,NPBS,OTHERS,PIO,PIS
      PARAMETER (MAXMEC=30)
      COMMON/GCTRAK/VECT(7),GETOT,GEKIN,VOUT(7),NMEC,LMEC(MAXMEC)
     + ,NAMEC(MAXMEC),NSTEP ,MAXNST,DESTEP,DESTEL,SAFETY,SLENG
     + ,STEP  ,SNEXT ,SFIELD,TOFG  ,GEKRAT,UPWGHT,IGNEXT,INWVOL
     + ,ISTOP ,IGAUTO,IEKBIN, ILOSL, IMULL,INGOTO,NLDOWN,NLEVIN
     + ,NLVSAV,ISTORY
C
C
            DIMENSION KEYSH(2)
            SAVE SHOLO2,SHOHI2
            DATA ICALL/0/

C
C **** Given a particle with energy EE, find the two nearest
C **** showers in the library, above and below. ( EE should be
C **** less than 1 Tev and greater than 0.5 Gev.
C **** Otherwise we call
C **** NUCER2 to measure the number of photons around the track.)
C **** Take the weighted average of the higher and lower library
C **** shower, and use it.
C
C           WRITE(6,*)'IN NUSLIB'
C
            IF ( ICALL .EQ. 0 ) THEN
              SHOLO2 = QN(JZSDIR+NWBIN+2*IDXSHO(1))**2
              SHOHI2 = QN(JZSDIR+NWBIN+2*IDXSHO(NSH))**2
            WRITE(6,*)'Highest, lowest shower = ',
     *          SQRT(SHOHI2),SQRT(SHOLO2)
C
              ICALL = 1
            ENDIF
            EE2 = VDOT(PVEC,PVEC,3)
            IF ( EE2 .LT. SHOLO2 ) THEN
              GOTO 997
            ELSEIF ( PARMAS  .LT. 0.00001 ) THEN
              GOTO 997
            ELSEIF ( PARMAS  .GT. 0.001 ) THEN
              GOTO 997
            ELSEIF ( EE2 .GT. SHOHI2 ) THEN
              GOTO 997
            ELSE
              CALL TIMED(TTT)
              EE = SQRT(EE2)
              IDXLO = -1
              IDXHI = -1
              EELO = -99999.9
              EEHI = +99999.9
CCC          NSH = IQN(JZSDIR+1) <- DEFINED IN COMMON
C          WRITE(6,*)'Trying to find good showers,'
C          WRITE(6,*)'among the ',NSH,' showers in the library,'
C          WRITE(6,*)'for a particle with energy ',EE
CCCC          DO 10 I=1,NSH
C            IF ( ICALL .EQ. 0 ) THEN
C              WRITE(6,*)'EE(',I,') = ',QN(JZSDIR+NWBIN+2*I)
C            ENDIF
C            IF ( QN(JZSDIR+NWBIN+2*I) .LT. EE ) THEN
C              IF ( QN(JZSDIR+NWBIN+2*I) .GT. EELO ) THEN
C                IDXLO = IQN(JZSDIR+NWBIN+2*I+1)
C                EELO = QN(JZSDIR+NWBIN+2*I)
C              ENDIF
C            ELSEIF ( QN(JZSDIR+NWBIN+2*I) .GT. EE ) THEN
C              IF ( QN(JZSDIR+NWBIN+2*I) .LT. EEHI ) THEN
C                IDXHI = IQN(JZSDIR+NWBIN+2*I+1)
C                EEHI = QN(JZSDIR+NWBIN+2*I)
C              ENDIF
C            ENDIF
CCCC   10     CONTINUE
C
C            Binary search:
C
              KK = LOCATF(EESHOW,NSH,EE)
              IF (( KK .EQ. 0 ).OR.(KK .LE. -NSH )) THEN
                IDXLO = -1
              ELSEIF ( KK .GT. 0 ) THEN
                IDXLO = IDXSHO(KK)
                IDXHI = IDXSHO(KK+1)
              ELSE
                IDXLO = IDXSHO(ABS(KK))
                IDXHI = IDXSHO(ABS(KK)+1)
              ENDIF
C
              ICALL = 1
C              WRITE(6,*)'IDXLO,IDXHI = ',IDXLO,IDXHI
              IF (( IDXLO .LE. 0 ) .OR. (IDXHI .LE. 0 )) THEN
                WRITE(6,*)'Cannot find two suitable showers!'
                GOTO 997
              ELSE
                INSLIB = .TRUE.
                CALL HFILL(4463,LOG(EE),0.,1.)
                CALL HFILL(4464,FLOAT(ABS(KK)),0.,1.)
                EELO = QN(JZSDIR+NWBIN+2*IDXLO)
                EEHI = QN(JZSDIR+NWBIN+2*IDXHI)
C              WRITE(6,*)'EEHI, EELO = ',EEHI,EELO
C
                IF ( EEHI .EQ. EELO ) THEN
                  FRAC = 1.0
                  FRACM1 = 0.0
                ELSE
                  FRAC = (EE-EELO)/(EEHI-EELO)
                  FRACM1 = 1.0-FRAC
                ENDIF
C
C       (** Primitive way of picking one shower or the other, **)
C       (** with probability according to FRAC                **)
C
C         write(6,*)'FRAC,FRACM1 = ',FRAC,FRACM1
                IF ( RNDM(FRAC) .LT. FRAC ) THEN
                  CALL HFILL(6679,FRAC,0.,1.)
                  FRAC = 1.0
                  FRACM1 = 0.0
                ELSE
                  CALL HFILL(6680,FRAC,0.,1.)
                  FRAC = 0.0
                  FRACM1 = 1.0
                ENDIF
C         write(6,*)'FRAC,FRACM1 = ',FRAC,FRACM1
C
C
C ****  RZNDIR & RZCDIR: define the naming directory and the current
C ****  working directory.'R'---> read. 'P'----> print
C ****  CHNPAH is the default name.
C
                CALL RZNDIR(CHNPAZ,' ')
                CALL RZCDIR(CHCPAZ,' ')
C
C **** Read a bank or a data structure ( CALL RZIN ):
C
                KEYSH(1) = IDXLO
                KEYSH(2) = 0
                CALL RZIN(IXSTOR,JZSHLO,2,KEYSH,9999,' ')
C
                KEYSH(1) = IDXHI
                KEYSH(2) = 0
                CALL RZIN(IXSTOR,JZSHHI,2,KEYSH,9999,' ')
C
C            WRITE(6,*)'JZSHLO,JZSHHI = ',JZSHLO,JZSHHI
C         WRITE(6,*)'NDHI,NDLO = ',IQN(JZSHLO-1),IQN(JZSHHI-1)
C
                CALL HFILL(4465,FRAC,0.,1.)
C
                CALL VZERO(SHPHOT,(MXZBIN+1)*(MXTBIN+1))
                NBZLO = IQN(JZSHLO+2)
                NBTLO = IQN(JZSHLO+3)
                NBZHI = IQN(JZSHHI+2)
                NBTHI = IQN(JZSHHI+3)
                NBZ = MAX(NBZLO,NBZHI)
                NBTH = MAX(NBTLO,NBTHI)
C              WRITE(6,*)'NBZLO,NBTLO = ',NBZLO,NBTLO
C              WRITE(6,*)'NBZHI,NBTHI = ',NBZHI,NBTHI
C              WRITE(6,*)'NBZ  ,NBTH   = ',NBZ,NBTH
C
                DO 15 IZ = 1,NBZLO
                  DO 20 IT = 1,NBTLO
C                  WRITE(6,*)'IZ,IT,index = ',IZ,IT,6+NBTLO*(IZ-1)+IT
C                  WRITE(6,*)'QN = ',QN(JZSHLO+6+NBTLO*(IZ-1)+IT)
                    SHPHOT(IZ,IT) = FRAC*QN(JZSHLO+6+NBTLO*(IZ-1)+IT)
   20             CONTINUE
   15           CONTINUE
CC            WRITE(6,*) 'each EE corisponding to the shower(LO):',
CC     &        SHPHOT
CC            WRITE(6,*) 'FRAC=', FRAC
                SHTOT = 0.0
                DO 30 IZ = 1,NBZHI
                  DO 40 IT = 1,NBTHI
                    SHPHOT(IZ,IT) = SHPHOT(IZ,IT) +
     +                     FRACM1*QN(JZSHHI+6+NBTHI*(IZ-1)+IT)
   40             CONTINUE
   30           CONTINUE
CC            WRITE(6,*) 'each EE corisponding to the shower(HI):',
CC     &        SHPHOT
C            WRITE(6,*) 'FRACM1=', FRACM1
                DO 42 IZ = 1,NBZ
                  DO 44 IT = 1,NBTH
                    SHTOT = SHTOT + SHPHOT(IZ,IT)
                    IF ( SHPHOT(IZ,IT) .LT. 0.0 ) THEN
                      write(6,*)'SHPHOT(',IZ,IT,') = ',SHPHOT(IZ,IT)
                      WRITE(6,*)'FRAC,FRACM1 = ',FRAC,FRACM1
                      WRITE(6,*)'SHLO,SHHI = '
     *                        ,QN(JZSHLO+6+NBTLO*(IZ-1)+IT)
     *                        ,QN(JZSHHI+6+NBTHI*(IZ-1)+IT)
                    ENDIF
   44             CONTINUE
   42           CONTINUE
CC        write(6,*) 'shphot(IZ,IT)=',shphot(IZ,IT)
C            WRITE(6,*)'SHTOT = ',SHTOT
                ZZMAX = ZBINN(NBZ)
                IF ( ZZMAX .LE. 0 ) THEN
C              WRITE(6,*)'ZZMAX = ',ZZMAX
                ELSE
                  PHPCM = SHTOT/ZZMAX
                  CALL HFILL(6676,PHPCM,0.,1.)
                ENDIF
C
C
C *** Take the shower, as stored in SHPHOT, and add to the PM tube hits.
C *** Then kill the particle.
C
                CALL NUCERL(PVEC,PARMAS,X,Y,Z,TLENG,T0,OK)
                ISTOP = 1
C
C ****  Drop a bank abd its depandents: CALL MZDROP
C
                CALL MZDROP(IXSTOR,JZSHLO,' ')
                CALL MZDROP(IXSTOR,JZSHHI,' ')
                CALL RZNDIR(CHNPAH,' ')
                CALL RZCDIR(CHCPAH,' ')
C
                CALL TIMED(TTT)
C            WRITE(6,*)'TTT = ',TTT
                CALL HFILL(4461,TTT,0.,1.)
                CALL HFILL(4462,EE,0.,1.)
C
              ENDIF
            ENDIF
C
            RETURN
C--------------------------
C
  997       CONTINUE
            IF ( CGEN ) THEN
              CALL NUCER2(PVEC,PARMAS,X,Y,Z,TLENG,T0,OK)
            ENDIF
C
            END
CDECK  ID>, NUCERL. 
            SUBROUTINE NUCERL(PVEC,PARMAS,X00,Y00,Z00,TLEN,T0,OK)
C----------------------------------------------------------------------
C Special routine for shower lib.
C
C Replace the current particle with a shower from the shower library,
C
C
C
C TLENG is the lenght of the track segment being considered (cm).
C_______________________________________________________________________
C
C

C
            INTEGER OK
            REAL PVEC(3)
            REAL AC
C
      COMMON/NUFLAG/FGEANT,KGEN,KFLUXE,KFLUXS,FWEI,
     *       IRNALG,IRNSAV,IEVTRS,IEVACT,NEVACT,
     *       KCERAL,FOPTIN,FDPST,FRSCAT,IVTXAL,FTISME,
     *       ZINVBS
      LOGICAL FGEANT,FWEI,FOPTIN,FDPST,FRSCAT,FTISME,ZINVBS
      PARAMETER ( MAXPM  = 250 )
      PARAMETER ( MXPMHI = 100*MAXPM )
      COMMON/NUCPM/PMINT(MAXPM),TFIRST(MAXPM),TAVE(MAXPM),NPMINT(MAXPM),
     +             NHITOT, IPMHIT(MXPMHI),
     +             TIMHIT(MXPMHI), CPHINT(MXPMHI), NPEINT(MXPMHI),
     +             PHPPMH(MXPMHI), ANGPMH(MXPMHI), BETPMH(MXPMHI),
     +             IPAHIT(MXPMHI)
      COMMON/NUPMPO/NNPM,X0PM(MAXPM),Y0PM(MAXPM),Z0PM(MAXPM),
     +              IORIPM(MAXPM), LAYSTR(MAXPM)
C     IORIPM specifies the orientation of the PM tubes:
C     IORIPM(i) = -1: PM tube number i is looking downward
C     IORIPM(i) = +1: PM tube number i is looking upward
C     LAYSTR(i) = the number of the layer or string that PM tube
C                 number i belongs to
      COMMON/NUGRID/ZSIZE,HGPOS,DISTX,DISTY,DISTZ,OFFLAY,OFFLX,OFFLY,
     *              OFFX,OFFY,NBINS,NLAY,NSIDE,NSTR,NPMSTR,IDCONF,
     *              ICONTY,PMDIAM,PMRAD,PMR2,PMAREA,XYFUZZ,
     *              VXPLSI, VXPLDI
      COMMON/NUMISC/START,  IPTYPE, PIN, EMU0, THEMU0, PHIMU0,
     +              ZTRMAX, IORIVX, IORITR, VOLD(7),
     +              KTMED,  REFIND, ATTLEN, CWLGT1, CWLGT2,
     +              CUTKC2, IPRHIS(3), TDIAM, INUTRA,
     +              ZCLAPP, DCLAPP
      INTEGER START
      PARAMETER ( MAXQQ  = 5000.)
      PARAMETER ( MAXTIM = 100.)
      PARAMETER ( MAXPS  = 5000)
      COMMON /MYSHOC/ XPRTRA(MAXPS), YPRTRA(MAXPS), ZPRTRA(MAXPS),
     +                DXPTRA(MAXPS), DYPTRA(MAXPS), DZPTRA(MAXPS),
     +                PEMAS1,PTOT1(MAXPS),PTOT2,IPS, IPSTRA,
     +                ZSHMAX,ZDIST(MAXPS),ETOT1(MAXPS),TOTPHO,
     +                DEGTHE(MAXQQ),
     +                GARMS, ELES,NEUS,MUS,NPBS,OTHERS,PIO,PIS
      PARAMETER ( MXZBIN = 50 )
      PARAMETER ( MXTBIN = 40 )
      PARAMETER ( MXSHOW = 5000 )
      COMMON/NUSHLC/SHPHOT(0:MXZBIN,0:MXTBIN)
     *             ,SHT(0:MXZBIN,0:MXTBIN),SHT2(0:MXZBIN,0:MXTBIN)
     *             ,SHLATE(0:MXZBIN,0:MXTBIN)
     *             ,NBZ,NBTH,ZBINN(0:MXZBIN),THBINN(0:MXTBIN)
     *             ,EESHOW(MXSHOW),IDXSHO(MXSHOW)
     *             ,USESHL,NWBIN,INSLIB,NSH,CGEN,NZBMAX,NTBMAX
     *             ,QQCER2,QQSLIB,NPCER2,NPSLIB
     *             ,SHOFIL,CHDIR,CHNPAH,CHCPAH,CHNPAZ,CHCPAZ
      CHARACTER*48 SHOFIL
      CHARACTER*5 CHDIR
      CHARACTER*20 CHNPAH,CHCPAH,CHNPAZ,CHCPAZ
      LOGICAL USESHL,INSLIB,CGEN
C
C Pick up GEANT constants - pi and c, for example:
      COMMON/GCONST/PI,TWOPI,PIBY2,DEGRAD,RADDEG,CLIGHT,BIG,EMASS
      COMMON/GCONSX/EMMU,PMASS,AVO
C
C
C Pick up some other GEANT variables:
      COMMON/GCKINE/IKINE,PKINE(10),ITRA,ISTAK,IVERT,IPART,ITRTYP
     +      ,NAPART(5),AMASS,CHARGE,TLIFE,VERT(3),PVERT(4),IPAOLD
C
      PARAMETER (MAXMEC=30)
      COMMON/GCTRAK/VECT(7),GETOT,GEKIN,VOUT(7),NMEC,LMEC(MAXMEC)
     + ,NAMEC(MAXMEC),NSTEP ,MAXNST,DESTEP,DESTEL,SAFETY,SLENG
     + ,STEP  ,SNEXT ,SFIELD,TOFG  ,GEKRAT,UPWGHT,IGNEXT,INWVOL
     + ,ISTOP ,IGAUTO,IEKBIN, ILOSL, IMULL,INGOTO,NLDOWN,NLEVIN
     + ,NLVSAV,ISTORY
C
C
C
C--- (X0PM,Y0PM,Z0PM) in /NUPMPO/ (in +KEEP,NUCPM) are the coordinates of
C--- each PMT in the Lab system.
C--- (X, Y, Z) are the origin of the new coordinates which are parallell
C--- moved to the A-th random point at the muon ray.
C--- (X1,Y1,Z1) are the coordinates of each PMT in the rotated system of
C--- the incoming muon.
C--- (X00,Y00,Z00) is a random point on the muon ray.
C--- PVEC(3) is the momentum vector of the muon; the unit is GeV/c.
C--- PARMAS  is the particle mass (GeV/c**2).
C--- ATTLEN  in /NUMISC/ is the attenuation length.
C--- TLEN is the lenght of the track segment being considered (cm).
C--- T0   is the starting time (s).
C
C
C Fine structure constant:
            PARAMETER (ALPHA = 1./137.036)
C
            DIMENSION X(MAXPM),   Y(MAXPM),   Z(MAXPM),
     +          X1(MAXPM),  Y1(MAXPM),  Z1(MAXPM),  D(MAXPM),
     +          X1SQ(MAXPM),Y1SQ(MAXPM),Z1SQ(MAXPM),DSQ(MAXPM)
            DIMENSION TPM(MAXPM), PMI(MAXPM), TPM1(MAXPM)
C
            DATA ICALL/0/
C
            SAVE
C
            IF (ICALL .EQ. 0) THEN
              BETAC  = 1./REFIND
              CINMED = CLIGHT/REFIND
              TWOPIA = TWOPI*ALPHA
              CWTERM = 1./CWLGT1 - 1./CWLGT2
              ICALL = 1
            ENDIF
C
C      WRITE(6,*) 'IN NUCERL, PVEC = ',PVEC
C      WRITE(6,*) '         , XYZ  = ',X00,Y00,Z00
C
            CALL VZERO(TPM, NNPM)
C
C---- Transform the coordinate.
C
C---- First translation:
C
            DO 60 I = 1,NNPM
              X(I) = X0PM(I) - X00
              Y(I) = Y0PM(I) - Y00
              Z(I) = Z0PM(I) - Z00
   60       CONTINUE
C
C---- Then rotation:
C
C----   The momentum of the muon is P, and its total energy is ETOT
C----   The Cherenkov angle is ANG = ACOS(1/BETA*REFIND)
C
            PSQUAR = PVEC(1)*PVEC(1) + PVEC(2)*PVEC(2) + PVEC(3)*PVEC(3)
            P      = SQRT(PSQUAR)
            IF (P .LT. CUTKC2) GOTO 999
C
            ETOT  = SQRT(PSQUAR + PARMAS*PARMAS)
            BETA  = P/ETOT
            IF (BETA .LE. BETAC) GOTO 999
C
            CXBETA = CLIGHT*BETA
C
            COSANG = 1.0/(BETA*REFIND)
            ANG    = ACOS(COSANG)
            SINANG = SIN(ANG)
            TANANG = TAN(ANG)

C
            PT = SQRT(PVEC(1)*PVEC(1) + PVEC(2)*PVEC(2))
            THETA = ATAN2(PT,PVEC(3))
C          WRITE(6,*) 'PVEC(3) SHOULD .GT. 0,IN NUCERL:',PVEC(3)
            IF (PT .GT. 1.0E-3) THEN
              ALFA = ATAN2(PVEC(2),PVEC(1))
            ELSE
              ALFA = 0.
            ENDIF
C
            COSTHE = COS(THETA)
            SINTHE = SIN(THETA)
            COSALF = COS(ALFA)
            SINALF = SIN(ALFA)
C
C      WRITE(6,*)'AFTER TRANSFORM'
C
C********************************************************************
C From now on, all PM tubes has new coordinates, in MUON's system:
C*********************************************************************
C
            DO 70 I = 1,NNPM
              X1(I) =  COSTHE*COSALF*X(I) +
     +           COSTHE*SINALF*Y(I) -
     -           SINTHE*Z(I)
              Y1(I) = -SINALF*X(I) + COSALF*Y(I)
              Z1(I) =  SINTHE*COSALF*X(I) +
     +           SINTHE*SINALF*Y(I) +
     +           COSTHE*Z(I)
              X1SQ(I) = X1(I)*X1(I)
              Y1SQ(I) = Y1(I)*Y1(I)
              Z1SQ(I) = Z1(I)*Z1(I)
C----   D(I) is the distance from each PMT to the muon ray
CCC        D(I) = SQRT(X1(I)*X1(I) + Y1(I)*Y1(I))
              DSQ(I) = X1SQ(I) + Y1SQ(I)
              D(I)   = SQRT(DSQ(I))
C
C            WRITE(6,*) ' NEW AND OLD COORDINATES OF PM TUBES:'
CCC  100   FORMAT(1X,'X1,Y1,Z1 = ',3F9.1,'  X0,Y0,Z0 = ',3F9.1)
CCC        WRITE(6,*) 'D(I)=       ',D(I)
   70       CONTINUE
C
C----   CALCULATE THE NUMBER OF PHOTONS RECEIVED BY EACH PMT
C
C
C
C ***** Loop over all PM tubes, and for each PM loop over
C ***** all bins in Z along the shower.
C
C      WRITE(6,*)'NNPM, NBZ = ',NNPM,NBZ
C
C
C            write(6,*) 'X1(I),Y1(I),Z1(I),X(I),Y(I),Z(I)'
C            WRITE(6,*) ,X1(I),Y1(I),Z1(I),X(I),Y(I),Z(I)
C
            DO 110 I = 1,NNPM
              PMI(I) = 0.0
              TPM(I) = 0.0
              TPM1(I) = 9999999.9
              RR = SQRT(X1(I)*X1(I)+Y1(I)*Y1(I))
              IF ( RR .GE. PMDIAM/2.) THEN
                RATIO = PMDIAM/(8.0*RR)
              ELSE
                RATIO = 0.25*PMDIAM*PMDIAM/(RR+0.5*PMDIAM)**2
              ENDIF
CC        RATIO = 1.0/(8.0*RR)
              DDOLD = 99999.9
              DDSUM = 0.0
              THSUM = 0.0
              QQQMAX = 0.0
              DD2 = SQRT(RR*RR+(Z1(I)-ZBINN(0))**2)
              DO 120 JZ = 1, NBZ
C
C ** For each Z bin, find the angular range covered by the tube.
C ** and calculate which fraction of the bin angular range this is.
C!!!!!!!!!!!!!!!!!!!
C Rough approximation, assuming PM far from shower,
C and much smaller than the bin size.
C!!!!!!!!!!!!!!!!!!!
C *** THC1(2) is the angle between Z's and the line which connecting
C *** ZBINN(JZ-1) and the centry of the PMT. the same as THC2( from
C *** ZBINN(JZ) point)
C ***DD1 and DD2 is the line which connecting ZBINN(JZ-1) [or zbinn(jz)]
C *** to the centry of PMT.  DD1 for this tube was DD2 of last tube.
C
                IF( (Z1(I)-ZBINN(JZ-1)) .LT. 1E-6) THEN
                  THC1 = PI/2.
                  THC2 = THC1
                ELSEIF ( RR .NE. 0.) THEN
                  THC1 = ATAN2(RR,Z1(I)-ZBINN(JZ-1))
                  THC2 = ATAN2(RR,Z1(I)-ZBINN(JZ))
                ENDIF
                DD1 = DD2
                DD2 = SQRT(RR*RR+(Z1(I)-ZBINN(JZ))**2)
                IF ( DD2 .GT. 15000.0 ) THEN
                  IF ( DD1 .GT. 15000.0 ) THEN
C              WRITE(6,*)'Too far from PM: DD1,DD2 = ',DD1,DD2
                    IF ( DD1+DD2 .GT. DDOLD ) THEN
C                WRITE(6,*)'and getting worse; skip this PM'
                      GOTO 110
                    ELSE
                      GOTO 120
                    ENDIF
                  ENDIF
                ENDIF
                DDOLD = DD1+DD2
C
C *** ANGSZ1(2) is the angle from JZ-1(JZ) to two edges of PMT.
C *** THL1(THH1) is the angle between Z's to one of PMT's edge,from
C *** point (JZ-1).
C
                ANGSZ1 = 2.0*ATAN2(0.5*PMDIAM,DD1)
                ANGSZ2 = 2.0*ATAN2(0.5*PMDIAM,DD2)
                THL1 = RADDEG*ABS(THC1-0.5*ANGSZ1)
                THH1 = RADDEG*ABS(THC1+0.5*ANGSZ1)
                THL2 = RADDEG*ABS(THC2-0.5*ANGSZ2)
                THH2 = RADDEG*ABS(THC2+0.5*ANGSZ2)
                IF ( THH1 .GT. 180. ) THEN
                  THH1 = 180.
                ENDIF
                IF ( THH2 .GT. 180. ) THEN
                  THH2 = 180.
                ENDIF
C
                NTL1 = NBINTH(THL1)
                NTH1 = NBINTH(THH1)
                NTL2 = NBINTH(THL2)
                NTH2 = NBINTH(THH2)
C
                IF (  ( SHPHOT(JZ,NTL1) .LT. 0.0 )
     *            .OR.( SHPHOT(JZ,NTH1) .LT. 0.0 )
     *            .OR.( SHPHOT(JZ,NTL2) .LT. 0.0 )
     *            .OR.( SHPHOT(JZ,NTH2) .LT. 0.0 )) THEN
       write(6,*) 'jz=', jz
        write(6,*)'thl1=', thl1
        write(6,*)'thl2=', thl2
        write(6,*)'thH1=', thH1
        write(6,*)'thh2=', thh2
        write(6,*)'NTL1=', NTL1
        write(6,*)'NTH1=', NTH1
        write(6,*)'NTL2=', NTL2
        write(6,*)'NTH2=', NTH2
            WRITE(6,*)'SHPHOT(NTL1) = ',SHPHOT(JZ,NTL1)
            WRITE(6,*)'SHPHOT(NTH1) = ',SHPHOT(JZ,NTH1)
            WRITE(6,*)'SHPHOT(NTL2) = ',SHPHOT(JZ,NTL2)
            WRITE(6,*)'SHPHOT(NTH2) = ',SHPHOT(JZ,NTH2)
          ENDIF

C
C
C ***** BHH is the slope of up-edge line of the triangle.(BLL is lower one)
C ***** ZHH is uper-point coordinate in Z's axis.(ZLL is a lower ONE)
C
                CALL HFILL(6678,FLOAT(NTH1-NTL1),0.,1.)
C
C If angle of THL1 .gt. 90 degree, there are no photons, we skip those cases:
C
                IF (THL1 .GT. 90.0 )  GOTO 110
C
                QQQ = 0.0
                IDUM = NBINZ(QQQ)
                ZWID = ZBINN(JZ)-ZBINN(JZ-1)
                JCASE = -1
C
C*** Find PM coverage in theta-zeta plane:
C*** How many theta bins are involved in this zeta bin between NTL1 and NTH2?
          IF ( NTL1 .EQ. NTH1 ) THEN
            IF ( NTL2 .EQ. NTH2 ) THEN
              IF ( NTL1 .EQ. NTL2 ) THEN
C
C***   All in one theta-zeta bin, no problem.
C
                JCASE = 1
                IF ( SHPHOT(JZ,NTL1) .EQ. 0.0 ) THEN
                  GOTO 120
                ENDIF
                THWID = THBINN(NTL1)-THBINN(NTL1-1)
C
C *** APM is the areas to be calculated:
C
                APM = ZWID*0.5*((THH1-THL1)+(THH2-THL2))
                FRAC = APM/(ZWID*THWID)
C
C *** Frac is the theta-zeta-fraction covered by the tube.
C *** The additional RATIO factor is for phi coverage
C *** and circular-not-square PM cross section.
C ***** QQQ is the number of photons in APM.
C
                QQQ = FRAC*RATIO*SHPHOT(JZ,NTL1)
              ELSEIF ( NTL2 .EQ. NTL1+1 ) THEN
C
C ***  Two theta-zeta bins, triangle calculations.
C ***** Case 2: THL1,THH1 are in the same bin (NTL1), and THL2,THH2 are in
C ***** next bin ( means: NTL2 = NTL1+1)
C
                JCASE = 2
                IF ( ( SHPHOT(JZ,NTL1) .LE. 0.0 )
     *          .AND.( SHPHOT(JZ,NTL2) .LE. 0.0 )) THEN
                   GOTO 120
                ENDIF
C
C ***** BHH is the slope of up-edge line of the triangle.(BLL is lower one)
C ***** ZHH is uper-point coordinate in Z's axis.ZLL is a lower coordinate.
C
                BHH = (ZBINN(JZ)-ZBINN(JZ-1))/(THL2-THL1)
                ZHH = ZBINN(JZ-1) + BHH*(THBINN(NTL1)-THL1)
                BLL = (ZBINN(JZ)-ZBINN(JZ-1))/(THH2-THH1)
                ZLL = ZBINN(JZ-1) + BLL*(THBINN(NTL1)-THH1)
                TRIHH = 0.5*(THBINN(NTL1)-THL1)*(ZHH-ZBINN(JZ-1))
                TRILL = 0.5*(THBINN(NTL1)-THH1)*(ZLL-ZBINN(JZ-1))
                APM = TRIHH-TRILL
                THWID = THBINN(NTL1)-THBINN(NTL1-1)
                FRAC = APM/(ZWID*THWID)
                QQQ = FRAC*RATIO*SHPHOT(JZ,NTL1)
                TRIHH = 0.5*(THL2-THBINN(NTL1))*(ZBINN(JZ)-ZHH)
                TRILL = 0.5*(THH2-THBINN(NTL1))*(ZBINN(JZ)-ZLL)
                APM = TRILL-TRIHH
                FRAC = APM/(ZWID*THBINN(NTL2)-THBINN(NTL2-1))
C                (*QQQ = QQQ + FRAC*RATIO*SHPHOT(JZ,NTL1+1)*)
                 QQQ = QQQ + FRAC*RATIO*SHPHOT(JZ,NTL2)
              ELSE
C ***  Several theta-zeta bins. Triangles at each end,
C ***  as above, and loop over the bins between.
C *****  Case 3: THL1,THH1 are in the same bin (NTL1), and THL2,THH2 are in
C *****  the same bin ( NTL2 ). There are some bins between NTL1 and NTL2.
C
                JCASE = 3
                BHH = (ZBINN(JZ)-ZBINN(JZ-1))/(THL2-THL1)
                ZHH = ZBINN(JZ-1) + BHH*(THBINN(NTL1)-THL1)
                BLL = (ZBINN(JZ)-ZBINN(JZ-1))/(THH2-THH1)
                ZLL = ZBINN(JZ-1) + BLL*(THBINN(NTL1)-THH1)
                TRIHH = 0.5*(THBINN(NTL1)-THL1)*(ZHH-ZBINN(JZ-1))
                TRILL = 0.5*(THBINN(NTL1)-THH1)*(ZLL-ZBINN(JZ-1))
                APM = TRIHH-TRILL
                THWID = THBINN(NTL1)-THBINN(NTL1-1)
                FRAC = APM/(ZWID*THWID)
                QQQ = FRAC*RATIO*SHPHOT(JZ,NTL1)
C          write(6,*) ' case 3,part 1:'
C	  write(6,*) 'frac',frac
C          write(6,*) ' ratio', ratio
C          write(6,*) ' shphot', shphot(jz,ntl1)
C	 WRITE(6,*)'CASE 3,QQQ1',QQQ
                DO 210 ITL = NTL1+1,NTL2-1
                  ZHH0 = ZHH
                  ZLL0 = ZLL
                  THWID = THBINN(ITL)-THBINN(ITL-1)
                  ZHH = ZHH0 + BHH*THWID
                  ZLL = ZLL0 + BLL*THWID
                  APM = 0.5*((ZHH-ZLL)+(ZHH0-ZLL0))*THWID
                  FRAC = APM/((ZBINN(JZ)-ZBINN(JZ-1))*THWID)
                  QQQ = QQQ + FRAC*RATIO*SHPHOT(JZ,ITL)
C                  QQQ = QQQ + FRAC*RATIO*SHPHOT(JZ,NTL2)
 210            CONTINUE
C          write(6,*) ' case 3,part 2:'
C	  write(6,*) 'frac',frac
C          write(6,*) ' ratio', ratio
C          write(6,*) ' shphot', shphot(jz,ntl1)
C	 WRITE(6,*)'CASE 3,QQQ2',QQQ

                TRIHH = 0.5*(THL2-THBINN(NTL2-1))
     *                  *(ZBINN(JZ)-ZHH)
                TRILL = 0.5*(THH2-THBINN(NTL2-1))
     *                  *(ZBINN(JZ)-ZLL)
                APM = TRILL-TRIHH
                FRAC = APM/((ZBINN(JZ)-ZBINN(JZ-1))
     *                 *(THBINN(NTL2)-THBINN(NTL2-1)))
                QQQ = QQQ + FRAC*RATIO*SHPHOT(JZ,NTL2)
C          write(6,*) ' case 3,part 3:'
C	  write(6,*) 'frac',frac
C          write(6,*) ' ratio', ratio
C          write(6,*) ' shphot', shphot(jz,ntl1)
C	 WRITE(6,*)'CASE 3,QQQ3',QQQ
              ENDIF
            ELSEIF ( NTL2 .EQ. NTL1 ) THEN
C
C ***  Two theta-zeta bins, complicated triangle calculations.
C *****  Case 4: THL1,THH1,THL2 in the same bin (NTL1), and THH2 in next
C *****  bin (NTH2 = NTL1+1).
C
              IF ( NTH2 .NE. NTL1 + 1 ) THEN
                WRITE(6,*)'CASE 4.5! NTL1,NTH2 = ',NTL1,NTH2
              ENDIF
              JCASE = 4
              BHH = ZWID/(THL2-THL1)
              ZHH = ZBINN(JZ-1) + BHH*(THBINN(NTL1)-THL1)
              THHH = ZWID/BHH + THL1
              BLL = ZWID/(THH2-THH1)
              ZLL = ZBINN(JZ-1) + BLL*(THBINN(NTL1)-THH1)
              TRIHH = 0.5*(THBINN(NTL1)-THL1)*(ZHH-ZBINN(JZ-1))
              TRIHHH = 0.5*(THBINN(NTL1)-THHH)*(ZHH-ZBINN(JZ))
              TRILL = 0.5*(THBINN(NTL1)-THH1)*(ZLL-ZBINN(JZ-1))
              APM = TRIHH-TRILL-TRIHHH
              THWID = THBINN(NTL1)-THBINN(NTL1-1)
              FRAC = APM/(ZWID*THWID)
              QQQ = FRAC*RATIO*SHPHOT(JZ,NTL1)
C              write(6,*) ' case 4,part 1:'
C	      write(6,*) 'frac',frac
C              write(6,*) ' shphot', shphot(jz,ntl1)
C	      WRITE(6,*)'CASE 4,QQQ1',QQQ
              TRILL = 0.5*(THH2-THBINN(NTL1))*(ZBINN(JZ)-ZLL)
              APM = TRILL
              FRAC = APM/((ZBINN(JZ)-ZBINN(JZ-1))
     *                *(THBINN(NTL1+1)-THBINN(NTL1)))
              QQQ = QQQ + FRAC*RATIO*SHPHOT(JZ,NTL1+1)
C              QQQ = QQQ + FRAC*RATIO*SHPHOT(JZ,NTL2)
C              write(6,*) ' case 4,part 2:'
C              write(6,*) 'frac',frac
C              write(6,*) ' ratio', ratio
C              write(6,*) ' shphot(NTLE+1)', shphot(jz,ntl1+1)
C              WRITE(6,*)'CASE 4,QQQ2',QQQ
            ELSE
C
C *****  Case 5: THL1,THH1 in the same bin (NTL1), THH2 in the bin(NTH2) which
C *****  next to THL2 [bin (NTH2=NTL2+1)], There are some bins between NTL1 and
C *****  THL2.
C
              WRITE(6,*)'CASE 5: L1,H1,L2,H2 = ',NTL1,NTH1,NTL2,NTH2
              JCASE = 5
              BHH = (ZBINN(JZ)-ZBINN(JZ-1))/(THL2-THL1)
              ZHH = ZBINN(JZ-1) + BHH*(THBINN(NTL1)-THL1)
              BLL = (ZBINN(JZ)-ZBINN(JZ-1))/(THH2-THH1)
              ZLL = ZBINN(JZ-1) + BLL*(THBINN(NTL1)-THH1)
              TRIHH = 0.5*(THBINN(NTL1)-THL1)*(ZHH-ZBINN(JZ-1))
              TRILL = 0.5*(THBINN(NTL1)-THH1)*(ZLL-ZBINN(JZ-1))
              APM = TRIHH-TRILL
              THWID = THBINN(NTL1)-THBINN(NTL1-1)
              FRAC = APM/(ZWID*THWID)
              QQQ = FRAC*RATIO*SHPHOT(JZ,NTL1)
C          write(6,*) ' case 5,part 1:'
C	  write(6,*) 'frac',frac
C          write(6,*) ' ratio', ratio
C          write(6,*) ' shphot', shphot(jz,ntl1)
C	 WRITE(6,*)'CASE 5,QQQ1',QQQ
              DO 220 ITL = NTL1+1,NTL2
                ZHH0 = ZHH
                ZLL0 = ZLL
                THWID = THBINN(ITL)-THBINN(ITL-1)
                ZHH = ZHH0 + BHH*THWID
                ZLL = ZLL0 + BLL*THWID
                APM = 0.5*((ZHH-ZLL)+(ZHH0-ZLL0))*THWID
                IF ( ZHH .GT. ZBINN(JZ) ) THEN
                  APM = APM-0.5*(ZHH-ZBINN(JZ))*
     *                 (THBINN(ITL)-THL2)
                ENDIF
                FRAC = APM/(ZWID*THWID)
                QQQ = QQQ + FRAC*RATIO*SHPHOT(JZ,ITL)
C          write(6,*) ' case 5,part 2:'
C	  write(6,*) 'frac',frac
C          write(6,*) ' ratio', ratio
C          write(6,*) ' shphot', shphot(jz,ntl1)
C	 WRITE(6,*)'CASE 5,QQQ2',QQQ
 220          CONTINUE
              TRILL = 0.5*(THH2-THBINN(NTL2))*(ZBINN(JZ)-ZLL)
              APM = TRILL
              FRAC = APM/((ZBINN(JZ)-ZBINN(JZ-1))
     *               *(THBINN(NTL2+1)-THBINN(NTL2)))
C              QQQ = QQQ + FRAC*RATIO*SHPHOT(JZ,NTL2+1)
              QQQ = QQQ + FRAC*RATIO*SHPHOT(JZ,NTL2)
C          write(6,*) ' case 5,part 3:'
C	  write(6,*) 'frac',frac
C          write(6,*) ' ratio', ratio
C          write(6,*) ' shphot', shphot(jz,ntl1)
C	 WRITE(6,*)'CASE 5,QQQ3',QQQ
            ENDIF
          ELSEIF ( NTL2 .EQ. NTH2 ) THEN
C
C *****  Case 6: THL1 is in a different bin, but THH1, THL2 and THH2 in the
C *****  same bin:
C
            IF ( NTH1 .EQ. NTL2 ) THEN
              JCASE = 6
C             (* First the triangle in the THL1 bin: *)
              BHH = ZWID/(THL2-THL1)
              ZHH = ZBINN(JZ-1) + BHH*(THBINN(NTL1)-THL1)
              BLL = (ZBINN(JZ)-ZBINN(JZ-1))/(THH2-THH1)
              TRIHH = 0.5*(THBINN(NTL1)-THL1)*(ZHH-ZBINN(JZ-1))
              APM = TRIHH
              FRAC = APM/((ZBINN(JZ)-ZBINN(JZ-1))
     *               *(THBINN(NTL1)-THBINN(NTL1-1)))
              QQQ = FRAC*RATIO*SHPHOT(JZ,NTL1)
C             (* Intermediate bins, if any: *)
              IF ( NTH1-NTL1 .GT. 1 ) THEN
                DO 610 ITL = NTL1+1,NTH1-1
                  ZHHL = ZHH
                  ZHH  = ZBINN(JZ-1) + BHH*(THBINN(ITL)-THL1)
                  APM = 0.5*(ZHHL+ZHH)*(THBINN(ITL)-THBINN(ITL-1))
                  FRAC = APM/((ZBINN(JZ)-ZBINN(JZ-1))
     *               *(THBINN(ITL)-THBINN(ITL-1)))
                  QQQ = QQQ + FRAC*RATIO*SHPHOT(JZ,ITL)
 610            CONTINUE
              ENDIF
C             (* Then the more complicated shape in the last bin:*)
              TRIHH = 0.5*(THL2-THBINN(NTL2-1))*(ZBINN(JZ)-ZHH)
              APM = (0.5*(THH2+THH1)-THBINN(NTH2-1))*ZWID - TRIHH
              FRAC = APM/((ZBINN(JZ)-ZBINN(JZ-1))
     *               *(THBINN(NTL2)-THBINN(NTL2-1)))
              QQQ = FRAC*RATIO*SHPHOT(JZ,NTL2)
            ELSEIF ( NTH1 .EQ. NTL2-1 ) THEN
              JCASE = 8
C             (* Small triangle in first bin: *)
              BHH = ZWID/(THL2-THL1)
              ZHH = ZBINN(JZ-1) + BHH*(THBINN(NTL1)-THL1)
              BLL = (ZBINN(JZ)-ZBINN(JZ-1))/(THH2-THH1)
              TRIHH = 0.5*(THBINN(NTL1)-THL1)*(ZHH-ZBINN(JZ-1))
              APM = TRIHH
              FRAC = APM/((ZBINN(JZ)-ZBINN(JZ-1))
     *               *(THBINN(NTL1)-THBINN(NTL1-1)))
              QQQ = FRAC*RATIO*SHPHOT(JZ,NTL1)
C             (* THH1-bin: *)
              TRIHHA = TRIHH
              ZHH = ZBINN(JZ-1) + BHH*(THBINN(NTH1)-THL1)
              TRIHH = 0.5*(THBINN(NTH1)-THL1)*(ZHH-ZBINN(JZ-1))
              ZLL = ZBINN(JZ-1) + BLL*(THBINN(NTH1)-THH1)
              TRILL = 0.5*(THBINN(NTH1)-THH1)*(ZLL-ZBINN(JZ-1))
              APM = TRIHH - TRIHHA - TRILL
              FRAC = APM/((ZBINN(JZ)-ZBINN(JZ-1))
     *               *(THBINN(NTH1)-THBINN(NTH1-1)))
              QQQ = FRAC*RATIO*SHPHOT(JZ,NTH1)
C             (* Last bin, with NTL2 and NTH2 *)
              TRIHH = 0.5*(THL2-THBINN(NTL2-1))*(ZBINN(JZ)-ZHH)
              TRILL = 0.5*(THH2-THBINN(NTL2-1))*(ZBINN(JZ)-ZLL)
              APM = TRILL-TRIHH
              FRAC = APM/(ZWID*THBINN(NTL2)-THBINN(NTL2-1))
              QQQ = QQQ + FRAC*RATIO*SHPHOT(JZ,NTL2)
            ELSE
              WRITE(6,*)'CASE 6.5 : L1,H1,L2,H2 = ',NTL1,NTH1,NTL2,NTH2
            ENDIF
          ELSE
C
C ***** The last case, THL1,THH1 in different bins(next) and THL2,THH2 also
C ***** in different bins:( some bins between THL1 and THL2)
C
            IF ( NTL1 .EQ. NTL2 ) THEN
              IF ( NTH1 .EQ. NTH2 ) THEN
                JCASE = 7
C               (* First the trapezoid in the L1L2 bin *)
                APM = (THBINN(NTL1)-0.5*(THL2+THL1))*ZWID
                FRAC = APM/(ZWID*(THBINN(NTL1)-THBINN(NTL1-1)))
                QQQ = FRAC*RATIO*SHPHOT(JZ,NTL1)
C               (* Then any intermediate bins, WHOLE bins! *)
                IF ( NTH1-NTL1 .GT. 1 ) THEN
                  DO 710 ITL = NTL1+1,NTH1-1
                    QQQ = QQQ + RATIO*SHPHOT(JZ,ITL)
 710              CONTINUE
                ENDIF
C               (* Last H1H2 trapezoid *)
                APM = (0.5*(THH1+THH2)-THBINN(NTH1-1))*ZWID
                FRAC = APM/(ZWID*(THBINN(NTH1)-THBINN(NTH1-1)))
                QQQ = QQQ + FRAC*RATIO*SHPHOT(JZ,NTH1)
              ELSE
                WRITE(6,*)'JCASE = 7.1, L1,H1,L2,H2 = '
     *                           ,NTL1,NTH1,NTL2,NTH2
              ENDIF
            ELSEIF ( NTH1 .EQ. NTL2 ) THEN
              JCASE = 9
C             (* Small triangle in first bin: *)
              BHH = ZWID/(THL2-THL1)
              ZHH = ZBINN(JZ-1) + BHH*(THBINN(NTL1)-THL1)
              BLL = (ZBINN(JZ)-ZBINN(JZ-1))/(THH2-THH1)
              TRIHH = 0.5*(THBINN(NTL1)-THL1)*(ZHH-ZBINN(JZ-1))
              APM = TRIHH
              FRAC = APM/((ZBINN(JZ)-ZBINN(JZ-1))
     *               *(THBINN(NTL1)-THBINN(NTL1-1)))
              QQQ = FRAC*RATIO*SHPHOT(JZ,NTL1)
C             (* NTH1+NTL2-bin *)
              TRIHH = 0.5*(THL2-THBINN(NTL2-1))*(ZBINN(JZ)-ZHH)
              ZLL = ZBINN(JZ-1) + BLL*(THBINN(NTH1)-THH1)
              TRILL = 0.5*(THBINN(NTH1)-THH1)*(ZLL-ZBINN(JZ-1))
              BNAREA = (ZBINN(JZ)-ZBINN(JZ-1))
     *                *(THBINN(NTH1)*(THBINN(NTH1-1)))
              APM = BNAREA - TRIHH - TRILL
              FRAC = APM/BNAREA
              QQQ = QQQ + FRAC*RATIO*SHPHOT(JZ,NTH1)
C             (* Last bin, with NTH2: *)
              TRIHH = 0.5*(THH2-THBINN(NTH2-1))*(ZBINN(JZ)-ZLL)
              APM = TRIHH
              FRAC = APM/((ZBINN(JZ)-ZBINN(JZ-1))
     *               *(THBINN(NTH2)-THBINN(NTH2-1)))
              QQQ = FRAC*RATIO*SHPHOT(JZ,NTH2)
            ELSEIF ( NTH1 .LT. NTL2 ) THEN
              JCASE = 10
C             (* Small triangle in first bin: *)
              BHH = ZWID/(THL2-THL1)
              ZHH = ZBINN(JZ-1) + BHH*(THBINN(NTL1)-THL1)
              BLL = (ZBINN(JZ)-ZBINN(JZ-1))/(THH2-THH1)
              TRIHH = 0.5*(THBINN(NTL1)-THL1)*(ZHH-ZBINN(JZ-1))
              APM = TRIHH
              FRAC = APM/((ZBINN(JZ)-ZBINN(JZ-1))
     *               *(THBINN(NTL1)-THBINN(NTL1-1)))
              QQQ = FRAC*RATIO*SHPHOT(JZ,NTL1)
C             (* Bins between NTL1 and NTH1 *)
              DO 102 ITL=NTL1+1,NTH1-1
                ZHHOLD = ZHH
                ZHH = ZBINN(JZ-1) + BHH*(THBINN(ITL)-THL1)
                APM = (0.5*(ZHH+ZHHOLD)-ZBINN(JZ-1))
     *                      *(THBINN(ITL)-THBINN(ITL-1))
                FRAC = APM/((ZBINN(JZ)-ZBINN(JZ-1))
     *               *(THBINN(ITL)-THBINN(ITL-1)))
                QQQ = QQQ + FRAC*RATIO*SHPHOT(JZ,ITL)
 102          CONTINUE
C             (* NTH1 bin; as above except lower right triangle *)
              ZHHOLD = ZHH
              ZHH = ZBINN(JZ-1) + BHH*(THBINN(NTH1)-THL1)
              ATRAP = (0.5*(ZHH+ZHHOLD)-ZBINN(JZ-1))
     *               *(THBINN(NTH1)-THBINN(NTH1-1))
              ZLL = ZBINN(JZ-1) + BLL*(THBINN(NTH1)-THH1)
              TRILL = 0.5*(ZLL-ZBINN(JZ-1))*(THBINN(NTH1)-THH1)
              APM = ATRAP - TRILL
              FRAC = APM/((ZBINN(JZ)-ZBINN(JZ-1))
     *               *(THBINN(NTH1)-THBINN(NTH1-1)))
              QQQ = QQQ + FRAC*RATIO*SHPHOT(JZ,NTH1)
C             (* Bins between NTH1 and NTL2 *)
              DO 103 ITL = NTH1+1,NTL2-1
                ZHHOLD = ZHH
                ZLLOLD = ZLL
                ZHH = ZBINN(JZ-1) + BHH*(THBINN(ITL)-THL1)
                ZLL = ZBINN(JZ-1) + BLL*(THBINN(ITL)-THH1)
                APM = 0.5*((ZHH+ZHHOLD)-(ZLL+ZLLOLD))
     *                    *(THBINN(ITL)-THBINN(ITL-1))
                FRAC = APM/((ZBINN(JZ)-ZBINN(JZ-1))
     *               *(THBINN(ITL)-THBINN(ITL-1)))
                QQQ = QQQ + FRAC*RATIO*SHPHOT(JZ,I)
 103          CONTINUE
C             (* NTL2 bin *)
              ZLLOLD = ZLL
              ZLL = ZBINN(JZ-1) + BLL*(THBINN(NTL2)-THH1)
              TRIHH = 0.5*(ZBINN(JZ)-ZHH)*(THL2-THBINN(NTL2-1))
              ATRAP = (ZBINN(JZ)-0.5*(ZLL+ZLLOLD))
     *               *(THBINN(NTL2)-THBINN(NTL2-1))
              APM = ATRAP-TRIHH
              FRAC = APM/((ZBINN(JZ)-ZBINN(JZ-1))
     *               *(THBINN(NTL2)-THBINN(NTL2-1)))
              QQQ = QQQ + FRAC*RATIO*SHPHOT(JZ,NTL2)
            ELSE
                WRITE(6,*)'JCASE = 7.2, L1,H1,L2,H2 = '
     *                           ,NTL1,NTH1,NTL2,NTH2
            ENDIF
          ENDIF
C
                CALL HFILL(6671,FLOAT(JCASE),0.,1.)
C
C                WRITE(6,*)'JCASE = ',JCASE
                CALL HFILL(6672,QQQ,0.,1.)
                CALL HFILL(6675,FLOAT(JCASE),QQQ,1.)
                CALL HFILL(6677,FLOAT(JCASE),DD1,1.)
                IF ( QQQ .GT. QQQMAX ) THEN
C                  WRITE(6,*)'QQQ,QQQMAX = ',QQQ,QQQMAX
                  QQQMAX = QQQ
                  JCMAX = JCASE
                ENDIF
C
                PMI(I) = PMI(I) + QQQ
                DDSUM = DDSUM + 0.5*(DD1+DD2)*QQQ
                THSUM = THSUM + 0.5*(THC1+THC2)*QQQ
                TT = T0 + 0.5*(ZBINN(JZ)+ZBINN(JZ-1))
     +            + 0.5*(DD1+DD2)*REFIND/CLIGHT
                TPM(I) = TPM(I) + PMI(I)*TT
                IF ( TT .LT. TPM1(I) ) THEN
                  TPM1(I) = TT
                ENDIF
C
  120         CONTINUE
C
              IF (PMI(I) .LE. 0.) GOTO 110
C
              QQSLIB = QQSLIB + PMI(I)
C
C--- AVERAGE PHOTON PATH
C
              PHPATH = DDSUM/PMI(I)
              CALL HFILL(6673,PHPATH,0.,PMI(I))
              RLAMBD = THSUM/PMI(I)
              CALL HFILL(6674,RADDEG*RLAMBD,0.,PMI(I))
C
C--- Average the time:
C
              TPM(I) = TPM(I)/PMI(I)
C
C       Now find a reasonable INTEGER number for the # Cherenkov
C       photons reaching the PM tube:
              CALL POISSN(PMI(I),NPHTMP,IERR)
C
              IF (NPHTMP .LE. 0) GOTO 110
C
C----   CALCULATE THE EFFICIENCY OF THE PHOTON IN EACH PMT
C
C       Calculate the angle between the photon path and the Z axis in lab.
C !!!!!!!!!!!!!!!!!!!!!
C  This assumes Cerenkov angle, and should be fixed!
C!!!!!!!!!!!!!!!!!!!!!!
C
              IF (PHPATH .GT. PMRAD) THEN
                COSLAM = COS(RLAMBD)
              ELSE
                COSLAM = 1.
                RLAMBD = 0.
              ENDIF
              IF (IORIPM(I) .GT. 0) RLAMBD = PI - RLAMBD
C
C       Now for the calculation of the number of photo-electrons...
C
              NPELEC = 0
C
C       Take into account the efficiency of the PM tube, including the
C       wavelength-dependent attenuation of the Cherenkov light, the
C       wavelength-dependent quantum efficiency, and the geometrical
C       acceptance of the tube:
              CALL NUOMEF(I,RLAMBD,PHPATH,BETA,NPHTMP,NPELEC)
C
              NPSLIB = NPSLIB + NPELEC
C
              IF ( NPELEC .GT. 100 ) THEN
                WRITE(6,*)'Too many photons!?'
                WRITE(6,*)'NPELEC,PMI,PHPATH = ',NPELEC,PMI(I),PHPATH
                WRITE(6,*)'Shower XYZ = ',X00,Y00,Z00
                WRITE(6,*)'PM-XYZ     = ',X1(I),Y1(I),Z1(I)
                WRITE(6,*)'QQQMAX,JCMAX = ',QQQMAX,JCMAX
              ENDIF
C
C      (The output from the above subroutine call is NPELEC.)
C
C     - Note: The jump goes to statement label 111, not 110.
              IF (NPELEC .LE. 0) GOTO 111
C
              IF (FTISME) THEN
C          Find the smearing of the time response of the PM tube to
C          the signal giving rise to these NPELEC photo-electrons:
                CALL NUTSME(I,NPELEC,TSMEAR)
                TPM(I) = TPM(I) + TSMEAR
              ENDIF
C
C       Calculate the 'average' arrival time of the Cherenkov photons:
              IF (TAVE(I) .EQ. 0.) THEN
                TAVE(I) = TPM(I)
              ELSE
                TAVE(I) = (TAVE(I)*FLOAT(NPMINT(I)) + TPM(I)
     &            *FLOAT(NPELEC))
     /            /FLOAT(NPMINT(I) + NPELEC)
              ENDIF
C
              IF (TPM1(I) .LT. TFIRST(I)) THEN
                TFIRST(I) = TPM1(I)
              ENDIF
C
              NPMINT(I) = NPMINT(I) + NPELEC
C
C     - Continue here, even if NPELEC = 0.
  111         CONTINUE
C
C     - Use this condition, instead of simply NPELEC = 0, for skipping
C       the 'hit by hit' stuff. The reason for doing this is that we
C       want to be able to play around with the efficiecies 'off-line'.
C       (Note: For 'hits' with 0 photo-electrons there will be no smea-
C       ring of the time response in the 'hit-by-hit' vector TIMHIT,
C       even if FTISME is set = TRUE.)
              IF ((PMI(I) .LE. 0.1) .AND. (NPELEC .LE. 0)) GOTO 110
C
C       Fill the 'hit by hit' vectors (note that zero skipping is
C       effective, because of the condition on PMI(I) above):
              IF (NHITOT .LT. MXPMHI) THEN
                NHITOT = NHITOT + 1
                IPMHIT(NHITOT) = I
                TIMHIT(NHITOT) = TPM(I)
                CPHINT(NHITOT) = PMI(I)
                NPEINT(NHITOT) = NPELEC
                PHPPMH(NHITOT) = PHPATH
                ANGPMH(NHITOT) = RLAMBD
                BETPMH(NHITOT) = BETA
                IPAHIT(NHITOT) = IPART
              ELSE
                WRITE(6,*) ' '
                WRITE(6,*) '*** NUCERL: More than ',MXPMHI,
     &            ' PM-tube hits!'
              ENDIF
C
  110       CONTINUE
C
  830       FORMAT(' '/' ','*** ERROR in NUCERL:'/' ','*** COSLAM = ',
     &        F7.3/
     +        ' ','*** (Adjusted to ',F4.1,')')
C
C
            RETURN
C
  999       CONTINUE
C      WRITE(6,*)'*** ERROR in NUCERL: *** Jump to 999 ***'
C
            RETURN
            END

CDECK  ID>, NBINZ.  
            INTEGER FUNCTION NBINZ(D)
C
      PARAMETER ( MXZBIN = 50 )
      PARAMETER ( MXTBIN = 40 )
      PARAMETER ( MXSHOW = 5000 )
      COMMON/NUSHLC/SHPHOT(0:MXZBIN,0:MXTBIN)
     *             ,SHT(0:MXZBIN,0:MXTBIN),SHT2(0:MXZBIN,0:MXTBIN)
     *             ,SHLATE(0:MXZBIN,0:MXTBIN)
     *             ,NBZ,NBTH,ZBINN(0:MXZBIN),THBINN(0:MXTBIN)
     *             ,EESHOW(MXSHOW),IDXSHO(MXSHOW)
     *             ,USESHL,NWBIN,INSLIB,NSH,CGEN,NZBMAX,NTBMAX
     *             ,QQCER2,QQSLIB,NPCER2,NPSLIB
     *             ,SHOFIL,CHDIR,CHNPAH,CHCPAH,CHNPAZ,CHCPAZ
      CHARACTER*48 SHOFIL
      CHARACTER*5 CHDIR
      CHARACTER*20 CHNPAH,CHCPAH,CHNPAZ,CHCPAZ
      LOGICAL USESHL,INSLIB,CGEN
C
            DATA ICALL/0/
C
            IF ( ICALL .EQ. 0 ) THEN
              ICALL = 1
              ZBINN(0) = 0.0
              ZBINN(1) = 50.0
              ZBINN(2) = 100.0
              ZBINN(3) = 150.0
              ZBINN(4) = 200.0
              ZBINN(5) = 300.0
              ZBINN(6) = 450.0
              ZBINN(7) = 600.0
              ZBINN(8) = 800.0
              ZBINN(9) = 1100.0
              ZBINN(10) = 1500.0
              ZBINN(11) = 2000.0
              ZBINN(12) = 99999.9
              NZBMAX = 11
            ENDIF
            NBINZ = NZBMAX + 1
C
            DO 10 IBINZ = 1,NZBMAX
              IF ( D .LT. ZBINN(IBINZ) ) THEN
                NBINZ = IBINZ
                RETURN
              ENDIF
   10       CONTINUE
C
            RETURN
            END
CDECK  ID>, NBINTH. 
            INTEGER FUNCTION NBINTH(THEG)
      PARAMETER ( MXZBIN = 50 )
      PARAMETER ( MXTBIN = 40 )
      PARAMETER ( MXSHOW = 5000 )
      COMMON/NUSHLC/SHPHOT(0:MXZBIN,0:MXTBIN)
     *             ,SHT(0:MXZBIN,0:MXTBIN),SHT2(0:MXZBIN,0:MXTBIN)
     *             ,SHLATE(0:MXZBIN,0:MXTBIN)
     *             ,NBZ,NBTH,ZBINN(0:MXZBIN),THBINN(0:MXTBIN)
     *             ,EESHOW(MXSHOW),IDXSHO(MXSHOW)
     *             ,USESHL,NWBIN,INSLIB,NSH,CGEN,NZBMAX,NTBMAX
     *             ,QQCER2,QQSLIB,NPCER2,NPSLIB
     *             ,SHOFIL,CHDIR,CHNPAH,CHCPAH,CHNPAZ,CHCPAZ
      CHARACTER*48 SHOFIL
      CHARACTER*5 CHDIR
      CHARACTER*20 CHNPAH,CHCPAH,CHNPAZ,CHCPAZ
      LOGICAL USESHL,INSLIB,CGEN
            REAL THEG
C
            IF ( ICALL .EQ. 0 ) THEN
              ICALL = 1
              THBINN(0) = 0.0
              THBINN(1) = 10.0
              THBINN(2) = 20.0
              THBINN(3) = 28.0
              THBINN(4) = 33.0
              THBINN(5) = 36.0
              THBINN(6) = 39.0
              THBINN(7) = 41.0
              THBINN(8) = 42.0
              THBINN(9) = 44.0
              THBINN(10) = 47.0
              THBINN(11) = 50.0
              THBINN(12) = 55.0
              THBINN(13) = 62.0
              THBINN(14) = 70.0
              THBINN(15) = 85.0
              THBINN(16) = 180.0
              NTBMAX = 15
            ENDIF
            NBINTH = NTBMAX + 1
C
            IF ( THEG .GT. 181.0 ) THEN
              WRITE(6,*)'ERROR IN NBINTH: THEG = ',THEG
            ENDIF
C
            DO 10 IBINTH = 1,NTBMAX
              IF ( THEG .LT. THBINN(IBINTH) ) THEN
                NBINTH = IBINTH
                RETURN
              ENDIF
   10       CONTINUE
C
            RETURN
            END
C+DECK,NUSLIG. *************************************************************
            SUBROUTINE NUSLIG(PVEC,PARMAS,X,Y,Z,TLENG,T0,OK)
C----------------------------------------------------------------------
C Special routine for shower lib.
C
C Start from the first Electron's track with GEANT
C Called from GUSTEP to generate the secondary Cerenkov light, with
C parametrization algorithm.
C
C
C Transforms the coordinates of one point of the first electron's
C track (we call it as a 'primary track') from the Lab to the IT's
C own track system ( parallel moving + rotation )
C
C
C TLENG is the lenght of the track segment being considered (cm).
C                                         AUTHOR  QIN SUN 931106
C_______________________________________________________________________
C
C
            INTEGER OK
            REAL PVEC(3)
            REAL AC
      COMMON/NUFLAG/FGEANT,KGEN,KFLUXE,KFLUXS,FWEI,
     *       IRNALG,IRNSAV,IEVTRS,IEVACT,NEVACT,
     *       KCERAL,FOPTIN,FDPST,FRSCAT,IVTXAL,FTISME,
     *       ZINVBS
      LOGICAL FGEANT,FWEI,FOPTIN,FDPST,FRSCAT,FTISME,ZINVBS
      PARAMETER ( MAXPM  = 250 )
      PARAMETER ( MXPMHI = 100*MAXPM )
      COMMON/NUCPM/PMINT(MAXPM),TFIRST(MAXPM),TAVE(MAXPM),NPMINT(MAXPM),
     +             NHITOT, IPMHIT(MXPMHI),
     +             TIMHIT(MXPMHI), CPHINT(MXPMHI), NPEINT(MXPMHI),
     +             PHPPMH(MXPMHI), ANGPMH(MXPMHI), BETPMH(MXPMHI),
     +             IPAHIT(MXPMHI)
      COMMON/NUPMPO/NNPM,X0PM(MAXPM),Y0PM(MAXPM),Z0PM(MAXPM),
     +              IORIPM(MAXPM), LAYSTR(MAXPM)
C     IORIPM specifies the orientation of the PM tubes:
C     IORIPM(i) = -1: PM tube number i is looking downward
C     IORIPM(i) = +1: PM tube number i is looking upward
C     LAYSTR(i) = the number of the layer or string that PM tube
C                 number i belongs to
      COMMON/NUGRID/ZSIZE,HGPOS,DISTX,DISTY,DISTZ,OFFLAY,OFFLX,OFFLY,
     *              OFFX,OFFY,NBINS,NLAY,NSIDE,NSTR,NPMSTR,IDCONF,
     *              ICONTY,PMDIAM,PMRAD,PMR2,PMAREA,XYFUZZ,
     *              VXPLSI, VXPLDI
      COMMON/NUMISC/START,  IPTYPE, PIN, EMU0, THEMU0, PHIMU0,
     +              ZTRMAX, IORIVX, IORITR, VOLD(7),
     +              KTMED,  REFIND, ATTLEN, CWLGT1, CWLGT2,
     +              CUTKC2, IPRHIS(3), TDIAM, INUTRA,
     +              ZCLAPP, DCLAPP
      INTEGER START
      PARAMETER ( MAXQQ  = 5000.)
      PARAMETER ( MAXTIM = 100.)
      PARAMETER ( MAXPS  = 5000)
      COMMON /MYSHOC/ XPRTRA(MAXPS), YPRTRA(MAXPS), ZPRTRA(MAXPS),
     +                DXPTRA(MAXPS), DYPTRA(MAXPS), DZPTRA(MAXPS),
     +                PEMAS1,PTOT1(MAXPS),PTOT2,IPS, IPSTRA,
     +                ZSHMAX,ZDIST(MAXPS),ETOT1(MAXPS),TOTPHO,
     +                DEGTHE(MAXQQ),
     +                GARMS, ELES,NEUS,MUS,NPBS,OTHERS,PIO,PIS
      PARAMETER ( MXZBIN = 50 )
      PARAMETER ( MXTBIN = 40 )
      PARAMETER ( MXSHOW = 5000 )
      COMMON/NUSHLC/SHPHOT(0:MXZBIN,0:MXTBIN)
     *             ,SHT(0:MXZBIN,0:MXTBIN),SHT2(0:MXZBIN,0:MXTBIN)
     *             ,SHLATE(0:MXZBIN,0:MXTBIN)
     *             ,NBZ,NBTH,ZBINN(0:MXZBIN),THBINN(0:MXTBIN)
     *             ,EESHOW(MXSHOW),IDXSHO(MXSHOW)
     *             ,USESHL,NWBIN,INSLIB,NSH,CGEN,NZBMAX,NTBMAX
     *             ,QQCER2,QQSLIB,NPCER2,NPSLIB
     *             ,SHOFIL,CHDIR,CHNPAH,CHCPAH,CHNPAZ,CHCPAZ
      CHARACTER*48 SHOFIL
      CHARACTER*5 CHDIR
      CHARACTER*20 CHNPAH,CHCPAH,CHNPAZ,CHCPAZ
      LOGICAL USESHL,INSLIB,CGEN
C
C Pick up GEANT constants - pi and c, for example:
      COMMON/GCONST/PI,TWOPI,PIBY2,DEGRAD,RADDEG,CLIGHT,BIG,EMASS
      COMMON/GCONSX/EMMU,PMASS,AVO
C
C
C Pick up some other GEANT variables:
      COMMON/GCKINE/IKINE,PKINE(10),ITRA,ISTAK,IVERT,IPART,ITRTYP
     +      ,NAPART(5),AMASS,CHARGE,TLIFE,VERT(3),PVERT(4),IPAOLD
C
      PARAMETER (MAXMEC=30)
      COMMON/GCTRAK/VECT(7),GETOT,GEKIN,VOUT(7),NMEC,LMEC(MAXMEC)
     + ,NAMEC(MAXMEC),NSTEP ,MAXNST,DESTEP,DESTEL,SAFETY,SLENG
     + ,STEP  ,SNEXT ,SFIELD,TOFG  ,GEKRAT,UPWGHT,IGNEXT,INWVOL
     + ,ISTOP ,IGAUTO,IEKBIN, ILOSL, IMULL,INGOTO,NLDOWN,NLEVIN
     + ,NLVSAV,ISTORY
C
C
C Fine structure constant:
            PARAMETER (ALPHA = 1./137.036)
            PARAMETER (AC = 1.0)
C Z1 is a length of the primary track, in centimeter.
            DIMENSION  COSTHE(MAXQQ),THE(MAXQQ),GX1(MAXQQ),
     +           GY1(MAXQQ),GZ1(MAXQQ),
     +           ACOPHI(MAXQQ),TOTQQ(MAXQQ),
     +           DEGPHI(MAXQQ),COSPHI(MAXQQ)
            INTEGER   NUMPHO,NUPHI
C

            DATA ICALL/0/
C
            SAVE

C
C                   IF (ICALL .GE. 0) THEN
            BETAC  = 1./REFIND
            CINMED = CLIGHT/REFIND
            TWOPIA = TWOPI*ALPHA
            CWTERM = 1./CWLGT1 - 1./CWLGT2

C                   ENDIF

C
C PTOT2  is total momentum of the CURRENT track
C Calculate  momentum components of primary and current track:
C Current track: total momentum: PTOT2. The monmentum component in X axis
C is as we calculated in the GUSTEP: PVEC(1)(-->Px) = VECT(4)*PTOT
C For the current  track:
            PVEX   = PVEC(1)
            PVEY   = PVEC(2)
            PVEZ   = PVEC(3)
            SPVEC1 = PVEX*PVEX
            SPVEC2 = PVEY*PVEY
            SPVEC3 = PVEZ*PVEZ
            SPTOT2 = SPVEC1 + SPVEC2 + SPVEC3
            PTOT2  = SQRT(SPTOT2)
            ETOT2  = SQRT(SPTOT2 + PARMAS*PARMAS)
            BETA2  = PTOT2/ETOT2
            IF ( PTOT2 .LT. 0. ) THEN
              WRITE(6,*) 'PTOT2 .LT. 0.00!!'
            ENDIF

            COSCDX = PVEX/PTOT2
            COSCDY = PVEY/PTOT2
            COSCDZ = PVEZ/PTOT2
            IF( COSCDZ .GT. 1.0) THEN
              COSCDZ = 0.999999
            ENDIF
C For each secondary point: ( as input values of X,Y,Z)
            RADIUS = SQRT(X*X+Y*Y)
C          CALL HFILL(6662,RADIUS,Z,1.)
C          CALL HFILL (6663,RADIUS,0.,1.)
C          CALL HFILL (6664,X,Z,1.)
C          call hfill (6665, ptot2,0.,1.)
C
C The primary track:
C PTOT1 is total momentum of primary track
C DXPTRA, DYPTRA, DZPTRA are DIRECTION COSINE components of the primary track
C
C        COSAEX = DXPTRA(IPSTRA)
C        COSAEY = DYPTRA(IPSTRA)
C        COSAEZ = DZPTRA(IPSTRA)

            COSAEX = DXPTRA(1)
            COSAEY = DYPTRA(1)
            COSAEZ = DZPTRA(1)

            ETOT1(1) = SQRT(PTOT1(1)*PTOT1(1)
     +                          + PEMAS1*PEMAS1)
            BETA1  = PTOT1(1)/ETOT1(1)

C
C----   The Cherenkov angle : (from the current tracks)
            COSANG = 1.0/(BETA2*REFIND)
            IF( ABS(COSANG ).GE. 1.) THEN
              COSANG = SIGN(1.,COSANG)
            ENDIF

            ANG = ACOS(COSANG)
C transforms CHERNK to degree:

            DECHNK = RADDEG*ANG
C      WRITE(6,*) 'COSANG,DECHNK=',COSANG,DECHNK
            IF (BETA1 .LE. BETAC) GOTO 999
            IF (BETA2 .LE. BETAC) GOTO 999
C
C    Number of Cherenkov photons, with wavelength between CWLGT1
C    and CWLGT2, emitted per centimetre along the current track
C    (the factor 1.E7 is a conversion from nm to cm):
C   (Only one cone is calculated:)

            REFBET     = REFIND*BETA2
            BNTERM     = (1. - 1./REFBET)*(1.+ 1./REFBET)
            QQ         = TWOPIA*BNTERM*CWTERM
            QQ         = QQ*1.E7
            NUMPHO     = INT(QQ)
            CHPHEV     = CHPHEV + QQ*TLENG
C        CALL HFILL(5554,CHPHEV,0.,1.)
            IF (QQ .LT. 0.) THEN
              WRITE(6,*) ' '
              WRITE(6,*) 'NUSLIB: Negative # photos per cm: QQ =',
     +                        QQ
              WRITE(6,*) '*** (Setting QQ = 0.)'
              QQ = 0.0
              NUMPHO = 0

            ENDIF
C

            IF (NUMPHO .GT. MAXQQ) THEN
              WRITE(6,*) ' '
              WRITE(6,*) '*** NUSLIB:'
              WRITE(6,*) '*** NUMPHO > MAXQQ: ',NUMPHO,' ,',MAXQQ
              WRITE(6,*) '*** SETTING NUMPHO = MAXQQ'
              NUMPHO = MAXQQ
            ENDIF

C Define direction cosine of new coordinates:( in the current track system)
C NEW X1's axis presented in old coordinate: We just consider the positive
C consines --- chose the simplest way for calculate the angular distribution
C in the positive Z's axis.
            IF( COSCDZ .GT. 0.99999) THEN
              CX1 = 0.
              CY1 = 0.
              CZ1 = 1.
            ELSE
              AA     = COSCDY*COSAEZ - COSCDZ*COSAEY
              BB     = COSCDZ*COSAEX - COSCDX*COSAEZ
              CC     = COSCDX*COSAEY - COSCDY*COSAEX
              SQRABC = SQRT(AA*AA + BB*BB + CC*CC)
C If the angle between primary and the secondary tracks is CLOSE to eachother,
C then the new X axis is defined as the cross of the Z's axis in Lab.system and
C the primary track is:
              IF (SQRABC .LE. 1.0E-6) THEN
                SQRTXY = SQRT(COSAEY*COSAEY+COSAEX*COSAEX)
                IF (SQRTXY .LE. 1.0E-6) THEN
                  COSBEX = 1.
                  COSBEY = 0.
                  COSBEZ = 0.
                ELSE
C              WRITE(6,*) 'NUSLIG:SQRABC .LE. 0:'
C              WRITE(6,*) 'COSAEX,COSAEY,COSAEZ=',COSAEX,COSAEY,COSAEZ

C              WRITE (6,*) ' SQRABC =', SQRABC
                  COSBEX =  COSAEY/SQRTXY
                  COSBEY = -COSAEX/SQRTXY
                  COSBEZ =  0.
                ENDIF
              ELSE
                COSBEX = AA/SQRABC
                COSBEY = BB/SQRABC
                COSBEZ = CC/SQRABC
              ENDIF

C NEW Y's axis presented in old coordinate, three components and direction
C cosins are:

              DD     = COSAEY*COSBEZ - COSAEZ*COSBEY
              EE     = COSAEZ*COSBEX - COSAEX*COSBEZ
              FF     = COSAEX*COSBEY - COSAEY*COSBEX
              SQRDEF = SQRT(DD*DD + EE*EE + FF*FF)

C        IF (SQRDEF .LE. 1.0E-6) THEN
C          WRITE(6,*) 'SQRDEF=', SQRDEF
C        ENDIF


              COSGAX = DD/SQRDEF
              COSGAY = EE/SQRDEF
              COSGAZ = FF/SQRDEF

C
C Calculate the position of point C, use the direction of the current
C track and selecte the length of AC is chosed as 1.0 centimeter.
C COSACY is the angle between AC and the new Y's axis is:
              COSACY = (COSCDX*COSGAX + COSCDY*COSGAY +
     +                         COSCDZ *COSGAZ)
              COSACZ = (COSCDX*COSAEX + COSCDY*COSAEY +
     +                         COSCDZ *COSAEZ)
              IF(COSACY .GE. 1.0) THEN
                COSACY = 1.
              ELSE
                CX1 = 0.
                CY1 = COSACY
                CZ1 = COSACZ
              ENDIF
            ENDIF
C  Calculate the constants :
            NUM  = 100

            XNUM = NUM
            RRTAN   = TAN(ANG)*AC
            PCQQ = QQ/XNUM
C
C       (UNIT:NANO.SECOND= SECONDS*1.0E+9)
C
            DO 11 I = 1,NUM+1
              PHI = (I-1)*3.141592654/XNUM
              CPHI = COS(PHI)
              CPHI = COS(PHI)
              CTHE = (CZ1+CY1*RRTAN*CPHI)*COSANG
              IF (ABS(CTHE).GT. 1.0 ) THEN
                DEGTHE(I) = 0.0
              ELSEIF (CTHE .EQ. 1.0) THEN
                DEGTHE(I) = 0.
              ELSEIF (CTHE .LE. -1.) THEN
                DEGTHE(I) = 180.
              ELSE
                THE(I) = ACOS(CTHE)
                DEGTHE(I)  = THE(I)*RADDEG
              ENDIF
C
C
C            CALL HFILL (5550,DEGTHE(I),0.,PCQQ)
C
              IF ( ABS((Z-ZPRTRA(1))-500.0) .LT. 50. ) THEN
                CALL HFILL(5553,DEGTHE(I),TOFG*1E9,PCQQ)
              ENDIF
              CALL HFILL(6661,DEGTHE(I),Z-ZPRTRA(1),PCQQ)
              CALL HFILL(6660,Z-ZPRTRA(1),0.,PCQQ)
   11       CONTINUE
C
CC        STOP
C        IF ( COSCDZ .LT. 0.1 ) THEN
C          CALL HFILL(5551,Z-ZPRTRA(1),0.,QQ*TLENG)
C        ELSE
C          TPROJ = TLENG*COSCDZ
C          IF ( TPROJ .LT. 5. ) THEN
C            CALL HFILL(5557,Z-ZPRTRA(1),0.,QQ*TLENG)
C          ELSE
C            NTLENG = NINT(TPROJ)
C            DO 73 JJ = 1,NTLENG
C              CALL HFILL(5557,Z-ZPRTRA(1)+FLOAT(JJ)
C     *                          ,0.,QQ/COSCDZ)
C   73       CONTINUE
C          ENDIF
C        ENDIF

C          The factpr 1.e9 is converte seconds to nano.seconds
C        STIME = STIME + (TLENG/CLIGHT*BETA2)*1.E9
            CALL HFILL (5552,Z-ZPRTRA(1),TOFG*1E9,1.)
C
C                     WRITE(6,*)'TOFG = ',TOFG


  999       CONTINUE
C
            RETURN
            END









