      PROGRAM TRAIN
 
      PARAMETER(MAXPAT=200)
      PARAMETER(MAXPAR=400)
 
 
      COMMON /JNDAT1/ MSTJN(20),PARJN(20),OIN(1000),OUT(1000)
      COMMON /PATS/ PAT(MAXPAT,MAXPAR),DESO(MAXPAT,5),NPAT,NPAR

      DIMENSION EREC(MAXPAT)
 
      CALL GETPAT
 
C...set number of layers:
      MSTJN(1)=4
C...set update frequency to number of training patterns
      MSTJN(2) = NPAT
C...set shape of threshold function to g(x) = x
      MSTJN(3) = 4
C...set number of nodes in input layer:
      MSTJN(10)=NPAR
C...set number of nodes in hidden layer:
      DO 45 J=1,MSTJN(1)-2
        MSTJN(10+J)=30
 45   CONTINUE
C...set number of nodes in output layer:
      MSTJN(9+MSTJN(1))=1
C...set learning rate
      PARJN(1) = 0.001
C...all other parameters and switches keeps their default value
C
C...initialize the net:
      CALL JNINIT
      open(unit=8,file='junk')
      call jndump(8)
      close(8)
C
      OPEN(UNIT=7,FILE='net4.weights',STATUS='OLD',ERR=38)
      CALL JNREAD(7)
      close(7)
 38   CONTINUE
C
      ERRMIN = 99999.9
      NCYCLE = 40000
C...loop over training cycles:
      DO 100 ICYCLE=1,NCYCLE
 
       NRIGHT = 0
       SUMSQR = 0.0
       SUMSQ2 = 0.0
C...loop over all patterns:
       DO 200 IPAT=1,NPAT
 
C...put pattern IPAT in OIN
        DO 300 I=1,MSTJN(10)
          OIN(I)=PAT(IPAT,I)
300     CONTINUE
 
C...put desired output pattern in OUT
        do 302 i=1,5
          OUT(I)=DESO(IPAT,I)
 302	CONTINUE
 
C...train the net
        CALL JNTRAL
 
C...check the performance of the net

CC        DEE = (OUT(1)-DESO(IPAT,1))/DESO(IPAT,1)
        DEE = (OUT(1)-DESO(IPAT,1))
        DESO(IPAT,5) = ABS(DEE)
        dee2= (EXP(OUT(1))-EXP(DESO(IPAT,1)))/EXP(DESO(IPAT,1))
        IF ( ABS(DEE) .LT. 0.2 ) THEN
          NRIGHT = NRIGHT+1
        ENDIF
        SUMSQR = SUMSQR + DEE*DEE
        SUMSQ2 = SUMSQ2 + DEE2*DEE2
        EREC(IPAT) = OUT(1)
        IF ( MOD(ICYCLE,100) .EQ. 1 ) THEN
          WRITE(6,227)OUT(1),DESO(IPAT,1),OUT(2),OUT(3),OUT(4)
     *           ,DEE,OUT(5)
CCC,OUT(1)/DESO(IPAT)
 227      FORMAT(1X,'OUT,DEE,E,E/E0 : ',7F7.2)
       ENDIF
       IF ( ICYCLE .EQ. NCYCLE ) THEN
          WRITE(6,227)EXP(OUT(1)),EXP(DESO(IPAT,1))
     *           ,EXP(OUT(1))/EXP(DESO(IPAT,1))-1.0
       ENDIF
 
200    CONTINUE
       IF ( SUMSQR .LT. ERRMIN ) THEN
         ERRMIN = SUMSQR
         WRITE(6,*)ICYCLE,'; ',NRIGHT,' CORRECT out of'
     *          ,NPAT,'; '
         WRITE(6,*)'Sigma(dlogE) = ',SQRT(SUMSQR/NPAT)
         WRITE(6,*)'Sigma(dE/E)  = ',SQRT(SUMSQ2/NPAT)
CCCC         CALL FLUSH(6)
         DO 220 IPAT = 1,NPAT
           WRITE(6,228)EXP(DESO(IPAT,1)+1.0),EXP(EREC(IPAT)+1.0)
 228       FORMAT(1X,2F10.3)
 220     CONTINUE
       ENDIF
       IF ( MOD(ICYCLE,10) .EQ. 1 ) THEN
C...write out success rate to keep track of the learning process
        WRITE(6,*)ICYCLE,'; ',NRIGHT,' correct out of'
     *          ,NPAT,'; Sigma(dE/E) = ',SQRT(SUMSQR/NPAT)
CCCCC        CALL FLUSH(6)
        IF ( MOD(ICYCLE,10000) .EQ. 1 ) THEN
          OPEN(8,FILE='net4.weights')
          call jndump(8)
          close(8)
        ENDIF
       ENDIF
       
100   CONTINUE
    
C...dump all information of the net
      CALL JNDUMP(6)
 
      END
C
C===========================================================================
C
      SUBROUTINE GETPAT
C
      PARAMETER(MAXPAT=200)
      PARAMETER(MAXPAR=400)
      COMMON /PATS/ PAT(MAXPAT,MAXPAR),DESO(MAXPAT,5),NPAT,NPAR
C
      SAVE S
      DIMENSION S(MAXPAR,MAXPAT),E(MAXPAT),SS(MAXPAR)
      DIMENSION NPM(4,9,9),ELAY(4),NPHLAY(4),IKLAY(4)
C       (* Initialize *)
        OPEN(10,FILE='rc_emu.nofit.dat',STATUS='OLD')
        READ(10,*)
        NPAT = 0
        NLAY = 4
        DO 10 I=1,MAXPAT
          NPAT = NPAT + 1
          KPAR = 0
          DO 20 J=1,NLAY
            READ(10,*)
            READ(10,911,END=11)((NPM(J,II,JJ),JJ=1,9),II=1,9)
 911        FORMAT(9I8)
            READ(10,*,END=11)E(NPAT),NPHOT,OFFX,OFFY,DIST,ATTLEN
            ELAY(J) = E(NPAT)
            DESO(NPAT,5) = 0.0
            IF ( E(NPAT) .GT. 5000.0 ) THEN
              IF ( J .EQ. 4 ) THEN
                WRITE(6,*)'E = ',E(NPAT),'; ACCEPTED'
              ENDIF
              DESO(NPAT,1) = LOG(0.001*E(NPAT))-1.0
              DESO(NPAT,2) = -1.0
              IF ( E(NPAT) .GT. 10000.0 ) THEN
                DESO(NPAT,3) =  1.0
                DESO(NPAT,4) = -1.0
              ELSE
                DESO(NPAT,3) = -1.0
                DESO(NPAT,4) =  1.0
              ENDIF
            ELSE
C
C             !!!!!!!!!!!
              IF ( J .EQ. 4 ) THEN
                WRITE(6,*)'E = ',E(NPAT),'; REJECTED'
                NPAT = NPAT-1
                GOTO 10
              ENDIF
C             !!!!!!!!!!!!
C
              DESO(NPAT,1) = -1.0
              DESO(NPAT,2) =  1.0
              DESO(NPAT,3) = -1.0
              DESO(NPAT,4) = -1.0
            ENDIF
            WRITE(6,*)'E,LOG(E) = ',E(NPAT),(DESO(NPAT,jj),JJ=1,4)
            NPHLAY(J) = NPHOT
  20      CONTINUE
          WRITE(6,*)ELAY
          DO 30 K=1,4
            NPHMAX = 0
            DO 35 M=1,4
              IF ( NPHLAY(M) .GT. NPHMAX ) THEN
                IKLAY(K) = M
                NPHMAX = NPHLAY(M)
              ENDIF
  35        CONTINUE
            NPHLAY(IKLAY(K)) = -NPHLAY(IKLAY(K))
  30      CONTINUE
          DO 40 K=1,4
            DO 26 KK = 1,4
              NSUM = 0
              DO 22 II= 5-KK,5+KK
                DO 24 JJ= 5-KK,5+KK
                  NSUM = NSUM + NPM(IKLAY(K),II,JJ)
  24            CONTINUE
  22          CONTINUE
              KPAR = KPAR+1
              PAT(NPAT,KPAR) = 0.0001*NSUM
  26        CONTINUE
            KPAR = KPAR+1
            PAT(NPAT,KPAR) = -0.0001*NPHLAY(IKLAY(K))
 40       CONTINUE
          NPAR = KPAR
  10    CONTINUE
  11    CONTINUE  
        NPAT = NPAT-1
        WRITE(6,*)'NPAT,NPAR = ',NPAT,NPAR
C
C        DO 50 I=1,NPAT
C          WRITE(6,*)(PAT(I,JJ),JJ=1,NPAR)
C 50     CONTINUE
C
      RETURN
      END



      SUBROUTINE MY_JNTRAL

C...JetNet subroutine TRaining ALgorithm

C...Uses a back-propagation algorithm to train the net

      PARAMETER(MAXV=2000,MAXM=50000,MAXI=1000,MAXO=1000)

      COMMON /JNDAT1/ MSTJN(20),PARJN(20),OIN(MAXI),OUT(MAXO)
      COMMON /JNDAT2/ TINV(10),IGFN(10),ETAL(10)
      COMMON /JNINT1/ O(MAXV),A(MAXV),D(MAXV),T(MAXV),DT(MAXV),
     &                W(MAXM),DW(MAXM),NSELF(MAXM)
      COMMON /JNINT2/ M(0:10),MV0(10),MM0(10),NG(10),NL,IPOTT,ER1,ER2


      MSTJN(7)=MSTJN(7)+1

      CALL JNFEED
      CALL JNDELT

      ERR=0.0

      DO 100 I=1,M(NL)
C       (*** Modified LSJ 910424 **)
        IF ( OUT(I) .GT. 1.0 ) THEN
          ERR=ERR + ((OUT(I)-O(JNINDX(NL,I,0)))/OUT(I))**2
        ELSE
          ERR=ERR + OUT(I)-O(JNINDX(NL,I,0))
        ENDIF
        OUT(I)=O(JNINDX(NL,I,0))
100   CONTINUE

      ERR=ERR/FLOAT(M(NL))
      PARJN(8)=ERR
      ER1=ER1+ERR
      ER2=ER2+ERR

      IF(MOD(MSTJN(7),MSTJN(2)).NE.0) RETURN

C...update only every MSTJN(2) calls

      PARJN(9)=ER1/FLOAT(MSTJN(2))
      ER1=0.0
      IF(MOD(MSTJN(7),MSTJN(2)*MSTJN(9)).EQ.0) THEN
       PARJN(10)=ER2/FLOAT(MSTJN(2)*MSTJN(9))
       ER2=0.0
      ENDIF

      IF(PARJN(5).GT.0.0) THEN

C...Include pruning factors

       DO 110 I=1,MM0(NL+1)
        DW(I)=DW(I)-2.0*PARJN(5)*PARJN(1)*W(I)/((1.0+W(I)**2)**2)
110    CONTINUE

       DO 120 I=1,MV0(NL+1)
        DT(I)=DT(I)-2.0*PARJN(5)*PARJN(1)*T(I)/((1.0+T(I)**2)**2)
120    CONTINUE

      ENDIF

      IF(MSTJN(5).EQ.0) THEN

C...Normal updating

      DO 200 I=1,MM0(NL+1)
       W(I)=W(I)+DW(I)*FLOAT(NSELF(I))
       DW(I)=DW(I)*PARJN(2)
200   CONTINUE

      DO 210 I=1,MV0(NL+1)
       T(I)=T(I)+DT(I)
       DT(I)=DT(I)*PARJN(2)
210   CONTINUE

      ELSE

C...Manhattan updating

      DO 300 I=1,MM0(NL+1)
       W(I)=W(I)+SIGN(PARJN(1),DW(I))*FLOAT(NSELF(I))
       DW(I)=DW(I)*PARJN(2)
300   CONTINUE

      DO 310 I=1,MV0(NL+1)
       T(I)=T(I)+SIGN(PARJN(1),DT(I))
       DT(I)=DT(I)*PARJN(2)
310   CONTINUE

      ENDIF

      RETURN
      END

C**** END OF JNTRAL ****************************************************
C***********************************************************************
