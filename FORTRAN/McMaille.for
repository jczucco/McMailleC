      PROGRAM McMAILLE
*
*     Version 4.00 parallelized with OpenMP for multi-core processors
*       but this version is slightly modified for monoprocessors
*
*     MAILLE in french = CELL in english 
*     Mc for Monte Carlo
*     Pronounce : MacMy
*
************************************************************************ 
*
*     A Monte Carlo and grid search code for indexing powder patterns
*
*     For more details see the documentation at 
*                   http://www.cristal.org/McMaille/ 
*             or    http://sdpd.univ-lemans.fr/McMaille/ 
*
*              by A. Le Bail - September 2002 for version 0.9
*                              as well as for versions 1.0, 2.0 and 3.0 
*                              October 2006 for version 4.00            
*                        alb@cristal.org 
*                        http://www.cristal.org/ 
*
*                        Résidence Cristal - Appt 213 
*                        2, rue de Gasperi 
*                        72100 Le Mans 
*                        FRANCE
*
*   Versions 0.9 : cubic only
*            1.0 : hexagonal/trigonal/rhombohedral, tetragonal, 
*                  orthorhombic added, plus .ckm and .prf files
*            2.0 : monoclinic and triclinic added in MC
*                  but not in grid search (too long)
*            3.0 : columnar peak shapes instead of Gaussian
*                  in versions 0.9-2.0
*                  no Le Bail fit contrarily to versions 0.9-2.0
*                  only fit by percentage of inclusion of the
*                  calculated column into the observed one
*            3.02: black box mode  
*            3.03: improved Monte Carlo
*            3.04: two-phases mode
*            4.00: automatisation improved : more chances to identify
*                   the correct cell in "black box" mode
*                  Identification of the Bravais lattice
*                  Parallelization by using OpenMP directives
*                  improving the speed with multicore processors
*                  speed x1.7 to 1.8 with "dual core" or "core duo"
*                  speed x3.6 expected with the quad core in 2007
*                  speed x79 expected with the 80-core in 2012...;-)
*                  
*                      
*
************************************************************************ 
*
*    Copyright (C) 2002-2006 Armel Le Bail 
*
* This program is free software; you can redistribute it and/or 
* modify it under the terms of the GNU General Public License 
* as published by the Free Software Foundation.
*
* This program is distributed in the hope that it will be useful, 
* but WITHOUT ANY WARRANTY; without even the implied warranty of 
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the 
* GNU General Public License for more details.
*
* You should have received a copy of the GNU General Public License 
* along with this program; if not, write to the Free Software 
* Foundation, Inc.,
* 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
*
*
*********************************************************************** 
*
*
*
C
      USE OMP_LIB
C
      CHARACTER*79 BUFFER
      CHARACTER*80 NAM
      CHARACTER*80 FILE
      CHARACTER*80 TEMPO
      CHARACTER*1 FEND,BL
      CHARACTER*11 MORE
      CHARACTER*12 INDXPROG
      CHARACTER*7 DATENOW
      CHARACTER*8 TIMENOW
      CHARACTER SELECT1*80,VERSION*8
      REAL*8 NTRIED,NTRIEDT,IGT,NTRIEDB,NCYCLES,NTIMELIM(6)
      real  time_begin,time_end,totaltime
      LOGICAL QEX
      LOGICAL(4) PRESSEDK / .FALSE. /
      PARAMETER(VERSION='4.00')
      PARAMETER(N_HKL=10000,N_DAT=100)
C
      DIMENSION NSYS(6),PSTARTB(6),DELTA(3),PSTART(3)
      DIMENSION CELPRE(6),CELOLD(6),W1(N_HKL),FM20(N_HKL),FF20(N_HKL)
      DIMENSION BPAR(6),CEL(6,N_HKL),RP(N_HKL),VGC(N_HKL),D(N_HKL)
      DIMENSION IFI(N_HKL),TEXT(20),LL(N_HKL),QO(N_HKL),IB(N_HKL)
      DIMENSION AFI(8),BB(8),KM(N_HKL),IHKL(3,N_HKL),TH3(N_HKL)
      DIMENSION PMI(6),PMA(6),NSOL(N_HKL),CNCALC(N_HKL),XFOM(N_HKL)
      DIMENSION DELTCT(3),ASTARTT(3),IMN(10),IM(10)
      DIMENSION HW(N_HKL),HW4(N_HKL),BBB(N_HKL),FCAL(N_HKL),HH(3)
      DIMENSION POS(16000),YOBS(16000),YCALC(16000)
      DIMENSION DUMP(N_HKL),SOMEGA(N_HKL),THETA(N_HKL),RMAX0(6)
      DIMENSION NHA(16000),NHB(16000),JHH(3,N_HKL),IREFS(N_HKL)
      DIMENSION ISYST(7),LLL(N_HKL),QL(6),RP2(N_HKL),KM2(N_HKL)
      DIMENSION KM3(100000),LL2(100000),ID1(100000),ID2(100000)
C
      COMMON/CAL/NHKL0,LHKL,NDAT,DMIN,SLABDA2,IHH(3,N_HKL),AL(3,3),PI,
     1CRI(N_HKL),DIFP(N_HKL),DIFM(N_HKL),TH2(N_HKL),FOBS(N_HKL),SUM_F,
     2NIND,W2(N_HKL),NMX,NDAT10
      COMMON/CAL2/IND(N_DAT,N_HKL)
C
C$OMP THREADPRIVATE(/cal/,/cal2/) 
C
C      CALL KMP_GET_STACKSIZE(ISIZE)
C      print *,ISIZE
C      ISIZE=200000000
C      CALL KMP_SET_STACKSIZE(ISIZE)
C      CALL KMP_GET_STACKSIZE(ISIZE)
C      print *,ISIZE
C
C  Search for the number of processors available
C
      IPROCS=1
C      IPROCS=OMP_GET_NUM_PROCS()
      PRINT *
      PRINT *,'Number of used processors :',IPROCS
      PRINT *
      PROCS=IPROCS  
      IF(PROCS.LT.2.)PROCS=1.
C
      PI=114.59156  ! 360./3.1415926
      N1=1
      N2=2
      X=0.
C
C Open files using name from command line and standard extensions.
C The OPEN statements may need to be changed for some computers.
C Subroutine SXNM gets the generic filename from the command line.
C If nothing is found, then the user is prompted for the filename.
C
      CALL MCMNAM(LN,NAM)
C	PRINT *,LN,NAM
      IF(NAM(LN-3:LN).EQ.'.exe')GOTO 334
      FILE=NAM(1:LN)
      LFILE=LN
      GO TO 335
C
334   PRINT 1
1     FORMAT('  Entry file (no extension) ??',$)
      READ 2,FILE
      WRITE(*,*)FILE
2     FORMAT(A20)
      LFILE=LEN(FILE)
      DO WHILE(FILE(LFILE:LFILE).EQ.' ')
      LFILE=LFILE-1
      ENDDO
C
335   LC0=19
      LCA=21
      TEMPO=FILE(1:LFILE)//'.inp'
      INQUIRE(FILE=TEMPO,EXIST=QEX)
      IF(.NOT.QEX) GO TO 3
      CALL FILEDEL(21,TEMPO)
3     CALL OPEN_WRITE1(21,TEMPO)
      TEMPO=FILE(1:LFILE)//'.dat'
      INQUIRE(FILE=TEMPO,EXIST=QEX)
      IF(QEX) GO TO 333
      PRINT *,'That file does not exist, try again...'
      GO TO 334
333   CALL OPEN_READ1(19,TEMPO)
4     READ(19,5,END=6)SELECT1
      IF(SELECT1(1:1).NE.'!')WRITE(21,5)SELECT1
      GO TO 4
5     FORMAT(A80)
6     CLOSE(LC0)
      CLOSE(LCA)
      TEMPO=FILE(1:LFILE)//'.inp'
      CALL OPEN_READ1(21,TEMPO)
      LPR=20
      TEMPO=FILE(1:LFILE)//'.imp'
      INQUIRE(FILE=TEMPO,EXIST=QEX)
      IF(.NOT.QEX) GO TO 7
      CALL FILEDEL(20,TEMPO)
7     CALL OPEN_WRITE1(20,TEMPO)
      WRITE(BUFFER,*)'McMaille version ',VERSION
      CALL PROGRESSVIEW (BUFFER)
      WRITE(BUFFER,*)'Data file : ',FILE(1:LFILE)
      CALL PROGRESSVIEW (BUFFER)
      WRITE (20,3000) VERSION,FILE(1:LFILE)
      WRITE(20,*)
      WRITE(20,*)'  Number of Processors :',IPROCS
      WRITE(20,*)
      CALL DATN(DATENOW,TIMENOW)
      CALL CPU_TIME(time_begin)
C
C.....READ PROBLEM IDENTIFICATION
C
      READ(LCA,8,ERR=3500)(TEXT(I),I=1,20)
8     FORMAT(20A4)
      WRITE(20,8)(TEXT(I),I=1,20)
C
C.....READ wavelength and type of calculation
C     SLABDA = wavelength
C     ZERO   = zeropoint to be added at thr beginning
C     NGRID  if = 1  : grid cell generation
C            if = 0  : Monte Carlo cell generation
C            if = 2  : both
C            if = 3  : black box Monte Carlo only...
C            if = -3 : black box Monte Carlo only, without triclinic...
C            if = 4  : black box Monte Carlo + grid search...
C
      READ(LCA,*,ERR=3501)SLABDA,ZERO,NGRID 
      IVERB=0
      IF(SLABDA.LT.0.)THEN
      IVERB=1
      SLABDA=-SLABDA
      ENDIF
      BB(2)=SLABDA
      BB(1)=0.  ! zeropoint after correction... = 0.
      AFI(2)=0. ! code for wavelength refinement
      AFI(1)=1. ! code for zeropoint refinement
      NOTRIC=0
      IF(NGRID.EQ.-3)THEN
      NOTRIC=1
      NGRID=3
      ENDIF
C
C.....READ codes for search in crystalline systems
C
C    NSYS(n)
C     n
C     1   Cubic
C     2   Hexagonal/trigonal/Rhombohedral
C     3   Tetragonal
C     4   Orthorhombic
C     5   Monoclinic
C     6   Triclinic
C
C     if NSYS(n)=0 : no search
C     if NSYS(n)=1 : search
C     if NSYS(2)=2 : search in rhombohedral
C
      NBLACK=0
      IF(NGRID.EQ.4)NBLACK=1
      IF(NGRID.EQ.4)NGRID=3
      IF(NGRID.EQ.3)THEN
      TEMPO=FILE(1:LFILE)//'-new.dat'
      INQUIRE(FILE=TEMPO,EXIST=QEX)
      IF(.NOT.QEX) GO TO 9325
      CALL FILEDEL(28,TEMPO)
9325  CALL OPEN_WRITE1(28,TEMPO)
      WRITE(28,8)(TEXT(I),I=1,20)
      WRITE(28,9660)
9660  FORMAT('! Wavelength, zeropoint, Ngrid')
      WRITE(28,9661)SLABDA,ZERO
9661  FORMAT(F9.6,2X,F7.4,' 0')
      WRITE(28,9662)
9662  FORMAT('! Codes for symmetry')
      WRITE(28,*)'1 0 0 0 0 0'
      DO 9300 I=1,6
9300  NSYS(I)=1
      IF(NOTRIC.EQ.1)NSYS(6)=0
      ELSE
      READ(LCA,*,ERR=3502)(NSYS(I),I=1,6)
      ENDIF
C
C.....Read the tolerated error on 2-theta values W
C          which is also the column width of the
C          columnar profile shape
C          and how many non-indexed reflections NIND
C
C     If W is given as negative, then individual W
C     values will be read later (triplets : 2-theta, I, Width) 
C          moreover, Width will be multiplied by -W
C          (use W = -1 for no change...)
C
      IF(NGRID.EQ.3)THEN
      NIND=3
      W=0.30*SLABDA/1.54056
      WRITE(28,9663)
9663  FORMAT('! W, Nind')
      WRITE(28,9664)W,NIND
9664  FORMAT(F5.3,I3)
      ELSE   
      READ(LCA,*,ERR=3503)W,NIND
      ENDIF
C
C.....Some printing
C
      WRITE(20,*)
      WRITE(20,9)SLABDA,ZERO
9     FORMAT(' Wavelength : ',F9.5,' Zeropoint : ',F8.4/)
      WRITE(20,*)
      WRITE(20,*)'====================================================
     1==========================='
      IF(NGRID.EQ.0)WRITE(20,*)'    Monte Carlo cell generation'
      IF(NGRID.EQ.1)WRITE(20,*)'        Grid cell generation'
      IF(NGRID.EQ.2)WRITE(20,*)'Both searches - Monte Carlo AND grid'
      IF(NGRID.EQ.3)WRITE(20,*)' -- Black box Monte Carlo search --'
      IF(NBLACK.EQ.1)WRITE(20,*)' Black box Monte Carlo + grid search '
      WRITE(20,*)'====================================================
     1==========================='
      WRITE(20,*)
      WRITE(20,10)W
10    FORMAT(' Width of the columnar profile shape, W  = ',F9.4) 
      WRITE(20,*)
      WRITE(20,331)NIND
331   FORMAT(' Max non-indexed reflections, NIND  = ',I4) 
C
C.....READ Min/Max parameters, volume, Rp
C...  RMI : If Rp < RMI, stop, should be the good cell
C...  RMAX : Keep a refined cell if Rp < Rmax
C...  RMAXREF : if Rp < RMAXREF, refine that cell by Monte Carlo
C
C
      IF(NGRID.EQ.3)THEN
      PMIN=2.
      PMAX=30.
      VMIN=8.
      VMAX=27000.
      RMI=0.02
      RMAX=0.15
      RMAXREF=0.40
      WRITE(28,9668)
9668  FORMAT('!Pmin, Pmax, Vmin, Vmax, Rmin, Rmax, Rmaxref')
      WRITE(28,*)' 2. 50. 8. 125000. 0.05 0.15 0.50'
      ELSE
      READ(LCA,*,ERR=3504)PMIN,PMAX,VMIN,VMAX,RMI,RMAX,RMAXREF
      ENDIF
      IF(PMIN.LT.0.)READ(LCA,*,ERR=3504)(PMI(I),PMA(I),I=1,6)
      WRITE(20,*)
      IF(PMIN.GT.0.)THEN
C	WRITE(20,*)' Min/Max a,b,c, V ',PMIN,PMAX,VMIN,VMAX
      DO 30 I=1,3
      J=I+3
      PMI(I)=PMIN
      PMA(I)=PMAX
      PMI(J)=60.
      PMA(J)=120.
30    CONTINUE
      ELSE
      WRITE(20,*)' Min/Max a cell parameter ',PMI(1),PMA(1)
      WRITE(20,*)' Min/Max b cell parameter ',PMI(2),PMA(2)
      WRITE(20,*)' Min/Max c cell parameter ',PMI(3),PMA(3)
      WRITE(20,*)' Min/Max alpha cell parameter ',PMI(4),PMA(4)
      WRITE(20,*)' Min/Max beta  cell parameter ',PMI(5),PMA(5)
      WRITE(20,*)' Min/Max gamma cell parameter ',PMI(6),PMA(6)
      ENDIF
C
C  NR is test for automatic Rmax decrease
C
      NR=0
      IF(RMAX.LT.0.)THEN
      NR=1
      RMAX=-RMAX
      ENDIF
C
C	WRITE(20,*)
C	WRITE(20,*)' Min/Max volumes ',VMIN,VMAX
      WRITE(20,*)
      WRITE(20,*)' Min/Max Rp, Rmaxref ',RMI,RMAX,RMAXREF
C
C.....According to NGRID, read either grid steps
C                              or Monte Carlo parameters 
C                              or both
      IF(NGRID.EQ.3)THEN
      SPAR=0.02
      SANG=0.05
C      WRITE(28,9665)
C9665  FORMAT('! Spar, Sang')
C      WRITE(28,*)'0.02 0.2'
      WRITE(20,*)
      WRITE(20,*)' Steps on (a,b,c) and angles ',SPAR,SANG
      WRITE(28,9666)
9666  FORMAT('! Ntests, Nruns')
      WRITE(28,*)'-100 20'
      WRITE(28,9667)
9667  FORMAT('!  2-theta   Intensity')
      ELSE
      IF(NGRID.EQ.0)GO TO 12
      IF(NGRID.EQ.1)GO TO 11
      IF(NGRID.EQ.2)GO TO 11
      WRITE(20,*)
      WRITE(20,*)' UNKNOWN NGRID PARAMETER : STOP'
      STOP
11    CONTINUE
C
C.....READ grid steps
C     SPAR = step on cell parameters
C     SANG = step on angles
C
      READ(LCA,*,ERR=3505)SPAR,SANG
      WRITE(20,*)
      WRITE(20,*)' Steps on (a,b,c) and angles ',SPAR,SANG
      IF(NGRID.EQ.2)GO TO 12
      GO TO 13
12    CONTINUE
C
C  Continue up to NTIMELIM tests
C  Save parameters if Rp < Rmax
C  Make NRUNS times those NTIMELIM tests
C  If NTIMELIM is negative, |NTIMELIM| will apply to cubic
C             and |NTIMELIM|*50. for tetragonal, hexagonal, 
C                 etc
C
      READ (LCA,*,ERR=3506)TIMLIM,NRUNS
      IF(TIMLIM.LT.0.)THEN
      TIMLIM=-TIMLIM
      WRITE(20,*)
      WRITE(20,*)' N of runs ',NRUNS
      NTIMELIM(1)=TIMLIM
      WRITE(20,*)' N of MC events in cubic        ',NTIMELIM(1)
      NTIMELIM(2)=NTIMELIM(1)*20.
      WRITE(20,*)' N of MC events in tetra/hexa   ',NTIMELIM(2)
      NTIMELIM(3)=NTIMELIM(2)
      NTIMELIM(4)=NTIMELIM(3)*20.
      WRITE(20,*)' N of MC events in orthorhombic ',NTIMELIM(4)
      NTIMELIM(5)=NTIMELIM(4)*20.
      WRITE(20,*)' N of MC events in monoclinic   ',NTIMELIM(5)
      NTIMELIM(6)=NTIMELIM(5)*20.
      WRITE(20,*)' N of MC events in triclinic    ',NTIMELIM(6)
      ELSE
      DO 8221 I=1,6
      NTIMELIM(I)=TIMLIM
8221  CONTINUE
      WRITE(20,*)
      WRITE(20,*)' N of MC events, N of runs ',TIMLIM,NRUNS
      ENDIF
C
13    CONTINUE
      ENDIF
C
C... Make a WARNING 
C
      IF(NGRID.EQ.3)THEN
      PRINT *
      PRINT *,'      This is the black box mode, can be long...'
      IF(NOTRIC.EQ.1)THEN
      PRINT *
      PRINT *,'               No triclinic search.'
      ENDIF
      ENDIF
      PRINT *
      PRINT *,'  To cancel and save, type K (capital letter) anytime'
      PRINT *
C
C.....And now, read couples of 2-theta and I values
C
      SUM_F=0.
      NDAT=1
14    IF(W.GT.0.)THEN
      READ(LCA,*,END=15)TH2(NDAT),FOBS(NDAT)
      ELSE
      READ(LCA,*,END=15)TH2(NDAT),FOBS(NDAT),W1(NDAT)
      ENDIF
      IF(TH2(NDAT).GE.180.)GO TO 3508
      NDAT=NDAT+1
      IF(NDAT.GT.100)WRITE(20,*)'Max data = 100 !'
      IF(NDAT.GT.100)STOP
      GO TO 14
15    CONTINUE
      NDAT=NDAT-1
C
C.... Verify if these are d values or 2-theta  
C
      IF(TH2(2).LT.TH2(1))THEN
      WRITE(20,*)
      WRITE(20,*)'    WARNING : DATA were given as d(A) values'
      WRITE(20,*)
      WRITE(*,*)
      WRITE(*,*)'    WARNING : DATA were given as d(A) values'
      WRITE(*,*)
      DO 1412 NDA=1,NDAT
1412  TH2(NDA)=ASIN(SLABDA/(2.*TH2(NDA)))*PI
      ENDIF
C
      DO 1413 NDA=1,NDAT
      IF(W.GT.0.)THEN
      W1(NDA)=W
      IF(NGRID.EQ.3)THEN
      WRITE(28,*)TH2(NDA),FOBS(NDA)
      ENDIF
      ELSE
      W1(NDA)=W1(NDA)*(-W)
      ENDIF
      IF(TH2(NDAT).GE.180.)GO TO 3508
C
C     Addition of the Zeropoint
C
      TH2(NDA)=TH2(NDA)+ZERO
      D(NDA)=SLABDA/(2.*SIN(TH2(NDA)/PI))
      QO(NDA)=1/D(NDA)**2
1413  CONTINUE
      IF(NGRID.EQ.3.AND.NDAT.GE.20)NDAT=20
      IF(NGRID.EQ.3)CLOSE(28)
      NHKL=NDAT
C
C... NDAT10 is the max limit for the number of calculated 
C           peak positions = 10 times the number of
C           observed peak positions          
C
      NDAT10=NDAT*10
      NDAT2=NDAT*2
C
      DO 9400 I=1,NHKL
9400  SUM_F=SUM_F+FOBS(I)
      NMAX=NDAT-NIND
C
C.....END OF DATA READING
C
      CLOSE (LCA,STATUS='DELETE')
C
C  Output of some data
C
      WRITE(20,*)
      WRITE(20,*)'    2-THETA     d(A)    Iobs       W'
      DO 21 NM=1,NHKL
      WRITE(20,22)TH2(NM),D(NM),FOBS(NM),W1(NM)
21    CONTINUE
22    FORMAT(2X,F10.3,F9.4,2F10.3)
      WRITE(20,*)
C
C...  Various starting values initialized
C
      SLABDA2=SLABDA**2/4.
      IF(PMIN.GT.0.)THEN
      DO 23 I=1,3
      DELTCT(I)=(120.-60.)/2.
      ASTARTT(I)=60.
      DELTA(I)=(PMAX-PMIN)/2.
23    PSTART(I)=PMIN
      DELTC=(120.-90.)/2.
      ASTART=90.
      ELSE
      DO 24 I=1,3
      J=I+3
      DELTCT(I)=(PMA(J)-PMI(J))/2.
      ASTARTT(I)=PMI(J)
      DELTA(I)=(PMA(I)-PMI(I))/2.
24    PSTART(I)=PMI(I)
      DELTC=(PMA(5)-PMI(5))/2.
      ASTART=PMI(5)
      ENDIF
      DO 1925 I=1,6
1925  RMAX0(I)=RMAX
      DO 25 J=1,NDAT
      W2(J)=W1(J)/2.
      DIFP(J)=TH2(J)+W2(J)
      DIFM(J)=TH2(J)-W2(J)
25    CONTINUE
C   DELTAB = zone explored in angstrom around a good cell
C   DMIN acts as a lower d limit for keeping reflections
      DELTAB=0.02
C      DMIN=SLABDA/(2.*SIN(TH2(NHKL)/PI))-DELTAB
      DMIN=D(NHKL)-DELTAB
C  Dmax values will help to determine max cell parameters
C    in black box mode
      DMAX1=D(1)+2.*DELTAB
      DMAX2=D(2)+2.*DELTAB
      DMAX3=D(3)+2.*DELTAB
C
C  Warning on the wavelength...
C
      WRITE(20,*)
      WRITE(20,*)' dmax = ',D(1)
      WRITE(20,*)' Be sure that your choice of max cell parameters'
      WRITE(20,*)'  ensures the exploration of all possibilities.'
      WRITE(20,*)
      WRITE(*,*)
      WRITE(*,*)' dmax = ',D(1)
      WRITE(*,*)' Be sure that your choice of max cell parameters'
      WRITE(*,*)'  ensures the exploration of all possibilities.'
      WRITE(*,*)
      IF(DMAX1.GT.30.)THEN
      WRITE(20,*)
      WRITE(20,*)' WARNING : dmax > 30 A.'
      WRITE(20,*)' Divide the wavelength by 10 and try again...'
      WRITE(20,*)' and then, multiply the cell parameters by 10.'
      WRITE(20,*)
      WRITE(*,*)
      WRITE(*,*)' WARNING : dmax > 30 A.'
      WRITE(*,*)' Divide the wavelength by 10 and try again...'
      WRITE(*,*)' and then, multiply the cell parameters by 10.'
      WRITE(*,*)
      ENDIF
C
      WRITE(20,*)
      WRITE(20,*)'   Dmin = ',DMIN
      DMIN=1/DMIN**2
C   DELTAD = zone explored in 1/100 of degrees around a good cell
      DELTAD=0.20
C  IGC = number of retained cells, IGM = max IGC
C  IREF = code for refining the best cell if it had
C         not Rp < Rmin
C  IGT = total number of cells satisfying to Rp < Rmax
C        including multiple identical cells
      IGC=0
      IGT=0.
      IREF=0
      IGM=10000
      IHR=0
      NRUNS2=1
      RPSMALL=1.
C
C  ESCAPE : a value of 0.15 means that in 15% of the tests,
C            a parameter change may be accepted even if that
C            change does not lead to any Rp or number of
C            indexed reflections improvement
C
      ESCAPE=0.15
C
      CALL DATN(DATENOW,TIMENOW)
      WRITE(20,*)
      WRITE(20,*)'====================================================
     1==========================='
      WRITE(20,*)'                RESULTS - RESULTS - RESULTS - RESULTS'
      WRITE(20,*)'====================================================
     1==========================='
      WRITE(20,*)'                EXPLORED CELL PARAMETERS AND VOLUMES:'
      WRITE(20,*)'====================================================
     1==========================='
      WRITE(20,*)
*------------------------------------------------------------------------- 
*     Initialisation
*
      CALL ESP_INIT(ISEED)
*
*-------------------------------------------------------------------------
C
C
C.....AND NOW : Generate the cell, either by Monte Carlo
C               or by grid search
C
      IF(NGRID.EQ.1)GO TO 700
C
C...  Cell generation by Monte Carlo
C
      PRINT *
      PRINT *,'Monte Carlo search :'
C
      IF(NSYS(1).EQ.0)GO TO 200
C
C    Cubic case
C
      PRINT *,'Cubic:        Rp     a        V     Nind'
C
      IFILE=1
      NCYCLES=200.
      CY=NCYCLES*1.1
      CELPRE(4)=90.
      CELPRE(5)=90.
      CELPRE(6)=90.
C
      IF(NGRID.EQ.3)THEN
      NRUNS=1
      PMIN=DMAX1*0.9
      PMAX=DMAX1*3.1
      VMIN=PMIN*PMIN*PMIN
      VMAX=PMAX*PMAX*PMAX
      NTIMELIM(1)=(VMAX-VMIN)*0.5
      IF(NTIMELIM(1).GT.10000.)NTIMELIM(1)=10000.
      DELTA(1)=(PMAX-PMIN)/2.
      PSTART(1)=PMIN
      ENDIF
C
      WRITE(20,*)
      WRITE(20,*)'Cubic Monte Carlo search :'
      WRITE(20,*)' Max a, V ',PMAX,VMAX
      WRITE(20,*)
C
      IF(IVERB.EQ.1)THEN
      WRITE(20,180)NRUN,NTIMELIM(1)
      WRITE(20,*)'====================================================
     1==========================='
      WRITE(20,*)' Rp  Trial number   a        V    Nind Icod'
      WRITE(20,*)
      ENDIF
180   FORMAT('  Results in cubic, run, tests :',I3,F12.0)
C
C     READ hkl Miller indices in cub.hkl
C
      TEMPO='cub.hkl'
      CALL OPEN_READ1(35,TEMPO)
      READ(35,*)NHKL0
      NHKL0=NDAT*6
      IF(NHKL0.GT.400)NHKL0=400
      DO 101 I=1,NHKL0
101   READ(35,*)(IHH(KK,I),KK=1,3)
      CLOSE(35)
C
      DO 199 NRUN=1,NRUNS
*------------------------------------------------------------------------- 
*     Initialisation
*
C      CALL ESP_INIT(ISEED)
*
*-------------------------------------------------------------------------
      RMAX=RMAXREF
      RMIN=RMAX
C
C...  here starts the loop
C
      INTEREST=0
      TMAX=NTIMELIM(1)/PROCS
      TTMAX=10.*NTIMELIM(1)
      NCELLS=NTIMELIM(1)
      IISEED=0
      NTRIED=0.
      NTRIEDT=0.
      NOUT=0
C
C$OMP PARALLEL DEFAULT(SHARED) COPYIN(/CAL/,/CAL2/)
C$OMP& PRIVATE(NCEL,NTRIEDB,DEL,V1,ICODE,LLHKL,IHKL,TH3,NCALC,
C$OMP& RMAX2,A,V2,BPAR,V3,PSTARTB,IPEN,ISEE,INDIC,DDT,DDQ,
C$OMP& DIFF,DIFF2)
C$OMP& FIRSTPRIVATE(iseed,iiseed,rmax0,ntried,ntriedt,nout,
C$OMP& celpre,rmin,rmax,bb,afi)
C$OMP DO
C
      DO 196 NCEL=1,NCELLS
      IF(NOUT.GE.1)GO TO 196
      IF(INTEREST.GE.1)GO TO 196
      IISEED=IISEED+1
      IF(IISEED.EQ.1)ISEED=((ISEED-NCEL*NRUN)/2)*2+1
C
102   CONTINUE
C
      NTRIEDB=0.
      CELPRE(1)=PSTART(1)+2.*DELTA(1)*RANDI(ISEED)
      NTRIED=NTRIED+1.
      GO TO 104
103   DEL=DELTAB*(1.-NTRIEDB/CY)
      CELPRE(1)=PSTARTB(1)+(DEL*(RANDI(ISEED)-0.5)*2.)
      NTRIEDB=NTRIEDB+1.
104   CONTINUE
      CELPRE(2)=CELPRE(1)
      CELPRE(3)=CELPRE(1)
      DO 105 I=1,3
      DO 105 J=1,3
105   AL(I,J)=0.0
      CALL DCELL(CELPRE,AL,V1)
      IF(NTRIED.GT.TMAX)THEN
      NOUT=NOUT+1
      GO TO 196
      ENDIF
      IF(NTRIEDB.NE.0.)GO TO 106
      IF(V1.GT.VMAX.OR.V1.LT.VMIN) THEN
      NTRIED=NTRIED-1.
      NTRIEDT=NTRIEDT+1.
      IF(NTRIEDT.GT.TTMAX)THEN
      NOUT=NOUT+1
      GO TO 196
      ENDIF
      GO TO 102
      ENDIF
C
106   CALL CALCUL1(DIFF,DIFF2)
      IF(NMX.GT.NDAT10)THEN
      NTRIED=NTRIED-1
      GO TO 102
      ENDIF
      IF(NTRIEDB.NE.0.)GO TO 114
C
C... Here are the 2 criteria for selecting a cell to be
C...      "refined" by Monte Carlo (NCYCLES events)  
C... Rp value satisfying (<0.5) ??? or enough hkl explained ???
C
      IF(LHKL.GE.NMAX)THEN
      RMAX=DIFF
      ICODE=2
      IF(DIFF.LE.RMAXREF)ICODE=1
      ELSE
      ICODE=1
      ENDIF
      IF(DIFF.GT.RMAX) GO TO 117
      IF(LHKL.LT.NMAX) GO TO 117
114   IF(DIFF.LE.RMAX)THEN
      LLHKL=LHKL
      RMAX=DIFF
      RMAX2=DIFF2
      A=CELPRE(1)
      V2=V1
      IF(DIFF.LT.RMIN)THEN
      RMIN=DIFF
      BPAR(1)=A
      V3=V1
      ENDIF
C
C... "Refine" that cell (by Monte Carlo too...)
C
      PSTARTB(1)=CELPRE(1)
      ENDIF
      IF(NTRIEDB.LE.NCYCLES)GO TO 103
      NTRIEDB=0.
      IF(RMAX.GE.RMAX0(1))GO TO 117
      IF(RMAX2.GE.0.15)GO TO 117
      IPEN=NDAT-LLHKL
      IF(IPEN.GT.NIND)GO TO 117
C
C$OMP CRITICAL(STORE1)
C
      IGC=IGC+1
C
C  Test if too much proposals, if yes decrease Rmax by 5%
C
      IGT=IGT+1.
      IF(NR.EQ.1)THEN
      IF(IGT.GT.50.)THEN
      IF((NTRIED/IGT).LT.100.)THEN
      IF(RMAX0(1).GT.0.2)THEN
      RMAX0(1)=RMAX0(1)-RMAX0(1)*0.05
      WRITE(20,*)'  Rmax reduced by 5%, now Rmax = ',RMAX0(1)
      WRITE(*,*)'  Rmax reduced by 5%, now Rmax = ',RMAX0(1)
      ENDIF
      ENDIF
      ENDIF
      ENDIF
C
      IF(IGC.GT.10000)THEN
      WRITE(20,*)'   More than 10000 good cells = STOP'
      PRINT *,'   More than 10000 good cells = STOP'
      IGC=IGC-1
      INTEREST=INTEREST+1
      ENDIF
      CEL(1,IGC)=A
      CEL(2,IGC)=A
      CEL(3,IGC)=A
      CEL(4,IGC)=90.
      CEL(5,IGC)=90.
      CEL(6,IGC)=90.
C
C$OMP END CRITICAL(STORE1)
C
C... Check for supercell
C
      CELPRE(1)=A
      CELPRE(2)=A
      CELPRE(3)=A
      DO 140 I=1,3
      DO 140 J=1,3
140   AL(I,J)=0.0
      CALL DCELL(CELPRE,AL,V1)
C
C$OMP CRITICAL(STORE2)
C
      CALL CALCUL2(DIFF,IHKL,TH3,NCALC,IGC)
      KM(IGC)=LLHKL
      KM2(IGC)=LHKL
      IFI(IGC)=IFILE
      NSOL(IGC)=1
      VGC(IGC)=V1
      RP(IGC)=RMAX
      RP2(IGC)=DIFF
      IF(RP(IGC).LT.RPSMALL)THEN
      RPSMALL=RP(IGC)
      ISEE=1
      ELSE
      ISEE=0
      ENDIF
      CALL SUPCEL(LHKL,IHKL,CEL,IGC,VGC,3)
      CALL BRAV(LHKL,IHKL,IBR)
      IB(IGC)=IBR
      A=CEL(1,IGC)
      CEL(2,IGC)=A
      CEL(3,IGC)=A
      V2=VGC(IGC)
C
C$OMP END CRITICAL(STORE2)
C
C... Check for interesting result
C
C      IF(INTEREST.GE.1)GO TO 196
      INDIC=1
      BB(3)=A
      BB(4)=A
      BB(5)=A
      BB(6)=90.
      BB(7)=90.
      BB(8)=90.
      AFI(3)=1
      AFI(4)=1
      AFI(5)=1
      AFI(6)=0
      AFI(7)=0
      AFI(8)=0
      CELPRE(1)=A
      CELPRE(2)=A
      CELPRE(3)=A
      DO 110 I=1,3
      DO 110 J=1,3
110   AL(I,J)=0.0
C
C$OMP CRITICAL(FOUND)
C
      IF(RP(IGC).LT.RMI)THEN
      INTEREST=INTEREST+1
      WRITE(*,1115)RMAX,A,V2,IPEN
      WRITE(20,*)
      WRITE(20,*)'====================================================
     1==========================='
      WRITE(20,*)
      WRITE(20,*)' YOU HAVE FOUND AN INTERESTING RESULT : Rp < Rmin !'
      WRITE(20,*)
      WRITE(*,*)
      WRITE(*,*)' YOU HAVE FOUND AN INTERESTING RESULT : Rp < Rmin !'
      WRITE(*,*)
C
C... Refine that cell
C
      CALL DCELL(CELPRE,AL,V1)
      CALL CALCUL2(DIFF,IHKL,TH3,NCALC,IGC)
      CALL CELREF2(INDIC,BB,AFI,LHKL,TH3,IHKL,DDT,DDQ)
      CNCALC(IGC)=NCALC
      IF(NDAT.GE.20)THEN
      FM20(IGC)=QO(20)/(2.*CNCALC(IGC)*DDQ)
      FF20(IGC)=20./(CNCALC(IGC)*DDT)
      ELSE
      PNDAT=NDAT
      FM20(IGC)=QO(NDAT)/(2.*CNCALC(IGC)*DDQ)
      FF20(IGC)=PNDAT/(CNCALC(IGC)*DDT)
      ENDIF
C      WRITE(20,7000)FM20(IGC)
C      WRITE(20,7001)FF20(IGC),DDT,NCALC
C	WRITE(20,*)
C      PRINT 7000,FM20(IGC)
C      PRINT 7001,FF20(IGC),DDT,NCALC
C	PRINT *
      IREF=1
      GO TO 197
      ELSE
C
C  Anyway, calculate the M20 and F20 values
C
      CALL DCELL(CELPRE,AL,V1)
      CALL CALCUL2(DIFF,IHKL,TH3,NCALC,IGC)
      CALL CELREF2(INDIC,BB,AFI,LHKL,TH3,IHKL,DDT,DDQ)
      CNCALC(IGC)=NCALC
      IF(NDAT.GE.20)THEN
      FM20(IGC)=QO(20)/(2.*CNCALC(IGC)*DDQ)
      FF20(IGC)=20./(CNCALC(IGC)*DDT)
      ELSE
      PNDAT=NDAT
      FM20(IGC)=QO(NDAT)/(2.*CNCALC(IGC)*DDQ)
      FF20(IGC)=PNDAT/(CNCALC(IGC)*DDT)
      ENDIF
      ENDIF
C
C Test if cell already found
C
      IF(IGC.GT.1)THEN
      DO 118 I=1,IGC-1
      IF(IFI(I).NE.IFILE)GO TO 118
      VDELT=VGC(IGC)/300.
      VP=VGC(IGC)+VDELT
      VM=VGC(IGC)-VDELT
      IF(VGC(I).GT.VP.OR.VGC(I).LT.VM)GO TO 118
      NSOL(I)=NSOL(I)+1
      IF(RP(IGC).LT.RP(I))THEN
      IF(ISEE.EQ.1)WRITE(*,1115)RMAX,A,V2,IPEN
      KM(I)=KM(IGC)
      VGC(I)=VGC(IGC)
      RP(I)=RP(IGC)
      CEL(1,I)=CEL(1,IGC)
      CEL(2,I)=CEL(2,IGC)
      CEL(3,I)=CEL(3,IGC)
      ENDIF
      IGC=IGC-1
      IF(NSOL(I).GT.5)THEN
      NTRIED=TMAX+1.
      NOUT=NOUT+1
      ENDIF
      GO TO 119
118   CONTINUE
      IF(IVERB.EQ.1)WRITE(20,115)RMAX,NTRIED,A,V2,IPEN,ICODE
      IF(ISEE.EQ.1)WRITE(*,1115)RMAX,A,V2,IPEN
119   CONTINUE
      ELSE
      IF(IVERB.EQ.1)WRITE(20,115)RMAX,NTRIED,A,V2,IPEN,ICODE
      IF(ISEE.EQ.1)WRITE(*,1115)RMAX,A,V2,IPEN
      ENDIF
197   CONTINUE
C
C$OMP END CRITICAL(FOUND)
C
C
C
C
C First criterium reinitialized to Rmaxref
C
117   RMAX=RMAXREF
C
C... Stop if max limit of Monte Carlo tests outpassed
C         or if K is pressed (tested every 30000 MC event)
C
C   END ON MC tests
C
196   CONTINUE
C
C$OMP END DO NOWAIT
C$OMP END PARALLEL
C
      IF(INTEREST.GE.1)GO TO 5000
      CALL KILLK(PRESSEDK)
      IF(RMIN.EQ.RMAX)GO TO 198
      IF(IVERB.EQ.1)THEN
      WRITE(20,*)
      WRITE(20,*)'Best result : a=',BPAR(1),'V=',V3,' Rp=',RMIN
      WRITE(20,*)
      CALL DATN(DATENOW,TIMENOW)
      ENDIF
      RMIN=RMAX
198   IF(PRESSEDK)GO TO 5000
C
C  END ON NRUNS
C
199    CONTINUE
C
200   IF(NSYS(2).EQ.0)GO TO 300
C
C    Hexagonal case
C
C
      IF(NGRID.EQ.3)IHR=1
290   CONTINUE
      IF(IHR.EQ.2)NSYS(2)=2
      IF(NSYS(2).EQ.1)THEN
      PRINT *,'Hexagonal:    Rp     a       c        V     Nind'
      RPSMALL=1.
      ELSE
      PRINT *,'Rhombohedral: Rp     a       c        V     Nind'
      RPSMALL=1.
      ENDIF
C
      IFILE=2
      NCYCLES=500.
      CY=NCYCLES*1.1
      CELPRE(4)=90.
      CELPRE(5)=90.
      CELPRE(6)=120.
C
      IF(NGRID.EQ.3)THEN
      NRUNS=10
      PMIN=2.
      PMAX=30.
      IF(NSYS(2).EQ.2)PMAX=60.
      PMA(3)=DMAX1*3.1
      IF(NSYS(2).EQ.2)PMA(3)=DMAX1*6.1
      IF(PMA(3).GT.PMAX)PMA(3)=PMAX
      PMA(1)=DMAX1*2.1
      IF(PMA(1).GT.PMAX)PMA(1)=PMAX
      PMA(2)=PMA(1)
      VMIN=8.
      DO 223 I=1,3
      PMI(I)=PMIN
      DELTA(I)=(PMA(I)-PMI(I))/2.
223   PSTART(I)=PMI(I)
      VMAX=PMA(1)*PMA(2)*PMA(3)
      IF(VMAX.GT.4000.)VMAX=4000.
      NTIMELIM(2)=VMAX*5.
      ENDIF
C
      WRITE(20,*)
      WRITE(20,*)'Hexagonal/Trigonal/Rhomboedral Monte Carlo search :'
      WRITE(20,*)' Max(a,c), V ',PMA(1),PMA(3),VMAX
      WRITE(20,*)
C
      IF(IVERB.EQ.1)THEN
      IF(NSYS(2).EQ.1)WRITE(20,280)NRUN,NTIMELIM(2)
      IF(NSYS(2).EQ.2)WRITE(20,281)NRUN,NTIMELIM(2)
      WRITE(20,*)'====================================================
     1==========================='
      WRITE(20,*)' Rp  Trial number    a      c        V  Nind Icod'
      WRITE(20,*)
      ENDIF
280   FORMAT('  Results in hexagonal, run, tests :',I3,F12.0)
281   FORMAT('  Results in rhombohedral, run, tests :',I3,F12.0)
C
C     READ hkl Miller indices in hex.hkl
C
      IF(NSYS(2).EQ.2)GO TO 260
      TEMPO='hex.hkl'
      GO TO 261
260   TEMPO='rho.hkl'
      IFILE=7
261   CALL OPEN_READ1(35,TEMPO)
      READ(35,*)NHKL0
      IF(NSYS(2).EQ.2)GO TO 262
      NHKL0=NDAT*12
      IF(NHKL0.GT.800)NHKL0=800
      GO TO 263
262   NHKL0=NDAT*12
      IF(NHKL0.GT.600)NHKL0=600
263   DO 201 I=1,NHKL0
201   READ(35,*)(IHH(KK,I),KK=1,3)
      CLOSE(35)
C
      DO 299 NRUN=1,NRUNS
*------------------------------------------------------------------------- 
*     Initialisation
*
C      CALL ESP_INIT(ISEED)
*
*-------------------------------------------------------------------------
      RMAX=RMAXREF
      RMIN=RMAX
C
C...  here starts the loop
C
      INTEREST=0
      TMAX=NTIMELIM(2)/PROCS
      TTMAX=10.*NTIMELIM(2)
      NCELLS=NTIMELIM(2)
      IISEED=0
      NTRIED=0
      NTRIEDT=0
      NOUT=0
      CELPRE(1)=PSTART(1)+2.*DELTA(1)*RANDI(ISEED)
      CELPRE(2)=CELPRE(1)
      CELPRE(3)=PSTART(3)+2.*DELTA(3)*RANDI(ISEED)
      CELOLD(1)=CELPRE(1)
      CELOLD(3)=CELPRE(3)
      RGLOB=1.
      NGLOB=0
C
C$OMP PARALLEL DEFAULT(SHARED) COPYIN(/CAL/,/CAL2/)
C$OMP& PRIVATE(NCEL,NTRIEDB,DEL,V1,ICODE,LLHKL,IHKL,TH3,NCALC,
C$OMP& RMAX2,A,C,V2,BPAR,V3,PSTARTB,IPEN,ISEE,INDIC,IP,
C$OMP& DIFF,DIFF2,DDT,DDQ)
C$OMP& FIRSTPRIVATE(iseed,iiseed,rmax0,ntried,ntriedt,nout,
C$OMP& celpre,celold,rglob,nglob,rmin,rmax,bb,afi)
C$OMP DO
C
      DO 296 NCEL=1,NCELLS
      IF(NOUT.GE.1)GO TO 296
      IF(INTEREST.GE.1)GO TO 296
      IISEED=IISEED+1
      IF(IISEED.EQ.1)ISEED=((ISEED-NCEL*NRUN)/2)*2+1
C
202   CONTINUE
C
C     Which parameter to vary ? a or c ?
C
      NTRIEDB=0.
      IP=3
      IF(RANDI(ISEED).GT.0.5)IP=1
      CELPRE(IP)=PSTART(IP)+2.*DELTA(IP)*RANDI(ISEED)
      NTRIED=NTRIED+1.
      GO TO 204
203   DEL=DELTAB*(1.-NTRIEDB/CY)
      I=3
      IF(RANDI(ISEED).GT.0.5)I=1
      CELPRE(I)=PSTARTB(I)+(DEL*(RANDI(ISEED)-0.5)*2.)
      NTRIEDB=NTRIEDB+1.
204   CONTINUE
      CELPRE(2)=CELPRE(1)
      DO 205 I=1,3
      DO 205 J=1,3
205   AL(I,J)=0.0
      CALL DCELL(CELPRE,AL,V1)
      IF(NTRIED.GT.TMAX)THEN
      NOUT=NOUT+1
      GO TO 296
      ENDIF
      IF(NTRIEDB.NE.0.)GO TO 206
      IF(V1.GT.VMAX.OR.V1.LT.VMIN) THEN
      NTRIED=NTRIED-1.
      NTRIEDT=NTRIEDT+1.
      IF(NTRIEDT.GT.TTMAX)THEN
      NOUT=NOUT+1
      GO TO 296
      ENDIF
      GO TO 202
      ENDIF
C
206   CALL CALCUL1(DIFF,DIFF2)
      IF(NMX.GT.NDAT10)THEN
      NTRIED=NTRIED-1
      GO TO 202
      ENDIF
      IF(NTRIEDB.NE.0.)GO TO 214
C
C... Rp value satisfying ???
C
      IF(DIFF.LT.RGLOB.OR.LHKL.GT.NGLOB)THEN
      RGLOB=DIFF
      NGLOB=LHKL
      CELOLD(IP)=CELPRE(IP)
      ENDIF
      IF(LHKL.GE.NMAX)THEN
      RMAX=DIFF
      ICODE=2
      IF(DIFF.LE.RMAXREF)ICODE=1
      ELSE
      ICODE=1
      ENDIF
      IF(DIFF.GT.RMAX) GO TO 217
      IF(LHKL.LT.NMAX) GO TO 217
214   IF(DIFF.LE.RMAX)THEN
      LLHKL=LHKL
      RMAX=DIFF
      RMAX2=DIFF2
      A=CELPRE(1)
      C=CELPRE(3)
      V2=V1
      IF(DIFF.LT.RMIN)THEN
      RMIN=DIFF
      BPAR(1)=A
      BPAR(3)=C
      V3=V1
      ENDIF
C
C... "Refine" that cell (by Monte Carlo too...)
C
      PSTARTB(1)=CELPRE(1)
      PSTARTB(3)=CELPRE(3)
      ENDIF
      IF(NTRIEDB.LE.NCYCLES)GO TO 203
      RGLOB=0.5
      NGLOB=NDAT2
      IF(IP.EQ.1)CELOLD(IP)=A
      IF(IP.EQ.3)CELOLD(IP)=C
      NTRIEDB=0.
      IF(RMAX.GE.RMAX0(2))GO TO 217
      IF(RMAX2.GE.0.15)GO TO 217
      IPEN=NDAT-LLHKL
      IF(IPEN.GT.NIND)GO TO 217
C
C$OMP CRITICAL (STORE1)
C
      IGC=IGC+1
C
C  Test if too much proposals, if yes decrease Rmax by 5%
C
      IGT=IGT+1.
      IF(NR.EQ.1)THEN
      IF(IGT.GT.50.)THEN
      IF((NTRIED/IGT).LT.1000.)THEN
      IF(RMAX0(2).GT.0.2)THEN
      RMAX0(2)=RMAX0(2)-RMAX0(2)*0.05
      WRITE(20,*)'  Rmax reduced by 5%, now Rmax = ',RMAX0(2)
      WRITE(*,*)'  Rmax reduced by 5%, now Rmax = ',RMAX0(2)
      ENDIF
      ENDIF
      ENDIF
      ENDIF
C
      IF(IGC.GT.10000)THEN
      WRITE(20,*)'   More than 10000 good cells = STOP'
      PRINT *,'   More than 10000 good cells = STOP'
      IGC=IGC-1
      INTEREST=INTEREST+1
C      GO TO 5000
      ENDIF
      CEL(1,IGC)=A
      CEL(2,IGC)=A
      CEL(3,IGC)=C
      CEL(4,IGC)=90.
      CEL(5,IGC)=90.
      CEL(6,IGC)=120.
C
C$OMP END CRITICAL(STORE1)
C
C... Check for supercell
C
      CELPRE(1)=A
      CELPRE(2)=A
      CELPRE(3)=C
      DO 240 I=1,3
      DO 240 J=1,3
240   AL(I,J)=0.0
      CALL DCELL(CELPRE,AL,V1)
C
C$OMP CRITICAL(STORE2)
C
      CALL CALCUL2(DIFF,IHKL,TH3,NCALC,IGC)
      KM(IGC)=LLHKL
      KM2(IGC)=LHKL
      IFI(IGC)=IFILE
      NSOL(IGC)=1
      VGC(IGC)=V1
      RP(IGC)=RMAX
      RP2(IGC)=DIFF
      IF(RP(IGC).LT.RPSMALL)THEN
      RPSMALL=RP(IGC)
      ISEE=1
      ELSE
      ISEE=0
      ENDIF
      CALL SUPCEL(LHKL,IHKL,CEL,IGC,VGC,2)
      CALL BRAV(LHKL,IHKL,IBR)
      IB(IGC)=IBR
      A=CEL(1,IGC)
      CEL(2,IGC)=A
      C=CEL(3,IGC)
      V2=VGC(IGC)
C
C$OMP END CRITICAL(STORE2)
C
C... Check for interesting result
C
C      IF(INTEREST.GE.1)GO TO 296
      INDIC=2
      BB(3)=A
      BB(4)=A
      BB(5)=C
      BB(6)=90.
      BB(7)=90.
      BB(8)=120.
      AFI(3)=1
      AFI(4)=1
      AFI(5)=1
      AFI(6)=0
      AFI(7)=0
      AFI(8)=0
      CELPRE(1)=A
      CELPRE(2)=A
      CELPRE(3)=C
      DO 210 I=1,3
      DO 210 J=1,3
210   AL(I,J)=0.0
C
C$OMP CRITICAL(FOUND)
C
      IF(RP(IGC).LT.RMI)THEN
      INTEREST=INTEREST+1
      WRITE(*,1215)RMAX,A,C,V2,IPEN
      WRITE(20,*)
      WRITE(20,*)'====================================================
     1==========================='
      WRITE(20,*)
      WRITE(20,*)' YOU HAVE FOUND AN INTERESTING RESULT : Rp < Rmin!'
      WRITE(*,*)
      WRITE(*,*)' YOU HAVE FOUND AN INTERESTING RESULT : Rp < Rmin!'
C
C... Refine that cell
C
      CALL DCELL(CELPRE,AL,V1)
      CALL CALCUL2(DIFF,IHKL,TH3,NCALC,IGC)
      CALL CELREF2(INDIC,BB,AFI,LHKL,TH3,IHKL,DDT,DDQ)
      CNCALC(IGC)=NCALC
      IF(NDAT.GE.20)THEN
      FM20(IGC)=QO(20)/(2.*CNCALC(IGC)*DDQ)
      FF20(IGC)=20./(CNCALC(IGC)*DDT)
      ELSE
      PNDAT=NDAT
      FM20(IGC)=QO(NDAT)/(2.*CNCALC(IGC)*DDQ)
      FF20(IGC)=PNDAT/(CNCALC(IGC)*DDT)
      ENDIF
C      WRITE(20,7000)FM20(IGC)
C      WRITE(20,7001)FF20(IGC),DDT,NCALC
C	WRITE(20,*)
C      PRINT 7000,FM20(IGC)
C      PRINT 7001,FF20(IGC),DDT,NCALC
C	PRINT *
      IREF=1
      GO TO 297
      ELSE
C
C  Anyway, calculate the M20 and F20 values
C
      CALL DCELL(CELPRE,AL,V1)
      CALL CALCUL2(DIFF,IHKL,TH3,NCALC,IGC)
      CALL CELREF2(INDIC,BB,AFI,LHKL,TH3,IHKL,DDT,DDQ)
      CNCALC(IGC)=NCALC
      IF(NDAT.GE.20)THEN
      FM20(IGC)=QO(20)/(2.*CNCALC(IGC)*DDQ)
      FF20(IGC)=20./(CNCALC(IGC)*DDT)
      ELSE
      PNDAT=NDAT
      FM20(IGC)=QO(NDAT)/(2.*CNCALC(IGC)*DDQ)
      FF20(IGC)=PNDAT/(CNCALC(IGC)*DDT)
      ENDIF
      ENDIF
C
C Test if cell already found
C
      IF(IGC.GT.1)THEN
      DO 218 I=1,IGC-1
      IF(IFI(I).NE.IFILE)GO TO 218
      VDELT=VGC(IGC)/300.
      VP=VGC(IGC)+VDELT
      VM=VGC(IGC)-VDELT
      IF(VGC(I).GT.VP.OR.VGC(I).LT.VM)GO TO 218
      ADELT=CEL(1,IGC)/500.
      AP=CEL(1,IGC)+ADELT
      AM=CEL(1,IGC)-ADELT
      IF(CEL(1,I).GT.AP.OR.CEL(1,I).LT.AM)GO TO 218
      NSOL(I)=NSOL(I)+1
      IF(RP(IGC).LT.RP(I))THEN
      IF(ISEE.EQ.1)WRITE(*,1215)RMAX,A,C,V2,IPEN
      KM(I)=KM(IGC)
      VGC(I)=VGC(IGC)
      RP(I)=RP(IGC)
      CEL(1,I)=CEL(1,IGC)
      CEL(2,I)=CEL(2,IGC)
      CEL(3,I)=CEL(3,IGC)
      ENDIF
      IGC=IGC-1
      IF(NSOL(I).GT.5)THEN
      NTRIED=TMAX+1.
      NOUT=NOUT+1
      ENDIF
      GO TO 219
218   CONTINUE
      IF(IVERB.EQ.1)WRITE(20,215)RMAX,NTRIED,A,C,V2,IPEN,ICODE
      IF(ISEE.EQ.1)WRITE(*,1215)RMAX,A,C,V2,IPEN
219   CONTINUE
      ELSE
      IF(IVERB.EQ.1)WRITE(20,215)RMAX,NTRIED,A,C,V2,IPEN,ICODE
      IF(ISEE.EQ.1)WRITE(*,1215)RMAX,A,C,V2,IPEN
      ENDIF
297   CONTINUE
C
C$OMP END CRITICAL(FOUND)
C
217   RMAX=RMAXREF
      IF(RANDI(ISEED).GT.ESCAPE)CELPRE(IP)=CELOLD(IP)
C
C  END ON MC TESTS
C
296   CONTINUE
C
C$OMP END DO NOWAIT
C$OMP END PARALLEL
C
C
C... Stop if max limit of Monte Carlo tests outpassed
C         or if K is pressed (tested every 30000 MC event)
C
      IF(INTEREST.GE.1)GO TO 5000
      CALL KILLK(PRESSEDK)
      IF(RMIN.EQ.RMAX)GO TO 298
      IF(IVERB.EQ.1)THEN
      WRITE(20,*)
      WRITE(20,*)'Best result : a = ',BPAR(1),' Rp = ',RMIN
      WRITE(20,*)'Best result : c = ',BPAR(3),'V = ',V3
      WRITE(20,*)
      CALL DATN(DATENOW,TIMENOW)
      ENDIF
      RMIN=RMAX
298   IF(PRESSEDK)GO TO 5000
C
C  END ON NRUNS
C
299   CONTINUE
C
C
      IHR=IHR+1
      IF(IHR.EQ.2)GO TO 290
300   IF(NSYS(3).EQ.0)GO TO 400
C
C    Tetragonal case
C
C
      RPSMALL=1.
      PRINT *,'Tetragonal:   Rp     a       c        V     Nind'
C
      IFILE=3
      NCYCLES=500.
      CY=NCYCLES*1.1
      CELPRE(4)=90.
      CELPRE(5)=90.
      CELPRE(6)=90.
C
      IF(NGRID.EQ.3)THEN
      NRUNS=10
      PMIN=2.
      PMAX=30.
      PMA(3)=DMAX1*4.1
      IF(PMA(3).GT.PMAX)PMA(3)=PMAX
      PMA(1)=DMAX1*2.1
      IF(PMA(1).GT.PMAX)PMA(1)=PMAX
      PMA(2)=PMA(1)
      VMIN=8.
      DO 323 I=1,3
      PMI(I)=PMIN
      DELTA(I)=(PMA(I)-PMI(I))/2.
323   PSTART(I)=PMI(I)
      VMAX=PMA(1)*PMA(2)*PMA(3)
      IF(VMAX.GT.4000.)VMAX=4000.
      NTIMELIM(3)=VMAX*5.
      ENDIF
C
      WRITE(20,*)
      WRITE(20,*)'Tetragonal Monte Carlo search :'
      WRITE(20,*)' Max(a,c), V ',PMA(1),PMA(3),VMAX
      WRITE(20,*)
C
      IF(IVERB.EQ.1)THEN
      WRITE(20,380)NRUN,NTIMELIM(3)
      WRITE(20,*)'====================================================
     1==========================='
      WRITE(20,*)' Rp  Trial number    a      c        V  Nind Icod'
      WRITE(20,*)
      ENDIF
380   FORMAT('  Results in tetragonal, run, tests :',I3,F12.0)
C
C     READ hkl Miller indices in tet.hkl
C
      TEMPO='tet.hkl'
      CALL OPEN_READ1(35,TEMPO)
      READ(35,*)NHKL0
      NHKL0=NDAT*12
      IF(NHKL0.GT.800)NHKL0=800
      DO 301 I=1,NHKL0
301   READ(35,*)(IHH(KK,I),KK=1,3)
      CLOSE(35)
C
      DO 399 NRUN=1,NRUNS
*------------------------------------------------------------------------- 
*     Initialisation
*
C      CALL ESP_INIT(ISEED)
*
*-------------------------------------------------------------------------
      RMAX=RMAXREF
      RMIN=RMAX
C
C...  here starts the loop
C
      INTEREST=0
      TMAX=NTIMELIM(3)/PROCS
      TTMAX=10.*NTIMELIM(3)
      NCELLS=NTIMELIM(3)
      IISEED=0
      NTRIED=0.
      NTRIEDT=0.
      NOUT=0
      CELPRE(1)=PSTART(1)+2.*DELTA(1)*RANDI(ISEED)
      CELPRE(2)=CELPRE(1)
      CELPRE(3)=PSTART(3)+2.*DELTA(3)*RANDI(ISEED)
      CELOLD(1)=CELPRE(1)
      CELOLD(3)=CELPRE(3)
      RGLOB=1.
      NGLOB=0
C
C$OMP PARALLEL DEFAULT(SHARED) COPYIN(/CAL/,/CAL2/)
C$OMP& PRIVATE(NCEL,NTRIEDB,DEL,V1,ICODE,LLHKL,IHKL,TH3,NCALC,
C$OMP& RMAX2,A,C,V2,BPAR,V3,PSTARTB,IPEN,ISEE,INDIC,IP,
C$OMP& DIFF,DIFF2,DDT,DDQ)
C$OMP& FIRSTPRIVATE(iseed,iiseed,rmax0,ntried,ntriedt,nout,
C$OMP& celpre,celold,rglob,nglob,rmin,rmax,bb,afi)
C$OMP DO
C
      DO 396 NCEL=1,NCELLS
      IF(NOUT.GE.1)GO TO 396
      IF(INTEREST.GE.1)GO TO 396
      IISEED=IISEED+1
      IF(IISEED.EQ.1)ISEED=((ISEED-NCEL*NRUN)/2)*2+1
C
302   CONTINUE
C
C     Which parameter to vary ? a or c ?
C
      NTRIEDB=0.
      IP=3
      IF(RANDI(ISEED).GT.0.5)IP=1
      CELPRE(IP)=PSTART(IP)+2.*DELTA(IP)*RANDI(ISEED)
      NTRIED=NTRIED+1.
      GO TO 304
303   DEL=DELTAB*(1.-NTRIEDB/CY)
      I=3
      IF(RANDI(ISEED).GT.0.5)I=1
      CELPRE(I)=PSTARTB(I)+(DEL*(RANDI(ISEED)-0.5)*2.)
      NTRIEDB=NTRIEDB+1.
304   CONTINUE
      CELPRE(2)=CELPRE(1)
      DO 305 I=1,3
      DO 305 J=1,3
305   AL(I,J)=0.0
      CALL DCELL(CELPRE,AL,V1)
      IF(NTRIED.GT.TMAX)THEN
      NOUT=NOUT+1
      GO TO 396
      ENDIF
      IF(NTRIEDB.NE.0.)GO TO 306
      IF(V1.GT.VMAX.OR.V1.LT.VMIN) THEN
      NTRIED=NTRIED-1.
      NTRIEDT=NTRIEDT+1.
      IF(NTRIEDT.GT.TTMAX)THEN
      NOUT=NOUT+1
      GO TO 396
      ENDIF
      GO TO 302
      ENDIF
C
306   CALL CALCUL1(DIFF,DIFF2)
      IF(NMX.GT.NDAT10)THEN
      NTRIED=NTRIED-1
      GO TO 302
      ENDIF
      IF(NTRIEDB.NE.0.)GO TO 314
C
C... Rp value satisfying ???
C
      IF(DIFF.LT.RGLOB.OR.LHKL.GT.NGLOB)THEN
      RGLOB=DIFF
      NGLOB=LHKL
      CELOLD(IP)=CELPRE(IP)
      ENDIF
      IF(LHKL.GE.NMAX)THEN
      RMAX=DIFF
      ICODE=2
      IF(DIFF.LE.RMAXREF)ICODE=1
      ELSE
      ICODE=1
      ENDIF
      IF(DIFF.GT.RMAX) GO TO 317
      IF(LHKL.LT.NMAX) GO TO 317
314   IF(DIFF.LE.RMAX)THEN
      LLHKL=LHKL
      RMAX=DIFF
      RMAX2=DIFF2
      A=CELPRE(1)
      C=CELPRE(3)
      V2=V1
      IF(DIFF.LT.RMIN)THEN
      RMIN=DIFF
      BPAR(1)=A
      BPAR(3)=C
      V3=V1
      ENDIF
C
C... "Refine" that cell (by Monte Carlo too...)
C
      PSTARTB(1)=CELPRE(1)
      PSTARTB(3)=CELPRE(3)
      ENDIF
      IF(NTRIEDB.LE.NCYCLES)GO TO 303
      RGLOB=0.5
      NGLOB=NDAT2
      IF(IP.EQ.1)CELOLD(IP)=A
      IF(IP.EQ.3)CELOLD(IP)=C
      NTRIEDB=0.
      IF(RMAX.GE.RMAX0(3))GO TO 317
      IF(RMAX2.GE.0.15)GO TO 317
      IPEN=NDAT-LLHKL
      IF(IPEN.GT.NIND)GO TO 317
C
C$OMP CRITICAL(STORE1)
C
      IGC=IGC+1
C
C  Test if too much proposals, if yes decrease Rmax by 5%
C
      IGT=IGT+1.
      IF(NR.EQ.1)THEN
      IF(IGT.GT.50.)THEN
      IF((NTRIED/IGT).LT.1000.)THEN
      IF(RMAX0(3).GT.0.2)THEN
      RMAX0(3)=RMAX0(3)-RMAX0(3)*0.05
      WRITE(20,*)'  Rmax reduced by 5%, now Rmax = ',RMAX0(3)
      WRITE(*,*)'  Rmax reduced by 5%, now Rmax = ',RMAX0(3)
      ENDIF
      ENDIF
      ENDIF
      ENDIF
C
      IF(IGC.GT.10000)THEN
      WRITE(20,*)'   More than 10000 good cells = STOP'
      PRINT *,'   More than 10000 good cells = STOP'
      IGC=IGC-1
      INTEREST=INTEREST+1
C      GO TO 5000
      ENDIF
      CEL(1,IGC)=A
      CEL(2,IGC)=A
      CEL(3,IGC)=C
      CEL(4,IGC)=90.
      CEL(5,IGC)=90.
      CEL(6,IGC)=90.
C
C$OMP END CRITICAL(STORE1)
C
C... Check for supercell
C
      CELPRE(1)=A
      CELPRE(2)=A
      CELPRE(3)=C
      DO 340 I=1,3
      DO 340 J=1,3
340   AL(I,J)=0.0
      CALL DCELL(CELPRE,AL,V1)
C
C$OMP CRITICAL(STORE2)
C
      CALL CALCUL2(DIFF,IHKL,TH3,NCALC,IGC)
      KM(IGC)=LLHKL
      KM2(IGC)=LHKL
      IFI(IGC)=IFILE
      NSOL(IGC)=1
      VGC(IGC)=V1
      RP(IGC)=RMAX
      RP2(IGC)=DIFF
      IF(RP(IGC).LT.RPSMALL)THEN
      RPSMALL=RP(IGC)
      ISEE=1
      ELSE
      ISEE=0
      ENDIF
      CALL SUPCEL(LHKL,IHKL,CEL,IGC,VGC,4)
      CALL BRAV(LHKL,IHKL,IBR)
      IB(IGC)=IBR
      A=CEL(1,IGC)
      CEL(2,IGC)=A
      C=CEL(3,IGC)
      V2=VGC(IGC)
C
C$OMP END CRITICAL(STORE2)
C
C... Check for interesting result
C
C      IF(INTEREST.GE.1)GO TO 396
      INDIC=2
      BB(3)=A
      BB(4)=A
      BB(5)=C
      BB(6)=90.
      BB(7)=90.
      BB(8)=90.
      AFI(3)=1
      AFI(4)=1
      AFI(5)=1
      AFI(6)=0
      AFI(7)=0
      AFI(8)=0
      CELPRE(1)=A
      CELPRE(2)=A
      CELPRE(3)=C
      DO 310 I=1,3
      DO 310 J=1,3
310   AL(I,J)=0.0
C
C$OMP CRITICAL(FOUND)
C
      IF(RP(IGC).LT.RMI)THEN
      INTEREST=INTEREST+1
      WRITE(*,1215)RMAX,A,C,V2,IPEN
      WRITE(20,*)
      WRITE(20,*)'====================================================
     1==========================='
      WRITE(20,*)
      WRITE(20,*)' YOU HAVE FOUND AN INTERESTING RESULT : Rp < Rmin!'
      WRITE(*,*)
      WRITE(*,*)' YOU HAVE FOUND AN INTERESTING RESULT : Rp < Rmin!'
C
C... Refine that cell
C
      CALL DCELL(CELPRE,AL,V1)
      CALL CALCUL2(DIFF,IHKL,TH3,NCALC,IGC)
      CALL CELREF2(INDIC,BB,AFI,LHKL,TH3,IHKL,DDT,DDQ)
      CNCALC(IGC)=NCALC
      IF(NDAT.GE.20)THEN
      FM20(IGC)=QO(20)/(2.*CNCALC(IGC)*DDQ)
      FF20(IGC)=20./(CNCALC(IGC)*DDT)
      ELSE
      PNDAT=NDAT
      FM20(IGC)=QO(NDAT)/(2.*CNCALC(IGC)*DDQ)
      FF20(IGC)=PNDAT/(CNCALC(IGC)*DDT)
      ENDIF
C      WRITE(20,7000)FM20(IGC)
C      WRITE(20,7001)FF20(IGC),DDT,NCALC
C	WRITE(20,*)
C      PRINT 7000,FM20(IGC)
C      PRINT 7001,FF20(IGC),DDT,NCALC
C	PRINT *
      IREF=1
      GO TO 397
      ELSE
C
C  Anyway, calculate the M20 and F20 values
C
      CALL DCELL(CELPRE,AL,V1)
      CALL CALCUL2(DIFF,IHKL,TH3,NCALC,IGC)
      CALL CELREF2(INDIC,BB,AFI,LHKL,TH3,IHKL,DDT,DDQ)
      CNCALC(IGC)=NCALC
      IF(NDAT.GE.20)THEN
      FM20(IGC)=QO(20)/(2.*CNCALC(IGC)*DDQ)
      FF20(IGC)=20./(CNCALC(IGC)*DDT)
      ELSE
      PNDAT=NDAT
      FM20(IGC)=QO(NDAT)/(2.*CNCALC(IGC)*DDQ)
      FF20(IGC)=PNDAT/(CNCALC(IGC)*DDT)
      ENDIF
      ENDIF
C
C Test if cell already found
C
      IF(IGC.GT.1)THEN
      DO 318 I=1,IGC-1
      IF(IFI(I).NE.IFILE)GO TO 318
      VDELT=VGC(IGC)/300.
      VP=VGC(IGC)+VDELT
      VM=VGC(IGC)-VDELT
      IF(VGC(I).GT.VP.OR.VGC(I).LT.VM)GO TO 318
      ADELT=CEL(1,IGC)/500.
      AP=CEL(1,IGC)+ADELT
      AM=CEL(1,IGC)-ADELT
      IF(CEL(1,I).GT.AP.OR.CEL(1,I).LT.AM)GO TO 318
      NSOL(I)=NSOL(I)+1
      IF(RP(IGC).LT.RP(I))THEN
      IF(ISEE.EQ.1)WRITE(*,1215)RMAX,A,C,V2,IPEN
      KM(I)=KM(IGC)
      VGC(I)=VGC(IGC)
      RP(I)=RP(IGC)
      CEL(1,I)=CEL(1,IGC)
      CEL(2,I)=CEL(2,IGC)
      CEL(3,I)=CEL(3,IGC)
      ENDIF
      IGC=IGC-1
      IF(NSOL(I).GT.5)THEN
      NTRIED=TMAX+1.
      NOUT=NOUT+1
      ENDIF
      GO TO 319
318   CONTINUE
      IF(IVERB.EQ.1)WRITE(20,215)RMAX,NTRIED,A,C,V2,IPEN,ICODE
      IF(ISEE.EQ.1)WRITE(*,1215)RMAX,A,C,V2,IPEN
319   CONTINUE
      ELSE
      IF(IVERB.EQ.1)WRITE(20,215)RMAX,NTRIED,A,C,V2,IPEN,ICODE
      IF(ISEE.EQ.1)WRITE(*,1215)RMAX,A,C,V2,IPEN
      ENDIF
397   CONTINUE
C
C$OMP END CRITICAL(FOUND)
C
C
317   RMAX=RMAXREF
      IF(RANDI(ISEED).GT.ESCAPE)CELPRE(IP)=CELOLD(IP)
C
C   END ON 4000 tests
C
396   CONTINUE
C
C$OMP END DO NOWAIT
C$OMP END PARALLEL
C
C
C... Stop if max limit of Monte Carlo tests outpassed
C         or if K is pressed (tested every 30000 MC event)
C
      IF(INTEREST.GE.1)GO TO 5000
      CALL KILLK(PRESSEDK)
      IF(RMIN.EQ.RMAX)GO TO 398
      IF(IVERB.EQ.1)THEN
      WRITE(20,*)
      WRITE(20,*)'Best result : a = ',BPAR(1),' Rp = ',RMIN
      WRITE(20,*)'Best result : c = ',BPAR(3),'V = ',V3
      WRITE(20,*)
      CALL DATN(DATENOW,TIMENOW)
      ENDIF
      RMIN=RMAX
398   IF(PRESSEDK)GO TO 5000
C
C  END ON NRUNS
C
399   CONTINUE
C
C
400   IF(NSYS(4).EQ.0)GO TO 500
C
C    Orthorhombic case
C
C
      RPSMALL=1.
      PRINT *,'Orthorhombic: Rp     a       b       c        V     Nind'
      IFILE=4
      NCYCLES=1000.
      CY=NCYCLES*1.1
      CELPRE(4)=90.
      CELPRE(5)=90.
      CELPRE(6)=90.
C
      IF(NGRID.EQ.3)THEN
      NRUNS2=6
      NRUNS=10
      PMIN=2.
      PMAX=20.
      PMA(1)=DMAX1*2.1
      PMA(2)=DMAX2*2.1
      PMA(3)=DMAX3*2.1
      IF(PMA(1).GT.PMAX)PMA(1)=PMAX
      IF(PMA(2).GT.PMAX)PMA(2)=PMAX
      IF(PMA(3).GT.PMAX)PMA(3)=PMAX
      VORTH=PMA(1)*PMA(2)*PMA(3)
      IF(VORTH.GT.3000.)VORTH=3000.
      VMAX=VORTH
      DO 423 I=1,3
      PMI(I)=PMIN
      DELTA(I)=(PMA(I)-PMI(I))/2.
423   PSTART(I)=PMI(I)
      ENDIF
C
      WRITE(20,*)
      WRITE(20,*)'Orthorhombic Monte Carlo search :'
      WRITE(20,*)' Max(a,b,c), V ',PMA(1),PMA(2),PMA(3),VMAX
      WRITE(20,*)
C
      IF(IVERB.EQ.1)THEN
      WRITE(20,480)NRUN,NTIMELIM(4)
      WRITE(20,*)'====================================================
     1==========================='
      WRITE(20,*)' Rp  Trial number    a   b   c       V  Nind icod'
      WRITE(20,*)
      ENDIF
480   FORMAT('  Results in orthorhombic, run, tests :',I3,F12.0)
C
C     READ hkl Miller indices in ort.hkl
C
      TEMPO='ort.hkl'
      CALL OPEN_READ1(35,TEMPO)
      READ(35,*)NHKL0
      NHKL0=NDAT*20
      IF(NHKL0.GT.1000)NHKL0=1000
      DO 401 I=1,NHKL0
401   READ(35,*)(IHH(KK,I),KK=1,3)
      CLOSE(35)
C
      DO 499 NRUN2=1,NRUNS2
C
      IF(NGRID.EQ.3)THEN
      IF(NRUN2.EQ.1)THEN
      VMIN=8.
      VMAX=500.
      IF(VORTH.LT.500.)VMAX=VORTH
      NTIMELIM(4)=(VMAX-VMIN)*20.
      ENDIF
      IF(NRUN2.EQ.2)THEN
      IF(VORTH.LT.500.)GO TO 500
      VMIN=500.
      VMAX=1000.
      IF(VORTH.LT.1000.)VMAX=VORTH
      NTIMELIM(4)=(VMAX-VMIN)*20.
      ENDIF
      IF(NRUN2.EQ.3)THEN
      IF(VORTH.LT.1000.)GO TO 500
      VMIN=1000.
      VMAX=1500.
      IF(VORTH.LT.1500.)VMAX=VORTH
      NTIMELIM(4)=(VMAX-VMIN)*20.
      ENDIF
      IF(NRUN2.EQ.4)THEN
      IF(VORTH.LT.1500.)GO TO 500
      VMIN=1500.
      VMAX=2000.
      IF(VORTH.LT.2000.)VMAX=VORTH
      NTIMELIM(4)=(VMAX-VMIN)*20.
      ENDIF
      IF(NRUN2.EQ.5)THEN
      IF(VORTH.LT.2000.)GO TO 500
      VMIN=2000.
      VMAX=2500.
      IF(VORTH.LT.2500.)VMAX=VORTH
      NTIMELIM(4)=(VMAX-VMIN)*20.
      ENDIF
      IF(NRUN2.EQ.6)THEN
      IF(VORTH.LT.2500.)GO TO 500
      VMIN=2500.
      VMAX=3000.
      IF(VORTH.LT.3000.)VMAX=VORTH
      NTIMELIM(4)=(VMAX-VMIN)*20.
      ENDIF
      ENDIF
C
      DO 499 NRUN=1,NRUNS
*------------------------------------------------------------------------- 
*     Initialisation
*
C      CALL ESP_INIT(ISEED)
*
*-------------------------------------------------------------------------
      RMAX=RMAXREF
      RMIN=RMAX
C
C...  here starts the loop
C
      INTEREST=0
      TMAX=NTIMELIM(4)/PROCS
      TTMAX=10.*NTIMELIM(4)
      NCELLS=NTIMELIM(4)
      IISEED=0
      NTRIED=0.
      NTRIEDT=0.
      NOUT=0
      CELPRE(1)=PSTART(1)+2.*DELTA(1)*RANDI(ISEED)
      CELPRE(2)=PSTART(2)+2.*DELTA(2)*RANDI(ISEED)
      CELPRE(3)=PSTART(3)+2.*DELTA(3)*RANDI(ISEED)
      CELOLD(1)=CELPRE(1)
      CELOLD(2)=CELPRE(2)
      CELOLD(3)=CELPRE(3)
      RGLOB=1.
      NGLOB=0
C
C$OMP PARALLEL DEFAULT(SHARED) COPYIN(/CAL/,/CAL2/)
C$OMP& PRIVATE(NCEL,NTRIEDB,DEL,V1,ICODE,LLHKL,IHKL,TH3,NCALC,
C$OMP& RMAX2,A,B,C,V2,BPAR,V3,PSTARTB,IPEN,ISEE,INDIC,IP,X,
C$OMP& DIFF,DIFF2,DDT,DDQ)
C$OMP& FIRSTPRIVATE(iseed,iiseed,rmax0,ntried,ntriedt,nout,
C$OMP& celpre,celold,rglob,nglob,rmin,rmax,bb,afi)
C$OMP DO
C
      DO 496 NCEL=1,NCELLS
      IF(NOUT.GE.1)GO TO 496
      IF(INTEREST.GE.1)GO TO 496
      IISEED=IISEED+1
      IF(IISEED.EQ.1)ISEED=((ISEED-NCEL*NRUN)/2)*2+1
C
402   CONTINUE
C
C     Which parameter to vary ? a or b or c ?
C
      NTRIEDB=0.
      X=RANDI(ISEED)
      IF(X.GE.0.AND.X.LT.0.33333)IP=1
      IF(X.GE.0.33333.AND.X.LT.0.66666)IP=2
      IF(X.GE.0.66666.AND.X.LE.1.)IP=3
      CELPRE(IP)=PSTART(IP)+2*DELTA(IP)*RANDI(ISEED)
      NTRIED=NTRIED+1.
      GO TO 404
403   DEL=DELTAB*(1.-NTRIEDB/CY)
      X=RANDI(ISEED)
      IF(X.GE.0.AND.X.LT.0.33333)I=1
      IF(X.GE.0.33333.AND.X.LT.0.66666)I=2
      IF(X.GE.0.66666.AND.X.LT.1.)I=3
      CELPRE(I)=PSTARTB(I)+(DEL*(RANDI(ISEED)-0.5)*2.)
      NTRIEDB=NTRIEDB+1.
404   CONTINUE
      DO 405 I=1,3
      DO 405 J=1,3
405   AL(I,J)=0.0
      CALL DCELL(CELPRE,AL,V1)
      IF(NTRIED.GT.TMAX)THEN
      NOUT=NOUT+1
      GO TO 496
      ENDIF
      IF(NTRIEDB.NE.0.)GO TO 406
      IF(V1.GT.VMAX.OR.V1.LT.VMIN) THEN
      NTRIED=NTRIED-1.
      NTRIEDT=NTRIEDT+1.
      IF(NTRIEDT.GT.TTMAX)THEN
      NOUT=NOUT+1
      GO TO 496
      ENDIF
      GO TO 402
      ENDIF
C
406   CALL CALCUL1(DIFF,DIFF2)
      IF(NMX.GT.NDAT10)THEN
      NTRIED=NTRIED-1
      GO TO 402
      ENDIF
      IF(NTRIEDB.NE.0.)GO TO 414
C
C... Rp value satisfying ???
C
      IF(DIFF.LT.RGLOB.OR.LHKL.GT.NGLOB)THEN
      RGLOB=DIFF
      NGLOB=LHKL
      CELOLD(IP)=CELPRE(IP)
      ENDIF
      IF(LHKL.GE.NMAX)THEN
      RMAX=DIFF
      ICODE=2
      IF(DIFF.LE.RMAXREF)ICODE=1
      ELSE
      ICODE=1
      ENDIF
      IF(DIFF.GT.RMAX) GO TO 417
      IF(LHKL.LT.NMAX) GO TO 417
414   IF(DIFF.LE.RMAX)THEN
      LLHKL=LHKL
      RMAX=DIFF
      RMAX2=DIFF2
      A=CELPRE(1)
      B=CELPRE(2)
      C=CELPRE(3)
      V2=V1
      IF(DIFF.LT.RMIN)THEN
      RMIN=DIFF
      BPAR(1)=A
      BPAR(2)=B
      BPAR(3)=C
      V3=V1
      ENDIF
C
C... "Refine" that cell (by Monte Carlo too...)
C
      PSTARTB(1)=CELPRE(1)
      PSTARTB(2)=CELPRE(2)
      PSTARTB(3)=CELPRE(3)
      ENDIF
      IF(NTRIEDB.LE.NCYCLES)GO TO 403
      RGLOB=0.5
      NGLOB=NDAT2
      IF(IP.EQ.1)CELOLD(IP)=A
      IF(IP.EQ.2)CELOLD(IP)=B
      IF(IP.EQ.3)CELOLD(IP)=C
      NTRIEDB=0.
      IF(RMAX.GE.RMAX0(4))GO TO 417
      IF(RMAX2.GE.0.15)GO TO 417
      IPEN=NDAT-LLHKL
      IF(IPEN.GT.NIND)GO TO 417
C
C$OMP CRITICAL(STORE1)
C
      IGC=IGC+1
C
C  Test if too much proposals, if yes decrease Rmax by 5%
C
      IGT=IGT+1.
      IF(NR.EQ.1)THEN
      IF(IGT.GT.50.)THEN
      IF((NTRIED/IGT).LT.10000.)THEN
      IF(RMAX0(4).GT.0.2)THEN
      RMAX0(4)=RMAX0(4)-RMAX0(4)*0.05
      WRITE(20,*)'  Rmax reduced by 5%, now Rmax = ',RMAX0(4)
      WRITE(*,*)'  Rmax reduced by 5%, now Rmax = ',RMAX0(4)
      ENDIF
      ENDIF
      ENDIF
      ENDIF
C
      IF(IGC.GT.10000)THEN
      WRITE(20,*)'   More than 10000 good cells = STOP'
      PRINT *,'   More than 10000 good cells = STOP'
      IGC=IGC-1
      INTEREST=INTEREST+1
      ENDIF
      CEL(1,IGC)=A
      CEL(2,IGC)=B
      CEL(3,IGC)=C
      CEL(4,IGC)=90.
      CEL(5,IGC)=90.
      CEL(6,IGC)=90.
C
C$OMP END CRITICAL(STORE1)
C
C... Check for supercell
C
      CELPRE(1)=A
      CELPRE(2)=B
      CELPRE(3)=C
      DO 440 I=1,3
      DO 440 J=1,3
440   AL(I,J)=0.0
      CALL DCELL(CELPRE,AL,V1)
C
C$OMP CRITICAL(STORE2)
C
      CALL CALCUL2(DIFF,IHKL,TH3,NCALC,IGC)
      KM(IGC)=LLHKL
      KM2(IGC)=LHKL
      IFI(IGC)=IFILE
      NSOL(IGC)=1
      VGC(IGC)=V1
      RP(IGC)=RMAX
      RP2(IGC)=DIFF
      IF(RP(IGC).LT.RPSMALL)THEN
      RPSMALL=RP(IGC)
      ISEE=1
      ELSE
      ISEE=0
      ENDIF
      CALL SUPCEL(LHKL,IHKL,CEL,IGC,VGC,1)
      CALL BRAV(LHKL,IHKL,IBR)
      IB(IGC)=IBR
      A=CEL(1,IGC)
      B=CEL(2,IGC)
      C=CEL(3,IGC)
      V2=VGC(IGC)
C
C$OMP END CRITICAL(STORE2)
C
C... Check for interesting result
C
C      IF(INTEREST.GE.1)GO TO 496
      INDIC=0
      BB(3)=A
      BB(4)=B
      BB(5)=C
      BB(6)=90.
      BB(7)=90.
      BB(8)=90.
      AFI(3)=1
      AFI(4)=1
      AFI(5)=1
      AFI(6)=0
      AFI(7)=0
      AFI(8)=0
      CELPRE(1)=A
      CELPRE(2)=B
      CELPRE(3)=C
      DO 410 I=1,3
      DO 410 J=1,3
410   AL(I,J)=0.0
C
C$OMP CRITICAL(FOUND)
C
      IF(RP(IGC).LT.RMI)THEN
      INTEREST=INTEREST+1
      WRITE(*,1415)RMAX,A,B,C,V2,IPEN
      WRITE(20,*)
      WRITE(20,*)'====================================================
     1==========================='
      WRITE(20,*)
      WRITE(20,*)' YOU HAVE FOUND AN INTERESTING RESULT : Rp < Rmin!'
      WRITE(*,*)
      WRITE(*,*)' YOU HAVE FOUND AN INTERESTING RESULT : Rp < Rmin!'
C
C... Refine that cell
C
      CALL DCELL(CELPRE,AL,V1)
      CALL CALCUL2(DIFF,IHKL,TH3,NCALC,IGC)
      CALL CELREF2(INDIC,BB,AFI,LHKL,TH3,IHKL,DDT,DDQ)
      CNCALC(IGC)=NCALC
      IF(NDAT.GE.20)THEN
      FM20(IGC)=QO(20)/(2.*CNCALC(IGC)*DDQ)
      FF20(IGC)=20./(CNCALC(IGC)*DDT)
      ELSE
      PNDAT=NDAT
      FM20(IGC)=QO(NDAT)/(2.*CNCALC(IGC)*DDQ)
      FF20(IGC)=PNDAT/(CNCALC(IGC)*DDT)
      ENDIF
C      WRITE(20,7000)FM20(IGC)
C      WRITE(20,7001)FF20(IGC),DDT,NCALC
C	WRITE(20,*)
C      PRINT 7000,FM20(IGC)
C      PRINT 7001,FF20(IGC),DDT,NCALC
C	PRINT *
      IREF=1
      GO TO 497
      ELSE
C
C  Anyway, calculate the M20 and F20 values
C
      CALL DCELL(CELPRE,AL,V1)
      CALL CALCUL2(DIFF,IHKL,TH3,NCALC,IGC)
      CALL CELREF2(INDIC,BB,AFI,LHKL,TH3,IHKL,DDT,DDQ)
      CNCALC(IGC)=NCALC
      IF(NDAT.GE.20)THEN
      FM20(IGC)=QO(20)/(2.*CNCALC(IGC)*DDQ)
      FF20(IGC)=20./(CNCALC(IGC)*DDT)
      ELSE
      PNDAT=NDAT
      FM20(IGC)=QO(NDAT)/(2.*CNCALC(IGC)*DDQ)
      FF20(IGC)=PNDAT/(CNCALC(IGC)*DDT)
      ENDIF
      ENDIF
C
C Test if cell already found
C
      IF(IGC.GT.1)THEN
      DO 418 I=1,IGC-1
      IF(IFI(I).NE.IFILE)GO TO 418
      VDELT=VGC(IGC)/300.
      VP=VGC(IGC)+VDELT
      VM=VGC(IGC)-VDELT
      IF(VGC(I).GT.VP.OR.VGC(I).LT.VM)GO TO 418
      ADELT=CEL(1,IGC)/500.
      AP=CEL(1,IGC)+ADELT
      AM=CEL(1,IGC)-ADELT
      NA=0
      IF(CEL(1,I).GT.AP.OR.CEL(1,I).LT.AM)NA=1
      BDELT=CEL(2,IGC)/500.
      BP=CEL(2,IGC)+BDELT
      BM=CEL(2,IGC)-BDELT
      NB=0
      IF(CEL(1,I).GT.BP.OR.CEL(1,I).LT.BM)NB=1
      CDELT=CEL(3,IGC)/500.
      CP=CEL(3,IGC)+CDELT
      CM=CEL(3,IGC)-CDELT
      NC=0
      IF(CEL(1,I).GT.CP.OR.CEL(1,I).LT.CM)NC=1
      IF(NA.EQ.1.AND.NB.EQ.1.AND.NC.EQ.1)GO TO 418
      NSOL(I)=NSOL(I)+1
      IF(RP(IGC).LT.RP(I))THEN
      IF(ISEE.EQ.1)WRITE(*,1415)RMAX,A,B,C,V2,IPEN
      KM(I)=KM(IGC)
      VGC(I)=VGC(IGC)
      RP(I)=RP(IGC)
      CEL(1,I)=CEL(1,IGC)
      CEL(2,I)=CEL(2,IGC)
      CEL(3,I)=CEL(3,IGC)
      ENDIF
      IGC=IGC-1
      IF(NSOL(I).GT.5)THEN
      NTRIED=TMAX+1.
      NOUT=NOUT+1
      ENDIF
      GO TO 419
418   CONTINUE
      IF(IVERB.EQ.1)WRITE(20,415)RMAX,NTRIED,A,B,C,V2,IPEN,ICODE
      IF(ISEE.EQ.1)WRITE(*,1415)RMAX,A,B,C,V2,IPEN
419   CONTINUE
      ELSE
      IF(IVERB.EQ.1)WRITE(20,415)RMAX,NTRIED,A,B,C,V2,IPEN,ICODE
      IF(ISEE.EQ.1)WRITE(*,1415)RMAX,A,B,C,V2,IPEN
      ENDIF
497   CONTINUE
C
C$OMP END CRITICAL(FOUND)
C
417   RMAX=RMAXREF
      IF(RANDI(ISEED).GT.ESCAPE)CELPRE(IP)=CELOLD(IP)
C
C  END ON MC tests
C
496   CONTINUE
C
C$OMP END DO NOWAIT
C$OMP END PARALLEL
C
C... Stop if max limit of Monte Carlo tests outpassed
C         or if K is pressed (tested every 30000 MC event)
C
      IF(INTEREST.GE.1)GO TO 5000
      CALL KILLK(PRESSEDK)
      IF(RMIN.EQ.RMAX)GO TO 498
      IF(IVERB.EQ.1)THEN
      WRITE(20,*)
      WRITE(20,*)'Best result : a = ',BPAR(1),' Rp = ',RMIN
      WRITE(20,*)'Best result : b = ',BPAR(2)
      WRITE(20,*)'Best result : c = ',BPAR(3),'V = ',V3
      WRITE(20,*)
      CALL DATN(DATENOW,TIMENOW)
      ENDIF
      RMIN=RMAX
498   IF(PRESSEDK)GO TO 5000
C
C  END ON NRUNS
C
499   CONTINUE
C
500   IF(NSYS(5).EQ.0)GO TO 600
C
C    Monoclinic case
C
C
      RPSMALL=1.
      PRINT *,
     1'Monoclinic:   Rp     a       b       c       bet     V     Nind'
      IFILE=5
      NCYCLES=2000.
      CY=NCYCLES*1.1
      CELPRE(4)=90.
      CELPRE(6)=90.
C
      IF(NGRID.EQ.3)THEN
      NRUNS=20
      NRUNS2=6
      PMIN=2.
      PMAX=20.
      PMA(1)=DMAX1*2.1
      PMA(2)=DMAX1*2.1
      PMA(3)=DMAX2*2.1
      IF(PMA(1).GT.PMAX)PMA(1)=PMAX
      IF(PMA(2).GT.PMAX)PMA(2)=PMAX
      IF(PMA(3).GT.PMAX)PMA(3)=PMAX
      VMON=PMA(1)*PMA(2)*PMA(3)
      IF(VMON.GT.3000.)VMON=3000.
      DO 523 I=1,3
      PMI(I)=PMIN
      DELTA(I)=(PMA(I)-PMI(I))/2.
523   PSTART(I)=PMI(I)
      ENDIF
C
      WRITE(20,*)
      WRITE(20,*)'Monoclinic Monte Carlo search :'
      WRITE(20,*)' Max(a,b,c), V ',PMA(1),PMA(2),PMA(3),VMAX
      WRITE(20,*)
C
      IF(IVERB.EQ.1)THEN
      WRITE(20,580)NRUN,NTIMELIM(5)
      WRITE(20,*)'====================================================
     1==========================='
      WRITE(20,*)' Rp  Trial number    a   b   c    bet   V  Nind Icod'
      WRITE(20,*)
      ENDIF
580   FORMAT('  Results in monoclinic, run, tests :',I3,F12.0)
C
C     READ hkl Miller indices in mon.hkl
C
      TEMPO='mon.hkl'
      CALL OPEN_READ1(35,TEMPO)
      READ(35,*)NHKL0
      NHKL0=NDAT*20
      IF(NHKL0.GT.1000)NHKL0=1000
      DO 501 I=1,NHKL0
501   READ(35,*)(IHH(KK,I),KK=1,3)
      CLOSE(35)
C
      DO 599 NRUN2=1,NRUNS2
C
      IF(NGRID.EQ.3)THEN
      IF(NRUN2.EQ.1)THEN
      VMIN=8.
      VMAX=500.
      IF(VMON.LT.500.)VMAX=VMON
      NTIMELIM(5)=(VMAX-VMIN)*200.
      ENDIF
      IF(NRUN2.EQ.2)THEN
      IF(VMON.LT.500.)GO TO 600
      VMIN=500.
      VMAX=1000.
      IF(VMON.LT.1000.)VMAX=VMON
      NTIMELIM(5)=(VMAX-VMIN)*200.
      ENDIF
      IF(NRUN2.EQ.3)THEN
      IF(VMON.LT.1000.)GO TO 600
      VMIN=1000.
      VMAX=1500.
      IF(VMON.LT.1500.)VMAX=VMON
      NTIMELIM(5)=(VMAX-VMIN)*200.
      ENDIF
      IF(NRUN2.EQ.4)THEN
      IF(VMON.LT.1500.)GO TO 600
      VMIN=1500.
      VMAX=2000.
      IF(VMON.LT.2000.)VMAX=VMON
      NTIMELIM(5)=(VMAX-VMIN)*200.
      ENDIF
      IF(NRUN2.EQ.5)THEN
      IF(VMON.LT.2000.)GO TO 600
      VMIN=2000.
      VMAX=2500.
      IF(VMON.LT.2500.)VMAX=VMON
      NTIMELIM(5)=(VMAX-VMIN)*200.
      ENDIF
      IF(NRUN2.EQ.6)THEN
      IF(VMON.LT.2500.)GO TO 600
      VMIN=2500.
      VMAX=3000.
      IF(VMON.LT.3000.)VMAX=VMON
      NTIMELIM(5)=(VMAX-VMIN)*200.
      ENDIF
      ENDIF
C
      DO 599 NRUN=1,NRUNS
*------------------------------------------------------------------------- 
*     Initialisation
*
C      CALL ESP_INIT(ISEED)
*
*-------------------------------------------------------------------------
      RMAX=RMAXREF
      RMIN=RMAX
C
C...  here starts the loop
C
      INTEREST=0
      TMAX=NTIMELIM(5)/PROCS
      TTMAX=10.*NTIMELIM(5)
      NCELLS=NTIMELIM(5)
      IISEED=0
      NTRIED=0.
      NTRIEDT=0.
      NOUT=0
      CELPRE(1)=PSTART(1)+2.*DELTA(1)*RANDI(ISEED)
      CELPRE(2)=PSTART(2)+2.*DELTA(2)*RANDI(ISEED)
      CELPRE(3)=PSTART(3)+2.*DELTA(3)*RANDI(ISEED)
      CELPRE(5)=ASTART+2.*DELTC*RANDI(ISEED)
      CELOLD(1)=CELPRE(1)
      CELOLD(2)=CELPRE(2)
      CELOLD(3)=CELPRE(3)
      CELOLD(5)=CELPRE(5)
      RGLOB=1.
      NGLOB=0
C$OMP PARALLEL DEFAULT(SHARED) COPYIN(/CAL/,/CAL2/)
C$OMP& PRIVATE(NCEL,NTRIEDB,DEL,DELD,V1,ICODE,LLHKL,IHKL,TH3,
C$OMP& RMAX2,A,B,C,BET,V2,BPAR,V3,PSTARTB,IPEN,ISEE,INDIC,IP,X,
C$OMP& DIFF,DIFF2,NCALC,DDT,DDQ)
C$OMP& FIRSTPRIVATE(iseed,iiseed,rmax0,ntried,ntriedt,nout,
C$OMP& celpre,celold,rglob,nold,rmin,rmax,bb,afi)
C$OMP DO
C
      DO 596 NCEL=1,NCELLS
      IF(NOUT.GE.1)GO TO 596
      IF(INTEREST.GE.1)GO TO 596
      IISEED=IISEED+1
      IF(IISEED.EQ.1)ISEED=((ISEED-NCEL*NRUN)/2)*2+1
C
502   CONTINUE
C
C     Which parameter to vary ? a or b or c or beta ?
C
      NTRIEDB=0.
      X=RANDI(ISEED)
      IF(X.GE.0.AND.X.LT.0.25)IP=1
      IF(X.GE.0.25.AND.X.LT.0.5)IP=2
      IF(X.GE.0.5.AND.X.LT.0.75)IP=3
      IF(X.GE.0.75.AND.X.LE.1.0)IP=5
      IF(IP.NE.5)THEN
      CELPRE(IP)=PSTART(IP)+2.*DELTA(IP)*RANDI(ISEED)
      ELSE
      CELPRE(IP)=ASTART+2.*DELTC*RANDI(ISEED)
      ENDIF
      NTRIED=NTRIED+1.
      GO TO 504
503   DEL=DELTAB*(1.-NTRIEDB/CY)
      DELD=DELTAD*(1.-NTRIEDB/CY)
      X=RANDI(ISEED)
      IF(X.GE.0.AND.X.LT.0.25)I=1
      IF(X.GE.0.25.AND.X.LT.0.5)I=2
      IF(X.GE.0.5.AND.X.LT.0.75)I=3
      IF(X.GE.0.75.AND.X.LE.1.0)I=5
      IF(I.NE.5)THEN
      CELPRE(I)=PSTARTB(I)+(DEL*(RANDI(ISEED)-0.5)*2.)
      ELSE
      CELPRE(I)=PSTARTB(I)+(DELD*(RANDI(ISEED)-0.5)*2.)
      ENDIF
      NTRIEDB=NTRIEDB+1.
504   CONTINUE
      DO 505 I=1,3
      DO 505 J=1,3
505   AL(I,J)=0.0
      CALL DCELL(CELPRE,AL,V1)
      IF(NTRIED.GT.TMAX)THEN
      NOUT=NOUT+1
      GO TO 596
      ENDIF
      IF(NTRIEDB.NE.0.)GO TO 506
      IF(V1.GT.VMAX.OR.V1.LT.VMIN) THEN
      NTRIED=NTRIED-1.
      NTRIEDT=NTRIEDT+1.
      IF(NTRIEDT.GT.TTMAX)THEN
      NOUT=NOUT+1
      GO TO 596
      ENDIF
      GO TO 502
      ENDIF
C
506   CALL CALCUL1(DIFF,DIFF2)
      IF(NMX.GT.NDAT10)THEN
      NTRIED=NTRIED-1
      GO TO 502
      ENDIF
      IF(NTRIEDB.NE.0.)GO TO 514
C
C... Rp value satisfying ???
C
      IF(DIFF.LT.RGLOB.OR.LHKL.GT.NGLOB)THEN
      RGLOB=DIFF
      NGLOB=LHKL
      CELOLD(IP)=CELPRE(IP)
      ENDIF
      IF(LHKL.GE.NMAX)THEN
      RMAX=DIFF
      ICODE=2
      IF(DIFF.LE.RMAXREF)ICODE=1
      ELSE
      ICODE=1
      ENDIF
      IF(DIFF.GT.RMAX) GO TO 517
      IF(LHKL.LT.NMAX) GO TO 517
514   IF(DIFF.LE.RMAX)THEN
      LLHKL=LHKL
      RMAX=DIFF
      RMAX2=DIFF2
      A=CELPRE(1)
      B=CELPRE(2)
      C=CELPRE(3)
      BET=CELPRE(5)
      V2=V1
      IF(DIFF.LT.RMIN)THEN
      RMIN=DIFF
      BPAR(1)=A
      BPAR(2)=B
      BPAR(3)=C
      BPAR(5)=BET
      V3=V1
      ENDIF
C
C... "Refine" that cell (by Monte Carlo too...)
C
      PSTARTB(1)=CELPRE(1)
      PSTARTB(2)=CELPRE(2)
      PSTARTB(3)=CELPRE(3)
      PSTARTB(5)=CELPRE(5)
      ENDIF
      IF(NTRIEDB.LE.NCYCLES)GO TO 503
      RGLOB=0.5
      NGLOB=NDAT2
      IF(IP.EQ.1)CELOLD(IP)=A
      IF(IP.EQ.2)CELOLD(IP)=B
      IF(IP.EQ.3)CELOLD(IP)=C
      IF(IP.EQ.5)CELOLD(IP)=BET
      NTRIEDB=0.
      IF(RMAX.GE.RMAX0(5))GO TO 517
      IF(RMAX2.GE.0.15)GO TO 517
      IPEN=NDAT-LLHKL
      IF(IPEN.GT.NIND)GO TO 517
C
C$OMP CRITICAL(STORE1)
C
      IGC=IGC+1
C
C  Test if too much proposals, if yes decrease Rmax by 5%
C
      IGT=IGT+1.
      IF(NR.EQ.1)THEN
      IF(IGT.GT.50.)THEN
      IF((NTRIED/IGT).LT.100000.)THEN
      IF(RMAX0(5).GT.0.2)THEN
      RMAX0(5)=RMAX0(5)-RMAX0(5)*0.05
      WRITE(20,*)'  Rmax reduced by 5%, now Rmax = ',RMAX0(5)
      WRITE(*,*)'  Rmax reduced by 5%, now Rmax = ',RMAX0(5)
      ENDIF
      ENDIF
      ENDIF
      ENDIF
C
      IF(IGC.GT.10000)THEN
      WRITE(20,*)'   More than 10000 good cells = STOP'
      PRINT *,'   More than 10000 good cells = STOP'
      IGC=IGC-1
      INTEREST=INTEREST+1
C      GO TO 5000
      ENDIF
      CEL(1,IGC)=A
      CEL(2,IGC)=B
      CEL(3,IGC)=C
      CEL(4,IGC)=90.
      CEL(5,IGC)=BET
      CEL(6,IGC)=90.
C
C$OMP END CRITICAL(STORE1)
C
C... Check for supercell
C
      CELPRE(1)=A
      CELPRE(2)=B
      CELPRE(3)=C
      CELPRE(5)=BET
      DO 540 I=1,3
      DO 540 J=1,3
540   AL(I,J)=0.0
      CALL DCELL(CELPRE,AL,V1)
C
C$OMP CRITICAL(STORE2)
C
      CALL CALCUL2(DIFF,IHKL,TH3,NCALC,IGC)
      KM(IGC)=LLHKL
      KM2(IGC)=LHKL
      IFI(IGC)=IFILE
      NSOL(IGC)=1
      VGC(IGC)=V1
      RP(IGC)=RMAX
      RP2(IGC)=DIFF
      IF(RP(IGC).LT.RPSMALL)THEN
      RPSMALL=RP(IGC)
      ISEE=1
      ELSE
      ISEE=0
      ENDIF
      CALL SUPCEL(LHKL,IHKL,CEL,IGC,VGC,1)
      CALL BRAV(LHKL,IHKL,IBR)
      IB(IGC)=IBR
      A=CEL(1,IGC)
      B=CEL(2,IGC)
      C=CEL(3,IGC)
      V2=VGC(IGC)
C
C$OMP END CRITICAL(STORE2)
C
C... Check for interesting result
C
C      IF(INTEREST.GE.1)GO TO 596
      INDIC=0
      BB(3)=A
      BB(4)=B
      BB(5)=C
      BB(6)=90.
      BB(7)=BET
      BB(8)=90.
      AFI(3)=1
      AFI(4)=1
      AFI(5)=1
      AFI(6)=0
      AFI(7)=1
      AFI(8)=0
      CELPRE(1)=A
      CELPRE(2)=B
      CELPRE(3)=C
      CELPRE(5)=BET
      DO 510 I=1,3
      DO 510 J=1,3
510   AL(I,J)=0.0
C
C$OMP CRITICAL(FOUND)
C
      IF(RP(IGC).LT.RMI)THEN
      INTEREST=INTEREST+1
      WRITE(*,1515)RMAX,A,B,C,BET,V2,IPEN
      WRITE(20,*)
      WRITE(20,*)'====================================================
     1==========================='
      WRITE(20,*)
      WRITE(20,*)' YOU HAVE FOUND AN INTERESTING RESULT : Rp < Rmin!'
      WRITE(*,*)
      WRITE(*,*)' YOU HAVE FOUND AN INTERESTING RESULT : Rp < Rmin!'
C
C... Refine that cell
C
      CALL DCELL(CELPRE,AL,V1)
      CALL CALCUL2(DIFF,IHKL,TH3,NCALC,IGC)
      CALL CELREF2(INDIC,BB,AFI,LHKL,TH3,IHKL,DDT,DDQ)
      CNCALC(IGC)=NCALC
      IF(NDAT.GE.20)THEN
      FM20(IGC)=QO(20)/(2.*CNCALC(IGC)*DDQ)
      FF20(IGC)=20./(CNCALC(IGC)*DDT)
      ELSE
      PNDAT=NDAT
      FM20(IGC)=QO(NDAT)/(2.*CNCALC(IGC)*DDQ)
      FF20(IGC)=PNDAT/(CNCALC(IGC)*DDT)
      ENDIF
C      WRITE(20,7000)FM20(IGC)
C      WRITE(20,7001)FF20(IGC),DDT,NCALC
C	WRITE(20,*)
C      PRINT 7000,FM20(IGC)
C      PRINT 7001,FF20(IGC),DDT,NCALC
C	PRINT *
      IREF=1
      GO TO 597
      ELSE
C
C  Anyway, calculate the M20 and F20 values
C
      CALL DCELL(CELPRE,AL,V1)
      CALL CALCUL2(DIFF,IHKL,TH3,NCALC,IGC)
      CALL CELREF2(INDIC,BB,AFI,LHKL,TH3,IHKL,DDT,DDQ)
      CNCALC(IGC)=NCALC
      IF(NDAT.GE.20)THEN
      FM20(IGC)=QO(20)/(2.*CNCALC(IGC)*DDQ)
      FF20(IGC)=20./(CNCALC(IGC)*DDT)
      ELSE
      PNDAT=NDAT
      FM20(IGC)=QO(NDAT)/(2.*CNCALC(IGC)*DDQ)
      FF20(IGC)=PNDAT/(CNCALC(IGC)*DDT)
      ENDIF
      ENDIF
C
C Test if cell already found
C
C
      IF(IGC.GT.1)THEN
      DO 518 I=1,IGC-1
      IF(IFI(I).NE.IFILE)GO TO 518
      VDELT=VGC(IGC)/300.
      VP=VGC(IGC)+VDELT
      VM=VGC(IGC)-VDELT
      IF(VGC(I).GT.VP.OR.VGC(I).LT.VM)GO TO 518
      BDELT=CEL(2,IGC)/500.
      BP=CEL(2,IGC)+BDELT
      BM=CEL(2,IGC)-BDELT
      IF(CEL(2,I).GT.BP.OR.CEL(2,I).LT.BM)GO TO 518
      BETDELT=CEL(5,IGC)/500.
      BETP=CEL(5,IGC)+BETDELT
      BETM=CEL(5,IGC)-BETDELT
      IF(CEL(5,I).GT.BETP.OR.CEL(5,I).LT.BETM)GO TO 518
      ADELT=CEL(1,IGC)/500.
      AP=CEL(1,IGC)+ADELT
      AM=CEL(1,IGC)-ADELT
      NA=0
      IF(CEL(1,I).GT.AP.OR.CEL(1,I).LT.AM)NA=1
      CDELT=CEL(3,IGC)/500.
      CP=CEL(3,IGC)+CDELT
      CM=CEL(3,IGC)-CDELT
      NC=0
      IF(CEL(1,I).GT.CP.OR.CEL(1,I).LT.CM)NC=1
      IF(NA.EQ.1.AND.NC.EQ.1)GO TO 518
      NSOL(I)=NSOL(I)+1
      IF(RP(IGC).LT.RP(I))THEN
      IF(ISEE.EQ.1)WRITE(*,1515)RMAX,A,B,C,BET,V2,IPEN
      KM(I)=KM(IGC)
      VGC(I)=VGC(IGC)
      RP(I)=RP(IGC)
      CEL(1,I)=CEL(1,IGC)
      CEL(2,I)=CEL(2,IGC)
      CEL(3,I)=CEL(3,IGC)
      CEL(5,I)=CEL(5,IGC)
      ENDIF
      IGC=IGC-1
      IF(NSOL(I).GT.5)THEN
      NTRIED=TMAX+1.
      NOUT=NOUT+1
      ENDIF
      GO TO 519
518   CONTINUE
      IF(IVERB.EQ.1)WRITE(20,515)RMAX,NTRIED,A,B,C,BET,V2,IPEN,ICODE
      IF(ISEE.EQ.1)WRITE(*,1515)RMAX,A,B,C,BET,V2,IPEN
519   CONTINUE
      ELSE
      IF(IVERB.EQ.1)WRITE(20,515)RMAX,NTRIED,A,B,C,BET,V2,IPEN,ICODE
      IF(ISEE.EQ.1)WRITE(*,1515)RMAX,A,B,C,BET,V2,IPEN
      ENDIF
597   CONTINUE
C
C$OMP END CRITICAL(FOUND)
C
517   RMAX=RMAXREF
      IF(RANDI(ISEED).GT.ESCAPE)CELPRE(IP)=CELOLD(IP)
C  
C  END ON MC tests
C
596   CONTINUE
C
C$OMP END DO NOWAIT
C$OMP END PARALLEL
C
C... Stop if max limit of Monte Carlo tests outpassed
C         or if K is pressed (tested every 30000 MC event)
C
      IF(INTEREST.GE.1)GO TO 5000
      CALL KILLK(PRESSEDK)
      IF(RMIN.EQ.RMAX)GO TO 598
      IF(IVERB.EQ.1)THEN
      WRITE(20,*)
      WRITE(20,*)'Best result : a =    ',BPAR(1),' Rp = ',RMIN
      WRITE(20,*)'Best result : b =    ',BPAR(2)
      WRITE(20,*)'Best result : c =    ',BPAR(3)
      WRITE(20,*)'Best result : beta = ',BPAR(5),'  V = ',V3
      WRITE(20,*)
      CALL DATN(DATENOW,TIMENOW)
      ENDIF
      RMIN=RMAX
598   IF(PRESSEDK)GO TO 5000
C
C  END ON NRUNS
C
599   CONTINUE
C
600   IF(NSYS(6).EQ.0)GO TO 700
C
C    Triclinic case
C
C
      RPSMALL=1.
      PRINT *,'Triclinic:    Rp
     1     a       b       c       alp    bet    gam     V     Nind'
      IFILE=6
      NCYCLES=5000.
      CY=NCYCLES*1.1
C
      IF(NGRID.EQ.3)THEN
      NRUNS=20
      NRUNS2=8
      PMIN=2.
      PMAX=20.
      PMA(1)=DMAX1*1.5
      PMA(2)=DMAX2*1.5
      PMA(3)=DMAX3*1.5
      IF(PMA(1).GT.PMAX)PMA(1)=PMAX
      IF(PMA(2).GT.PMAX)PMA(2)=PMAX
      IF(PMA(3).GT.PMAX)PMA(3)=PMAX
      VTRIC=PMA(1)*PMA(2)*PMA(3)
      IF(VTRIC.GT.2000.)VTRIC=2000.
      VMAX=VTRIC
      DO 623 I=1,3
      PMI(I)=PMIN
      DELTA(I)=(PMA(I)-PMI(I))/2.
623   PSTART(I)=PMI(I)
      ENDIF
C
      WRITE(20,*)
      WRITE(20,*)'Triclinic Monte Carlo search :'
      WRITE(20,*)' Max(a,b,c), V ',PMA(1),PMA(2),PMA(3),VTRIC
      WRITE(20,*)
C
      IF(IVERB.EQ.1)THEN
      WRITE(20,680)NRUN,NTIMELIM(6)
      WRITE(20,*)'====================================================
     1==========================='
      WRITE(20,*)' Rp Trial number  a  b  c  alp bet gam  V  Nind icod'
      WRITE(20,*)
      ENDIF
680   FORMAT('  Results in triclinic, run, tests :',I3,F12.0)
C
C     READ hkl Miller indices in tri.hkl
C
      TEMPO='tri.hkl'
      CALL OPEN_READ1(35,TEMPO)
      READ(35,*)NHKL0
      NHKL0=NDAT*20
      IF(NHKL0.GT.1000)NHKL0=1000
      DO 601 I=1,NHKL0
601   READ(35,*)(IHH(KK,I),KK=1,3)
      CLOSE(35)
C
      DO 699 NRUN2=1,NRUNS2
C
      IF(NGRID.EQ.3)THEN
      IF(NRUN2.EQ.1)THEN
      VMIN=8.
      VMAX=250.
      IF(VTRIC.LT.250.)VMAX=VTRIC
      NTIMELIM(6)=(VMAX-VMIN)*4000.
      ENDIF
      IF(NRUN2.EQ.2)THEN
      IF(VTRIC.LT.250.)GO TO 700
      VMIN=250.
      VMAX=500.
      IF(VTRIC.LT.500.)VMAX=VTRIC
      NTIMELIM(6)=(VMAX-VMIN)*4000.
      ENDIF
      IF(NRUN2.EQ.3)THEN
      IF(VTRIC.LT.500.)GO TO 700
      VMIN=500.
      VMAX=750.
      IF(VTRIC.LT.750.)VMAX=VTRIC
      NTIMELIM(6)=(VMAX-VMIN)*4000.
      ENDIF
      IF(NRUN2.EQ.4)THEN
      IF(VTRIC.LT.750.)GO TO 700
      VMIN=750.
      VMAX=1000.
      IF(VTRIC.LT.1000.)VMAX=VTRIC
      NTIMELIM(6)=(VMAX-VMIN)*4000.
      ENDIF
      IF(NRUN2.EQ.5)THEN
      IF(VTRIC.LT.1000.)GO TO 700
      VMIN=1000.
      VMAX=1250.
      IF(VTRIC.LT.1250.)VMAX=VTRIC
      NTIMELIM(6)=(VMAX-VMIN)*4000.
      ENDIF
      IF(NRUN2.EQ.6)THEN
      IF(VTRIC.LT.1250.)GO TO 700
      VMIN=1250.
      VMAX=1500.
      IF(VTRIC.LT.1500.)VMAX=VTRIC
      NTIMELIM(6)=(VMAX-VMIN)*4000.
      ENDIF
      IF(NRUN2.EQ.7)THEN
      IF(VTRIC.LT.1500.)GO TO 700
      VMIN=1500.
      VMAX=1750.
      IF(VTRIC.LT.1750.)VMAX=VTRIC
      NTIMELIM(6)=(VMAX-VMIN)*4000.
      ENDIF
      IF(NRUN2.EQ.8)THEN
      IF(VTRIC.LT.1750.)GO TO 700
      VMIN=1750.
      VMAX=2000.
      IF(VTRIC.LT.2000.)VMAX=VTRIC
      NTIMELIM(6)=(VMAX-VMIN)*4000.
      ENDIF
      ENDIF
C
      DO 699 NRUN=1,NRUNS
*------------------------------------------------------------------------- 
*     Initialisation
*
C      CALL ESP_INIT(ISEED)
*
*-------------------------------------------------------------------------
      RMAX=RMAXREF
      RMIN=RMAX
C
C...  here starts the loop
C
      INTEREST=0
      TMAX=NTIMELIM(6)/PROCS
      TTMAX=10.*NTIMELIM(6)
      NCELLS=NTIMELIM(6)
      IISEED=0
      NTRIED=0.
      NTRIEDT=0.
      NOUT=0
C
C...  here starts the loop
C
      CELPRE(1)=PSTART(1)+2.*DELTA(1)*RANDI(ISEED)
      CELPRE(2)=PSTART(2)+2.*DELTA(2)*RANDI(ISEED)
      CELPRE(3)=PSTART(3)+2.*DELTA(3)*RANDI(ISEED)
      CELPRE(4)=ASTARTT(1)+2.*DELTCT(1)*RANDI(ISEED)
      CELPRE(5)=ASTARTT(2)+2.*DELTCT(2)*RANDI(ISEED)
      CELPRE(6)=ASTARTT(3)+2.*DELTCT(3)*RANDI(ISEED)
      CELOLD(1)=CELPRE(1)
      CELOLD(2)=CELPRE(2)
      CELOLD(3)=CELPRE(3)
      CELOLD(4)=CELPRE(4)
      CELOLD(5)=CELPRE(5)
      CELOLD(6)=CELPRE(6)
      RGLOB=1.
      NGLOB=0
C
C$OMP PARALLEL DEFAULT(SHARED) COPYIN(/CAL/,/CAL2/)
C$OMP& PRIVATE(NCEL,NTRIEDB,DEL,DELD,V1,ICODE,LLHKL,IHKL,TH3,
C$OMP& RMAX2,A,B,C,ALP,BET,GAM,V2,BPAR,V3,PSTARTB,IPEN,ISEE,INDIC,IP,X,
C$OMP& DIFF,DIFF2,IP2,ANG,NCALC,DDT,DDQ)
C$OMP& FIRSTPRIVATE(iseed,iiseed,rmax0,ntried,ntriedt,nout,
C$OMP& celpre,celold,rglob,nglob,rmin,rmax,bb,afi)
C$OMP DO
C
      DO 696 NCEL=1,NCELLS
      IF(NOUT.GE.1)GO TO 696
      IF(INTEREST.GE.1)GO TO 696
      IISEED=IISEED+1
      IF(IISEED.EQ.1)ISEED=((ISEED-NCEL*NRUN)/2)*2+1
C
602   CONTINUE
C
C     Which parameter to vary ? a or b or c or alpha or beta or gamma ?
C
      NTRIEDB=0.
      X=RANDI(ISEED)
      IF(X.GE.0.AND.X.LT.0.16666)IP=1
      IF(X.GE.0.16666.AND.X.LT.0.33333)IP=2
      IF(X.GE.0.33333.AND.X.LT.0.5)IP=3
      IF(X.GE.0.5.AND.X.LT.0.66666)IP=4
      IF(X.GE.0.66666.AND.X.LT.0.83333)IP=5
      IF(X.GE.0.83333.AND.X.LE.1.0)IP=6
      IF(IP.NE.4.AND.IP.NE.5.AND.IP.NE.6)THEN
      CELPRE(IP)=PSTART(IP)+2.*DELTA(IP)*RANDI(ISEED)
      ELSE
      CELPRE(IP)=ASTARTT(IP-3)+2.*DELTCT(IP-3)*RANDI(ISEED)
      ENDIF
      ANG=CELPRE(4)+CELPRE(5)+CELPRE(6)
      IF(ANG.GE.360..AND.ANG.LE.180.)GO TO 696
      NTRIED=NTRIED+1.
      GO TO 604
603   DEL=DELTAB*(1.-NTRIEDB/CY)
      DELD=DELTAD*(1.-NTRIEDB/CY)
      X=RANDI(ISEED)
      IF(X.GE.0.AND.X.LT.0.16666)IP2=1
      IF(X.GE.0.16666.AND.X.LT.0.33333)IP2=2
      IF(X.GE.0.33333.AND.X.LT.0.5)IP2=3
      IF(X.GE.0.5.AND.X.LT.0.66666)IP2=4
      IF(X.GE.0.66666.AND.X.LT.0.83333)IP2=5
      IF(X.GE.0.83333.AND.X.LE.1.0)IP2=6
      IF(IP2.NE.4.AND.IP2.NE.5.AND.IP2.NE.6)THEN
      CELPRE(IP2)=PSTARTB(IP2)+(DEL*(RANDI(ISEED)-0.5)*2.)
      ELSE
      CELPRE(IP2)=PSTARTB(IP2)+(DELD*(RANDI(ISEED)-0.5)*2.)
      ENDIF
      ANG=CELPRE(4)+CELPRE(5)+CELPRE(6)
      IF(ANG.GE.360..AND.ANG.LE.180.)GO TO 603
      NTRIEDB=NTRIEDB+1.
604   CONTINUE
      DO 605 I=1,3
      DO 605 J=1,3
605   AL(I,J)=0.0
      CALL DCELL(CELPRE,AL,V1)
      IF(NTRIED.GT.TMAX)THEN
      NOUT=NOUT+1
      GO TO 696
      ENDIF
      IF(NTRIEDB.NE.0.)GO TO 606
      IF(V1.GT.VMAX.OR.V1.LT.VMIN) THEN
      NTRIED=NTRIED-1.
      NTRIEDT=NTRIEDT+1.
      IF(NTRIEDT.GT.TTMAX)THEN
      NOUT=NOUT+1
      GO TO 696
      ENDIF
      GO TO 602
      ENDIF
C
606   CALL CALCUL1(DIFF,DIFF2)
      IF(NMX.GT.NDAT10)THEN
      NTRIED=NTRIED-1
      GO TO 602
      ENDIF
      IF(NTRIEDB.NE.0.)GO TO 614
C
C... Rp value satisfying ???
C
      IF(DIFF.LT.RGLOB.OR.LHKL.GT.NGLOB)THEN
      RGLOB=DIFF
      NGLOB=LHKL
      CELOLD(IP)=CELPRE(IP)
      ENDIF
      IF(LHKL.GE.NMAX)THEN
      RMAX=DIFF
      ICODE=2
      IF(DIFF.LE.RMAXREF)ICODE=1
      ELSE
      ICODE=1
      ENDIF
      IF(DIFF.GT.RMAX) GO TO 617
      IF(LHKL.LT.NMAX) GO TO 617
614   IF(DIFF.LE.RMAX)THEN
      LLHKL=LHKL
      RMAX=DIFF
      RMAX2=DIFF2
      A=CELPRE(1)
      B=CELPRE(2)
      C=CELPRE(3)
      ALP=CELPRE(4)
      BET=CELPRE(5)
      GAM=CELPRE(6)
      V2=V1
      IF(DIFF.LT.RMIN)THEN
      RMIN=DIFF
      BPAR(1)=A
      BPAR(2)=B
      BPAR(3)=C
      BPAR(4)=ALP
      BPAR(5)=BET
      BPAR(6)=GAM
      V3=V1
      ENDIF
C
C... "Refine" that cell (by Monte Carlo too...)
C
      PSTARTB(1)=CELPRE(1)
      PSTARTB(2)=CELPRE(2)
      PSTARTB(3)=CELPRE(3)
      PSTARTB(4)=CELPRE(4)
      PSTARTB(5)=CELPRE(5)
      PSTARTB(6)=CELPRE(6)
      ENDIF
      IF(NTRIEDB.LE.NCYCLES)GO TO 603
      RGLOB=0.5
      NGLOB=NDAT2
      IF(IP.EQ.1)CELOLD(IP)=A
      IF(IP.EQ.2)CELOLD(IP)=B
      IF(IP.EQ.3)CELOLD(IP)=C
      IF(IP.EQ.4)CELOLD(IP)=ALP
      IF(IP.EQ.5)CELOLD(IP)=BET
      IF(IP.EQ.6)CELOLD(IP)=GAM
      NTRIEDB=0.
      IF(RMAX.GE.RMAX0(6))GO TO 617
      IF(RMAX2.GE.0.15)GO TO 617
      IPEN=NDAT-LLHKL
      IF(IPEN.GT.NIND)GO TO 617
C
C$OMP CRITICAL(STORE1)
C
      IGC=IGC+1
C
C  Test if too much proposals, if yes decrease Rmax by 5%
C
      IGT=IGT+1.
      IF(NR.EQ.1)THEN
      IF(IGT.GT.50.)THEN
      IF((NTRIED/IGT).LT.100000.)THEN
      IF(RMAX0(6).GT.0.2)THEN
      RMAX0(6)=RMAX0(6)-RMAX0(6)*0.05
      WRITE(20,*)'  Rmax reduced by 5%, now Rmax = ',RMAX0(6)
      WRITE(*,*)'  Rmax reduced by 5%, now Rmax = ',RMAX0(6)
      ENDIF
      ENDIF
      ENDIF
      ENDIF
C
      IF(IGC.GT.10000)THEN
      WRITE(20,*)'   More than 10000 good cells = STOP'
      PRINT *,'   More than 10000 good cells = STOP'
      IGC=IGC-1
      INTEREST=INTEREST+1
C      GO TO 5000
      ENDIF
      CEL(1,IGC)=A
      CEL(2,IGC)=B
      CEL(3,IGC)=C
      CEL(4,IGC)=ALP
      CEL(5,IGC)=BET
      CEL(6,IGC)=GAM
C
C$OMP END CRITICAL(STORE1)
C
C... Check for supercell
C
      CELPRE(1)=A
      CELPRE(2)=B
      CELPRE(3)=C
      CELPRE(4)=ALP
      CELPRE(5)=BET
      CELPRE(6)=GAM
      DO 640 I=1,3
      DO 640 J=1,3
640   AL(I,J)=0.0
      CALL DCELL(CELPRE,AL,V1)
C
C$OMP CRITICAL(STORE2)
C
      CALL CALCUL2(DIFF,IHKL,TH3,NCALC,IGC)

      KM(IGC)=LLHKL
      KM2(IGC)=LHKL
      IFI(IGC)=IFILE
      NSOL(IGC)=1
      VGC(IGC)=V1
      RP(IGC)=RMAX
      RP2(IGC)=DIFF
      IF(RP(IGC).LT.RPSMALL)THEN
      RPSMALL=RP(IGC)
      ISEE=1
      ELSE
      ISEE=0
      ENDIF
      CALL SUPCEL(LHKL,IHKL,CEL,IGC,VGC,1)
      CALL BRAV(LHKL,IHKL,IBR)
      IB(IGC)=IBR
      A=CEL(1,IGC)
      B=CEL(2,IGC)
      C=CEL(3,IGC)
      V2=VGC(IGC)
C
C$OMP END CRITICAL(STORE2)
C
C... Check for interesting result
C
C      IF(INTEREST.GE.1)GO TO 696
      INDIC=0
      BB(3)=A
      BB(4)=B
      BB(5)=C
      BB(6)=ALP
      BB(7)=BET
      BB(8)=GAM
      AFI(3)=1
      AFI(4)=1
      AFI(5)=1
      AFI(6)=1
      AFI(7)=1
      AFI(8)=1
      CELPRE(1)=A
      CELPRE(2)=B
      CELPRE(3)=C
      CELPRE(4)=ALP
      CELPRE(5)=BET
      CELPRE(6)=GAM
      DO 610 I=1,3
      DO 610 J=1,3
610   AL(I,J)=0.0
C
C$OMP CRITICAL(FOUND)
C
      IF(RP(IGC).LT.RMI)THEN
      INTEREST=INTEREST+1
      WRITE(*,1615)RMAX,A,B,C,ALP,BET,GAM,V2,IPEN
      WRITE(20,*)
      WRITE(20,*)'====================================================
     1==========================='
      WRITE(20,*)
      WRITE(20,*)' YOU HAVE FOUND AN INTERESTING RESULT : Rp < Rmin!'
      WRITE(*,*)
      WRITE(*,*)' YOU HAVE FOUND AN INTERESTING RESULT : Rp < Rmin!'
C
C... Refine that cell
C
      CALL DCELL(CELPRE,AL,V1)
      CALL CALCUL2(DIFF,IHKL,TH3,NCALC,IGC)
      CALL CELREF2(INDIC,BB,AFI,LHKL,TH3,IHKL,DDT,DDQ)
      CNCALC(IGC)=NCALC
      IF(NDAT.GE.20)THEN
      FM20(IGC)=QO(20)/(2.*CNCALC(IGC)*DDQ)
      FF20(IGC)=20./(CNCALC(IGC)*DDT)
      ELSE
      PNDAT=NDAT
      FM20(IGC)=QO(NDAT)/(2.*CNCALC(IGC)*DDQ)
      FF20(IGC)=PNDAT/(CNCALC(IGC)*DDT)
      ENDIF
C      WRITE(20,7000)FM20(IGC)
C      WRITE(20,7001)FF20(IGC),DDT,NCALC
C	WRITE(20,*)
C      PRINT 7000,FM20(IGC)
C      PRINT 7001,FF20(IGC),DDT,NCALC
C	PRINT *
      IREF=1
      GO TO 697
      ELSE
C
C  Anyway, calculate the M20 and F20 values
C
      CALL DCELL(CELPRE,AL,V1)
      CALL CALCUL2(DIFF,IHKL,TH3,NCALC,IGC)
      CALL CELREF2(INDIC,BB,AFI,LHKL,TH3,IHKL,DDT,DDQ)
      CNCALC(IGC)=NCALC
      IF(NDAT.GE.20)THEN
      FM20(IGC)=QO(20)/(2.*CNCALC(IGC)*DDQ)
      FF20(IGC)=20./(CNCALC(IGC)*DDT)
      ELSE
      PNDAT=NDAT
      FM20(IGC)=QO(NDAT)/(2.*CNCALC(IGC)*DDQ)
      FF20(IGC)=PNDAT/(CNCALC(IGC)*DDT)
      ENDIF
      ENDIF
C
C Test if cell already found
C
C
      IF(IGC.GT.1)THEN
      DO 618 I=1,IGC-1
      IF(IFI(I).NE.IFILE)GO TO 618
      VDELT=VGC(IGC)/300.
      VP=VGC(IGC)+VDELT
      VM=VGC(IGC)-VDELT
      IF(VGC(I).GT.VP.OR.VGC(I).LT.VM)GO TO 618
      ADELT=CEL(1,IGC)/500.
      AP=CEL(1,IGC)+ADELT
      AM=CEL(1,IGC)-ADELT
      NA=0
      IF(CEL(1,I).GT.AP.OR.CEL(1,I).LT.AM)NA=1
      BDELT=CEL(2,IGC)/500.
      BP=CEL(2,IGC)+BDELT
      BM=CEL(2,IGC)-BDELT
      NB=0
      IF(CEL(1,I).GT.BP.OR.CEL(1,I).LT.BM)NB=1
      CDELT=CEL(3,IGC)/500.
      CP=CEL(3,IGC)+CDELT
      CM=CEL(3,IGC)-CDELT
      NC=0
      IF(CEL(1,I).GT.CP.OR.CEL(1,I).LT.CM)NC=1
      IF(NA.EQ.1.AND.NB.EQ.1.AND.NC.EQ.1)GO TO 618
      NA=0
      IF(CEL(2,I).GT.AP.OR.CEL(2,I).LT.AM)NA=1
      NB=0
      IF(CEL(2,I).GT.BP.OR.CEL(2,I).LT.BM)NB=1
      NC=0
      IF(CEL(2,I).GT.CP.OR.CEL(2,I).LT.CM)NC=1
      IF(NA.EQ.1.AND.NB.EQ.1.AND.NC.EQ.1)GO TO 618
      NSOL(I)=NSOL(I)+1
      IF(RP(IGC).LT.RP(I))THEN
      IF(ISEE.EQ.1)
C     1WRITE(*,1615)RMAX,A,B,C,ALP,BET,GAM,V2,IPEN
     1WRITE(*,1615)RMAX,A,B,C,ALP,BET,GAM,V2,IPEN
      KM(I)=KM(IGC)
      VGC(I)=VGC(IGC)
      RP(I)=RP(IGC)
      CEL(1,I)=CEL(1,IGC)
      CEL(2,I)=CEL(2,IGC)
      CEL(3,I)=CEL(3,IGC)
      CEL(4,I)=CEL(4,IGC)
      CEL(5,I)=CEL(5,IGC)
      CEL(6,I)=CEL(6,IGC)
      ENDIF
      IGC=IGC-1
      IF(NSOL(I).GT.5)THEN
      NTRIED=TMAX+1.
      NOUT=NOUT+1
      ENDIF
      GO TO 619
618   CONTINUE
      IF(IVERB.EQ.1)
     1WRITE(20,615)RMAX,NTRIED,A,B,C,ALP,BET,GAM,V2,IPEN,ICODE
      IF(ISEE.EQ.1)
     1WRITE(*,1615)RMAX,A,B,C,ALP,BET,GAM,V2,IPEN
619   CONTINUE
      ELSE
      IF(IVERB.EQ.1)
     1WRITE(20,615)RMAX,NTRIED,A,B,C,ALP,BET,GAM,V2,IPEN,ICODE
      IF(ISEE.EQ.1)
     1WRITE(*,1615)RMAX,A,B,C,ALP,BET,GAM,V2,IPEN
      ENDIF
697   CONTINUE
C
C$OMP END CRITICAL(FOUND)
C
617   RMAX=RMAXREF
      IF(RANDI(ISEED).GT.ESCAPE)CELPRE(IP)=CELOLD(IP)
C
C  END ON MC tests
C
696   CONTINUE
C
C$OMP END DO NOWAIT
C$OMP END PARALLEL
C
C... Stop if max limit of Monte Carlo tests outpassed
C         or if K is pressed (tested every 30000 MC event)
C
      IF(INTEREST.GE.1)GO TO 5000
      CALL KILLK(PRESSEDK)
616   CONTINUE
      IF(RMIN.EQ.RMAX)GO TO 698
      IF(IVERB.EQ.1)THEN
      WRITE(20,*)
      WRITE(20,*)'Best result : a =    ',BPAR(1),' Rp = ',RMIN
      WRITE(20,*)'Best result : b =    ',BPAR(2)
      WRITE(20,*)'Best result : c =    ',BPAR(3)
      WRITE(20,*)'Best result : alph = ',BPAR(4)
      WRITE(20,*)'Best result : beta = ',BPAR(5)
      WRITE(20,*)'Best result : gamm = ',BPAR(6),'  V = ',V3
      WRITE(20,*)
      CALL DATN(DATENOW,TIMENOW)
      ENDIF
      RMIN=RMAX
698   IF(PRESSEDK)GO TO 5000
C
C  END ON NRUNS
C
699   CONTINUE
C
C
C   SECOND PART OF McMAILLE : GRID SEARCH...
C
700   IF(NGRID.EQ.0)  GO TO 5000
      IF(NBLACK.EQ.0.AND.NGRID.EQ.3) GO TO 5000
      PRINT *
      PRINT *,'Grid search :'
C
C...  Cell generation by systematic grid
C
C
      IF(NSYS(1).EQ.0)GO TO 1200
C
C    Cubic case
C
      RPSMALL=1.
      PRINT *,'Cubic:        Rp     a       V     Nind'
*------------------------------------------------------------------------- 
*     Initialisation
*
C      CALL ESP_INIT(ISEED)
*
*-------------------------------------------------------------------------
      IFILE=1
      RMAX=RMAXREF
      RMIN=RMAX
      NTRIED=0.
      NCYCLES=200.
      CY=NCYCLES*1.1
      CELPRE(4)=90.
      CELPRE(5)=90.
      CELPRE(6)=90.
      IF(IVERB.EQ.1)THEN
      WRITE(20,*)
      WRITE(20,*)'====================================================
     1==========================='
      WRITE(20,*)'Grid search :'
      WRITE(20,*)'   Results in cubic :'
      WRITE(20,*)'====================================================
     1==========================='
      ENDIF
C
      IF(NGRID.EQ.3)THEN
      PMIN=2.
      PMAX=DMAX1*3.1
      PMI(1)=PMIN
      PMA(1)=PMAX
      VMIN=8.
      VMAX=PMAX*PMAX*PMAX
      IF(IVERB.EQ.1)WRITE(20,*)' Max a, V ',PMAX,VMAX
      ENDIF
C
      IF(IVERB.EQ.1)THEN
      WRITE(20,*)' Rp  Trial number    a         V  Nind Icod'
      WRITE(20,*)
      ENDIF
C
C     READ hkl Miller indices in cub.hkl
C
      TEMPO='cub.hkl'
      CALL OPEN_READ1(35,TEMPO)
      READ(35,*)NHKL0
      NHKL0=NDAT*6
      IF(NHKL0.GT.400)NHKL0=400
      DO 1101 I=1,NHKL0
1101  READ(35,*)(IHH(KK,I),KK=1,3)
      CLOSE(35)
C
C...  here starts the loop
C
      CELPRE(1)=PMI(1)-SPAR
1102  CONTINUE
      NTRIEDB=0.
      CELPRE(1)=CELPRE(1)+SPAR
      NTRIED=NTRIED+1.
      GO TO 1104
1103  DEL=DELTAB*(1.-NTRIEDB/CY)
      CELPRE(1)=PSTARTB(1)+(DEL*(RANDI(ISEED)-0.5)*2.)
      NTRIEDB=NTRIEDB+1.
1104  CONTINUE
      CELPRE(2)=CELPRE(1)
      CELPRE(3)=CELPRE(1)
      DO 1105 I=1,3
      DO 1105 J=1,3
1105  AL(I,J)=0.0
      CALL DCELL(CELPRE,AL,V1)
      IF(CELPRE(1).GT.PMA(1).AND.NTRIEDB.EQ.0.)GO TO 1116
      IF(NTRIEDB.NE.0.)GO TO 1106
      IF(V1.GT.VMAX.OR.V1.LT.VMIN) THEN
      NTRIED=NTRIED-1.
      GO TO 1102
      ENDIF
C
1106  CALL CALCUL1(DIFF,DIFF2)
      IF(NMX.GT.NDAT10)THEN
      NTRIED=NTRIED-1
      GO TO 1102
      ENDIF
      IF(NTRIEDB.NE.0.)GO TO 1114
C
C... Rp value satisfying ???
C
      IF(LHKL.GE.NMAX)THEN
      RMAX=DIFF
      ICODE=2
      IF(DIFF.LE.RMAXREF)ICODE=1
      ELSE
      ICODE=1
      ENDIF
      CELOLD(1)=CELPRE(1)
      IF(DIFF.GT.RMAX) GO TO 1117
      IF(LHKL.LT.NMAX) GO TO 1117
1114  IF(DIFF.LE.RMAX)THEN
      LLHKL=LHKL
      RMAX=DIFF
      RMAX2=DIFF2
      A=CELPRE(1)
      V2=V1
      IF(DIFF.LT.RMIN)THEN
      RMIN=DIFF
      BPAR(1)=A
      V3=V1
      ENDIF
C
C... "Refine" that cell (by Monte Carlo too...)
C
      PSTARTB(1)=CELPRE(1)
      ENDIF
      IF(NTRIEDB.LE.NCYCLES)GO TO 1103
      NTRIEDB=0.
      IF(RMAX.GE.RMAX0(1))GO TO 1117
      IF(RMAX2.GE.0.15)GO TO 1117
      IPEN=NDAT-LLHKL
      IF(IPEN.GT.NIND)GO TO 1117
      IGC=IGC+1
C
C  Test if too much proposals, if yes decrease Rmax by 5%
C
      IGT=IGT+1.
      IF(NR.EQ.1)THEN
      IF(IGT.GT.50.)THEN
      IF((NTRIED/IGT).LT.100.)THEN
      IF(RMAX0(1).GT.0.1)THEN
      RMAX0(1)=RMAX0(1)-RMAX0(1)*0.05
      WRITE(20,*)'  Rmax reduced by 5%, now Rmax = ',RMAX0(1)
      WRITE(*,*)'  Rmax reduced by 5%, now Rmax = ',RMAX0(1)
      ENDIF
      ENDIF
      ENDIF
      ENDIF
C
      IF(IGC.GT.10000)THEN
      WRITE(20,*)'   More than 10000 good cells = STOP'
      PRINT *,'   More than 10000 good cells = STOP'
      IGC=IGC-1
      GO TO 5000
      ENDIF
      CEL(1,IGC)=A
      CEL(2,IGC)=A
      CEL(3,IGC)=A
      CEL(4,IGC)=90.
      CEL(5,IGC)=90.
      CEL(6,IGC)=90.
C
C... Check for supercell
C
      CELPRE(1)=A
      CELPRE(2)=A
      CELPRE(3)=A
      DO 1140 I=1,3
      DO 1140 J=1,3
1140  AL(I,J)=0.0
      CALL DCELL(CELPRE,AL,V1)
      CALL CALCUL2(DIFF,IHKL,TH3,NCALC,IGC)
      KM(IGC)=LLHKL
      KM2(IGC)=LHKL
      IFI(IGC)=IFILE
      NSOL(IGC)=1
      VGC(IGC)=V1
      RP(IGC)=RMAX
      RP2(IGC)=DIFF
      IF(RP(IGC).LT.RPSMALL)THEN
      RPSMALL=RP(IGC)
      ISEE=1
      ELSE
      ISEE=0
      ENDIF
      CALL SUPCEL(LHKL,IHKL,CEL,IGC,VGC,3)
      CALL BRAV(LHKL,IHKL,IBR)
      IB(IGC)=IBR
      A=CEL(1,IGC)
      CEL(2,IGC)=A
      CEL(3,IGC)=A
      V2=VGC(IGC)
C
C... Check for interesting result
C
      IF(RP(IGC).LT.RMI)THEN
      WRITE(*,1115)RMAX,A,V2,IPEN
      WRITE(20,*)
      WRITE(20,*)'====================================================
     1==========================='
      WRITE(20,*)
      WRITE(20,*)' YOU HAVE FOUND AN INTERESTING RESULT : Rp < Rmin!'
      WRITE(*,*)
      WRITE(*,*)' YOU HAVE FOUND AN INTERESTING RESULT : Rp < Rmin!'
C
C... Refine that cell
C
      INDIC=1
      BB(3)=A
      BB(4)=A
      BB(5)=A
      BB(6)=90.
      BB(7)=90.
      BB(8)=90.
      AFI(3)=1
      AFI(4)=1
      AFI(5)=1
      AFI(6)=0
      AFI(7)=0
      AFI(8)=0
      CELPRE(1)=A
      CELPRE(2)=A
      CELPRE(3)=A
      DO 1110 I=1,3
      DO 1110 J=1,3
1110  AL(I,J)=0.0
      CALL DCELL(CELPRE,AL,V1)
      CALL CALCUL2(DIFF,IHKL,TH3,NCALC,IGC)
      CALL CELREF(INDIC,BB,AFI,LHKL,TH3,IHKL,DDT,DDQ)
      IF(NDAT.GE.20)THEN
      CNCALC(IGC)=NCALC
      FM20(IGC)=QO(20)/(2.*CNCALC(IGC)*DDQ)
      FF20(IGC)=20./(CNCALC(IGC)*DDT)
      WRITE(20,7000)FM20(IGC)
      WRITE(20,7001)FF20(IGC),DDT,NCALC
      WRITE(20,*)
      PRINT 7000,FM20(IGC)
      PRINT 7001,FF20(IGC),DDT,NCALC
      PRINT *
      ENDIF
      IREF=1
      GO TO 5000
      ENDIF
C
C Test if cell already found
C
      IF(IGC.GT.1)THEN
      DO 1118 I=1,IGC-1
      IF(IFI(I).NE.IFILE)GO TO 1118
      VDELT=VGC(IGC)/300.
      VP=VGC(IGC)+VDELT
      VM=VGC(IGC)-VDELT
      IF(VGC(I).GT.VP.OR.VGC(I).LT.VM)GO TO 1118
      NSOL(I)=NSOL(I)+1
      IF(RP(IGC).LT.RP(I))THEN
      IF(ISEE.EQ.1)WRITE(*,1115)RMAX,A,V2,IPEN
      KM(I)=KM(IGC)
      VGC(I)=VGC(IGC)
      RP(I)=RP(IGC)
      CEL(1,I)=CEL(1,IGC)
      CEL(2,I)=CEL(2,IGC)
      CEL(3,I)=CEL(3,IGC)
      ENDIF
      IGC=IGC-1
      GO TO 1119
1118  CONTINUE
      IF(IVERB.EQ.1)WRITE(20,115)RMAX,NTRIED,A,V2,IPEN,ICODE
      IF(ISEE.EQ.1)WRITE(*,1115)RMAX,A,V2,IPEN
1119  CONTINUE
      ELSE
      IF(IVERB.EQ.1)WRITE(20,115)RMAX,NTRIED,A,V2,IPEN,ICODE
      IF(ISEE.EQ.1)WRITE(*,1115)RMAX,A,V2,IPEN
      ENDIF
C
C
C
1117  RMAX=RMAXREF
      CELPRE(1)=CELOLD(1)
C
C... Stop if max limit of grid tests outpassed
C         or if K is pressed (tested every 30000 MC event)
C
      IF(CELPRE(1).GT.PMA(1))GO TO 1116
      TKILL=MOD(NTRIED,30000.)
      IF(TKILL.GE.0.)THEN
      CALL KILLK(PRESSEDK)
      IF(PRESSEDK)GO TO 1116
      ENDIF
      GO TO 1102
1116  CONTINUE
      IF(RMIN.EQ.RMAX)GO TO 1198
      IF(IVERB.EQ.1)THEN
      WRITE(20,*)
      WRITE(20,*)'Best result : a=',BPAR(1),'V=',V3,' Rp= ',RMIN
      WRITE(20,*)
      CALL DATN(DATENOW,TIMENOW)
      ENDIF
      RMIN=RMAX
1198  IF(PRESSEDK)GO TO 5000
C
C
1200  IF(NSYS(2).EQ.0)GO TO 1300
C
C
C    Hexagonal case
C
C
      IHR=0
      IF(NGRID.EQ.3)IHR=1
1290  CONTINUE
      IF(IHR.EQ.2)NSYS(2)=2
      IF(IVERB.EQ.1)THEN
      IF(NSYS(2).EQ.1)PRINT *,
     1'Hexagonal:    Rp     a      c       V     Nind'
      IF(NSYS(2).EQ.2)PRINT *,
     1'Rhombohedral: Rp     a      c       V     Nind'
      ENDIF
      RPSMALL=1.
*------------------------------------------------------------------------- 
*     Initialisation
*
C      CALL ESP_INIT(ISEED)
*
*-------------------------------------------------------------------------
      IFILE=2
      RMAX=RMAXREF
      RMIN=RMAX
      NTRIED=0.
      NCYCLES=500.
      CY=NCYCLES*1.1
      CELPRE(4)=90.
      CELPRE(5)=90.
      CELPRE(6)=120.
      IF(IVERB.EQ.1)THEN
      WRITE(20,*)
      WRITE(20,*)'====================================================
     1==========================='
      WRITE(20,*)'Grid search :'
      IF(NSYS(2).EQ.1)WRITE(20,*)'   Results in hexagonal'
      IF(NSYS(2).EQ.2)WRITE(20,*)'   Results in rhombohedral'
      WRITE(20,*)'====================================================
     1==========================='
      ENDIF
C
      IF(NGRID.EQ.3)THEN
      PMIN=2.
      PMAX=30.
      PMI(1)=PMIN
      PMI(3)=PMIN
      PMA(3)=DMAX1*6.1
      IF(PMA(3).GT.PMAX)PMA(3)=PMAX
      PMA(1)=DMAX1*2.1
      IF(PMA(1).GT.PMAX)PMA(1)=PMAX
      PMA(2)=PMA(1)
      VMIN=8.
      VMAX=PMA(1)*PMA(2)*PMA(3)
      IF(VMAX.GT.4000.)VMAX=4000.
      IF(IVERB.EQ.1)WRITE(20,*)' Max(a,c), V ',PMA(1),PMA(3),VMAX
      ENDIF
C
      IF(IVERB.EQ.1)THEN
      WRITE(20,*)' Rp  Trial number    a      c        V  Nind Icod'
      WRITE(20,*)
      ENDIF
C
C     READ hkl Miller indices in hex.hkl
C
      IF(NSYS(2).EQ.2)GO TO 1260
      TEMPO='hex.hkl'
      GO TO 1261
1260  TEMPO='rho.hkl'
      IFILE=7
1261  CALL OPEN_READ1(35,TEMPO)
      READ(35,*)NHKL0
      IF(NSYS(2).EQ.2)GO TO 1262
      NHKL0=NDAT*12
      IF(NHKL0.GT.800)NHKL0=800
      GO TO 1263
1262  NHKL0=NDAT*12
      IF(NHKL0.GT.600)NHKL0=600
1263  DO 1201 I=1,NHKL0
1201  READ(35,*)(IHH(KK,I),KK=1,3)
      CLOSE(35)
C
C...  here starts the loop
C
      CELPRE(1)=PMI(1)-SPAR
      CELPRE(2)=CELPRE(1)
      CELPRE(3)=PMI(3)-SPAR
      IFIN=1
C
1202  CONTINUE
C
C     Which parameter to vary ? a or c ?
C
      NTRIEDB=0.
      IF(IFIN.EQ.1)THEN
      CELPRE(1)=CELPRE(1)+SPAR
      IF(CELPRE(1).GT.PMA(1))GO TO 1216
      IFIN=0
      NTRIED=NTRIED+1.
      ENDIF
      CELPRE(3)=CELPRE(3)+SPAR
      IF(CELPRE(3).GT.PMA(3))THEN
      CELPRE(3)=PMI(3)-SPAR
      IFIN=1
      GO TO 1202
      ENDIF
      NTRIED=NTRIED+1.
      GO TO 1204
1203  DEL=DELTAB*(1.-NTRIEDB/CY)
      I=3
      IF(RANDI(ISEED).GT.0.5)I=1
      CELPRE(I)=PSTARTB(I)+(DEL*(RANDI(ISEED)-0.5)*2.)
      NTRIEDB=NTRIEDB+1.
1204  CONTINUE
      CELPRE(2)=CELPRE(1)
      DO 1205 I=1,3
      DO 1205 J=1,3
1205  AL(I,J)=0.0
      CALL DCELL(CELPRE,AL,V1)
      IF(CELPRE(1).GT.PMA(1).AND.NTRIEDB.EQ.0.)GO TO 1216
      IF(NTRIEDB.NE.0.)GO TO 1206
      IF(V1.GT.VMAX.OR.V1.LT.VMIN) THEN
      NTRIED=NTRIED-1.
      GO TO 1202
      ENDIF
C
1206  CALL CALCUL1(DIFF,DIFF2)
      IF(NMX.GT.NDAT10)THEN
      NTRIED=NTRIED-1
      GO TO 1202
      ENDIF
      IF(NTRIEDB.NE.0.)GO TO 1214
C
C... Rp value satisfying ???
C
      IF(LHKL.GE.NMAX)THEN
      RMAX=DIFF
      ICODE=2
      IF(DIFF.LE.RMAXREF)ICODE=1
      ELSE
      ICODE=1
      ENDIF
      CELOLD(1)=CELPRE(1)
      CELOLD(3)=CELPRE(3)
      IF(DIFF.GT.RMAX) GO TO 1217
      IF(LHKL.LT.NMAX) GO TO 1217
1214  IF(DIFF.LE.RMAX)THEN
      LLHKL=LHKL
      RMAX=DIFF
      RMAX2=DIFF2
      A=CELPRE(1)
      C=CELPRE(3)
      V2=V1
      IF(DIFF.LT.RMIN)THEN
      RMIN=DIFF
      BPAR(1)=A
      BPAR(3)=C
      V3=V1
      ENDIF
C
C... "Refine" that cell (by Monte Carlo too...)
C
      PSTARTB(1)=CELPRE(1)
      PSTARTB(3)=CELPRE(3)
      ENDIF
      IF(NTRIEDB.LE.NCYCLES)GO TO 1203
      NTRIEDB=0.
      IF(RMAX.GE.RMAX0(2))GO TO 1217
      IF(RMAX2.GE.0.15)GO TO 1217
      IPEN=NDAT-LLHKL
      IF(IPEN.GT.NIND)GO TO 1217
      IGC=IGC+1
C
C  Test if too much proposals, if yes decrease Rmax by 5%
C
      IGT=IGT+1.
      IF(NR.EQ.1)THEN
      IF(IGT.GT.50.)THEN
      IF((NTRIED/IGT).LT.1000.)THEN
      IF(RMAX0(2).GT.0.1)THEN
      RMAX0(2)=RMAX0(2)-RMAX0(2)*0.05
      WRITE(20,*)'  Rmax reduced by 5%, now Rmax = ',RMAX0(2)
      WRITE(*,*)'  Rmax reduced by 5%, now Rmax = ',RMAX0(2)
      ENDIF
      ENDIF
      ENDIF
      ENDIF
C
      IF(IGC.GT.10000)THEN
      WRITE(20,*)'   More than 10000 good cells = STOP'
      PRINT *,'   More than 10000 good cells = STOP'
      IGC=IGC-1
      GO TO 5000
      ENDIF
      CEL(1,IGC)=A
      CEL(2,IGC)=A
      CEL(3,IGC)=C
      CEL(4,IGC)=90.
      CEL(5,IGC)=90.
      CEL(6,IGC)=120.
C
C... Check for supercell
C
      CELPRE(1)=A
      CELPRE(2)=A
      CELPRE(3)=C
      DO 1240 I=1,3
      DO 1240 J=1,3
1240  AL(I,J)=0.0
      CALL DCELL(CELPRE,AL,V1)
      CALL CALCUL2(DIFF,IHKL,TH3,NCALC,IGC)
      KM(IGC)=LLHKL
      KM2(IGC)=LHKL
      IFI(IGC)=IFILE
      NSOL(IGC)=1
      VGC(IGC)=V1
      RP(IGC)=RMAX
      RP2(IGC)=DIFF
      IF(RP(IGC).LT.RPSMALL)THEN
      RPSMALL=RP(IGC)
      ISEE=1
      ELSE
      ISEE=0
      ENDIF
      CALL SUPCEL(LHKL,IHKL,CEL,IGC,VGC,2)
      CALL BRAV(LHKL,IHKL,IBR)
      IB(IGC)=IBR
      A=CEL(1,IGC)
      CEL(2,IGC)=A
      C=CEL(3,IGC)
      V2=VGC(IGC)
C
C... Check for interesting result
C
      IF(RP(IGC).LT.RMI)THEN
      WRITE(*,1215)RMAX,A,C,V2,IPEN
      WRITE(20,*)
      WRITE(20,*)'====================================================
     1==========================='
      WRITE(20,*)
      WRITE(20,*)' YOU HAVE FOUND AN INTERESTING RESULT : Rp < Rmin!'
      WRITE(*,*)
      WRITE(*,*)' YOU HAVE FOUND AN INTERESTING RESULT : Rp < Rmin!'
C
C... Refine that cell
C
      INDIC=2
      BB(3)=A
      BB(4)=A
      BB(5)=C
      BB(6)=90.
      BB(7)=90.
      BB(8)=120.
      AFI(3)=1
      AFI(4)=1
      AFI(5)=1
      AFI(6)=0
      AFI(7)=0
      AFI(8)=0
      CELPRE(1)=A
      CELPRE(2)=A
      CELPRE(3)=C
      DO 1210 I=1,3
      DO 1210 J=1,3
1210  AL(I,J)=0.0
      CALL DCELL(CELPRE,AL,V1)
      CALL CALCUL2(DIFF,IHKL,TH3,NCALC,IGC)
      CALL CELREF(INDIC,BB,AFI,LHKL,TH3,IHKL,DDT,DDQ)
      IF(NDAT.GE.20)THEN
      CNCALC(IGC)=NCALC
      FM20(IGC)=QO(20)/(2.*CNCALC(IGC)*DDQ)
      FF20(IGC)=20./(CNCALC(IGC)*DDT)
      WRITE(20,7000)FM20(IGC)
      WRITE(20,7001)FF20(IGC),DDT,NCALC
      WRITE(20,*)
      PRINT 7000,FM20(IGC)
      PRINT 7001,FF20(IGC),DDT,NCALC
      PRINT *
      ENDIF
      IREF=1
      GO TO 5000
      ENDIF
C
C Test if cell already found
C
      IF(IGC.GT.1)THEN
      DO 1218 I=1,IGC-1
      IF(IFI(I).NE.IFILE)GO TO 1218
      VDELT=VGC(IGC)/300.
      VP=VGC(IGC)+VDELT
      VM=VGC(IGC)-VDELT
      IF(VGC(I).GT.VP.OR.VGC(I).LT.VM)GO TO 1218
      ADELT=CEL(1,IGC)/500.
      AP=CEL(1,IGC)+ADELT
      AM=CEL(1,IGC)-ADELT
      IF(CEL(1,I).GT.AP.OR.CEL(1,I).LT.AM)GO TO 1218
      NSOL(I)=NSOL(I)+1
      IF(RP(IGC).LT.RP(I))THEN
      IF(ISEE.EQ.1)WRITE(*,1215)RMAX,A,C,V2,IPEN
      KM(I)=KM(IGC)
      VGC(I)=VGC(IGC)
      RP(I)=RP(IGC)
      CEL(1,I)=CEL(1,IGC)
      CEL(2,I)=CEL(2,IGC)
      CEL(3,I)=CEL(3,IGC)
      ENDIF
      IGC=IGC-1
      GO TO 1219
1218  CONTINUE
      IF(IVERB.EQ.1)WRITE(20,215)RMAX,NTRIED,A,C,V2,IPEN,ICODE
      IF(ISEE.EQ.1)WRITE(*,1215)RMAX,A,C,V2,IPEN
1219  CONTINUE
      ELSE
      IF(IVERB.EQ.1)WRITE(20,215)RMAX,NTRIED,A,C,V2,IPEN,ICODE
      IF(ISEE.EQ.1)WRITE(*,1215)RMAX,A,C,V2,IPEN
      ENDIF
C
C
1217  RMAX=RMAXREF
      CELPRE(1)=CELOLD(1)
      CELPRE(3)=CELOLD(3)
C
C... Stop if max limit of grid Carlo tests outpassed
C         or if K is pressed (tested every 30000 MC event)
C
      IF(CELPRE(1).GT.PMA(1))GO TO 1216
      TKILL=MOD(NTRIED,30000.)
      IF(TKILL.GE.0.)THEN
      CALL KILLK(PRESSEDK)
      IF(PRESSEDK)GO TO 1216
      ENDIF
      GO TO 1202
1216  CONTINUE
      IF(RMIN.EQ.RMAX)GO TO 1298
      IF(IVERB.EQ.1)THEN
      WRITE(20,*)
      WRITE(20,*)'Best result : a = ',BPAR(1),' Rp = ',RMIN
      WRITE(20,*)'Best result : c = ',BPAR(3),'V = ',V3
      WRITE(20,*)
      CALL DATN(DATENOW,TIMENOW)
      ENDIF
      RMIN=RMAX
1298  IF(PRESSEDK)GO TO 5000
C
C
      IHR=IHR+1
      IF(IHR.EQ.2)GO TO 1290
1300  IF(NSYS(3).EQ.0)GO TO 1400
C
C    Tetragonal case
C
C
      RPSMALL=1.
      PRINT *,'Tetragonal:   Rp     a       c        V     Nind'
*------------------------------------------------------------------------- 
*     Initialisation
*
C      CALL ESP_INIT(ISEED)
*
*-------------------------------------------------------------------------
      IFILE=3
      RMAX=RMAXREF
      RMIN=RMAX
      NTRIED=0.
      NCYCLES=500.
      CY=NCYCLES*1.1
      CELPRE(4)=90.
      CELPRE(5)=90.
      CELPRE(6)=90.
      IF(IVERB.EQ.1)THEN
      WRITE(20,*)
      WRITE(20,*)'====================================================
     1==========================='
      WRITE(20,*)'Grid search :'
      WRITE(20,*)'   Results in tetragonal :'
      WRITE(20,*)'====================================================
     1==========================='
      ENDIF
C
      IF(NGRID.EQ.3)THEN
      PMIN=2.
      PMAX=30.
      PMI(1)=PMIN
      PMI(3)=PMIN
      PMA(1)=DMAX1*2.1
      IF(PMA(1).GT.PMAX)PMA(1)=PMAX
      PMA(3)=DMAX1*4.
      IF(PMA(3).GT.PMAX)PMA(3)=PMAX
      VMIN=8.
      VMAX=PMA(1)*PMA(2)*PMA(3)
      IF(VMAX.GT.4000.)VMAX=4000.
      IF(IVERB.EQ.1)WRITE(20,*)' Max(a,c), V ',PMA(1),PMA(3),VMAX
      ENDIF
C
      IF(IVERB.EQ.1)THEN
      WRITE(20,*)' Rp  Trial number    a      c        V  Nind Icod'
      WRITE(20,*)
      ENDIF
C
C     READ hkl Miller indices in tet.hkl
C
      TEMPO='tet.hkl'
      CALL OPEN_READ1(35,TEMPO)
      READ(35,*)NHKL0
      NHKL0=NDAT*12
      IF(NHKL0.GT.800)NHKL0=800
      DO 1301 I=1,NHKL0
1301  READ(35,*)(IHH(KK,I),KK=1,3)
      CLOSE(35)
C
C...  here starts the loop
C
      CELPRE(1)=PMI(1)-SPAR
      CELPRE(2)=CELPRE(1)
      CELPRE(3)=PMI(3)-SPAR
      IFIN=1
C
1302  CONTINUE
C
C     Which parameter to vary ? a or c ?
C
      NTRIEDB=0.
      IF(IFIN.EQ.1)THEN
      CELPRE(1)=CELPRE(1)+SPAR
      IF(CELPRE(1).GT.PMA(1))GO TO 1316
      IFIN=0
      NTRIED=NTRIED+1.
      ENDIF
      CELPRE(3)=CELPRE(3)+SPAR
      IF(CELPRE(3).GT.PMA(3))THEN
      CELPRE(3)=PMI(3)-SPAR
      IFIN=1
      GO TO 1302
      ENDIF
      NTRIED=NTRIED+1.
      GO TO 1304
1303  DEL=DELTAB*(1.-NTRIEDB/CY)
      I=3
      IF(RANDI(ISEED).GT.0.5)I=1
      CELPRE(I)=PSTARTB(I)+(DEL*(RANDI(ISEED)-0.5)*2.)
      NTRIEDB=NTRIEDB+1.
1304  CONTINUE
      CELPRE(2)=CELPRE(1)
      DO 1305 I=1,3
      DO 1305 J=1,3
1305  AL(I,J)=0.0
      CALL DCELL(CELPRE,AL,V1)
      IF(CELPRE(1).GT.PMA(1).AND.NTRIEDB.EQ.0.)GO TO 1316
      IF(NTRIEDB.NE.0.)GO TO 1306
      IF(V1.GT.VMAX.OR.V1.LT.VMIN) THEN
      NTRIED=NTRIED-1
      GO TO 1302
      ENDIF
C
1306  CALL CALCUL1(DIFF,DIFF2)
      IF(NMX.GT.NDAT10)THEN
      NTRIED=NTRIED-1
      GO TO 1302
      ENDIF
      IF(NTRIEDB.NE.0.)GO TO 1314
C
C... Rp value satisfying ???
C
      IF(LHKL.GE.NMAX)THEN
      RMAX=DIFF
      ICODE=2
      IF(DIFF.LE.RMAXREF)ICODE=1
      ELSE
      ICODE=1
      ENDIF
      CELOLD(1)=CELPRE(1)
      CELOLD(3)=CELPRE(3)
      IF(DIFF.GT.RMAX) GO TO 1317
      IF(LHKL.LT.NMAX) GO TO 1317
1314  IF(DIFF.LE.RMAX)THEN
      LLHKL=LHKL
      RMAX=DIFF
      RMAX2=DIFF2
      A=CELPRE(1)
      C=CELPRE(3)
      V2=V1
      IF(DIFF.LT.RMIN)THEN
      RMIN=DIFF
      BPAR(1)=A
      BPAR(3)=C
      V3=V1
      ENDIF
C
C... "Refine" that cell (by Monte Carlo too...)
C
      PSTARTB(1)=CELPRE(1)
      PSTARTB(3)=CELPRE(3)
      ENDIF
      IF(NTRIEDB.LE.NCYCLES)GO TO 1303
      NTRIEDB=0.
      IF(RMAX.GE.RMAX0(3))GO TO 1317
      IF(RMAX2.GE.0.15)GO TO 1317
      IPEN=NDAT-LLHKL
      IF(IPEN.GT.NIND)GO TO 1317
      IGC=IGC+1
C
C  Test if too much proposals, if yes decrease Rmax by 5%
C
      IGT=IGT+1.
      IF(NR.EQ.1)THEN
      IF(IGT.GT.50.)THEN
      IF((NTRIED/IGT).LT.1000.)THEN
      IF(RMAX0(3).GT.0.1)THEN
      RMAX0(3)=RMAX0(3)-RMAX0(3)*0.05
      WRITE(20,*)'  Rmax reduced by 5%, now Rmax = ',RMAX0(3)
      WRITE(*,*)'  Rmax reduced by 5%, now Rmax = ',RMAX0(3)
      ENDIF
      ENDIF
      ENDIF
      ENDIF
C
      IF(IGC.GT.10000)THEN
      WRITE(20,*)'   More than 10000 good cells = STOP'
      PRINT *,'   More than 10000 good cells = STOP'
      IGC=IGC-1
      GO TO 5000
      ENDIF
      CEL(1,IGC)=A
      CEL(2,IGC)=A
      CEL(3,IGC)=C
      CEL(4,IGC)=90.
      CEL(5,IGC)=90.
      CEL(6,IGC)=90.
C
C... Check for supercell
C
      CELPRE(1)=A
      CELPRE(2)=A
      CELPRE(3)=C
      DO 1340 I=1,3
      DO 1340 J=1,3
1340  AL(I,J)=0.0
      CALL DCELL(CELPRE,AL,V1)
      CALL CALCUL2(DIFF,IHKL,TH3,NCALC,IGC)
      KM(IGC)=LLHKL
      KM2(IGC)=LHKL
      IFI(IGC)=IFILE
      NSOL(IGC)=1
      VGC(IGC)=V1
      RP(IGC)=RMAX
      RP2(IGC)=DIFF
      IF(RP(IGC).LT.RPSMALL)THEN
      RPSMALL=RP(IGC)
      ISEE=1
      ELSE
      ISEE=0
      ENDIF
      CALL SUPCEL(LHKL,IHKL,CEL,IGC,VGC,4)
      CALL BRAV(LHKL,IHKL,IBR)
      IB(IGC)=IBR
      A=CEL(1,IGC)
      CEL(2,IGC)=A
      C=CEL(3,IGC)
      V2=VGC(IGC)
C
C... Check for interesting result
C
      IF(RP(IGC).LT.RMI)THEN
      WRITE(*,1215)RMAX,A,C,V2,IPEN
      WRITE(20,*)
      WRITE(20,*)'====================================================
     1==========================='
      WRITE(20,*)
      WRITE(20,*)' YOU HAVE FOUND AN INTERESTING RESULT : Rp < Rmin!'
      WRITE(*,*)
      WRITE(*,*)' YOU HAVE FOUND AN INTERESTING RESULT : Rp < Rmin!'
C
C... Refine that cell
C
      INDIC=2
      BB(3)=A
      BB(4)=A
      BB(5)=C
      BB(6)=90.
      BB(7)=90.
      BB(8)=90.
      AFI(3)=1
      AFI(4)=1
      AFI(5)=1
      AFI(6)=0
      AFI(7)=0
      AFI(8)=0
      CELPRE(1)=A
      CELPRE(2)=A
      CELPRE(3)=C
      DO 1310 I=1,3
      DO 1310 J=1,3
1310  AL(I,J)=0.0
      CALL DCELL(CELPRE,AL,V1)
      CALL CALCUL2(DIFF,IHKL,TH3,NCALC,IGC)
      CALL CELREF(INDIC,BB,AFI,LHKL,TH3,IHKL,DDT,DDQ)
      IF(NDAT.GE.20)THEN
      CNCALC(IGC)=NCALC
      FM20(IGC)=QO(20)/(2.*CNCALC(IGC)*DDQ)
      FF20(IGC)=20./(CNCALC(IGC)*DDT)
      WRITE(20,7000)FM20(IGC)
      WRITE(20,7001)FF20(IGC),DDT,NCALC
      WRITE(20,*)
      PRINT 7000,FM20(IGC)
      PRINT 7001,FF20(IGC),DDT,NCALC
      PRINT *
      ENDIF
      IREF=1
      GO TO 5000
      ENDIF
C
C Test if cell already found
C
      IF(IGC.GT.1)THEN
      DO 1318 I=1,IGC-1
      IF(IFI(I).NE.IFILE)GO TO 1318
      VDELT=VGC(IGC)/300.
      VP=VGC(IGC)+VDELT
      VM=VGC(IGC)-VDELT
      IF(VGC(I).GT.VP.OR.VGC(I).LT.VM)GO TO 1318
      ADELT=CEL(1,IGC)/500.
      AP=CEL(1,IGC)+ADELT
      AM=CEL(1,IGC)-ADELT
      IF(CEL(1,I).GT.AP.OR.CEL(1,I).LT.AM)GO TO 1318
      NSOL(I)=NSOL(I)+1
      IF(RP(IGC).LT.RP(I))THEN
      IF(ISEE.EQ.1)WRITE(*,1215)RMAX,A,C,V2,IPEN
      KM(I)=KM(IGC)
      VGC(I)=VGC(IGC)
      RP(I)=RP(IGC)
      CEL(1,I)=CEL(1,IGC)
      CEL(2,I)=CEL(2,IGC)
      CEL(3,I)=CEL(3,IGC)
      ENDIF
      IGC=IGC-1
      GO TO 1319
1318  CONTINUE
      IF(IVERB.EQ.1)WRITE(20,215)RMAX,NTRIED,A,C,V2,IPEN,ICODE
      IF(ISEE.EQ.1)WRITE(*,1215)RMAX,A,C,V2,IPEN
1319  CONTINUE
      ELSE
      IF(IVERB.EQ.1)WRITE(20,215)RMAX,NTRIED,A,C,V2,IPEN,ICODE
      IF(ISEE.EQ.1)WRITE(*,1215)RMAX,A,C,V2,IPEN
      ENDIF
C
C
1317  RMAX=RMAXREF
      CELPRE(1)=CELOLD(1)
      CELPRE(3)=CELOLD(3)
C
C... Stop if max limit of grid tests outpassed
C         or if K is pressed (tested every 30000 MC event)
C
      IF(CELPRE(1).GT.PMA(1))GO TO 1316
      TKILL=MOD(NTRIED,30000.)
      IF(TKILL.GE.0.)THEN
      CALL KILLK(PRESSEDK)
      IF(PRESSEDK)GO TO 1316
      ENDIF
      GO TO 1302
1316  CONTINUE
      IF(RMIN.EQ.RMAX)GO TO 1398
      IF(IVERB.EQ.1)THEN
      WRITE(20,*)
      WRITE(20,*)'Best result : a = ',BPAR(1),' Rp = ',RMIN
      WRITE(20,*)'Best result : c = ',BPAR(3),'V = ',V3
      WRITE(20,*)
      CALL DATN(DATENOW,TIMENOW)
      ENDIF
      RMIN=RMAX
1398  IF(PRESSEDK)GO TO 5000
C
C
1400  IF(NSYS(4).EQ.0)GO TO 1500
C
C    Orthorhombic case
C
C
      RPSMALL=1.
      PRINT *,'Orthorhombic: Rp     a       b       c        V     Nind'
*------------------------------------------------------------------------- 
*     Initialisation
*
C      CALL ESP_INIT(ISEED)
*
*-------------------------------------------------------------------------
      IFILE=4
      RMAX=RMAXREF
      RMIN=RMAX
      NTRIED=0.
      NCYCLES=1000.
      CY=NCYCLES*1.1
      CELPRE(4)=90.
      CELPRE(5)=90.
      CELPRE(6)=90.
      IF(IVERB.EQ.1)THEN
      WRITE(20,*)
      WRITE(20,*)'====================================================
     1==========================='
      WRITE(20,*)'Grid search :'
      WRITE(20,*)'   Results in orthorhombic :'
      WRITE(20,*)'====================================================
     1==========================='
      ENDIF
C
      IF(NGRID.EQ.3)THEN
      PMIN=2.
      PMAX=20.
      PMI(1)=PMIN
      PMA(1)=DMAX1*2.1
      IF(PMA(1).GT.PMAX)PMA(1)=PMAX
      PMI(2)=PMIN
      PMA(2)=DMAX2*2.1
      IF(PMA(2).GT.PMAX)PMA(2)=PMAX
      PMI(3)=PMIN
      PMA(3)=DMAX3*2.1
      IF(PMA(3).GT.PMAX)PMA(3)=PMAX
      VMIN=8.
      VMAX=PMA(1)*PMA(2)*PMA(3)
      IF(VMAX.GT.2000.)VMAX=2000.
      IF(IVERB.EQ.1)
     1WRITE(20,*)' Max (a,b,c), V ',PMA(1),PMA(2),PMA(3),VMAX
      ENDIF
C
      IF(IVERB.EQ.1)THEN
      WRITE(20,*)' Rp  Trial number    a   b   c        V  Nind Icod'
      WRITE(20,*)
      ENDIF
C
C     READ hkl Miller indices in ort.hkl
C
      TEMPO='ort.hkl'
      CALL OPEN_READ1(35,TEMPO)
      READ(35,*)NHKL0
      NHKL0=NDAT*20
      IF(NHKL0.GT.1000)NHKL0=1000
      DO 1401 I=1,NHKL0
1401  READ(35,*)(IHH(KK,I),KK=1,3)
      CLOSE(35)
C
C...  here starts the loop
C
      CELPRE(1)=PMI(1)-SPAR
      CELPRE(2)=PMI(2)-SPAR
      CELPRE(3)=PMI(3)-SPAR
      IFIN=1
      IFIN2=1
C
1402  CONTINUE
C
C     Which parameter to vary ? a or b or c ?
C
      NTRIEDB=0.
      IF(IFIN.EQ.1)THEN
      CELPRE(1)=CELPRE(1)+SPAR
      PRINT *,'  a = ',CELPRE(1)
      IF(CELPRE(1).GT.PMA(1))GO TO 1416
      IFIN=0
      NTRIED=NTRIED+1.
      ENDIF
      IF(IFIN2.EQ.1)THEN
      CELPRE(2)=CELPRE(2)+SPAR
      IFIN2=0
      ENDIF
      IF(CELPRE(2).GT.PMA(2))THEN
      CELPRE(2)=PMI(2)-SPAR
      IFIN=1
      GO TO 1402
      ENDIF
      NTRIED=NTRIED+1.
      CELPRE(3)=CELPRE(3)+SPAR
      IF(CELPRE(3).GT.PMA(3))THEN
      CELPRE(3)=PMI(3)-SPAR
      IFIN2=1
      GO TO 1402
      ENDIF
      NTRIED=NTRIED+1.
      GO TO 1404
1403  DEL=DELTAB*(1.-NTRIEDB/CY)
      X=RANDI(ISEED)
      IF(X.GE.0.AND.X.LT.0.33333)I=1
      IF(X.GE.0.33333.AND.X.LT.0.66666)I=2
      IF(X.GE.0.66666.AND.X.LE.1.)I=3
      CELPRE(I)=PSTARTB(I)+(DEL*(RANDI(ISEED)-0.5)*2.)
      NTRIEDB=NTRIEDB+1.
1404  CONTINUE
      DO 1405 I=1,3
      DO 1405 J=1,3
1405  AL(I,J)=0.0
      CALL DCELL(CELPRE,AL,V1)
      IF(CELPRE(1).GT.PMA(1).AND.NTRIEDB.EQ.0.)GO TO 1416
      IF(NTRIEDB.NE.0.)GO TO 1406
      IF(V1.GT.VMAX.OR.V1.LT.VMIN) THEN
      NTRIED=NTRIED-1.
      GO TO 1402
      ENDIF
C
1406  CALL CALCUL1(DIFF,DIFF2)
      IF(NMX.GT.NDAT10)THEN
      NTRIED=NTRIED-1
      GO TO 1402
      ENDIF
      IF(NTRIEDB.NE.0.)GO TO 1414
C
C... Rp value satisfying ???
C
      IF(LHKL.GE.NMAX)THEN
      RMAX=DIFF
      ICODE=2
      IF(DIFF.LE.RMAXREF)ICODE=1
      ELSE
      ICODE=1
      ENDIF
      CELOLD(1)=CELPRE(1)
      CELOLD(2)=CELPRE(2)
      CELOLD(3)=CELPRE(3)
      IF(DIFF.GT.RMAX) GO TO 1417
      IF(LHKL.LT.NMAX) GO TO 1417
1414  IF(DIFF.LE.RMAX)THEN
      LLHKL=LHKL
      RMAX=DIFF
      RMAX2=DIFF2
      A=CELPRE(1)
      B=CELPRE(2)
      C=CELPRE(3)
      V2=V1
      IF(DIFF.LT.RMIN)THEN
      RMIN=DIFF
      BPAR(1)=A
      BPAR(2)=B
      BPAR(3)=C
      V3=V1
      ENDIF
C
C... "Refine" that cell (by Monte Carlo too...)
C
      PSTARTB(1)=CELPRE(1)
      PSTARTB(2)=CELPRE(2)
      PSTARTB(3)=CELPRE(3)
      ENDIF
      IF(NTRIEDB.LE.NCYCLES)GO TO 1403
      NTRIEDB=0.
      IF(RMAX.GE.RMAX0(4))GO TO 1417
      IF(RMAX2.GE.0.15)GO TO 1417
      IPEN=NDAT-LLHKL
      IF(IPEN.GT.NIND)GO TO 1417
      IGC=IGC+1
C
C  Test if too much proposals, if yes decrease Rmax by 5%
C
      IGT=IGT+1.
      IF(NR.EQ.1)THEN
      IF(IGT.GT.50.)THEN
      IF((NTRIED/IGT).LT.10000.)THEN
      IF(RMAX0(4).GT.0.1)THEN
      RMAX0(4)=RMAX0(4)-RMAX0(4)*0.05
      WRITE(20,*)'  Rmax reduced by 5%, now Rmax = ',RMAX0(4)
      WRITE(*,*)'  Rmax reduced by 5%, now Rmax = ',RMAX0(4)
      ENDIF
      ENDIF
      ENDIF
      ENDIF
C
      IF(IGC.GT.10000)THEN
      WRITE(20,*)'   More than 10000 good cells = STOP'
      PRINT *,'   More than 10000 good cells = STOP'
      IGC=IGC-1
      GO TO 5000
      ENDIF
      CEL(1,IGC)=A
      CEL(2,IGC)=B
      CEL(3,IGC)=C
      CEL(4,IGC)=90.
      CEL(5,IGC)=90.
      CEL(6,IGC)=90.
C
C... Check for supercell
C
      CELPRE(1)=A
      CELPRE(2)=B
      CELPRE(3)=C
      DO 1440 I=1,3
      DO 1440 J=1,3
1440  AL(I,J)=0.0
      CALL DCELL(CELPRE,AL,V1)
      CALL CALCUL2(DIFF,IHKL,TH3,NCALC,IGC)
      KM(IGC)=LLHKL
      KM2(IGC)=LHKL
      IFI(IGC)=IFILE
      NSOL(IGC)=1
      VGC(IGC)=V1
      RP(IGC)=RMAX
      RP2(IGC)=DIFF
      IF(RP(IGC).LT.RPSMALL)THEN
      RPSMALL=RP(IGC)
      ISEE=1
      ELSE
      ISEE=0
      ENDIF
      CALL SUPCEL(LHKL,IHKL,CEL,IGC,VGC,1)
      CALL BRAV(LHKL,IHKL,IBR)
      IB(IGC)=IBR
      A=CEL(1,IGC)
      B=CEL(2,IGC)
      C=CEL(3,IGC)
      V2=VGC(IGC)
C
C... Check for interesting result
C
      IF(RP(IGC).LT.RMI)THEN
      WRITE(*,1415)RMAX,A,B,C,V2,IPEN
      WRITE(20,*)
      WRITE(20,*)'====================================================
     1==========================='
      WRITE(20,*)
      WRITE(20,*)' YOU HAVE FOUND AN INTERESTING RESULT : Rp < Rmin!'
      WRITE(*,*)
      WRITE(*,*)' YOU HAVE FOUND AN INTERESTING RESULT : Rp < Rmin!'
C
C... Refine that cell
C
      INDIC=0
      BB(3)=A
      BB(4)=B
      BB(5)=C
      BB(6)=90.
      BB(7)=90.
      BB(8)=90.
      AFI(3)=1
      AFI(4)=1
      AFI(5)=1
      AFI(6)=0
      AFI(7)=0
      AFI(8)=0
      CELPRE(1)=A
      CELPRE(2)=B
      CELPRE(3)=C
      DO 1410 I=1,3
      DO 1410 J=1,3
1410  AL(I,J)=0.0
      CALL DCELL(CELPRE,AL,V1)
      CALL CALCUL2(DIFF,IHKL,TH3,NCALC,IGC)
      CALL CELREF(INDIC,BB,AFI,LHKL,TH3,IHKL,DDT,DDQ)
      IF(NDAT.GE.20)THEN
      CNCALC(IGC)=NCALC
      FM20(IGC)=QO(20)/(2.*CNCALC(IGC)*DDQ)
      FF20(IGC)=20./(CNCALC(IGC)*DDT)
      WRITE(20,7000)FM20(IGC)
      WRITE(20,7001)FF20(IGC),DDT,NCALC
      WRITE(20,*)
      PRINT 7000,FM20(IGC)
      PRINT 7001,FF20(IGC),DDT,NCALC
      PRINT *
      ENDIF
      IREF=1
      GO TO 5000
      ENDIF
C
C Test if cell already found
C
      IF(IGC.GT.1)THEN
      DO 1418 I=1,IGC-1
      IF(IFI(I).NE.IFILE)GO TO 1418
      VDELT=VGC(IGC)/300.
      VP=VGC(IGC)+VDELT
      VM=VGC(IGC)-VDELT
      IF(VGC(I).GT.VP.OR.VGC(I).LT.VM)GO TO 1418
      ADELT=CEL(1,IGC)/500.
      AP=CEL(1,IGC)+ADELT
      AM=CEL(1,IGC)-ADELT
      NA=0
      IF(CEL(1,I).GT.AP.OR.CEL(1,I).LT.AM)NA=1
      BDELT=CEL(2,IGC)/500.
      BP=CEL(2,IGC)+BDELT
      BM=CEL(2,IGC)-BDELT
      NB=0
      IF(CEL(1,I).GT.BP.OR.CEL(1,I).LT.BM)NB=1
      CDELT=CEL(3,IGC)/500.
      CP=CEL(3,IGC)+CDELT
      CM=CEL(3,IGC)-CDELT
      NC=0
      IF(CEL(1,I).GT.CP.OR.CEL(1,I).LT.CM)NC=1
      IF(NA.EQ.1.AND.NB.EQ.1.AND.NC.EQ.1)GO TO 1418
      NSOL(I)=NSOL(I)+1
      IF(RP(IGC).LT.RP(I))THEN
      IF(ISEE.EQ.1)WRITE(*,1415)RMAX,A,B,C,V2,IPEN
      KM(I)=KM(IGC)
      VGC(I)=VGC(IGC)
      RP(I)=RP(IGC)
      CEL(1,I)=CEL(1,IGC)
      CEL(2,I)=CEL(2,IGC)
      CEL(3,I)=CEL(3,IGC)
      ENDIF
      IGC=IGC-1
      GO TO 1419
1418  CONTINUE
      IF(IVERB.EQ.1)WRITE(20,415)RMAX,NTRIED,A,B,C,V2,IPEN,ICODE
      IF(ISEE.EQ.1)WRITE(*,1415)RMAX,A,B,C,V2,IPEN
1419  CONTINUE
      ELSE
      IF(IVERB.EQ.1)WRITE(20,415)RMAX,NTRIED,A,B,C,V2,IPEN,ICODE
      IF(ISEE.EQ.1)WRITE(*,1415)RMAX,A,B,C,V2,IPEN
      ENDIF
C
C
1417  RMAX=RMAXREF
      CELPRE(1)=CELOLD(1)
      CELPRE(2)=CELOLD(2)
      CELPRE(3)=CELOLD(3)
C
C... Stop if max limit of grid tests outpassed
C         or if K is pressed (tested every 30000 MC event)
C
      IF(CELPRE(1).GT.PMA(1))GO TO 1416
      TKILL=MOD(NTRIED,30000.)
      IF(TKILL.GE.0.)THEN
      CALL KILLK(PRESSEDK)
      IF(PRESSEDK)GO TO 1416
      ENDIF
      GO TO 1402
1416  CONTINUE
      IF(RMIN.EQ.RMAX)GO TO 1498
      IF(IVERB.EQ.1)THEN
      WRITE(20,*)
      WRITE(20,*)'Best result : a = ',BPAR(1),' Rp = ',RMIN
      WRITE(20,*)'Best result : b = ',BPAR(2)
      WRITE(20,*)'Best result : c = ',BPAR(3),'V = ',V3
      WRITE(20,*)
      CALL DATN(DATENOW,TIMENOW)
      ENDIF
      RMIN=RMAX
1498  IF(PRESSEDK)GO TO 5000
C
C
1500  IF(NSYS(5).EQ.0)GO TO 1600
C
C    Monoclinic case - would be too long in grid search, but...
C
C
      RPSMALL=1.
      PRINT *,
     1'Monoclinic:   Rp     a       b       c       bet     V     Nind'
*------------------------------------------------------------------------- 
*     Initialisation
*
C      CALL ESP_INIT(ISEED)
*
*-------------------------------------------------------------------------
      IFILE=5
      RMAX=RMAXREF
      RMIN=RMAX
      NTRIED=0.
      NCYCLES=2000.
      CY=NCYCLES*1.1
      CELPRE(4)=90.
      CELPRE(6)=90.
      IF(IVERB.EQ.1)THEN
      WRITE(20,*)
      WRITE(20,*)'====================================================
     1==========================='
      WRITE(20,*)'Grid search :'
      WRITE(20,*)'   Results in monoclinic :'
      WRITE(20,*)'====================================================
     1==========================='
      ENDIF
C
      IF(NGRID.EQ.3)THEN
      PMIN=2.
      PMAX=20.
      PMI(1)=PMIN
      PMA(1)=DMAX1*2.1
      IF(PMA(1).GT.PMAX)PMA(1)=PMAX
      PMI(2)=PMIN
      PMA(2)=DMAX1*2.1
      IF(PMA(2).GT.PMAX)PMA(2)=PMAX
      PMI(3)=PMIN
      PMA(3)=DMAX2*2.1
      IF(PMA(3).GT.PMAX)PMA(3)=PMAX
      VMIN=8.
      VMAX=PMA(1)*PMA(2)*PMA(3)
      IF(VMAX.GT.2000.)VMAX=2000.
      IF(IVERB.EQ.1)
     1WRITE(20,*)' Max (a,b,c), V ',PMA(1),PMA(2),PMA(3),VMAX
      ENDIF
C
      IF(IVERB.EQ.1)THEN
      WRITE(20,*)' Rp  Trial number    a   b   c     bet  V  Nind Icod'
      WRITE(20,*)
      ENDIF
C
C     READ hkl Miller indices in mon.hkl
C
      TEMPO='mon.hkl'
      CALL OPEN_READ1(35,TEMPO)
      READ(35,*)NHKL0
      NHKL0=NDAT*20
      IF(NHKL0.GT.1000)NHKL0=1000
      DO 1501 I=1,NHKL0
1501  READ(35,*)(IHH(KK,I),KK=1,3)
      CLOSE(35)
C
C...  here starts the loop
C
      CELPRE(1)=PMI(1)-SPAR
      CELPRE(2)=PMI(2)-SPAR
      CELPRE(3)=PMI(3)-SPAR
      CELPRE(5)=PMI(5)-SANG
      IFIN1=1
      IFIN2=1
      IFIN3=1
C
1502  CONTINUE
C
C     Which parameter to vary ? a or b or c or bet ?
C
      NTRIEDB=0.
      IF(IFIN1.EQ.1)THEN
      CELPRE(1)=CELPRE(1)+SPAR
      PRINT *,'  a = ',CELPRE(1)
      IF(CELPRE(1).GT.PMA(1))GO TO 1516
      IFIN1=0
      NTRIED=NTRIED+1.
      ENDIF
      IF(IFIN2.EQ.1)THEN
      CELPRE(2)=CELPRE(2)+SPAR
      IFIN2=0
      ENDIF
      IF(CELPRE(2).GT.PMA(2))THEN
      CELPRE(2)=PMI(2)-SPAR
      IFIN1=1
      GO TO 1502
      ENDIF
      NTRIED=NTRIED+1.
      IF(IFIN3.EQ.1)THEN
      CELPRE(3)=CELPRE(3)+SPAR
      IFIN3=0
      ENDIF
      IF(CELPRE(3).GT.PMA(3))THEN
      CELPRE(3)=PMI(3)-SPAR
      IFIN2=1
      GO TO 1502
      ENDIF
      NTRIED=NTRIED+1.
      CELPRE(5)=CELPRE(5)+SANG
      IF(CELPRE(5).GT.PMA(5))THEN
      CELPRE(5)=PMI(5)-SANG
      IFIN3=1
      GO TO 1502
      ENDIF
      NTRIED=NTRIED+1.
      GO TO 1504
1503  DEL=DELTAB*(1.-NTRIEDB/CY)
      DELD=DELTAD*(1.-NTRIEDB/CY)
      X=RANDI(ISEED)
      IF(X.GE.0.AND.X.LT.0.25)I=1
      IF(X.GE.0.25.AND.X.LT.0.5)I=2
      IF(X.GE.0.5.AND.X.LT.0.75)I=3
      IF(X.GE.0.75.AND.X.LE.1.0)I=5
      IF(I.NE.5)THEN
      CELPRE(I)=PSTARTB(I)+(DEL*(RANDI(ISEED)-0.5)*2.)
      ELSE
      CELPRE(I)=PSTARTB(I)+(DELD*(RANDI(ISEED)-0.5)*2.)
      ENDIF
      NTRIEDB=NTRIEDB+1.
1504  CONTINUE
      DO 1505 I=1,3
      DO 1505 J=1,3
1505  AL(I,J)=0.0
      CALL DCELL(CELPRE,AL,V1)
      IF(CELPRE(1).GT.PMA(1).AND.NTRIEDB.EQ.0.)GO TO 1516
      IF(NTRIEDB.NE.0.)GO TO 1506
      IF(V1.GT.VMAX.OR.V1.LT.VMIN) THEN
      NTRIED=NTRIED-1.
      GO TO 1502
      ENDIF
C
1506  CALL CALCUL1(DIFF,DIFF2)
      IF(NMX.GT.NDAT10)THEN
      NTRIED=NTRIED-1
      GO TO 1502
      ENDIF
      IF(NTRIEDB.NE.0.)GO TO 1514
C
C... Rp value satisfying ???
C
      IF(LHKL.GE.NMAX)THEN
      RMAX=DIFF
      ICODE=2
      IF(DIFF.LE.RMAXREF)ICODE=1
      ELSE
      ICODE=1
      ENDIF
      CELOLD(1)=CELPRE(1)
      CELOLD(2)=CELPRE(2)
      CELOLD(3)=CELPRE(3)
      CELOLD(5)=CELPRE(5)
      IF(DIFF.GT.RMAX) GO TO 1517
      IF(LHKL.LT.NMAX) GO TO 1517
1514  IF(DIFF.LE.RMAX)THEN
      LLHKL=LHKL
      RMAX=DIFF
      RMAX2=DIFF2
      A=CELPRE(1)
      B=CELPRE(2)
      C=CELPRE(3)
      BET=CELPRE(5)
      V2=V1
      IF(DIFF.LT.RMIN)THEN
      RMIN=DIFF
      BPAR(1)=A
      BPAR(2)=B
      BPAR(3)=C
      BPAR(5)=BET
      V3=V1
      ENDIF
C
C... "Refine" that cell (by Monte Carlo too...)
C
      PSTARTB(1)=CELPRE(1)
      PSTARTB(2)=CELPRE(2)
      PSTARTB(3)=CELPRE(3)
      PSTARTB(5)=CELPRE(5)
      ENDIF
      IF(NTRIEDB.LE.NCYCLES)GO TO 1503
      NTRIEDB=0.
      IF(RMAX.GE.RMAX0(5))GO TO 1517
      IF(RMAX2.GE.0.15)GO TO 1517
      IPEN=NDAT-LLHKL
      IF(IPEN.GT.NIND)GO TO 1517
      IGC=IGC+1
C
C  Test if too much proposals, if yes decrease Rmax by 5%
C
      IGT=IGT+1.
      IF(NR.EQ.1)THEN
      IF(IGT.GT.50.)THEN
      IF((NTRIED/IGT).LT.100000.)THEN
      IF(RMAX0(5).GT.0.1)THEN
      RMAX0(5)=RMAX0(5)-RMAX0(5)*0.05
      WRITE(20,*)'  Rmax reduced by 5%, now Rmax = ',RMAX0(5)
      WRITE(*,*)'  Rmax reduced by 5%, now Rmax = ',RMAX0(5)
      ENDIF
      ENDIF
      ENDIF
      ENDIF
C
      IF(IGC.GT.10000)THEN
      WRITE(20,*)'   More than 10000 good cells = STOP'
      PRINT *,'   More than 10000 good cells = STOP'
      IGC=IGC-1
      GO TO 5000
      ENDIF
      CEL(1,IGC)=A
      CEL(2,IGC)=B
      CEL(3,IGC)=C
      CEL(4,IGC)=90.
      CEL(5,IGC)=BET
      CEL(6,IGC)=90.
C
C... Check for supercell
C
      CELPRE(1)=A
      CELPRE(2)=B
      CELPRE(3)=C
      CELPRE(5)=BET
      DO 1540 I=1,3
      DO 1540 J=1,3
1540  AL(I,J)=0.0
      CALL DCELL(CELPRE,AL,V1)
      CALL CALCUL2(DIFF,IHKL,TH3,NCALC,IGC)
      KM(IGC)=LLHKL
      KM2(IGC)=LHKL
      IFI(IGC)=IFILE
      NSOL(IGC)=1
      VGC(IGC)=V1
      RP(IGC)=RMAX
      RP2(IGC)=DIFF
      IF(RP(IGC).LT.RPSMALL)THEN
      RPSMALL=RP(IGC)
      ISEE=1
      ELSE
      ISEE=0
      ENDIF
      CALL SUPCEL(LHKL,IHKL,CEL,IGC,VGC,1)
      A=CEL(1,IGC)
      B=CEL(2,IGC)
      C=CEL(3,IGC)
      V2=VGC(IGC)
C
C... Check for interesting result
C
      IF(RP(IGC).LT.RMI)THEN
      WRITE(*,1515)RMAX,A,B,C,BET,V2,IPEN
      WRITE(20,*)
      WRITE(20,*)'====================================================
     1==========================='
      WRITE(20,*)
      WRITE(20,*)' YOU HAVE FOUND AN INTERESTING RESULT : Rp < Rmin!'
      WRITE(*,*)
      WRITE(*,*)' YOU HAVE FOUND AN INTERESTING RESULT : Rp < Rmin!'
C
C... Refine that cell
C
      INDIC=0
      BB(3)=A
      BB(4)=B
      BB(5)=C
      BB(6)=90.
      BB(7)=BET
      BB(8)=90.
      AFI(3)=1
      AFI(4)=1
      AFI(5)=1
      AFI(6)=0
      AFI(7)=1
      AFI(8)=0
      CELPRE(1)=A
      CELPRE(2)=B
      CELPRE(3)=C
      CELPRE(5)=BET
      DO 1510 I=1,3
      DO 1510 J=1,3
1510  AL(I,J)=0.0
      CALL DCELL(CELPRE,AL,V1)
      CALL CALCUL2(DIFF,IHKL,TH3,NCALC,IGC)
      CALL CELREF(INDIC,BB,AFI,LHKL,TH3,IHKL,DDT,DDQ)
      IF(NDAT.GE.20)THEN
      CNCALC(IGC)=NCALC
      FM20(IGC)=QO(20)/(2.*CNCALC(IGC)*DDQ)
      FF20(IGC)=20./(CNCALC(IGC)*DDT)
      WRITE(20,7000)FM20(IGC)
      WRITE(20,7001)FF20(IGC),DDT,NCALC
      WRITE(20,*)
      PRINT 7000,FM20(IGC)
      PRINT 7001,FF20(IGC),DDT,NCALC
      PRINT *
      ENDIF
      IREF=1
      GO TO 5000
      ENDIF
C
C Test if cell already found
C
      IF(IGC.GT.1)THEN
      DO 1518 I=1,IGC-1
      IF(IFI(I).NE.IFILE)GO TO 1518
      VDELT=VGC(IGC)/300.
      VP=VGC(IGC)+VDELT
      VM=VGC(IGC)-VDELT
      IF(VGC(I).GT.VP.OR.VGC(I).LT.VM)GO TO 1518
      BDELT=CEL(2,IGC)/500.
      BP=CEL(2,IGC)+BDELT
      BM=CEL(2,IGC)-BDELT
      IF(CEL(2,I).GT.BP.OR.CEL(2,I).LT.BM)GO TO 1518
      BETDELT=CEL(5,IGC)/500.
      BETP=CEL(5,IGC)+BETDELT
      BETM=CEL(5,IGC)-BETDELT
      IF(CEL(5,I).GT.BETP.OR.CEL(5,I).LT.BETM)GO TO 1518
      ADELT=CEL(1,IGC)/500.
      AP=CEL(1,IGC)+ADELT
      AM=CEL(1,IGC)-ADELT
      NA=0
      IF(CEL(1,I).GT.AP.OR.CEL(1,I).LT.AM)NA=1
      CDELT=CEL(3,IGC)/500.
      CP=CEL(3,IGC)+CDELT
      CM=CEL(3,IGC)-CDELT
      NC=0
      IF(CEL(1,I).GT.CP.OR.CEL(1,I).LT.CM)NC=1
      IF(NA.EQ.1.AND.NC.EQ.1)GO TO 1518
      NSOL(I)=NSOL(I)+1
      IF(RP(IGC).LT.RP(I))THEN
      IF(ISEE.EQ.1)WRITE(*,1515)RMAX,A,B,C,BET,V2,IPEN
      KM(I)=KM(IGC)
      VGC(I)=VGC(IGC)
      RP(I)=RP(IGC)
      CEL(1,I)=CEL(1,IGC)
      CEL(2,I)=CEL(2,IGC)
      CEL(3,I)=CEL(3,IGC)
      CEL(5,I)=CEL(5,IGC)
      ENDIF
      IGC=IGC-1
      GO TO 1519
1518  CONTINUE
      IF(IVERB.EQ.1)
     1WRITE(20,515)RMAX,NTRIED,A,B,C,BET,V2,IPEN,ICODE
      IF(ISEE.EQ.1)WRITE(*,1515)RMAX,A,B,C,BET,V2,IPEN
1519  CONTINUE
      ELSE
      IF(IVERB.EQ.1)
     1WRITE(20,515)RMAX,NTRIED,A,B,C,BET,V2,IPEN,ICODE
      IF(ISEE.EQ.1)WRITE(*,1515)RMAX,A,B,C,BET,V2,IPEN
      ENDIF
C
C
1517  RMAX=RMAXREF
      CELPRE(1)=CELOLD(1)
      CELPRE(2)=CELOLD(2)
      CELPRE(3)=CELOLD(3)
      CELPRE(5)=CELOLD(5)
C
C... Stop if max limit of grid tests outpassed
C         or if K is pressed (tested every 30000 MC event)
C
      IF(CELPRE(1).GT.PMA(1))GO TO 1516
      TKILL=MOD(NTRIED,30000.)
      IF(TKILL.GE.0.)THEN
      CALL KILLK(PRESSEDK)
      IF(PRESSEDK)GO TO 1516
      ENDIF
      GO TO 1502
1516  CONTINUE
      IF(RMIN.EQ.RMAX)GO TO 1598
      IF(IVERB.EQ.1)THEN
      WRITE(20,*)
      WRITE(20,*)'Best result : a =   ',BPAR(1),'  Rp = ',RMIN
      WRITE(20,*)'Best result : b =   ',BPAR(2)
      WRITE(20,*)'Best result : c =   ',BPAR(3)
      WRITE(20,*)'Best result : bet = ',BPAR(5),'  V = ',V3
      WRITE(20,*)
      CALL DATN(DATENOW,TIMENOW)
      ENDIF
      RMIN=RMAX
1598  IF(PRESSEDK)GO TO 5000
C
C
1600  IF(NSYS(6).EQ.0)GO TO 5000
C
C    Triclinic case - would be too long in grid search...
C
C
5000  CONTINUE
      IF(IGC.EQ.0)GO TO 6000
      IMEM=0
C
C   Prepare sorted output for CHEKCELL
C
      WRITE(20,*)
      WRITE(20,*)'====================================================
     1==========================='
      WRITE(20,*)'  FINAL LIST OF CELL PROPOSALS, sorted by McM20 :'
      WRITE(20,*)'====================================================
     1==========================='
      WRITE(20,*)
      WRITE(20,*)'      Global list as produced in the .ckm file'
      WRITE(20,*)'            (IN=number of indexed lines)'
      WRITE(20,*)'  The correct cell has some chances to be just below'
      WRITE(20,*)
      WRITE(20,19993)
19993 FORMAT('IN  F.o.M.    Volume   V/V1      a        b        c      
     1 alpha   beta    gamma   Bravais lattice')
      WRITE(20,*)
      TEMPO=FILE(1:LFILE)//'.ckm'
      INQUIRE(FILE=TEMPO,EXIST=QEX)
      IF(.NOT.QEX) GO TO 1998
      CALL FILEDEL(24,TEMPO)
1998  CALL OPEN_WRITE1(24,TEMPO)
      TEMPO=FILE(1:LFILE)//'.mcm'
      INQUIRE(FILE=TEMPO,EXIST=QEX)
      IF(.NOT.QEX) GO TO 2066
      CALL FILEDEL(25,TEMPO)
2066  CALL OPEN_WRITE1(25,TEMPO)
C
C  Calculate new F.o.M
C      
      DO 18000 J=1,IGC
      IF(RP(J).LT.0.001)RP(J)=0.001
      X2=1.
      IF(IB(J).EQ.1)X2=4.
      IF(IB(J).EQ.2)X2=2.
      IF(IB(J).EQ.3)X2=2.
      IF(IB(J).EQ.4)X2=2.
      IF(IB(J).EQ.5)X2=6.
      IF(IFI(J).EQ.7)X2=6.
      X3=1.
      IF(IFI(J).EQ.1)X3=6.
      IF(IFI(J).EQ.2)X3=4.
      IF(IFI(J).EQ.3)X3=4.
      IF(IFI(J).EQ.4)X3=2.
      IF(IFI(J).EQ.7)X3=6.
      XFOM(J)=1./RP(J)*100./CNCALC(J)*X2*X3
18000 CONTINUE
C
      CALL SORT(IGC,XFOM,LL)
      IMEM=IMEM+1
      IM(IMEM)=LL(IGC)
      INDXPROG='McMaille4.00'
      IBEST=LL(IGC)
      DO 2000 I=1,IGC
      IF(I.GT.20)GO TO 20000
      J=LL(IGC+1-I)
      IPEDIG=J
      DO 2067 K=1,6
2067  CELPRE(K)=CEL(K,J)
      CALL DCELL(CELPRE,AL,V1)
      QL(1)=AL(1,1)*10000.
      QL(2)=AL(2,2)*10000.
      QL(3)=AL(3,3)*10000.
      QL(4)=AL(2,3)*10000.
      QL(5)=AL(1,3)*10000.
      QL(6)=AL(1,2)*10000.
      IF(IB(J).EQ.1)BL='I'
      IF(IB(J).EQ.2)BL='A'
      IF(IB(J).EQ.3)BL='B'
      IF(IB(J).EQ.4)BL='C'
      IF(IB(J).EQ.5)BL='F'
      IF(IB(J).EQ.6)BL='P'
      IF(IFI(J).EQ.7)BL='R'
      IF(IFI(J).EQ.1)MORE='Cubic *****'
      IF(IFI(J).EQ.2)MORE='Hexag **** '
      IF(IFI(J).EQ.3)MORE='Tetra **** '
      IF(IFI(J).EQ.4)MORE='Ortho ***  '
      IF(IFI(J).EQ.5)MORE='           '
      IF(IFI(J).EQ.6)MORE='           '
      IF(IFI(J).EQ.7)MORE='Rhomb **** '
      VR=VGC(J)/VGC(LL(IGC))
      WRITE(24,2002)KM(J),XFOM(J),VGC(J),VR,(CEL(K,J),K=1,6)
      WRITE(20,20031)KM(J),XFOM(J),VGC(J),VR,(CEL(K,J),K=1,6),BL,MORE
      WRITE(25,2069)KM(J),XFOM(J),VGC(J),VR,BL,INDXPROG,
     1DATENOW,TIMENOW,IPEDIG,(CEL(K,J),K=1,6),(QL(K),K=1,6)
2000  CONTINUE
20000 CONTINUE
2001  FORMAT(F8.3,F11.3,F6.2,2I4,3F9.4,3F8.3)
2002  FORMAT(I2,F8.2,F11.3,F6.2,39X,3F9.4,3F8.3)
20031 FORMAT(I2,F8.2,F11.3,F6.2,3X,3F9.4,3F8.3,4X,A1,2X,A11)
2069  FORMAT(I2,F8.2,F11.3,F6.2,1X,A1,1X,A12,1X,A7,1X,A8,I7,
     1 3F9.4,3F8.3,6F10.4)
      WRITE(20,*)
      WRITE(20,*)'====================================================
     1==========================='
      WRITE(20,*)'    FINAL LIST OF CELL PROPOSALS, sorted by F(20) :'
      WRITE(20,*)'====================================================
     1==========================='
      WRITE(20,*)
      WRITE(20,*)'      Global list (IN=number of indexed lines)'
      WRITE(20,*)'  The correct cell has some chances to be just below'
      WRITE(20,*)
      WRITE(20,19994)
19994 FORMAT('IN  F(20)     Volume   V/V1      a        b        c      
     1 alpha   beta    gamma   Bravais lattice')
      WRITE(20,*)
      CALL SORT(IGC,FF20,LL)
      IMEM=IMEM+1
      IM(IMEM)=LL(IGC)
      DO 20001 I=1,IGC
      IF(I.GT.20)GO TO 20002
      J=LL(IGC+1-I)
      IPEDIG=J
      DO 20671 K=1,6
20671 CELPRE(K)=CEL(K,J)
      CALL DCELL(CELPRE,AL,V1)
      IF(IB(J).EQ.1)BL='I'
      IF(IB(J).EQ.2)BL='A'
      IF(IB(J).EQ.3)BL='B'
      IF(IB(J).EQ.4)BL='C'
      IF(IB(J).EQ.5)BL='F'
      IF(IB(J).EQ.6)BL='P'
      IF(IFI(J).EQ.7)BL='R'
      IF(IFI(J).EQ.1)MORE='Cubic *****'
      IF(IFI(J).EQ.2)MORE='Hexag **** '
      IF(IFI(J).EQ.3)MORE='Tetra **** '
      IF(IFI(J).EQ.4)MORE='Ortho ***  '
      IF(IFI(J).EQ.5)MORE='           '
      IF(IFI(J).EQ.6)MORE='           '
      IF(IFI(J).EQ.7)MORE='Rhomb **** '
      VR=VGC(J)/VGC(LL(IGC))
      WRITE(20,20031)KM(J),FF20(J),VGC(J),VR,(CEL(K,J),K=1,6),BL,MORE
20001 CONTINUE
20002 CONTINUE
      WRITE(20,*)
      WRITE(20,*)'====================================================
     1==========================='
      WRITE(20,*)'    FINAL LIST OF CELL PROPOSALS, sorted by M(20) :'
      WRITE(20,*)'====================================================
     1==========================='
      WRITE(20,*)
      WRITE(20,*)'      Global list   (IN=number of indexed lines)'
      WRITE(20,*)'  The correct cell has some chances to be just below'
      WRITE(20,*)
      WRITE(20,19995)
19995 FORMAT('IN  M(20)     Volume   V/V1      a        b        c      
     1 alpha   beta    gamma   Bravais lattice')
      WRITE(20,*)
      CALL SORT(IGC,FM20,LL)
      IMEM=IMEM+1
      IM(IMEM)=LL(IGC)
      DO 20003 I=1,IGC
      IF(I.GT.20)GO TO 20004
      J=LL(IGC+1-I)
      IPEDIG=J
      DO 20672 K=1,6
20672 CELPRE(K)=CEL(K,J)
      CALL DCELL(CELPRE,AL,V1)
      IF(IB(J).EQ.1)BL='I'
      IF(IB(J).EQ.2)BL='A'
      IF(IB(J).EQ.3)BL='B'
      IF(IB(J).EQ.4)BL='C'
      IF(IB(J).EQ.5)BL='F'
      IF(IB(J).EQ.6)BL='P'
      IF(IFI(J).EQ.7)BL='R'
      IF(IFI(J).EQ.1)MORE='Cubic *****'
      IF(IFI(J).EQ.2)MORE='Hexag **** '
      IF(IFI(J).EQ.3)MORE='Tetra **** '
      IF(IFI(J).EQ.4)MORE='Ortho ***  '
      IF(IFI(J).EQ.5)MORE='           '
      IF(IFI(J).EQ.6)MORE='           '
      IF(IFI(J).EQ.7)MORE='Rhomb **** '
      VR=VGC(J)/VGC(LL(IGC))
      WRITE(20,20031)KM(J),FM20(J),VGC(J),VR,(CEL(K,J),K=1,6),BL,MORE
20003 CONTINUE
20004 CONTINUE
      WRITE(20,*)
      WRITE(20,*)'====================================================
     1==========================='
      WRITE(20,*)'  See for the highest F.o.M. above the cell(s) with 
     1highest symmetry, if any'
      WRITE(20,*)'  (Cubic, hexagonal, etc), they could correspond to 
     1the the right solution'
      WRITE(20,*)'====================================================
     1==========================='
      WRITE(20,*)
      CLOSE(24)
      CLOSE(25)
C
C   Make output for cells sorted by symmetry
C
      CALL SORT(IGC,RP,LL)
      WRITE(20,*)
      WRITE(20,*)'  Cells sorted by symmetry'
      WRITE(20,*)
      WRITE(20,1999)
1999  FORMAT('    Rp     Vol     Vol/V1 Ind Nsol    a        b     
     1    c      alpha  beta  gamma')
19992 FORMAT('    Rp2    Vol     Vol/V1 Ind Nsol    a        b     
     1    c      alpha  beta  gamma')
      DO 2018 I=1,7
2018  ISYST(I)=0
      DO 2019 I=1,IGC
      IF(IFI(I).EQ.1)ISYST(1)=1
      IF(IFI(I).EQ.2)ISYST(2)=1
      IF(IFI(I).EQ.3)ISYST(3)=1
      IF(IFI(I).EQ.4)ISYST(4)=1
      IF(IFI(I).EQ.5)ISYST(5)=1
      IF(IFI(I).EQ.6)ISYST(6)=1
      IF(IFI(I).EQ.7)ISYST(7)=1
2019  CONTINUE
      DO 2020 JIFI=1,7
      IF(ISYST(JIFI).EQ.0)GO TO 2020
      WRITE(20,*)
      GO TO (2021,2022,2023,2024,2025,2026,2027)JIFI
2021  WRITE(20,*)'   Cubic cells'
      TEMPO=FILE(1:LFILE)//'_cub.ckm'
      INQUIRE(FILE=TEMPO,EXIST=QEX)
      IF(.NOT.QEX) GO TO 2521
      CALL FILEDEL(24,TEMPO)
2521  CALL OPEN_WRITE1(24,TEMPO)
      GO TO 2028
2022  WRITE(20,*)'   Hexagonal/trigonal cells'
      TEMPO=FILE(1:LFILE)//'_hex.ckm'
      INQUIRE(FILE=TEMPO,EXIST=QEX)
      IF(.NOT.QEX) GO TO 2522
      CALL FILEDEL(24,TEMPO)
2522  CALL OPEN_WRITE1(24,TEMPO)
      GO TO 2028
2023  WRITE(20,*)'   Tetragonal cells'
      TEMPO=FILE(1:LFILE)//'_tet.ckm'
      INQUIRE(FILE=TEMPO,EXIST=QEX)
      IF(.NOT.QEX) GO TO 2523
      CALL FILEDEL(24,TEMPO)
2523  CALL OPEN_WRITE1(24,TEMPO)
      GO TO 2028
2024  WRITE(20,*)'   Orthorhombic cells'
      TEMPO=FILE(1:LFILE)//'_ort.ckm'
      INQUIRE(FILE=TEMPO,EXIST=QEX)
      IF(.NOT.QEX) GO TO 2524
      CALL FILEDEL(24,TEMPO)
2524  CALL OPEN_WRITE1(24,TEMPO)
      GO TO 2028
2025  WRITE(20,*)'   Monoclinic cells'
      TEMPO=FILE(1:LFILE)//'_mon.ckm'
      INQUIRE(FILE=TEMPO,EXIST=QEX)
      IF(.NOT.QEX) GO TO 2525
      CALL FILEDEL(24,TEMPO)
2525  CALL OPEN_WRITE1(24,TEMPO)
      GO TO 2028
2026  WRITE(20,*)'   Triclinic cells'
      TEMPO=FILE(1:LFILE)//'_tri.ckm'
      INQUIRE(FILE=TEMPO,EXIST=QEX)
      IF(.NOT.QEX) GO TO 2526
      CALL FILEDEL(24,TEMPO)
2526  CALL OPEN_WRITE1(24,TEMPO)
      GO TO 2028
2027  WRITE(20,*)'   Rhombohedral cells'
      TEMPO=FILE(1:LFILE)//'_rho.ckm'
      INQUIRE(FILE=TEMPO,EXIST=QEX)
      IF(.NOT.QEX) GO TO 2527
      CALL FILEDEL(24,TEMPO)
2527  CALL OPEN_WRITE1(24,TEMPO)
2028  CONTINUE
      WRITE(20,*)
      JJJ=0
      DO 2006 I=1,IGC
      IF(I.GT.20)GO TO 20060
      J=LL(I)
      IF(IFI(J).NE.JIFI)GO TO 2006
      JJJ=JJJ+1
      LLL(JJJ)=J
      IF(RP(J).LT.0.001)RP(J)=0.001
      X=1./RP(J)*5.
      VR=VGC(J)/VGC(LLL(1))
      WRITE(20,2001)RP(J),VGC(J),VR,KM(J),NSOL(J),(CEL(K,J),K=1,6)
      WRITE(24,2002)KM(J),X,VGC(J),VR,(CEL(K,J),K=1,6)
2006  CONTINUE
20060 NSOLMAX=NSOL(LLL(1))
      IF(NSOLMAX.GT.5)THEN
      WRITE(20,*)
      WRITE(20,*)'WARNING - WARNING - WARNING :'
      WRITE(20,*)'Same solution found Nsol = ',NSOLMAX,' times,'
      WRITE(20,*)'you should probably reduce the test numbers...'
      WRITE(20,*)
      ENDIF
      CLOSE(24)
2020  CONTINUE
C
C... Refine the "best" cell if this was not already done
C
C
C     READ hkl Miller indices in *.hkl
C
      IFILE=IFI(IBEST)
      GO TO (31,32,33,34,35,36,37)IFILE
31    TEMPO='cub.hkl'
      CALL OPEN_READ1(35,TEMPO)
      READ(35,*)NHKL0
      NHKL0=NDAT*6
      IF(NHKL0.GT.400)NHKL0=400
      DO 41 I=1,NHKL0
41    READ(35,*)(IHH(KK,I),KK=1,3)
      CLOSE(35)
      GO TO 50
32    TEMPO='hex.hkl'
      CALL OPEN_READ1(35,TEMPO)
      READ(35,*)NHKL0
      NHKL0=NDAT*12
      IF(NHKL0.GT.800)NHKL0=800
      DO 42 I=1,NHKL0
42    READ(35,*)(IHH(KK,I),KK=1,3)
      CLOSE(35)
      GO TO 50
33    TEMPO='tet.hkl'
      CALL OPEN_READ1(35,TEMPO)
      READ(35,*)NHKL0
      NHKL0=NDAT*12
      IF(NHKL0.GT.800)NHKL0=800
      DO 43 I=1,NHKL0
43    READ(35,*)(IHH(KK,I),KK=1,3)
      CLOSE(35)
      GO TO 50
34    TEMPO='ort.hkl'
      CALL OPEN_READ1(35,TEMPO)
      READ(35,*)NHKL0
      NHKL0=NDAT*20
      IF(NHKL0.GT.1000)NHKL0=1000
      DO 44 I=1,NHKL0
44    READ(35,*)(IHH(KK,I),KK=1,3)
      CLOSE(35)
      GO TO 50
35    TEMPO='mon.hkl'
      CALL OPEN_READ1(35,TEMPO)
      READ(35,*)NHKL0
      NHKL0=NDAT*20
      IF(NHKL0.GT.1000)NHKL0=1000
      DO 45 I=1,NHKL0
45    READ(35,*)(IHH(KK,I),KK=1,3)
      CLOSE(35)
      GO TO 50
36    TEMPO='tri.hkl'
      CALL OPEN_READ1(35,TEMPO)
      READ(35,*)NHKL0
      NHKL0=NDAT*20
      IF(NHKL0.GT.1000)NHKL0=1000
      DO 46 I=1,NHKL0
46    READ(35,*)(IHH(KK,I),KK=1,3)
      CLOSE(35)
      GO TO 50
37    TEMPO='rho.hkl'
      CALL OPEN_READ1(35,TEMPO)
      READ(35,*)NHKL0
      NHKL0=NDAT*12
      IF(NHKL0.GT.600)NHKL0=600
      DO 47 I=1,NHKL0
47    READ(35,*)(IHH(KK,I),KK=1,3)
      CLOSE(35)
C
50    CONTINUE
C
C      IF(IREF.EQ.1)GO TO 5900
      WRITE(20,*)
      WRITE(20,*)'    "Best" cell with largest McM20 :'
      WRITE(20,*)'    --------------------------------'
      WRITE(20,*)
      J=IBEST
      IFILE=IFI(J)
      INDIC=0
      IF(IFILE.EQ.1)INDIC=1
      IF(IFILE.EQ.2)INDIC=2
      IF(IFILE.EQ.3)INDIC=2
      IF(IFILE.EQ.7)INDIC=2
      BB(3)=CEL(1,J)
      BB(4)=CEL(2,J)
      BB(5)=CEL(3,J)
      BB(6)=CEL(4,J)
      BB(7)=CEL(5,J)
      BB(8)=CEL(6,J)
      AFI(3)=1.
      AFI(4)=1.
      AFI(5)=1.
      AFI(6)=1.
      AFI(7)=1.
      AFI(8)=1.
      IF(IFILE.EQ.1)AFI(6)=0.
      IF(IFILE.EQ.1)AFI(7)=0.
      IF(IFILE.EQ.1)AFI(8)=0.
      IF(IFILE.EQ.2)AFI(6)=0.
      IF(IFILE.EQ.2)AFI(7)=0.
      IF(IFILE.EQ.2)AFI(8)=0.
      IF(IFILE.EQ.3)AFI(6)=0.
      IF(IFILE.EQ.3)AFI(7)=0.
      IF(IFILE.EQ.3)AFI(8)=0.
      IF(IFILE.EQ.4)AFI(6)=0.
      IF(IFILE.EQ.4)AFI(7)=0.
      IF(IFILE.EQ.4)AFI(8)=0.
      IF(IFILE.EQ.5)AFI(6)=0.
      IF(IFILE.EQ.5)AFI(8)=0.
      IF(IFILE.EQ.7)AFI(6)=0.
      IF(IFILE.EQ.7)AFI(7)=0.
      IF(IFILE.EQ.7)AFI(8)=0.
      CELPRE(1)=BB(3)
      CELPRE(2)=BB(4)
      CELPRE(3)=BB(5)
      CELPRE(4)=BB(6)
      CELPRE(5)=BB(7)
      CELPRE(6)=BB(8)
      BB(1)=0.  ! zeropoint after correction... = 0.
      DO 6010 I=1,3
      DO 6010 JJ=1,3
6010  AL(I,JJ)=0.0
      CALL DCELL(CELPRE,AL,V1)
      CALL CALCUL2(DIFF,IHKL,TH3,NCALC,J)
      CALL CELREF(INDIC,BB,AFI,LHKL,TH3,IHKL,DDT,DDQ)
      IF(NDAT.GE.20)THEN
      CNCALC(IGC)=NCALC
      FM20(J)=QO(20)/(2.*CNCALC(J)*DDQ)
      FF20(J)=20./(CNCALC(J)*DDT)
      WRITE(20,7000)FM20(J)
      WRITE(20,7001)FF20(J),DDT,NCALC
      WRITE(20,*)
      PRINT 7000,FM20(J)
      PRINT 7001,FF20(J),DDT,NCALC
      PRINT *
      ENDIF
5900  CONTINUE
C
C   Make a .prf with the best solution
C
C
C.....START TO GENERATE THE "OBSERVED" PROFILE
C
C   Sum on observed peak positions
C
      XNHKL=NDAT
      HMIN=100.
      DO 1800 I=1,NDAT
      TH2(I)=TH2(I)+BB(1)
1800  FOBS(I)=FOBS(I)/SUM_F*XNHKL*50.
      SUM_F=XNHKL*50.
      FLOG =4.*ALOG(2.)
      SLOG =1./SQRT(FLOG)
      IF(W.LT.0.)THEN
      W=0.
      X=0.
      DO 1801 I=1,NHKL
      W=W+W1(I)
      X=X+1.
1801  CONTINUE
      W=W/X
      ENDIF
C
C  Use FWHM 1/2 the entered value... 
C
      W=(W/2.)**2
C
C
      DO 1802 NM=1,NHKL
      DD=SLABDA/(2.*SIN(TH2(NM)/PI))
      SSS=1/DD**2
C
C     Needs to calculate halfwidths and positions 
C
      D(NM)=SQRT(1./SSS)
      SINTH=SLABDA*SLABDA*SSS/4.
      COSTH=1.-SINTH
      TANTH=SQRT(SINTH/COSTH)
      THETA(NM)=ATAN(TANTH)*PI
      HW(NM)=SQRT(W)
      IF(HMIN.GT.HW(NM))HMIN=HW(NM)
      HW(NM)=HW(NM)*SLOG
      HW4(NM)=HW(NM)/SLOG*2.0
      BBB(NM)=HW(NM)*HW(NM)
1802  CONTINUE
      DMI=D(NHKL)-0.02
C
C  Step and positions
C  The minimum step is FWHMmin/STEPN
C
      STEP=HMIN/4.
      ASTEP=STEP*1000.
      ISTEP=ASTEP
      ASTEP=ISTEP
      STEP=ASTEP/1000.
      POS(1)=STEP
      PO=STEP
      DO 1803 I=2,16000
      PO=PO+STEP
      IF(PO.GT.160.)GO TO 1804
1803  POS(I)=PO
1804  NPTS=I-1
*  Sum at each point
      DO 1805 K=1,NPTS
      YOBS(K)=0.
*  SUM on all hkl
      DO 1806 J=1,NHKL
      DELTAP=POS(K)-THETA(J)
      IF(ABS(DELTAP).GT.HW4(J))GO TO 1806
      DELT=DELTAP*DELTAP
      OMEG=EXP(-DELT/BBB(J))/HW(J)
      YOBS(K)=YOBS(K)+OMEG*FOBS(J)
1806  CONTINUE
      IF(YOBS(K).LT.0.01)YOBS(K)=0.
1805  CONTINUE
      L=NPTS+1
      DO 1807 K=1,NPTS
      L=L-1
1807  IF(YOBS(L).NE.0.)GO TO 1808
1808  N2=L
      N1=1
C
      DO 1809 I=1,6
      K=I+2
1809  CELPRE(I)=BB(K)
      DO 55 I=1,3
      DO 55 J=1,3
55    AL(I,J)=0.0
      CALL DCELL(CELPRE,AL,V1)
C
C...  Keep only the hkl for d > dmin
C
      JH=0
      DO 59 I=1,NHKL0
      DO 56 KK=1,3
56    HH(KK)=IHH(KK,I)
C     HH IS INDICES OF REFLECTION
      X=0.
      DO 57 J=1,3
      DO 57 K=J,3
57    X=AL(J,K)*HH(J)*HH(K)+X
      X=1/SQRT(X)
      IF(X.LT.DMI)GO TO 59
      JH=JH+1
      DO 58 KK=1,3
58    JHH(KK,JH)=IHH(KK,I)
      FCAL(JH)=50.
      D(JH)=X
C     X IS D(hkl) FOR REFLECTION HH
59    CONTINUE
      NHKL=JH
C
C   Again, calculate 2-theta, etc.
C
      DO 60 NM=1,NHKL
C
      DD=D(NM)
      SSS=1/DD**2
C
C     Needs to calculate halfwidths and positions again
C
      SINTH=SLABDA*SLABDA*SSS/4.
      COSTH=1.-SINTH
      TANTH=SQRT(SINTH/COSTH)
      THETA(NM)=ATAN(TANTH)*PI
      HW(NM)=SQRT(W)
      IF(HMIN.GT.HW(NM))HMIN=HW(NM)
      HW(NM)=HW(NM)*SLOG
      HW4(NM)=HW(NM)/SLOG*2.0
      BBB(NM)=HW(NM)*HW(NM)
60    CONTINUE
C
C...  Calculate best Yobs
C
*  Sum at each point
*
      DO 1810 K=1,NHKL
      SOMEGA(K)=0.
1810  FOBS(K)=0.
      DO 1811 K=1,N2
      YCALC(K)=0.
*  SUM on all hkl
      KPOS=0
      OMEGT=0.
      DO 1812 J=1,NHKL
      DELTAP=POS(K)-THETA(J)
      IF(ABS(DELTAP).GT.HW4(J))GO TO 1812
      IF(KPOS.EQ.0)NHA(K)=J
      KPOS=1
      NHB(K)=J
      DELT=DELTAP*DELTAP
      OMEG=EXP(-DELT/BBB(J))/HW(J)
      SOMEGA(J)=SOMEGA(J)+OMEG
      DUMP(J)=FCAL(J)*OMEG
      YCALC(K)=YCALC(K)+DUMP(J)
1812  CONTINUE
      IF(YCALC(K).EQ.0.)GO TO 1813
      YOY=YOBS(K)/YCALC(K)
      DO 1814 J=NHA(K),NHB(K)
      FOBS(J)=FOBS(J)+DUMP(J)*YOY
1814  CONTINUE
1813  CONTINUE
1811  CONTINUE
      DO 1815 K=1,NHKL
      IF(SOMEGA(K).EQ.0.)GO TO 1816
      FOBS(K)=FOBS(K)/SOMEGA(K)
      GO TO 1815
1816  FOBS(K)=0.
1815  CONTINUE
      DO 61 J=1,NHKL
      FCAL(J)=FOBS(J)
61    CONTINUE
C
C...  2 more iterations by Le Bail fit
C
      DO 1850 KITER=1,2
*  Sum at each point
      DO 1817 K=1,NHKL
      SOMEGA(K)=0.
1817  FOBS(K)=0.
      DO 1818 K=1,N2
      YCALC(K)=0.
*  SUM on all hkl
      OMEGT=0.
      DO 1819 J=NHA(K),NHB(K)
      DELTAP=POS(K)-THETA(J)
      IF(ABS(DELTAP).GT.HW4(J))GO TO 1819
      DELT=DELTAP*DELTAP
      OMEG=EXP(-DELT/BBB(J))/HW(J)
      SOMEGA(J)=SOMEGA(J)+OMEG
      DUMP(J)=FCAL(J)*OMEG
      YCALC(K)=YCALC(K)+DUMP(J)
1819  CONTINUE
      IF(YCALC(K).EQ.0.)GO TO 1820
      YOY=YOBS(K)/YCALC(K)
      DO 1821 J=NHA(K),NHB(K)
      FOBS(J)=FOBS(J)+DUMP(J)*YOY
1821  CONTINUE
1820  CONTINUE
1818  CONTINUE
      DO 1822 K=1,NHKL
      IF(SOMEGA(K).EQ.0.)GO TO 1823
      FOBS(K)=FOBS(K)/SOMEGA(K)
      GO TO 1822
1823  FOBS(K)=0.
1822  CONTINUE
      DO 62 J=1,NHKL
      FCAL(J)=FOBS(J)
62    CONTINUE
1850  CONTINUE
C
C
C
      DIFF=0.
*  Sum at each point
      DO 1824 K=1,N2
      YCALC(K)=0.
      IF(YOBS(K).EQ.0.)GO TO 1824
*  SUM on all hkl
      OMEGT=0.
      DO 1825 J=NHA(K),NHB(K)
      DELTAP=POS(K)-THETA(J)
      IF(ABS(DELTAP).GT.HW4(J))GO TO 1825
      DELT=DELTAP*DELTAP
      OMEG=EXP(-DELT/BBB(J))/HW(J)
      YCALC(K)=YCALC(K)+FCAL(J)*OMEG
1825  CONTINUE
1824  CONTINUE
      SUM_Y=0.
      DO 1826 K=1,N2
      DIFF=DIFF+ABS(YOBS(K)-YCALC(K))
      SUM_Y=SUM_Y+YOBS(K)
1826  CONTINUE
      DIFF=DIFF/SUM_Y
      WRITE(20,*)
      WRITE(20,*)' Final Rp on the .prf = ',DIFF
      WRITE(20,*)
C
C     Make the .prf
C
      TEMPO=FILE(1:LFILE)//'.prf'
      INQUIRE(FILE=TEMPO,EXIST=QEX)
      IF(.not.QEX) go to 1827
      CALL FILEDEL(12,TEMPO)
1827  CALL OPEN_WRITE1(12,TEMPO)
      IVERS=8
      ZERO=0.
      WRITE(12,'(20A4)')(TEXT(I),I=1,20)
      WRITE(12,'(A,F9.4)') '3111   1.0000',ZERO
      THMAX=POS(N2)
      THMIN=POS(N1)
      ICN=NHKL
      DO 1828 I=1,NHKL
      I1=JHH(1,I)
      I2=JHH(2,I)
      I3=JHH(3,I)
1828  IREFS(I)=256*(256*(256*(8+1)+128+I1)+128+I2)+128+I3
      ICZ=NHKL
      NPTS=N2-N1+1
      WRITE(12,'(2F8.3,F8.5,I8)') THMAX,THMIN,STEP,IVERS
      AMDA1=SLABDA
      AMDA2=SLABDA
      NPAT1=1
      NVK=0
      WRITE(12,'(I3,I7,5f10.5)')NPAT1,NPTS,AMDA1,AMDA2,ZERO
      WRITE(12,'(8I5,8i3)')ICN,NVK
      WRITE(12,'(8(F7.0,1X))') (YOBS(J),J=N1,N2)
      WRITE(12,'(8(F7.0,1X))') (YCALC(J),J=N1,N2)
      WRITE(12,'(8I10)') (IREFS(J),J=1,ICZ)
      WRITE(12,'(10F8.3)')(THETA(J),J=1,ICZ)
      NEXCRG=0
      EXCRG=0
      WRITE(12,'(I8)') NEXCRG
      WRITE(12,'(2F8.2)') EXCRG,EXCRG
      CLOSE(12)
C
C
C
      IF(IGC.EQ.1)GO TO 6000
C
C   Sort cells by volume
C
      WRITE(20,*)
      WRITE(20,*)'====================================================
     1==========================='
      WRITE(20,*)'       CELL PROPOSALS sorted by increasing volume :'
      WRITE(20,*)'====================================================
     1==========================='
      WRITE(20,*)
      WRITE(20,*)
      CALL SORT(IGC,VGC,LL)
      IMEM=IMEM+1
      IM(IMEM)=LL(1)
      WRITE(20,1999)
      DO 2003 I=1,IGC
      IF(I.GT.20)GO TO 20030
      J=LL(I)
      IF(RP(J).LT.0.001)RP(J)=0.001
      VR=VGC(J)/VGC(LL(1))
      WRITE(20,2001)RP(J),VGC(J),VR,KM(J),NSOL(J),(CEL(K,J),K=1,6)
2003  CONTINUE
20030 NSOLMAX=NSOL(LL(1))
      IF(NSOLMAX.GT.5)THEN
      WRITE(20,*)
      WRITE(20,*)'WARNING - WARNING - WARNING :'
      WRITE(20,*)'Same solution found Nsol = ',NSOLMAX,' times,'
      WRITE(20,*)'you should probably reduce the test numbers...'
      WRITE(20,*)
      ENDIF
C
C... Refine the "best" cell if this was not already done
C          this time with smallest volume and if Rp < 10%
C
      IF(IREF.EQ.1)GO TO 5901
      J=LL(1)
      IF(RP(J).GT.0.10)GO TO 5901
C
C     READ hkl Miller indices in *.hkl
C
      IFILE=IFI(LL(1))
      GO TO (1731,1732,1733,1734,1735,1736,1737)IFILE
1731  TEMPO='cub.hkl'
      CALL OPEN_READ1(35,TEMPO)
      READ(35,*)NHKL0
      NHKL0=NDAT*6
      IF(NHKL0.GT.400)NHKL0=400
      DO 1741 I=1,NHKL0
1741  READ(35,*)(IHH(KK,I),KK=1,3)
      CLOSE(35)
      GO TO 1750
1732  TEMPO='hex.hkl'
      CALL OPEN_READ1(35,TEMPO)
      READ(35,*)NHKL0
      NHKL0=NDAT*12
      IF(NHKL0.GT.800)NHKL0=800
      DO 1742 I=1,NHKL0
1742  READ(35,*)(IHH(KK,I),KK=1,3)
      CLOSE(35)
      GO TO 1750
1733  TEMPO='tet.hkl'
      CALL OPEN_READ1(35,TEMPO)
      READ(35,*)NHKL0
      NHKL0=NDAT*12
      IF(NHKL0.GT.800)NHKL0=800
      DO 1743 I=1,NHKL0
1743  READ(35,*)(IHH(KK,I),KK=1,3)
      CLOSE(35)
      GO TO 1750
1734  TEMPO='ort.hkl'
      CALL OPEN_READ1(35,TEMPO)
      READ(35,*)NHKL0
      NHKL0=NDAT*20
      IF(NHKL0.GT.1000)NHKL0=1000
      DO 1744 I=1,NHKL0
1744  READ(35,*)(IHH(KK,I),KK=1,3)
      CLOSE(35)
      GO TO 1750
1735  TEMPO='mon.hkl'
      CALL OPEN_READ1(35,TEMPO)
      READ(35,*)NHKL0
      NHKL0=NDAT*20
      IF(NHKL0.GT.1000)NHKL0=1000
      DO 1745 I=1,NHKL0
1745  READ(35,*)(IHH(KK,I),KK=1,3)
      CLOSE(35)
      GO TO 1750
1736  TEMPO='tri.hkl'
      CALL OPEN_READ1(35,TEMPO)
      READ(35,*)NHKL0
      NHKL0=NDAT*20
      IF(NHKL0.GT.1000)NHKL0=1000
      DO 1746 I=1,NHKL0
1746  READ(35,*)(IHH(KK,I),KK=1,3)
      CLOSE(35)
      GO TO 1750
1737  TEMPO='rho.hkl'
      CALL OPEN_READ1(35,TEMPO)
      READ(35,*)NHKL0
      NHKL0=NDAT*12
      IF(NHKL0.GT.600)NHKL0=600
      DO 1747 I=1,NHKL0
1747  READ(35,*)(IHH(KK,I),KK=1,3)
      CLOSE(35)
C
C
1750  CONTINUE
      WRITE(20,*)
      WRITE(20,*)'    "Best" cell with smallest volume :'
      WRITE(20,*)'    ---------------------------------'
      WRITE(20,*)
      IFILE=IFI(J)
      INDIC=0
      IF(IFILE.EQ.1)INDIC=1
      IF(IFILE.EQ.2)INDIC=2
      IF(IFILE.EQ.3)INDIC=2
      IF(IFILE.EQ.7)INDIC=2
      BB(3)=CEL(1,J)
      BB(4)=CEL(2,J)
      BB(5)=CEL(3,J)
      BB(6)=CEL(4,J)
      BB(7)=CEL(5,J)
      BB(8)=CEL(6,J)
      AFI(3)=1.
      AFI(4)=1.
      AFI(5)=1.
      AFI(6)=1.
      AFI(7)=1.
      AFI(8)=1.
      IF(IFILE.EQ.1)AFI(6)=0.
      IF(IFILE.EQ.1)AFI(7)=0.
      IF(IFILE.EQ.1)AFI(8)=0.
      IF(IFILE.EQ.2)AFI(6)=0.
      IF(IFILE.EQ.2)AFI(7)=0.
      IF(IFILE.EQ.2)AFI(8)=0.
      IF(IFILE.EQ.3)AFI(6)=0.
      IF(IFILE.EQ.3)AFI(7)=0.
      IF(IFILE.EQ.3)AFI(8)=0.
      IF(IFILE.EQ.4)AFI(6)=0.
      IF(IFILE.EQ.4)AFI(7)=0.
      IF(IFILE.EQ.4)AFI(8)=0.
      IF(IFILE.EQ.5)AFI(6)=0.
      IF(IFILE.EQ.5)AFI(8)=0.
      IF(IFILE.EQ.7)AFI(6)=0.
      IF(IFILE.EQ.7)AFI(7)=0.
      IF(IFILE.EQ.7)AFI(8)=0.
      CELPRE(1)=BB(3)
      CELPRE(2)=BB(4)
      CELPRE(3)=BB(5)
      CELPRE(4)=BB(6)
      CELPRE(5)=BB(7)
      CELPRE(6)=BB(8)
      BB(1)=0.  ! zeropoint after correction... = 0.
      DO 6011 I=1,3
      DO 6011 JJ=1,3
6011  AL(I,JJ)=0.0
      CALL DCELL(CELPRE,AL,V1)
      CALL CALCUL2(DIFF,IHKL,TH3,NCALC,J)
      CALL CELREF(INDIC,BB,AFI,LHKL,TH3,IHKL,DDT,DDQ)
      IF(NDAT.GE.20)THEN
      CNCALC(J)=NCALC
      FM20(J)=QO(20)/(2.*CNCALC(J)*DDQ)
      FF20(J)=20./(CNCALC(J)*DDT)
      WRITE(20,7000)FM20(J)
      WRITE(20,7001)FF20(J),DDT,NCALC
      WRITE(20,*)
      PRINT 7000,FM20(J)
      PRINT 7001,FF20(J),DDT,NCALC
      PRINT *
      ENDIF
5901  CONTINUE
C
C
C
      IF(IGC.EQ.1)GO TO 6000
C
C   Sort cells, the most frequently found first
C
      WRITE(20,*)
      WRITE(20,*)'====================================================
     1==========================='
      WRITE(20,*)'       CELL PROPOSALS most frequently found :'
      WRITE(20,*)'====================================================
     1==========================='
      WRITE(20,*)
      WRITE(20,*)
      CALL SORT2(IGC,NSOL,LL)
      IMEM=IMEM+1
      IM(IMEM)=LL(IGC)
      WRITE(20,1999)
      DO 2004 I=1,IGC
      IF(I.GT.20)GO TO 20040
      J=LL(IGC+1-I)
      IF(NSOL(J).LT.2)GO TO 2004
      IF(RP(J).LT.0.001)RP(J)=0.001
      VR=VGC(J)/VGC(LL(IGC))
      WRITE(20,2001)RP(J),VGC(J),VR,KM(J),NSOL(J),(CEL(K,J),K=1,6)
2004  CONTINUE
20040 CONTINUE
C
C
C   Sort cells with largest number of peak indexed + small Rp2
C
      WRITE(20,*)
      WRITE(20,*)'====================================================
     1==========================='
      WRITE(20,*)'CELLS with small Rp2 + largest number of peak indexed'
      WRITE(20,*)'====================================================
     1==========================='
      WRITE(20,*)'Rp2 is on peak indexed only, and width divided by 2,'
      WRITE(20,*)'while Rp is on all peaks, and large width.'
      WRITE(20,*)
      CALL SORT2(IGC,RP2,LL)
      IMEM=IMEM+1
      IM(IMEM)=LL(I)
      WRITE(20,19992)
      DO 2005 I=1,IGC
      IF(I.GT.20)GO TO 20050
      J=LL(I)
      IF(RP2(J).LT.0.001)RP2(J)=0.001
      IF(RP2(J).GT.0.30) GO TO 2005
      IF(KM2(J).LT.NMAX) GO TO 2005
      VR=VGC(J)/VGC(LL(IGC))
      WRITE(20,2001)RP2(J),VGC(J),VR,KM2(J),NSOL(J),(CEL(K,J),K=1,6)
2005  CONTINUE
20050 CONTINUE
C
C
C   Sort associations of two cells with largest number of peak indexed
C
      IF(RMAX0(1).LT.0.5)GO TO 6000
      WRITE(20,*)
      WRITE(20,*)'====================================================
     1==========================='
      WRITE(20,*)'Double cells with largest number of peak indexed'
      WRITE(20,*)'====================================================
     1==========================='
      WRITE(20,*)
      WRITE(20,*)'WARNING - WARNING - WARNING - WARNING - WARNING'
      WRITE(20,*)'           This is the two-phase mode'
      WRITE(20,*)'    It could be better to go back to the lab' 
      WRITE(20,*)'         and try and make a pure sample'
      WRITE(20,*)
      II=0
      DO 2009 I=1,IGC-1
      IF(RP2(I).GT.0.30) GO TO 2009
      IF(KM2(I).LT.NMAX) GO TO 2009
      DO 2008 J=I+1,IGC
      IF(RP2(J).GT.0.30) GO TO 2008
      IF((RP2(I)+RP2(J)).GT.0.40) GO TO 2008
      IF(KM2(J).LT.NMAX) GO TO 2008
      II=II+1
      IF(II.GT.99999)GO TO 2011
      KM3(II)=0
      DO 2007 K=1,NDAT
      KM3(II)=KM3(II)+IND(K,I)+IND(K,J)
      IF((IND(K,I)*IND(K,J)).EQ.1)KM3(II)=KM3(II)-1
2007  CONTINUE
      IF(KM3(II).LT.(NDAT-10))THEN
      II=II-1
      GO TO 2008
      ENDIF
      ID1(II)=I
      ID2(II)=J
2008  CONTINUE
2009  CONTINUE
2011  IGC2=II
      CALL SORT3(IGC2,KM3,LL2)
      WRITE(20,19992)
      VR=1.
      TEMPO=FILE(1:LFILE)//'_two.ckm'
      INQUIRE(FILE=TEMPO,EXIST=QEX)
      IF(.NOT.QEX) GO TO 2530
      CALL FILEDEL(24,TEMPO)
2530  CALL OPEN_WRITE1(24,TEMPO)
      DO 2010 I=1,IGC2
      IF(I.GT.1000)GO TO 2012
      JJ=LL2(IGC2+1-I)
      J=ID1(LL2(IGC2+1-I))
      IF(RP2(J).LT.0.001)RP2(J)=0.001
      X=1./RP2(J)*5.
      WRITE(20,2001)RP2(J),VGC(J),VR,KM3(JJ),NSOL(J),(CEL(K,J),K=1,6)
      WRITE(24,2002)KM3(JJ),X,VGC(J),VR,(CEL(K,J),K=1,6)
      J=ID2(LL2(IGC2+1-I))
      IF(RP2(J).LT.0.001)RP2(J)=0.001
      X=1./RP2(J)*5.
      WRITE(20,2001)RP2(J),VGC(J),VR,KM2(J),NSOL(J),(CEL(K,J),K=1,6)
      WRITE(24,2002)KM2(J),X,VGC(J),VR,(CEL(K,J),K=1,6)
      WRITE(20,*)
2010  CONTINUE
2012  CONTINUE
      CLOSE(24)
C
C
6000  CONTINUE
      WRITE(20,*)
      WRITE(20,*)
      WRITE(20,*)
      WRITE(20,*)
      WRITE(20,*)
      WRITE(20,*)
      WRITE(20,*)
      WRITE(20,*)
      WRITE(20,*)
      WRITE(20,*)
      WRITE(20,*)
      WRITE(20,*)
      WRITE(20,*)
      WRITE(20,*)
      WRITE(20,*)
      WRITE(20,*)
      WRITE(20,*)
      WRITE(20,*)
      WRITE(20,*)
      WRITE(20,*)'====================================================
     1==========================='
      WRITE(20,*)'         THE SELECTION OF THE "BEST" CELL'
      WRITE(20,*)'based on McM20, Rp, F(20), M(20), V, high symmetry ?'
      WRITE(20,*)'           DEPENDS ON YOU, EXCLUSIVELY.'
      WRITE(20,*)
      WRITE(20,*)'                    However...'
      WRITE(20,*)'  It is suggested that the correct cell could be :'
      WRITE(20,*)'====================================================
     1==========================='
      WRITE(20,*)
      WRITE(20,19996)
19996 FORMAT('IN  F.o.M.    Volume             a        b        c      
     1 alpha   beta    gamma   Bravais lattice')
C
C  Ultimate analysis...
C
      DO 20010 I=1,IMEM
      IMN(I)=1
20010 CONTINUE
      IF(IMEM.GT.1)THEN
      DO 20008 I=1,IMEM-1
      IF(IMN(I).EQ.0)GO TO 20008
      DO 30008 K=I+1,IMEM
      IF(IMN(K).EQ.0)GO TO 30008
      IF(IM(K).EQ.IM(I))THEN
      IMN(I)=IMN(I)+1
      IMN(K)=0
      ENDIF
30008 CONTINUE
20008 CONTINUE
      ENDIF
      J=IM(1)
      IF(IB(J).EQ.1)BL='I'
      IF(IB(J).EQ.2)BL='A'
      IF(IB(J).EQ.3)BL='B'
      IF(IB(J).EQ.4)BL='C'
      IF(IB(J).EQ.5)BL='F'
      IF(IB(J).EQ.6)BL='P'
      IF(IFI(J).EQ.7)BL='R'
      IF(IFI(J).EQ.1)MORE='Cubic *****'
      IF(IFI(J).EQ.2)MORE='Hexag **** '
      IF(IFI(J).EQ.3)MORE='Tetra **** '
      IF(IFI(J).EQ.4)MORE='Ortho ***  '
      IF(IFI(J).EQ.5)MORE='           '
      IF(IFI(J).EQ.6)MORE='           '
      IF(IFI(J).EQ.7)MORE='Rhomb **** '
      WRITE(20,20032)KM(J),XFOM(J),VGC(J),(CEL(K,J),K=1,6),BL,MORE
      WRITE(20,*)'   Found ',IMN(1),' time(s) head of the best lists'
20032 FORMAT(I2,F8.2,F11.3,9X,3F9.4,3F8.3,4X,A1,2X,A11)
C
      IF(IMEM.GT.1)THEN
      XF1=XFOM(IM(1))/2.
      IMEMT=IMEM
      DO 20011 I=2,IMEM
      IF(IMN(I).EQ.0)THEN
      IMEMT=IMEMT-1
      GO TO 20011
      ENDIF
      J=IM(I)
      IF(XFOM(J).LT.XF1)THEN
      IMEMT=IMEMT-1
      IMN(I)=0
      ENDIF
20011 CONTINUE
      IF(IMEMT.GT.1)THEN
      WRITE(20,*)
      WRITE(20,*)
      WRITE(20,*)'   Other(s) having some chance :'
      WRITE(20,*)
      DO 20009 I=2,IMEM
      IF(IMN(I).EQ.0)GO TO 20009
      J=IM(I)
      IF(IB(J).EQ.1)BL='I'
      IF(IB(J).EQ.2)BL='A'
      IF(IB(J).EQ.3)BL='B'
      IF(IB(J).EQ.4)BL='C'
      IF(IB(J).EQ.5)BL='F'
      IF(IB(J).EQ.6)BL='P'
      IF(IFI(J).EQ.7)BL='R'
      IF(IFI(J).EQ.1)MORE='Cubic *****'
      IF(IFI(J).EQ.2)MORE='Hexag **** '
      IF(IFI(J).EQ.3)MORE='Tetra **** '
      IF(IFI(J).EQ.4)MORE='Ortho ***  '
      IF(IFI(J).EQ.5)MORE='           '
      IF(IFI(J).EQ.6)MORE='           '
      IF(IFI(J).EQ.7)MORE='Rhomb **** '
      WRITE(20,20032)KM(J),XFOM(J),VGC(J),(CEL(K,J),K=1,6),BL,MORE
      WRITE(20,*)'   Found ',IMN(I),' time(s) head of the best lists'
20009 CONTINUE
      ENDIF
      ENDIF
      WRITE(20,*)
      WRITE(20,*)'====================================================
     1==========================='
      WRITE(20,*)
      WRITE(20,*)
      WRITE(20,*)
      WRITE(20,*)
      WRITE(20,*)
      WRITE(20,*)
C
C
C
7000  FORMAT('   M(20) = ',F9.2)
7001  FORMAT('   F(20) = ',F9.2,' (',F8.4,',',I4,')')
10115 FORMAT(2X,2I10)
115   FORMAT(F5.3,F12.0,F8.4,F9.1,2I3)
215   FORMAT(F5.3,F12.0,2F8.4,F9.1,2I3)
415   FORMAT(F5.3,F12.0,3F8.4,F9.1,2I3)
515   FORMAT(F5.3,F12.0,3F8.4,F7.2,F9.1,2I3)
615   FORMAT(F5.3,F12.0,3F8.4,3F7.2,F9.1,2I3)
C1115  FORMAT(I3,F10.0,1X,F5.3,F8.4,F9.1,I3)
C1215  FORMAT(I3,F10.0,1X,F5.3,2F8.4,F9.1,I3)
C1415  FORMAT(I3,F10.0,1X,F5.3,3F8.4,F9.1,I3)
C1515  FORMAT(I3,F10.0,1X,F5.3,3F8.4,F7.2,F9.1,I3)
C1615  FORMAT(I3,F10.0,1X,F5.3,3F8.4,3F7.2,F9.1,I3)
1115  FORMAT(14X,F5.3,F8.4,F9.1,I3)
1215  FORMAT(14X,F5.3,2F8.4,F9.1,I3)
1415  FORMAT(14X,F5.3,3F8.4,F9.1,I3)
1515  FORMAT(14X,F5.3,3F8.4,F7.2,F9.1,I3)
1615  FORMAT(14X,F5.3,3F8.4,3F7.2,F9.1,I3)
3000  FORMAT(/1x,'=============================================
     -===================='/
     -        1x,'McMaille version ',a,' by A. Le Bail - 2006 -
     - alb@cristal.org'/
     -        1x,'=============================================
     -===================='//
     -        1x,'Using generic filename : ',A)
      GO TO 3600
3500  WRITE(20,*)'  Error reading the first line : TEXT'
      WRITE (BUFFER,*)'  Error reading the first line : TEXT'
      CALL PROGRESSVIEW (BUFFER)
      GO TO 3600
3501  WRITE(20,*)'  Error reading lambda, etc, line 2'
      WRITE (BUFFER,*)'  Error reading lambda, etc, line 2'
      CALL PROGRESSVIEW (BUFFER)
      GO TO 3600
3502  WRITE(20,*)'  Error reading symmetry code, line 3'
      WRITE (BUFFER,*)'  Error reading symmetry code, line 3'
      CALL PROGRESSVIEW (BUFFER)
      GO TO 3600
3503  WRITE(20,*)'  Error reading U, V, W, step'
      WRITE (BUFFER,*)'  Error reading U, V, W, Step'
      CALL PROGRESSVIEW (BUFFER)
      GO TO 3600
3504  WRITE(20,*)'  Error reading Pmin, Pmax, Vmin'
      WRITE (BUFFER,*)'  Error reading Pmin, Pmax, Vmin'
      CALL PROGRESSVIEW (BUFFER)
      GO TO 3600
3505  WRITE(20,*)'  Error reading grid steps'
      WRITE (BUFFER,*)'  Error reading grid steps'
      CALL PROGRESSVIEW (BUFFER)
      GO TO 3600
3506  WRITE(20,*)'  Error reading NTIMELIM'
      WRITE (BUFFER,*)'  Error reading NTIMELIM'
      CALL PROGRESSVIEW (BUFFER)
      GO TO 3600
3507  WRITE(20,*)'  Error reading nstart, rmax, test'
      WRITE (BUFFER,*)'  Error reading nstart, rmax, test'
      CALL PROGRESSVIEW (BUFFER)
      GO TO 3600
3508  WRITE(20,*)'  Error reading data angle > 180'
      WRITE (BUFFER,*)'  Error reading data angle > 180'
      CALL PROGRESSVIEW (BUFFER)
3600  CONTINUE
      CALL DATN(DATENOW,TIMENOW)
      CALL CPU_TIME(time_end)
      totaltime=time_end-time_begin
      write(20,*)' Total CPU time elapsed in seconds : ',totaltime
      WRITE(*,3955)
3955  FORMAT(//'  Type any character and a RETURN to continue : ',$)
      READ(*,*)FEND
      CLOSE(20)
      CLOSE(28)
      STOP
      END
*-----------------------------------------------------------------------
      SUBROUTINE PROGRESSVIEW(BUFFER)
      CHARACTER*79 BUFFER
      WRITE(*,1)BUFFER
1     FORMAT(A79)
      RETURN
      END
*-----------------------------------------------------------------------
      SUBROUTINE DATN(DATENOW,TIMENOW)
C
C Obviously, this subroutine is compiler specific...
C
C In case of problem, remove that subroutine and all 
C the CALL DATN inside the program.
C
C Date will be printed as follows :
C
C     10-Apr-2000   Hour:  0 Min: 29 Sec: 20 
C
C     dd-mmm-yyyy          hh     mm      ss 
C
C     DT_VALUES(3),MONTHS(DT_VALUES(2)),DT_VALUES(1), 
C     DT_VALUES(5),DT_VALUES(6),DT_VALUES(7) 
C
C
      CHARACTER*7 DATENOW
      CHARACTER*8 TIMENOW
      CHARACTER*8 DA
      CHARACTER*10 TIME
      CHARACTER*5 ZONE
      CHARACTER*2 YEAR
      INTEGER*4 DT_VALUES(8)    ! Return values from DATE_AND_TIME
      CHARACTER*3 MONTHS(12)    ! Month names
      DATA MONTHS /'Jan','Feb','Mar','Apr','May','Jun',
     1             'Jul','Aug','Sep','Oct','Nov','Dec'/

C
C Get current date values
C
       CALL DATE_AND_TIME (DA,TIME,ZONE,DT_VALUES) 
C
C Format date
C
      WRITE(20,*)
      WRITE(20,10) DT_VALUES(3),
     1 MONTHS(DT_VALUES(2)),DT_VALUES(1),DT_VALUES(5),
     2 DT_VALUES(6),DT_VALUES(7)
      WRITE(20,*)
10    FORMAT(1X,I2,'-',A3,'-',I4,5X,I2,' hour ',I2,' min ',I2,
     -' Sec ')
      IF(DT_VALUES(1).EQ.2006)YEAR='06'
      IF(DT_VALUES(1).EQ.2007)YEAR='07'
      IF(DT_VALUES(1).EQ.2008)YEAR='08'
      IF(DT_VALUES(1).EQ.2009)YEAR='09'
      IF(DT_VALUES(1).EQ.2010)YEAR='10'
      WRITE(DATENOW,'(I2,A3,A2)')
     1DT_VALUES(3),MONTHS(DT_VALUES(2)),YEAR
      WRITE(TIMENOW,'(I2,A,I2,A,I2)')
     1DT_VALUES(5),':',DT_VALUES(6),':',DT_VALUES(7)
C
C  End of compiler specific subroutine
C
      RETURN
      END
      SUBROUTINE OPEN_READ1(UNIT,FILE)
      integer UNIT
      CHARACTER*(*) FILE
      OPEN (UNIT,FILE=FILE,STATUS='old')
      RETURN
      END
*
      SUBROUTINE OPEN_WRITE1(UNIT,FILE)
      integer UNIT
      CHARACTER*(*) FILE
      OPEN (UNIT,FILE=FILE,STATUS='NEW')
      RETURN
      END
*
      SUBROUTINE FILEDEL(UNIT,FILE)
      integer UNIT
      CHARACTER*(*) FILE
      OPEN (UNIT,FILE=FILE,STATUS='old')
      close(UNIT,STATUS='delete')
      RETURN
      END
*
*----------------------------------------------------------------------- 
*
*     SUBROUTINE OPEN_WRITE OPENs a FILE for writing, filling in the 
*     supplied extension IF none is supplied with the FILE name.
*     Allows use of system specIFic facilities of OPEN statement.
*
      SUBROUTINE OPEN_WRITE(UNIT,FILE,extension)
      integer UNIT
      CHARACTER*(*) FILE,extension,temp*80
      temp=FILE
      i=index(FILE,'.')
      l=LEN(FILE)
      DO WHILE(FILE(l:l).EQ.' ')
         l=l-1
      ENDDO
      IF (i.EQ.0) temp=FILE(1:l)//extension
      OPEN (UNIT,FILE=temp,STATUS='old',err=10)
      close(UNIT,STATUS='delete')
10    OPEN (UNIT,FILE=temp,STATUS='NEW')
      RETURN
      END
*
*----------------------------------------------------------------------- 
*
*     SUBROUTINE ESP_INIT must set the seed for the random number 
*     generator and obtain the current cpu clock reading in seconds 
*
      SUBROUTINE ESP_INIT(ISEED) 
C
C PORTLIB/DFPORT is compiler spcecific part, introduced for 
C using the intrinsic function SECNDS(X) which returns the 
C  (time in seconds since midnight - X) 
C
C      USE PORTLIB
!      USE IFPORT   ! Disabled by Zucco
      X=0.
      ISEED=INT(SECNDS(X)*100.)+1
*
*  run the random number generator N times
*  for avoiding effects of starting value
*
      N=ISEED/2000
      IF(N.LE.0)N=100
      IF(N.GE.1000)N=1000
      DO 10 I=1,N
      BIDON=RANDI(ISEED)
10    ISEED=ISEED/3
      ISEED=ISEED*2+1
C      WRITE(20,*)' ISEED = ',ISEED
      RETURN
      END
*
*********************************************************************** 
*
      SUBROUTINE DCELL(CELLN,AL,V)
      DIMENSION CELLN(6),AL(3,3),RCELLN(6),CELL(6)
      DEGRAD=3.14159/180.0
      DO 30 I=1,6
30    CELL(I)=CELLN(I)
      DO 31 I=1,3
      L=I+3
      IF(CELL(L)-90.0)32,33,32
33    CELL(L)=0.0
      GOTO 31
32    CELL(L)=COS(DEGRAD*CELL(L))
31    CONTINUE
      CALL TRCL(CELL,RCELLN,V)
C     RCELL IS THE RECIPROCAL CELL CONSTANTS
      DO 34 I=1,3
      AL(I,I)=RCELLN(I)*RCELLN(I)
      CALL PERM(I,J,K)
      IF(J-K)35,36,36
35    AL(J,K)=2.0*RCELLN(J)*RCELLN(K)*RCELLN(I+3)
      GOTO 34
36    AL(K,J)=2.0*RCELLN(J)*RCELLN(K)*RCELLN(I+3) 
34    CONTINUE
C      DO 37 I=4,6
C37    RCELLN(I)=ACOS(RCELLN(I))/DEGRAD
      RETURN
      END
*
*----------------------------------------------------------------------- 
*
      SUBROUTINE TRCL(CELLN,RCELLN,V)
C     TRANSFORMS REAL CELL TO RECIPROCAL OR VICE VERSA 
C     INPUT CELL IS IN ARRAY CELL AS LENGTHS AND COSINES
      DIMENSION CELLN(6),RCELLN(6)
      DIMENSION SINA(3)
      ABC=1.0
      PROD=2.0
      V=-2.0
      DO 10 I=1,3
      L=I+3
      SINA(I)=1.0-CELLN(L)**2
      V=V+SINA(I)
      SINA(I)=SQRT(SINA(I))
      PROD=PROD*CELLN(L)
10    ABC=ABC*CELLN(I)
      V=ABC*SQRT(V+PROD)
C      V IS CELL VOLUME
C     PUT INVERTED CELL INTO RCELL
      DO 20 I=1,3
      CALL PERM(I,J,K)
      RCELLN(I)=CELLN(J)*CELLN(K)*SINA(I)/V
      L=I+3
20    RCELLN(L)=(CELLN(J+3)*CELLN(K+3)-CELLN(L))/(SINA(J)
     2 *SINA(K))
      RETURN
      END
*
*----------------------------------------------------------------------- 
*
      SUBROUTINE PERM(I,J,K)
C     PERMS USEFUL COMBINATIONS OF INTEGERS IN THE RANGE 1 TO 3
      IF(I-2)10,20,30
10    J=2
      K=3
      RETURN
20    J=3
      K=1
      RETURN
30    J=1
      K=2
      RETURN
      END
*
********************************************************* 
*
      SUBROUTINE SORT(N,A,L)
C     ********************************************************
C     THE SUBROUTINE SORT APPLIES TO THE ARRAY A WITH J ELEMENTS.
C     THE INDICES OF THE ORDERED ARRAY ARE PLACED IN ARRAY L
C     TO OBTAIN THE ORDERED ARRAY REPLACE THE INDICES I OF A(I)
C     WITH I=L(I)---REF.CACM 271
C     *********************************************************
C
      PARAMETER(N_HKL=10000)
      DIMENSION IUT(10),ILT(10)
      DIMENSION A(N_HKL),L(N_HKL)
      J=N
      I=1
      M=1
      DO 10 K=1,J
10    L(K)=K
20    IF(J-I-1)140,140,30
30    IP=(J+I)/2
      ITT=L(IP)
      T=A(ITT)
      L(IP)=L(I)
      IQ=J
      K=I+1
40    IF(K.GT.IQ)GOTO 90
      LK=L(K)
      IF(A(LK).LE.T)GOTO 80
      IQ=IQ
50    IF(IQ.LT.K)GOTO 70
      LQ=L(IQ)
      IF(A(LQ).GE.T)GOTO 60
      IX=L(K)
      L(K)=L(IQ)
      L(IQ)=IX
      IQ=IQ-1
      GOTO 80
60    IQ=IQ-1
      GOTO 50
70    CONTINUE
      IQ=K-1
      GOTO 100
80    K=K+1
      GOTO 40
90    CONTINUE
100   L(I)=L(IQ)
      L(IQ)=ITT
      IF(2*IQ-I-J)120,120,110
110   ILT(M)=I
      IUT(M)=IQ-1
      I=IQ+1
      GOTO 130
120   ILT(M)=IQ+1
      IUT(M)=J
      J=IQ-1
130   M=M+1
      GOTO 20
140   IF(I-J)150,170,170
150   LI=L(I)
      LJ=L(J)
      IF(A(LI)-A(LJ))170,170,160
160   IX=L(I)
      L(I)=L(J)
      L(J)=IX
170   M=M-1
      IF(M)190,190,180
180   I=ILT(M)
      J=IUT(M)
      GOTO 20
190   CONTINUE
      RETURN
      END
*
*
*
********************************************************* 
*
      SUBROUTINE SORT2(N,NA,L)
C     ********************************************************
C     THE SUBROUTINE SORT APPLIES TO THE ARRAY A WITH J ELEMENTS.
C     THE INDICES OF THE ORDERED ARRAY ARE PLACED IN ARRAY L
C     TO OBTAIN THE ORDERED ARRAY REPLACE THE INDICES I OF A(I)
C     WITH I=L(I)---REF.CACM 271
C     *********************************************************
C
      PARAMETER(N_HKL=10000)
      DIMENSION IUT(10),ILT(10)
      DIMENSION NA(N_HKL),L(N_HKL)
      J=N
      I=1
      M=1
      DO 10 K=1,J
10    L(K)=K
20    IF(J-I-1)140,140,30
30    IP=(J+I)/2
      ITT=L(IP)
      NT=NA(ITT)
      L(IP)=L(I)
      IQ=J
      K=I+1
40    IF(K.GT.IQ)GOTO 90
      LK=L(K)
      IF(NA(LK).LE.NT)GOTO 80
      IQ=IQ
50    IF(IQ.LT.K)GOTO 70
      LQ=L(IQ)
      IF(NA(LQ).GE.NT)GOTO 60
      IX=L(K)
      L(K)=L(IQ)
      L(IQ)=IX
      IQ=IQ-1
      GOTO 80
60    IQ=IQ-1
      GOTO 50
70    CONTINUE
      IQ=K-1
      GOTO 100
80    K=K+1
      GOTO 40
90    CONTINUE
100   L(I)=L(IQ)
      L(IQ)=ITT
      IF(2*IQ-I-J)120,120,110
110   ILT(M)=I
      IUT(M)=IQ-1
      I=IQ+1
      GOTO 130
120   ILT(M)=IQ+1
      IUT(M)=J
      J=IQ-1
130   M=M+1
      GOTO 20
140   IF(I-J)150,170,170
150   LI=L(I)
      LJ=L(J)
      IF(NA(LI)-NA(LJ))170,170,160
160   IX=L(I)
      L(I)=L(J)
      L(J)=IX
170   M=M-1
      IF(M)190,190,180
180   I=ILT(M)
      J=IUT(M)
      GOTO 20
190   CONTINUE
      RETURN
      END
*
*
*
********************************************************* 
*
      SUBROUTINE SORT3(N,NA,L)
C     ********************************************************
C     THE SUBROUTINE SORT APPLIES TO THE ARRAY A WITH J ELEMENTS.
C     THE INDICES OF THE ORDERED ARRAY ARE PLACED IN ARRAY L
C     TO OBTAIN THE ORDERED ARRAY REPLACE THE INDICES I OF A(I)
C     WITH I=L(I)---REF.CACM 271
C     *********************************************************
C
      PARAMETER(N_HKL=100000)
      DIMENSION IUT(10),ILT(10)
      DIMENSION NA(N_HKL),L(N_HKL)
      J=N
      I=1
      M=1
      DO 10 K=1,J
10    L(K)=K
20    IF(J-I-1)140,140,30
30    IP=(J+I)/2
      ITT=L(IP)
      NT=NA(ITT)
      L(IP)=L(I)
      IQ=J
      K=I+1
40    IF(K.GT.IQ)GOTO 90
      LK=L(K)
      IF(NA(LK).LE.NT)GOTO 80
      IQ=IQ
50    IF(IQ.LT.K)GOTO 70
      LQ=L(IQ)
      IF(NA(LQ).GE.NT)GOTO 60
      IX=L(K)
      L(K)=L(IQ)
      L(IQ)=IX
      IQ=IQ-1
      GOTO 80
60    IQ=IQ-1
      GOTO 50
70    CONTINUE
      IQ=K-1
      GOTO 100
80    K=K+1
      GOTO 40
90    CONTINUE
100   L(I)=L(IQ)
      L(IQ)=ITT
      IF(2*IQ-I-J)120,120,110
110   ILT(M)=I
      IUT(M)=IQ-1
      I=IQ+1
      GOTO 130
120   ILT(M)=IQ+1
      IUT(M)=J
      J=IQ-1
130   M=M+1
      GOTO 20
140   IF(I-J)150,170,170
150   LI=L(I)
      LJ=L(J)
      IF(NA(LI)-NA(LJ))170,170,160
160   IX=L(I)
      L(I)=L(J)
      L(J)=IX
170   M=M-1
      IF(M)190,190,180
180   I=ILT(M)
      J=IUT(M)
      GOTO 20
190   CONTINUE
      RETURN
      END
*
*
      SUBROUTINE CELREF(INDI,BBB,AFIN,NHKL,THETA,JHKL,DDT,DDQ)
*
C.....******************************************************************
C.....
C.....     PROGRAMME *** CELREF ***                         7/10/78
C.....
C.....     AUTEURS : JEAN LAUGIER & ALAIN FILHOL    20/10/78
C.....
C.....     AFFINEMENT LES PARAMETRES DE MAILLE
C.....                DU DECALAGE DE ZERO
C.....                DE LA LONGUEUR D'ONDE
C.....     A PARTIR DES ANGLES THETA DE BRAGG OBSERVES
C.....
C.....     METHODE : MOINDRES CARRES NON LINEAIRES
C.....
C.....     DONNEES :
C.....        1- CARTE COMMENTAIRE          FORMAT(16A5)
C.....        2- INDIC,IFIN                 FORMAT(2I)
C.....           INDIC  : CONTRAINTE D'AFFINEMENT
C.....                    0/1/2 POUR (A,B,C INDEPENDANTS)/(A=B=C)/(A=B)
C.....           IFIN   : NOMBRE MAXIMUM DE CYCLES D'AFFINEMENT
C.....
C.....        3- B(2),B(1)                       FORMAT(2F)
C.....           B(2)   : LONGUEUR D'ONDE
C.....           B(1)   : DECALAGE SYSTEMATIQUE DE ZERO
C.....
C.....        4- AFI(2),AFI(1)                   FORMAT(2F)
C.....           AFI(2) : 0/1 AFFINER (NON)/(OUI) LA LONGUEUR D'ONDE
C.....           AFI(1) :  "    "         "       LE DECALAGE DE ZERO
C.....
C.....        5- B(3 A 8)                        FORMAT(6F)
C.....           B(3)   : PARAMETRE A
C.....           B(4)   : PARAMETRE B
C.....           B(5)   : PARAMETRE C
C.....           B(6)   : ANGLE ALPHA
C.....           B(7)   : ANGLE BETA
C.....           B(8)   : ANGLE GAMMA
C.....
C.....        6- AFI(3 A 8)                      FORMAT(6F)
C.....           AFI(3 A 8) : 0/1 POUR AFFINER (NON)/(OUI)
C.....                        LES PARAMETRES DE MAILLE
C.....
C.....        7 A 7+NR- H,K,L,THETA              FORMAT(3I,F)
C.....          ("NR" NOMBRE D'OBSERVATIONS)
C.....          (DERNIERE CARTE : H,K,L=0 0 0)
C.....           H,K,L  : INDICES DE MILLER
C.....           THETA  : ANGLE DE BRAGG OBSERVE
C.....
C.....*****************************************************************
C
C
      INTEGER H(200),K(200),L(200)
      REAL*8 ICLE(3)
      DIMENSION THETA(10000),SIG(8),DUM(3),JHKL(3,10000)
      DIMENSION BBB(8),AFIN(8)
      EXTERNAL CALC
      COMMON/TROC/IWR,IRID
      COMMON/TRUC/QQ(200,10),BB(10),B(10),H,K,L,NPAF,AFI(10),NR,
     1            INDIC,PDS(200),NPAF2
C$OMP THREADPRIVATE(/TROC/,/TRUC/)
      DATA IWR/20/,ICLE/'A=B=C','A=B  ','NO   '/
      DATA SIG/8*0./
C
C----- A MODIFIER EN CAS DE CHANGEMENT DES DIMENSIONS
      NDMAX=200
C-----
      RD=180./ACOS(-1.)
C.....
C.....ENTREE DES DONNEES
   10 FORMAT(' PROGRAM *** CELREF ***  (J.LAUGIER & A.FILHOL 10/78)'/)
   20 FORMAT(20A4)
   30 FORMAT(10I)
   40 FORMAT(6F)
   50 FORMAT(3I,F)
   60 FORMAT(1H ,/
     1      ' OBSERVABLE NUMBER    : ',I5/
     1      ' ITERATION NUMBER : ',I5/
     1      ' REFINEMENT CONSTRAINTS : ',A5)
   70 FORMAT(' NUMBER OF INDEPENDENT PARAMETERS : ',I5/)
   80 FORMAT(' INITIAL VALUES :')
   90 FORMAT(' FINAL VALUES   : (STANDARD DEVIATIONS : 2nd LINE)')
  100 FORMAT(4X,'ZERO',4X,'LAMBDA',6X,'A',8X,'B',8X,'C',6X,
     1      'ALPHA',5X,'BETA',4X,'GAMMA')
  110 FORMAT(2X,F6.3,3X,F7.4,3(2X,F7.4),3(2X,F7.3))
  120 FORMAT(' RECIPROCAL CELL : ',3(2X,F7.5),3(2X,F7.3)/
     1       ' VOLUME (A**3)  : ',F12.3/)
  130 FORMAT(6X,F2.0,7(7X,F2.0))
  140 FORMAT(1H ,3X,1HH,5X,1HK,5X,1HL,2X,'TH(OBS)',4X,'TH-ZERO',4X,
     1      'TH(CALC)',5X,'DIFF.'/)
  150 FORMAT(2X,3(I3,3X),4(F7.3,4X))
  160 FORMAT(' SQRT(SUM(TH O-C)**2)/(NREF-NPAR))*1000 = ',F10.4/
     1 ' FACTEUR R : ',F10.4)
  170 FORMAT(' ##### ERROR REFLEXION : ',3I4,F8.3)
  180 FORMAT(' ##### DATA NUMBER GREATER THAN ',I4,' #####')
  190 FORMAT(' ##### IMPOSSIBLE TO REFINE ALL PARAMETERS TOGETHER'
     1 '##### THINK, PLEASE ! #####')
C
C.....ENTREE DES DONNEES
C      call open_read1(15,tmp)
C11    call open_write1(16,tmp)
11      WRITE(IWR,10)
C      READ(IRID,20)ITITR
C      READ(IRID,*)INDIC,IFIN
      INDIC=INDI
      DDT=0.
      DDQ=0.
      IFIN=10
      IF(INDIC.EQ.0.OR.INDIC.GT.3) INDIC=3
C      READ(IRID,*)B(2),B(1)
      B(1)=0.
C      READ(IRID,*)AFI(2),AFI(1)
      DO 5500 I=1,8
      AFI(I)=AFIN(I)
5500  B(I)=BBB(I)

C      READ(IRID,*)(B(I),I=3,8),(AFI(I),I=3,8)
      IFFI=IFIX(AFI(3)+AFI(4)+AFI(5)+0.1)
      IF(IFFI.EQ.0.OR.INDIC.EQ.3)GOTO 230
      IK=3-INDIC
      DO 220 I=1,IK
      IF(INDIC-2)200,210,210
  200 AFI(2+INDIC+I)=0.
      GOTO 220
  210 AFI(2+INDIC-1+I)=0.
  220 CONTINUE
      AFI(3)=1.
C
  230 NR=0
      DO 240 NR=1,NHKL
      IF(NR.GT.NDMAX)GOTO 380
      H(NR)=JHKL(1,NR)
      K(NR)=JHKL(2,NR)
      L(NR)=JHKL(3,NR)
      IHKL=IABS(H(NR))+IABS(K(NR))+IABS(L(NR))
      IF(IHKL.EQ.0)GOTO 260
      IF(THETA(NR).LE.0.)GOTO 370
  250 PDS(NR)=1.
  240 THETA(NR)=THETA(NR)/RD/2.   !!!!! 2*THETA EN THETA
C.....
  260 NR=NHKL
      WRITE(IWR,60)NR,IFIN,ICLE(INDIC)
      WRITE(IWR,80)
      WRITE(IWR,100)
      WRITE(IWR,130)(AFI(I),I=1,8)
      WRITE(IWR,110)(B(I),I=1,8)
      B(1)=B(1)/RD
      B(2)=B(2)*0.5
      DO 270 I=6,8
  270 B(I)=B(I)/RD
C   ...... CALCUL DES PARAMETRES MAILLE RECIPROQUE
      CALL INVER(B,DUM,VOLUM,0)
C   ......
      DO 280 I=1,3
  280 DUM(I)=B(5+I)*RD
      WRITE(IWR,120)(B(I),I=3,5),DUM,VOLUM
C
C....."NPAF" : NOMBRE DE PARAMETRES A AFFINER
C....."BB()" : TABLEAU DES PARAMETRES A AFFINER
      J=0
      DO 290 I=1,8
      IF(AFI(I).EQ.0.)GOTO 290
      J=J+1
      BB(J)=B(I)
  290 CONTINUE
      NPAF=J
      WRITE(IWR,70)NPAF
      IF(NPAF.EQ.8)GOTO 390
      NPAF2=NPAF+2
C   ......AFFINEMENT   (PDS() : POIDS (NON-UTILISE POUR CETTE VERSION))
      CALL MCRNL(QQ,NDMAX,THETA,BB,NPAF,NR,PDS,IFIN,CALC)
C   ......
C
C.....NOUVELLES VALEURS DES PARAMETRES
      J=0
      DO 300 I=1,8
      IF(AFI(I).EQ.0.)GOTO 300
      J=J+1
      B(I)=BB(J)
  300 CONTINUE
C   ......VALEURS DES ANGLES THETA CALCULES
      CALL FONC(THETA,R,RR)
C   ......
C
C.....CALCUL DES ECARTS TYPE
C....."SIG()" LES ECARTS TYPE DES PARAMETRES DE MAILLE QU'IL
C.....        CONTIENT SONT CEUX DES PARAMETRES RECIPROQUES.
      JJ=0
      DO 330 I=1,8
      IF(AFI(I))320,310,320
  310 SIG(I)=0.
      GOTO 330
  320 JJ=JJ+1
      SIG(I)=SQRT(QQ(JJ,JJ)*R)
  330 CONTINUE
      IF(INDIC.EQ.1.OR.INDIC.EQ.2) SIG(4)=SIG(3)
      IF(INDIC.EQ.1) SIG(5)=SIG(3)
      DO 340 I=1,3
      BB(I)=B(I+2)
  340 BB(I+3)=B(I+5)*RD
C
C   ......RETOUR A LA MAILLE DIRECTE (ET ECARTS TYPE CORRESPONDANTS)
      CALL INVER(B,SIG,VOLUM,1)
C   ......
      SIG(1)=SIG(1)*RD
      SIG(2)=SIG(2)*2.
C
C.....SORTIE DES RESULTATS
      VOLUM=1./VOLUM
C
C  Zeropoint in 2-theta to be added (same sense as TREOR, ITO, etc)
C
      B(1)=-B(1)*RD*2.
      SIG(1)=SIG(1)*2.
C
      B(2)=B(2)*2.
      DO 350 I=6,8
      SIG(I)=SIG(I)*RD
  350 B(I)=B(I)*RD
      WRITE(IWR,90)
      WRITE(IWR,*)
      WRITE(IWR,100)
      WRITE(IWR,110)(B(I),I=1,8)
      WRITE(IWR,110)(SIG(I),I=1,8)
      WRITE(IWR,120)(BB(I),I=1,6),VOLUM
      WRITE(IWR,140)
      WRITE(6,*)
      WRITE(6,90)
      WRITE(6,100)
      WRITE(6,110)(B(I),I=1,8)
      WRITE(6,110)(SIG(I),I=1,8)
      WRITE(6,120)(BB(I),I=1,6),VOLUM
      DO 351 I=1,8
351   BBB(I)=B(I)
      NPOUR=0
      DO 360 I=1,NR
      Y1=THETA(I)*RD
      Y2=Y1+B(1)/2.
      Y3=QQ(I,NPAF2)*RD + B(1)/2.
      Y4=Y2-Y3
      CONTINUE
      Y1=Y1*2.
      Y2=Y2*2.
      Y3=Y3*2.
      Y4=Y4*2.
      IF(I.LE.20)THEN
      DDT=DDT+ABS(Y4)
      QO=(2.*SIN(Y2/2.*3.141593/180.)/BBB(2))**2
      QC=(2.*SIN(Y3/2.*3.141593/180.)/BBB(2))**2
      DDQ=DDQ+ABS(QO-QC)
      ENDIF
360   WRITE(IWR,150)H(I),K(I),L(I),Y1,Y2,Y3,Y4
      WRITE(IWR,*)
      DDT=DDT/20.
      DDQ=DDQ/20.
      R=SQRT(R)*1000.
C      WRITE(IWR,160)R,RR
      IF(NPOUR.GT.0)GOTO 400
C      print 365
365   FORMAT('  Do you want a JCPDS-type output with hkl, dobs and
     2 dcalc ?'/
     3'  YES=0, NO=1 '$)
C      READ(5,*)NREP
      NREP=1
      IF(NREP.EQ.1)GOTO 400
      WRITE(IWR,366)
  366 FORMAT(1H ,3X,1HH,5X,1HK,5X,1HL,5X,'D(OBS)',4X,'D(CALC)'/)
      DO 367 I=1,NR
      Y0=B(2)/2.
      Y1=Y0/SIN(THETA(I)-B(1)/RD)
      Y2=Y0/SIN(QQ(I,NPAF2)-B(1)/RD)
  367 WRITE(IWR,368)H(I),K(I),L(I),Y1,Y2
  368 FORMAT(2X,3(I3,3X),2(F7.4,5X))
C      print 369
  369 FORMAT(/' The JCPDS output is at the end of the .imp file')
      GOTO 400
C
C.....MESSAGES D'ERREUR
  370 WRITE(IWR,170)H(NR),K(NR),L(NR),THETA(NR)
      NR=NR-1
C      GOTO 240
  380 WRITE(IWR,180)NDMAX
      GOTO 400
  390 WRITE(IWR,190)
  400 CONTINUE
      RETURN
      END
C
C
      SUBROUTINE CELREF2(INDI,BBB,AFIN,NHKL,THETA,JHKL,DDT,DDQ)
*
      INTEGER H(200),K(200),L(200)
      DIMENSION THETA(10000),SIG(8),DUM(3),JHKL(3,10000)
      DIMENSION BBB(8),AFIN(8)
      EXTERNAL CALC
      COMMON/TROC/IWR,IRID
      COMMON/TRUC/QQ(200,10),BB(10),B(10),H,K,L,NPAF,AFI(10),NR,
     1            INDIC,PDS(200),NPAF2
C$OMP THREADPRIVATE(/TROC/,/TRUC/)
      DATA IWR/20/
      DATA SIG/8*0./
      NDMAX=200
      RD=180./ACOS(-1.)
      PIP=3.141593/360.
      INDIC=INDI
      DDT=0.
      DDQ=0.
      IFIN=10
      IF(INDIC.EQ.0.OR.INDIC.GT.3) INDIC=3
      B(1)=0.
      DO 5500 I=1,8
      AFI(I)=AFIN(I)
5500  B(I)=BBB(I)
      IFFI=IFIX(AFI(3)+AFI(4)+AFI(5)+0.1)
      IF(IFFI.EQ.0.OR.INDIC.EQ.3)GOTO 230
      IK=3-INDIC
      DO 220 I=1,IK
      IF(INDIC-2)200,210,210
  200 AFI(2+INDIC+I)=0.
      GOTO 220
  210 AFI(2+INDIC-1+I)=0.
  220 CONTINUE
      AFI(3)=1.
C
  230 NR=0
      DO 240 NR=1,NHKL
      IF(NR.GT.NDMAX)GOTO 400
      H(NR)=JHKL(1,NR)
      K(NR)=JHKL(2,NR)
      L(NR)=JHKL(3,NR)
      IHKL=IABS(H(NR))+IABS(K(NR))+IABS(L(NR))
      IF(IHKL.EQ.0)GOTO 260
      IF(THETA(NR).LE.0.)GOTO 400
  250 PDS(NR)=1.
  240 THETA(NR)=THETA(NR)/RD/2.   
C.....
  260 NR=NHKL
      B(1)=B(1)/RD
      B(2)=B(2)*0.5
      DO 270 I=6,8
  270 B(I)=B(I)/RD
      CALL INVER(B,DUM,VOLUM,0)
      DO 280 I=1,3
  280 DUM(I)=B(5+I)*RD
      J=0
      DO 290 I=1,8
      IF(AFI(I).EQ.0.)GOTO 290
      J=J+1
      BB(J)=B(I)
  290 CONTINUE
      NPAF=J
      IF(NPAF.EQ.8)GOTO 400
      NPAF2=NPAF+2
      CALL MCRNL(QQ,NDMAX,THETA,BB,NPAF,NR,PDS,IFIN,CALC)
      J=0
      DO 300 I=1,8
      IF(AFI(I).EQ.0.)GOTO 300
      J=J+1
      B(I)=BB(J)
  300 CONTINUE
      CALL FONC(THETA,R,RR)
      JJ=0
      DO 330 I=1,8
      IF(AFI(I))320,310,320
  310 SIG(I)=0.
      GOTO 330
  320 JJ=JJ+1
      SIG(I)=SQRT(QQ(JJ,JJ)*R)
  330 CONTINUE
      IF(INDIC.EQ.1.OR.INDIC.EQ.2) SIG(4)=SIG(3)
      IF(INDIC.EQ.1) SIG(5)=SIG(3)
      DO 340 I=1,3
      BB(I)=B(I+2)
  340 BB(I+3)=B(I+5)*RD
      SIG(1)=SIG(1)*RD
      SIG(2)=SIG(2)*2.
      B(1)=-B(1)*RD*2.
      SIG(1)=SIG(1)*2.
      B(2)=B(2)*2.
      DO 350 I=6,8
      SIG(I)=SIG(I)*RD
  350 B(I)=B(I)*RD
      DO 351 I=1,8
351   BBB(I)=B(I)
      NPOUR=0
      DO 360 I=1,NR
      Y1=THETA(I)*RD
      Y2=Y1+B(1)/2.
      Y3=QQ(I,NPAF2)*RD + B(1)/2.
      Y4=Y2-Y3
      CONTINUE
      Y1=Y1*2.
      Y2=Y2*2.
      Y3=Y3*2.
      Y4=Y4*2.
      IF(I.LE.20)THEN
      DDT=DDT+ABS(Y4)
      QO=(2.*SIN(Y2*PIP)/BBB(2))**2
      QC=(2.*SIN(Y3*PIP)/BBB(2))**2
      DDQ=DDQ+ABS(QO-QC)
      ENDIF
360   CONTINUE
      DDT=DDT/20.
      DDQ=DDQ/20.
400   CONTINUE
      RETURN
      END
C
C
      SUBROUTINE CALC
C.....--------------
      INTEGER H(200),K(200),L(200)
      DIMENSION Q(200,10)
      COMMON/TRUC/QQ(200,10),BB(10),B(10),H,K,L,NPAF,AFI(10),NR,
     *            INDIC,PDS(200),NPAF2
      COMMON/TROC/IWR,IRID
C$OMP THREADPRIVATE(/TROC/,/TRUC/)
      J=0
      DO 1 I=1,8
      IF(AFI(I).EQ.0.)GO TO 1
      J=J+1
      B(I)=BB(J)
1     CONTINUE
      IF(INDIC.EQ.1.OR.INDIC.EQ.2)B(4)=B(3)
      IF(INDIC.EQ.1)B(5)=B(3)
      AE=B(3)
      BE=B(4)
      CE=B(5)
      CAE=COS(B(6))
      CBE=COS(B(7))
      CCE=COS(B(8))
      DO 2 I=1,NR
      DD=(AE*H(I))**2+(BE*K(I))**2+(CE*L(I))**2
     *+2.*(H(I)*K(I)*AE*BE*CCE+K(I)*L(I)*BE*CE*CAE
     *+L(I)*H(I)*CE*AE*CBE)
      D=1./SQRT(DD)
      RAD=SQRT(1.-(B(2)**2)*DD)
      F=B(2)*D/RAD
      Q(I,1)=1.
      Q(I,2)=1./(D*RAD)
      Q(I,3)=F*H(I)*(H(I)*AE+K(I)*BE*CCE+L(I)*CE*CBE)
      Q(I,4)=F*K(I)*(K(I)*BE+L(I)*CE*CAE+H(I)*AE*CCE)
      Q(I,5)=F*L(I)*(L(I)*CE+H(I)*AE*CBE+K(I)*BE*CAE)
      Q(I,6)=-F*K(I)*L(I)*BE*CE*SIN(B(6))
      Q(I,7)=-F*L(I)*H(I)*CE*AE*SIN(B(7))
      Q(I,8)=-F*H(I)*K(I)*AE*BE*SIN(B(8))
2     QQ(I,NPAF2)=B(1)+ASIN(B(2)/D)
      DO 3 IR=1,NR
      J=0
      DO 3 I=1,8
      IF(AFI(I).EQ.0.)GO TO 3
      J=J+1
      QQ(IR,J)=Q(IR,I)
3     CONTINUE
C      WRITE(IWR,5)(BB(I),I=1,NPAF)
5     FORMAT(2X,8(E13.7,2X),/)
      RETURN
      END
C
C
      SUBROUTINE FONC(THETA,R,RR)
C.....-----------------------
      DIMENSION THETA(200)
      INTEGER H(200),K(200),L(200)
      COMMON/TRUC/QQ(200,10),BB(10),B(10),H,K,L,NPAF,AFI(10),NR,
     *            INDIC,PDS(200),NPAF2
C$OMP THREADPRIVATE(/TRUC/)
      R1=0.
      R2=0.
      AE=B(3)
      BE=B(4)
      CE=B(5)
      CAE=COS(B(6))
      CBE=COS(B(7))
      CCE=COS(B(8))
      DO 2 I=1,NR
      DD=(AE*H(I))**2+(BE*K(I))**2+(CE*L(I))**2
     *+2.*(H(I)*K(I)*AE*BE*CCE+K(I)*L(I)*BE*CE*CAE
     *+L(I)*H(I)*CE*AE*CBE)
      D=1./SQRT(DD)
      YC=B(1)+ASIN(B(2)/D)
      R1=R1+(YC-THETA(I))**2
      R2=R2+THETA(I)**2
2     QQ(I,NPAF2)=YC
      R=R1/(NR-NPAF)
      RR=SQRT(R1/R2)
      RETURN
      END
C
C
      SUBROUTINE MCRNL(Q,ID,Y,B,M,N,P,IFIN,CALC)
C.....------------------------------------------
C
C.....PROGRAMME DE MOINDRES CARRES NON LINEAIRE
C.....UTILISANT L'INVERSION DE MATRICE TABLEAU A UNE DIMENSION
C.....LIMITE A 15 PARAMETRES
C..... M = NBRE DE PARAMETRES
C..... N = NBRE DE DONNEES
C..... ID = DIMENSION DU TABLEAU Q(ID,M+1) 
C..... A((M*(M-1)/2 + M) = ZONE DE TRAVAIL
C
      DIMENSION Q(ID,1),Y(1),B(1),P(1),A(60)
      MM=2*M
      M1=M+1
      M2=M+2
   10 IFIN=IFIN-1
C
C-----CALCUL DES DERIVEES PARTIELLES ET DES Y
      CALL CALC
C
C-----LES Y CALCULES SONT DANS LA COLONNE M+2
      DO 30 J=1,M
      R=0.
      DO 20 I=1,N
   20 R=R+(Y(I)-Q(I,M2))*Q(I,J)*P(I)
   30 Q(J,M1)=R
C
C-----CONSTRUCTION DE LA MATRICE SYMETRIQUE A=Q*QT
      NO=0
      DO 50 IQ=1,M
      DO 50 IL=IQ,M
      NO=NO+1
      R=0.
      DO 40 I=1,N
   40 R=R+Q(I,IQ)*Q(I,IL)*P(I)
   50 A(NO)=R
C
C-----INVERSION DE LA MATRICE A
      CALL MATINV(A,M,IER)
      IF(IFIN)110,110,60
C
C-----CALCUL DES NOUVEAUX PARAMETRES
   60 DO 100 I=1,M
      R=0.
      IMM=(I-1)*(MM-I)/2
      DO 90 J=1,M
      IF(J-I)80,70,70
   70 L=IMM+J
      GOTO 90
   80 L=(J-1)*(MM-J)/2+I
   90 R=R+A(L)*Q(J,M1)
  100 B(I)=B(I)+R
      GOTO 10
C
C-----REMISE DES ELEMENTS DIAGONAUX DANS Q(I,I),APRES LE DERNIER CYCLE
  110 DO 120 I=1,M
      L=(I-1)*(MM-I)/2+I
  120 Q(I,I)=A(L)
      RETURN
      END
C
C
      SUBROUTINE MATINV(AM,N,NFAIL)
C.....-----------------------------
      DIMENSION AM(1)
      DOUBLE PRECISION SUMA
C     ********** SEGMENT 1 OF CHOLESKI INVERSION **********
C     ***** FACTOR MATRIX INTO LOWER TRIANGLE X TRANSPOSE *****
      K=1
      IF(N-1)8,10,20
    8 NFAIL= K
      GOTO 210
   10 AM(1)=1.0/AM(1)
      GOTO 200
C     ***** LOOP M OF A(L,M) *****
   20 DO 110 M=1,N
      IMAX=M-1
C     ***** LOOP L OF A(L,M) *****
      DO 100 L=M,N
      SUMA=0.0
      KLI=L
      KMI=M
      IF(IMAX)50,50,30
C     *****SUM OVER I=1,M-1 A(L,I)*A(M,I) *****
   30 DO 40 I=1,IMAX
      SUMA=SUMA+AM(KLI)*AM(KMI)
      J=N-I
      KLI=KLI+J
   40 KMI=KMI+J
C     *****TERM=C(L,M)-SUM *****
   50 TERM=AM(K)-SUMA
      IF(L-M)60,60,90
   60 IF(TERM)80,80,70
C     ***** A(M,M)=SQRT(TERM) *****
   70 DENOM=SQRT(TERM)
      AM(K)=DENOM
      GOTO 100
   80 NFAIL= K
      GOTO 210
C     ***** A(L,M)=TERM/A(M,M) *****
   90 AM(K)=TERM/DENOM
  100 K=K+1
  110 CONTINUE
C     ********** SEGMENT 2 OF CHOLESKI INVERSION **********
C     *****INVERSION OF TRIANGULAR MATRIX*****
  120 AM(1)=1.0/AM(1)
      KDM=1
C     ***** STEP L OF B(L,M) *****
      DO 150 L=2,N
      KDM=KDM+N-L+2
C     ***** RECIPROCAL OF DIAGONAL TERM *****
      TERM = 1.0/AM(KDM)
      AM(KDM)=TERM
      KMI=0
      KLI=L
      IMAX=L-1
C     ***** STEP M OF B(L,M) *****
      DO 140 M=1,IMAX
      K=KLI
C     ***** SUM TERMS *****
      SUMA=0.0
      DO 130 I=M,IMAX
      II=KMI+I
      SUMA=SUMA-AM(KLI)*AM(II)
  130 KLI=KLI+N-I
C     ***** MULT SUM * RECIP OF DIAGONAL *****
      AM(K)=SUMA*TERM
      J=N-M
      KLI=K+J
  140 KMI=KMI+J
  150 CONTINUE
C     ********** SEGMENT 3 OF CHOLESKI INVERSION **********
C     *****PREMULTIPLY LOWER TRIANGLE BY TRANSPOSE*****
  160 K=1
      DO 190 M=1,N
      KLI=K
      DO 180 L=M,N
      KMI=K
      IMAX=N-L+1
      SUMA=0.0
      DO 170 I=1,IMAX
      SUMA=SUMA+AM(KLI)*AM(KMI)
      KLI=KLI+1
  170 KMI=KMI+1
      AM(K)=SUMA
  180 K=K+1
  190 CONTINUE
  200 NFAIL=0
  210 RETURN
      END
C
C
      SUBROUTINE INVER(B,DB,VOLUM,IV)
C.....
C.....  1-CALCUL LES PARAMETRES MAILLE INVERSE
C.....  2-CALCULE LES ECARTS TYPE DES PARAMETRES
C.....   IV=0 OPTION 1  IV=1 OPTIONS 1 & 2
C.....
      DIMENSION B(8),DB(8)
      DIMENSION AD(6),D(6),SINP(3),COSP(3),SP(3),SS(3),
     1          CC(3),DQD(3),SIG(6)
C
      CABC2=0.
      DO 10 I=1,3
      AD(I)=B(I+2)
      COSP(I)=COS(B(I+5))
      SINP(I)=SIN(B(I+5))
   10 CONTINUE
      DO 15 I=1,3
      J=MOD(I,3)+1
      K=MOD(I+1,3)+1
      DQD(I)=COSP(I)-COSP(J)*COSP(K)
      SS(I)=SINP(J)*SINP(K)
      CC(I)=-DQD(I)/SS(I)
      CABC2=CABC2+COSP(I)*COSP(I)
   15 CONTINUE
C
      Q2=1.-CABC2+2.*COSP(1)*COSP(2)*COSP(3)
      Q=SQRT(Q2)
      VOLUM=AD(1)*AD(2)*AD(3)*Q
C
      DO 20 I=1,3
      B(I+2)=SINP(I)/(AD(I)*Q)
      B(I+5)=ACOS(CC(I))
      SP(I)=SIN(B(I+5))
   20 CONTINUE
C
      IF(IV.EQ.0)GOTO 70
C
C.....DERIVEES DES PARAMETRES A , B , C
      DO 30 I=1,3
      J=MOD(I,3)+1
      K=MOD(I+1,3)+1
      D(I)=-SINP(I)/AD(I)
      D(J)=0.
      D(K)=0.
      D(I+3)=COSP(I)-SINP(I)*SINP(I)*DQD(I)/(Q*Q)
      D(J+3)=-SS(K)*DQD(J)/Q2
      D(K+3)=-SS(J)*DQD(K)/Q2
      CALL SIGMA(D,DB,R)
   30 SIG(I)=SQRT(R)/(Q*AD(I))
C
C.....DERIVEES DES ANGLES DE LA MAILLE
      DO 50 I=1,3
      DO 40 JJ=1,3
   40 D(JJ)=0.
      J=MOD(I,3)+1
      K=MOD(I+1,3)+1
      D(I+3)=SINP(I)/SS(I)
      D(J+3)=COSP(K)/SINP(K) + COSP(J)*CC(I)/SINP(J)
      D(K+3)=COSP(J)/SINP(J) + COSP(K)*CC(I)/SINP(K)
      CALL SIGMA(D,DB,R)
   50 SIG(I+3)=SQRT(R)/SP(I)
C
      DO 60 I=1,6
   60 DB(I+2)=SIG(I)
C
   70 RETURN
      END
C
C
      SUBROUTINE SIGMA(D,DB,R)
C.....
      DIMENSION DB(8),D(6)
      R=0.
      DO 10 I=1,6
   10 R=R + (D(I)*DB(I+2))**2
      RETURN
      END
C
C... Here is the program heart...
C
      SUBROUTINE CALCUL1(DIFF,DIFF2)
C
      PARAMETER(N_HKL=10000)   
C
      DIMENSION THETA(N_HKL),PERC(N_HKL),CR(N_HKL),DE(N_HKL)
C
      COMMON/CAL/NHKL0,LHKL,NDAT,DMIN,SLABDA2,IHH(3,N_HKL),AL(3,3),PI,
     1CRI(N_HKL),DIFP(N_HKL),DIFM(N_HKL),TH2(N_HKL),FOBS(N_HKL),SUM_F,
     2NIND,W2(N_HKL),NHKL,NDAT10
C
C$OMP THREADPRIVATE(/cal/) 
C
C...  Keep only the hkl for d > dmin
C
      JH=0
      DO 109 I=1,NHKL0
      X=0.
      DO 107 J=1,3
      DO 107 K=J,3
107   X=AL(J,K)*IHH(J,I)*IHH(K,I)+X
      IF(X.GT.DMIN)GO TO 109
      JH=JH+1
C     X IS 1/D(hkl)**2 FOR REFLECTION IHH
C
C     This should be optimized for speed :
C     working only on X, not calculating 2-theta...
C     - in fact, tests show that only 10-15% is gained -
C
      SINTH=SLABDA2*X
      SINTH=SQRT(SINTH)
      THETA(JH)=ASIN(SINTH)*PI
      CR(JH)=0.
109   CONTINUE
      NHKL=JH
      IF(NHKL.GT.NDAT10)RETURN
C
C...  Comparison with the data
C
      LHKL=0
      DO 113 J=1,NDAT
      CRI(J)=0.
      PERC(J)=0.
      DEMAX=W2(J)
      DO 111 K=1,NHKL
      IF(CR(K).EQ.1.)GO TO 111
      IF(THETA(K).LE.DIFP(J).AND.THETA(K).GE.DIFM(J))THEN
      DE(K)=ABS(THETA(K)-TH2(J))
      IF(DE(K).LE.DEMAX)THEN
      L=K
      DEMAX=DE(K)
      CRI(J)=1.
      ENDIF
      ENDIF
111   CONTINUE
C
C  PERC = percentage of columnar overlap for that peak
C
C  Potential problem here because only one reflection
C  overlapping the most closely with the column is 
C  included (if CRI =1)...
C
      IF(CRI(J).EQ.1) THEN
      PERC(J)=1.-DE(L)/W2(J)
      LHKL=LHKL+1
      CR(L)=1.
      ENDIF
113   CONTINUE
C
C...  Calculate "R"
C
      DIFF1=0.
      SUM_F2=0.
      DO 1122 K=1,NDAT
      SUM_F2=SUM_F2+FOBS(K)*CRI(K)
1122  DIFF1=DIFF1+CRI(K)*FOBS(K)*PERC(K)
      DIFF=1.-DIFF1/SUM_F
      DIFF2=1.-DIFF1/SUM_F2
      RETURN
      END
C
      SUBROUTINE CALCUL2(DIFF,IHKL,TH3,NCALC,IGC)
C
      PARAMETER(N_HKL=10000,N_DAT=100)   
C
      DIMENSION THETA(N_HKL),PERC(N_HKL),IHKL(3,N_HKL),DE(N_HKL)
      DIMENSION TH3(N_HKL),JHKL(3,N_HKL),CR(N_HKL),QC(N_HKL)
C
      COMMON/CAL/NHKL0,LHKL,NDAT,DMIN,SLABDA2,IHH(3,N_HKL),AL(3,3),PI,
     1CRI(N_HKL),DIFP(N_HKL),DIFM(N_HKL),TH2(N_HKL),FOBS(N_HKL),SUM_F,
     2NIND,W2(N_HKL),NHKL,NDAT10
      COMMON/CAL2/IND(N_DAT,N_HKL)
C
C$OMP THREADPRIVATE(/cal/,/cal2/) 
C
C...  Keep only the hkl for d > dmin
C
      JH=0
      DO 109 I=1,NHKL0
      X=0.
      DO 107 J=1,3
      DO 107 K=J,3
107   X=AL(J,K)*IHH(J,I)*IHH(K,I)+X
      IF(X.GT.DMIN)GO TO 109
      JH=JH+1
      DO 108 J=1,3
108   JHKL(J,JH)=IHH(J,I)
C     X IS 1/D(hkl)**2 FOR REFLECTION IHH
C
C     This should be optimized for speed :
C     working only on X, not calculating 2-theta...
C
      QC(JH)=X
      SINTH=SLABDA2*X
      SINTH=SQRT(SINTH)
      THETA(JH)=ASIN(SINTH)*PI
      CR(JH)=0.
109   CONTINUE
      NHKL=JH
C
C...  Comparison with the data
C
      LHKL=0
      NCALC=0
      JJ=0
      DO 113 J=1,NDAT
      CRI(J)=0.
      PERC(J)=0.
CCC
CCC  Eliminating too spurious peaks here ???
CCC    tolerance on width decreased by a factor 3
CCC
      DEMAX=W2(J)*0.33333
CCC
      DO 111 K=1,NHKL
      IF(CR(K).EQ.1.)GO TO 111
      IF(THETA(K).LE.DIFP(J).AND.THETA(K).GE.DIFM(J))THEN
      DE(K)=ABS(THETA(K)-TH2(J))
      IF(DE(K).LE.DEMAX)THEN
      L=K
      DEMAX=DE(K)
      CRI(J)=1.
      ENDIF
      ENDIF
111   CONTINUE
C
C  PERC = percentage of columnar overlap for that peak
C
C  Potential problem here because only one reflection
C  overlapping the most closely with the column is 
C  included (if CRI =1)...
C
      IF(CRI(J).EQ.1) THEN
      PERC(J)=1.-DE(L)/(W2(J)*0.33333)
      LHKL=LHKL+1
      CR(L)=1.
      DO 112 L2=1,3
112   IHKL(L2,LHKL)=JHKL(L2,L)
      TH3(LHKL)=TH2(J)
      JJ=JJ+CRI(J)
      IF(JJ.LE.20)JJJ=J
      ENDIF
113   CONTINUE
C  NCALC is for M(20) FoM
      DO 114 K=1,NHKL
      IF(TH2(JJJ).GE.THETA(K))NCALC=NCALC+1
114   CONTINUE
C
C...  Calculate "R"
C
      DIFF=0.
C
C  Change here with SUM_F2 being only on explained reflections...
C
      SUM_F2=0.
      DO 1122 K=1,NDAT
      IND(K,IGC)=CRI(K)
      SUM_F2=SUM_F2+FOBS(K)*CRI(K)
1122  DIFF=DIFF+CRI(K)*FOBS(K)*PERC(K)
      DIFF=1.-DIFF/SUM_F2
      RETURN
      END
C
C
C
      SUBROUTINE KILLK(PRESSEDK)
C
C    Checks if the 'K' keystroke has been pressed
C
!      USE DFLIB   ! disabled library by Zucco
      LOGICAL(4) PRESSED / .FALSE. /
      LOGICAL(4) PRESSEDK
      CHARACTER(1) KEY
      PRESSED = PEEKCHARQQ ( )
      IF(PRESSED)THEN
      KEY = GETCHARQQ()
      IF(KEY.EQ.'K')PRESSEDK=.TRUE.
      ENDIF
      RETURN
      END
      real function randi(ix)
      implicit integer*4 (a-z)
      real*8 x,s
*
*     A random number generator using the recursion IX=IX*A MOD P
*     where A=7**5 and P=2**31-1. The value returned is in the 
*     range 0.<= ran <1.
*
*     In this form fairly portable as does not require knowledge
*     of data storage, but does not adhere to the standard in two respects:
*     1) Assumes an integer word length of at least 32 bits 
*     2) Assumes that a positive integer less than 2**16 may be
*        floated without loss of digits.
*
*     This code is based on code published by Linus Schrage in
*     T.O.M.S. vol.5 no.2 June 1979 (pp 132-138)
*
*     The method employed is a multiplicative congruential one using a 
*     multiplier of 7**5 and taking the modulo to 2**31-1, i.e. the
*     generator number, x, is updated on each call to the value
*     x*7**5  modulo (2**31-1). The result returned is calculated as a 
*     real number having the value x/(2**31-1)  
*
      parameter (a=16807)
      parameter (b15=32768)
      parameter (b16=65536)
      parameter (p=2147483647)
      parameter (s=1d0/p)
*     a=7**5, b15=2**15, b16=2**16, p= 2**31-1
*
*     Get the 15 high order bits and 16 low order bits of ix
*
      xhi=ix/b16
      xlo=ix-xhi*b16
*
*     Multiply low order part
*
      xalo=xlo*a
*
*     Get high order part of product to carry
*
      leftlo=xalo/b16
*
*     And so obtain high order part of product
*
      xahi=xhi*a+leftlo
*
*     Obtain 32nd bit (overflow) of full product
*
      k=xahi/b15
*
*     Put all the bits of the product together and subtract p
*     (Must be in this form to prevent overflow. The ()'s are essential)
*
      ix=(((xalo-leftlo*b16 )-p) +(xahi-k*b15)*b16) +k
*
*     If <0 add p back again
*
      if (ix.lt.0) ix=ix+p
*
*     Finally multiply by 1/(2**31-1) to obtain number in range 0-1
*
      xhi=ix/b16
      x=dble(xhi)*65536.0d0+dble(ix-xhi*b16)
      x=x*s
      randi=sngl(x)
      return
      end
C
C... Check for supercell, and reduce to minimal cell
C
      SUBROUTINE SUPCEL(N,IHKL,CEL,L,VGC,JS)
      PARAMETER(N_HKL=10000)
      DIMENSION CEL(6,N_HKL),VGC(N_HKL),IHKL(3,N_HKL)
      DIMENSION ID2(3),ID3(3),ID4(3),ID5(3),ID6(3),IHMAX(3)
C
C  Sum on the h, k, l
      DO 10 J=1,3
C  At the end, if the following codes are still = 1
C     then there could be a common divider
      ID2(J)=1
      ID3(J)=1
      ID4(J)=1
      ID5(J)=1
      ID6(J)=1
      IHMAX(J)=0
C  Sum on the hkl data
      DO 9 I=1,N
C  Test on dividing h, k, l by 2, 3, 4, 5, 6
      K=IABS(IHKL(J,I))
      IF(K.GT.IHMAX(J))IHMAX(J)=K
      IF(K.EQ.0)GO TO 9
      IF(K.EQ.1.OR.K.EQ.7.OR.K.EQ.11.OR.K.EQ.13.OR.K.EQ.17.OR.K.EQ.19.OR
     1.K.EQ.21)THEN
      ID2(J)=0
      ID3(J)=0
      ID4(J)=0
      ID5(J)=0
      ID6(J)=0
      GO TO 10
      ENDIF
      IF(K.EQ.2)THEN
      ID3(J)=0
      ID4(J)=0
      ID5(J)=0
      ID6(J)=0
      ENDIF
      IF(K.EQ.3)THEN
      ID2(J)=0
      ID4(J)=0
      ID5(J)=0
      ID6(J)=0
      ENDIF
      IF(K.EQ.4)THEN
      ID3(J)=0
      ID5(J)=0
      ID6(J)=0
      ENDIF
      IF(K.EQ.5)THEN
      ID2(J)=0
      ID3(J)=0
      ID4(J)=0
      ID6(J)=0      
      ENDIF
      IF(K.EQ.6)THEN
      ID4(J)=0
      ID5(J)=0      
      ENDIF
      IF(K.EQ.8)THEN
      ID3(J)=0
      ID5(J)=0      
      ID6(J)=0      
      ENDIF
      IF(K.EQ.9)THEN
      ID2(J)=0
      ID4(J)=0
      ID5(J)=0
      ID6(J)=0      
      ENDIF
      IF(K.EQ.10)THEN
      ID3(J)=0
      ID4(J)=0
      ID6(J)=0
      ENDIF
      IF(K.EQ.12)ID5(J)=0
      IF(K.EQ.14)THEN
      ID3(J)=0
      ID4(J)=0
      ID5(J)=0
      ID6(J)=0      
      ENDIF
      IF(K.EQ.15)THEN
      ID2(J)=0
      ID4(J)=0
      ID6(J)=0
      ENDIF
      IF(K.EQ.16)THEN
      ID3(J)=0
      ID5(J)=0
      ID6(J)=0
      ENDIF
      IF(K.EQ.18)THEN
      ID4(J)=0
      ID5(J)=0
      ENDIF
      IF(K.EQ.20)THEN
      ID3(J)=0
      ID6(J)=0
      ENDIF
9     CONTINUE
10    CONTINUE
      IF(JS.EQ.3)THEN
      IF(ID2(1).NE.ID2(2))GO TO 20
      IF(ID2(1).NE.ID2(3))GO TO 20
      IF(ID3(1).NE.ID3(2))GO TO 20
      IF(ID3(1).NE.ID3(3))GO TO 20
      IF(ID4(1).NE.ID4(2))GO TO 20
      IF(ID4(1).NE.ID4(3))GO TO 20
      IF(ID5(1).NE.ID5(2))GO TO 20
      IF(ID5(1).NE.ID5(3))GO TO 20
      IF(ID6(1).NE.ID6(2))GO TO 20
      IF(ID6(1).NE.ID6(3))GO TO 20
      ENDIF
      IAB2=1
      IAB3=1
      IAB4=1
      IAB5=1
      IAB6=1
      IF(JS.EQ.2)THEN
      IF(ID2(1).NE.ID2(2))IAB2=0
      IF(ID3(1).NE.ID3(2))IAB3=0
      IF(ID4(1).NE.ID4(2))IAB4=0
      IF(ID5(1).NE.ID5(2))IAB5=0
      IF(ID6(1).NE.ID6(2))IAB6=0
      ENDIF
      IF(JS.EQ.4)THEN
      IF(ID2(1).NE.ID2(2))IAB2=0
      IF(ID3(1).NE.ID3(2))IAB3=0
      IF(ID4(1).NE.ID4(2))IAB4=0
      IF(ID5(1).NE.ID5(2))IAB5=0
      IF(ID6(1).NE.ID6(2))IAB6=0
      IHH=1
      DO 8 J=1,N
      K=IHKL(1,J)+IHKL(2,J)
      DO 8 I=1,21
      II=I*2-1 
      IF(K.EQ.II)IHH=0
8     CONTINUE
      IF(IHH.EQ.1)THEN
      XX=1.41421
      CEL(1,L)=CEL(1,L)/XX
      CEL(2,L)=CEL(2,L)/XX
      VGC(L)=VGC(L)/XX/XX
      ENDIF
      ENDIF
      DO 11 J=1,2
      IF(ID6(J).EQ.1)THEN
      ID3(J)=0
      ID2(J)=0
      ENDIF
      IF(ID4(J).EQ.1)ID2(J)=0
      IF(ID2(J).EQ.1.AND.IHMAX(J).GE.2.AND.IAB2.EQ.1)THEN
      XX=2.
      CEL(J,L)=CEL(J,L)/XX
      VGC(L)=VGC(L)/XX
      ENDIF
      IF(ID3(J).EQ.1.AND.IHMAX(J).GE.3.AND.IAB3.EQ.1)THEN
      XX=3.
      CEL(J,L)=CEL(J,L)/XX
      VGC(L)=VGC(L)/XX
      ENDIF
      IF(ID4(J).EQ.1.AND.IHMAX(J).GE.4.AND.IAB4.EQ.1)THEN
      XX=4.
      CEL(J,L)=CEL(J,L)/XX
      VGC(L)=VGC(L)/XX
      ENDIF
      IF(ID5(J).EQ.1.AND.IHMAX(J).GE.5.AND.IAB5.EQ.1)THEN
      XX=5.
      CEL(J,L)=CEL(J,L)/XX
      VGC(L)=VGC(L)/XX
      ENDIF
      IF(ID6(J).EQ.1.AND.IHMAX(J).GE.6.AND.IAB6.EQ.1)THEN
      XX=6.
      CEL(J,L)=CEL(J,L)/XX
      VGC(L)=VGC(L)/XX
      ENDIF
11    CONTINUE
      J=3
      IF(ID6(J).EQ.1)THEN
      ID3(J)=0
      ID2(J)=0
      ENDIF
      IF(ID4(J).EQ.1)ID2(J)=0
      IF(ID2(J).EQ.1.AND.IHMAX(J).GE.2)THEN
      XX=2.
      CEL(J,L)=CEL(J,L)/XX
      VGC(L)=VGC(L)/XX
      ENDIF
      IF(ID3(J).EQ.1.AND.IHMAX(J).GE.3)THEN
      XX=3.
      CEL(J,L)=CEL(J,L)/XX
      VGC(L)=VGC(L)/XX
      ENDIF
      IF(ID4(J).EQ.1.AND.IHMAX(J).GE.4)THEN
      XX=4.
      CEL(J,L)=CEL(J,L)/XX
      VGC(L)=VGC(L)/XX
      ENDIF
      IF(ID5(J).EQ.1.AND.IHMAX(J).GE.5)THEN
      XX=5.
      CEL(J,L)=CEL(J,L)/XX
      VGC(L)=VGC(L)/XX
      ENDIF
      IF(ID6(J).EQ.1.AND.IHMAX(J).GE.6)THEN
      XX=6.
      CEL(J,L)=CEL(J,L)/XX
      VGC(L)=VGC(L)/XX
      ENDIF
20    CONTINUE
      RETURN
      END
C
C  Test on Bravais lattice
C
      SUBROUTINE BRAV(N,IHKL,IBR)
      INTEGER HSUM,H,K,L
      PARAMETER(N_HKL=10000)
      DIMENSION IHKL(3,N_HKL)
      ibr=6
C
C     I lattice
C
      j=0
      do i=1,n
C      write(*,*)ihkl(1,i),ihkl(2,i),ihkl(3,i)
      hsum=ihkl(1,i)+ihkl(2,i)+ihkl(3,i)
      k=1
      if(2*(hsum/2).eq.hsum)k=0
      j=j+k
      enddo
      if(j.eq.0)then
      ibr=1
      return
      endif
C
C     F lattice
C
      j=0
      do i=1,n
      h=1
      if(2*(ihkl(1,i)/2).eq.ihkl(1,i))h=0
      k=1
      if(2*(ihkl(2,i)/2).eq.ihkl(2,i))k=0
      l=1
      if(2*(ihkl(3,i)/2).eq.ihkl(3,i))l=0
      jj=h+k+l
      if(jj.ne.0)then
      h=1
      if(2*(ihkl(1,i)/2).ne.ihkl(1,i))h=0
      k=1
      if(2*(ihkl(2,i)/2).ne.ihkl(2,i))k=0
      l=1
      if(2*(ihkl(3,i)/2).ne.ihkl(3,i))l=0
      jj=h+k+l
      endif
      j=j+jj
      enddo
      if(j.eq.0)then
      ibr=5
      return
      endif
C
C     A lattice
C
      j=0
      do i=1,n
      hsum=ihkl(2,i)+ihkl(3,i)
      k=1
      if(2*(hsum/2).eq.hsum)k=0
      j=j+k
      enddo
      if(j.eq.0)then
      ibr=2
      return
      endif
C
C     B lattice
C
      j=0
      do i=1,n
      hsum=ihkl(1,i)+ihkl(3,i)
      k=1
      if(2*(hsum/2).eq.hsum)k=0
      j=j+k
      enddo
      if(j.eq.0)then
      ibr=3
      return
      endif
C
C     C lattice
C
      j=0
      do i=1,n
      hsum=ihkl(1,i)+ihkl(2,i)
      k=1
      if(2*(hsum/2).eq.hsum)k=0
      j=j+k
      enddo
      if(j.eq.0)then
      ibr=4
      return
      endif
      RETURN
      END
C
C ------------------------------------------------------------
C
      SUBROUTINE MCMNAM(LN,NAM)
C
C Get generic filename (NAME) and its length (LN) from the command
C line. LN set to 0 if the file names are to be defined externally
C
      CHARACTER*1 KS
      CHARACTER*80 NAM,KR
      KR=' '
      CALL GETARG(IARGC(),KR)
      LN=0
      NAM=' '
        DO 2 I=1,80
        KS=KR(I:I)
        IF(KS.EQ.' ')GOTO 2
        LN=LN+1
        NAM(LN:LN)=KS
   2    CONTINUE
      RETURN
      END
C
C ------------------------------------------------------------
C
