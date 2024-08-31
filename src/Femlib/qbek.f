C     SPDX-FileCopyrightText: 2023 SAP SE
C
C     SPDX-License-Identifier: Apache-2.0
C
C     This file is part of FEDEM - https://openfedem.org
C
      SUBROUTINE QBEK31 (AEK,AELPAT,AX,AY,AALF,ADB,ADS,TH,RNY,
     +                   ISHEAR,IEL,IPRINT,LUNPRF,IERR)
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
C =======================-----------------------------------------------
C ===  FEDEM                                              MODULE: QBEK31
C =======================-----------------------------------------------
C
C     PURPOSE/
C     METHOD :
C
C     THE SUBROUTINE QBEK31 COMPUTES THE STIFFNESS MATRIX
C     OF A GENERAL QUADRILATERAL PLATE BENDING ELEMENT (WITH 4
C     NODES AND 12 NODAL DEGREES OF FREEDOM) AND THE AREA RATIOS
C     FOR LUMPING UNIFORM DISTRIBUTED LOAD TO NODAL LOAD
C
C     THE ELEMENT STIFFNESS MATRIX CREATED SATISFIES THE INDIVI-
C     DUAL ELEMENT TEST AND IS SYMMETRICAL. THE HIGHER MODES
C     (H-MODES) ONLY HAVE TO BE LINEARLY INDEPENDENT OF EACHOTHER
C     TO ENSURE MATRIX G TO BE REGULAR
C
C     ALL OF THE SHAPE FUNCTIONS (INCLUDING BOTH RC-MODES AND
C     H-MODES) ARE ASSUMED IN A LOCAL CARTESIAN COORDINATE SYSTEM
C     IN WHICH THE CENTROID OF THE ELEMENT IS TAKEN AS THE ORIGIN
C     AND THE TWO AXES ARE TAKEN PARALELL TO THE GLOBAL AXES
C
C     THE SHAPE FUNCTIONS EMPLOYED ARE:
C
C     1 , X , Y , X2 , XY , Y2 , X3 , X2Y , XY2 , Y3 , H1, H2
C
C       H1 = F1C*F2       (F1C = F1*F1*F1)
C       H2 = F2C*F1       (F2C = F2*F2*F2)
C
C     F1 AND F2 ARE FORMED BY THE EDGE EQUATIONS OF THE ELEMENT
C
C     SHEAR DEFORMATION MAY AND MAY NOT INCLUDED
C
C     THE EQUATION USED IS:
C
C       EK = (L*D*LT)/A + HHT*EKQH*HH
C
C     ARGUMENT INPUT:
C
C       AX,AY   - ARRAYS CONTAINING THE GLOBAL COORDINATE
C                 VALUES OF THE ELEMENT'S NODES
C       AALF    - THE ARRAY CONTAINNING ROTATION ANGLES
C                 AT ELEMENT NODES
C       ADB     - THE ELASTICITY MATRIX CORRESPONDING TO
C                 BENDING
C       ADS     - THE ELASTICITY MATRIX CORRESPONDING TO
C                 SHEAR
C       LUNPRF  - LOGICAL PRINT UNIT NUMBER
C
C     ARGUMENT OUTPUT:
C
C       AEK     - THE ELEMENT STIFFNESS MATRIX
C       AELPAT  - THE ARRAY CONTAINING THE RATIOS FOR LUMPING
C                 UNIFORM LOAD TO ELEMENT'S NODES
C       IERR    - ERROR FLAG
C
C     SWITCHES:
C
C       IPRINT  - PRINTING SWITCH
C        0 :      NO PRINTING
C        1 :      PRINTING
C
C       ISHEAR  - SHEAR DEFORMATION SWITCH
C        0 :      NO SHEAR DEFORMATION INCLUDED
C        1 :      SHEAR DEFORMATION INCLUDED
C
C Error KET:
C  If the the vector from node1 to node2 is orthogonal with the vector
C  from node3 to node 4 then the routine CENTPO returns NaN
C
C
C PROGRAMMED BY: XIUXI WANG
C DATE/VERSION : 82-04-06/1.0
C
C LAST UPDATE  : 86-04-30/1.1 BY: EGIL GIERTSEN - SINTEF DIV.71
C ======================================================================
C
C                                                - EXTERNAL VARIABLES
      DIMENSION AX(4),AY(4),AALF(4),ADB(3,3),ADS(2,2),AEK(12,12),
     +          AELPAT(12)
C                                                - INTERNAL VARIABLES
      DIMENSION AXY(2,4),AEKQH(6,6),AH(12,12),AHH(6,12),ALT(3,12),
     +          TRANM(12,12),TEM(12,12),TMP(144)
C                                                - COMMON
      COMMON / XPARAX / AREA,S,SI,SI2,SI3,SI4,A1,A2,B1,B2,C1,C2
C
C-----------------------------------------------------------
C     EXECUTABLE SECTION
C-----------------------------------------------------------
C
      IERR = 0
C
      IF (IPRINT .LE. 0)               GO TO   10
C
      WRITE(LUNPRF,6000) IEL
C
C----------------------------------------------------------
C     CALCULATE THE COORDINATE VALUES OF NODAL POINTS
C     REFERRED TO THE LOCAL COORDINATE SYSTEM
C----------------------------------------------------------
C
   10 CALL CENTPO (AX,AY,AXY)
C
      IF (IPRINT .LE. 0)               GO TO   30
C
      WRITE(LUNPRF,6100)
      CALL WRIME (AXY,2,4,LUNPRF)
C
   30 CALL DCOPY (12*12,0.0D0,0,AEK,1)
      CALL DCOPY (12,0.0D0,0,AELPAT,1)
C
C------------------------------------------------------------
C     CALCULATE AREA RATIOS FOR LUMPING UNIFORM LOAD
C------------------------------------------------------------
C
      AREA = 0.5*((AXY(1,1)-AXY(1,3))*(AXY(2,2)-AXY(2,4))-
     +            (AXY(1,2)-AXY(1,4))*(AXY(2,1)-AXY(2,3)))
C
      IF (AREA .LE. 0)                 GO TO 903
C
      AI   = 1.0/AREA
C
      AC1  = AXY(1,1)*AXY(2,2)-AXY(2,1)*AXY(1,2)
      AC2  = AXY(1,2)*AXY(2,3)-AXY(2,2)*AXY(1,3)
      AC3  = AXY(1,3)*AXY(2,4)-AXY(2,3)*AXY(1,4)
      AC4  = AXY(1,4)*AXY(2,1)-AXY(2,4)*AXY(1,1)
C
      AELPAT(1)  = 0.25*(AC4+AC1)
      AELPAT(4)  = 0.25*(AC1+AC2)
      AELPAT(7)  = 0.25*(AC2+AC3)
      AELPAT(10) = 0.25*(AC3+AC4)
C
      IF (IPRINT .LE. 0)               GO TO   60
C
      WRITE(LUNPRF,6900)
      CALL WRIME (AELPAT,1,12,LUNPRF)
C
C---------------------------------------------------------------
C     CALCULATE THE BASIC PART OF THE ELEMENT STIFFNESS
C     MATRIX CONTRIBUTED FROM RC-MODES
C---------------------------------------------------------------
C
   60 CALL FORMLT(ALT,AXY)
C
C     TMP(3,12) = ADB(3,3)*ALT(3,12)
      CALL DGEMM ('N','N',3,12,3,1.0D0,ADB,3,ALT,3,0.0D0,TMP,3)
C     AEK(12,12) = AEK(12,12) + AI*ALT(3,12)^t*TMP(3,12)
      CALL DGEMM ('T','N',12,12,3,AI,ALT,3,TMP,3,1.0D0,AEK,12)
C
      IF (IPRINT .LE. 0)               GO TO   80
C
      WRITE(LUNPRF,6200)
      CALL WRIME (AEK,12,12,LUNPRF)
C
C------------------------------------------------------------
C     CALCULATE THE PART OF THE ELEMENT STIFFNESS MATRIX
C     CONTRIBUTED FROM H-MODES
C------------------------------------------------------------
C
   80 CONTINUE
      S   = 0.5D0*SQRT(AREA)
      SI  = 1.0D0/S
      SI2 = SI*SI
      SI3 = SI2*SI
      SI4 = SI3*SI
C
      CALL DSCAL (8,SI,AXY,1)
C
      A1 = 0.25*(-AXY(2,1)-AXY(2,2)+AXY(2,3)+AXY(2,4))
      A2 = 0.25*(-AXY(2,1)+AXY(2,2)+AXY(2,3)-AXY(2,4))
      B1 = 0.25*(+AXY(1,1)+AXY(1,2)-AXY(1,3)-AXY(1,4))
      B2 = 0.25*(+AXY(1,1)-AXY(1,2)-AXY(1,3)+AXY(1,4))
      C1 = 0.25*(-AXY(1,1)*AXY(2,4)-AXY(1,2)*AXY(2,3)+AXY(1,3)*AXY(2,2)
     +           +AXY(1,4)*AXY(2,1))
      C2 = 0.25*(-AXY(1,1)*AXY(2,2)+AXY(1,2)*AXY(2,1)+AXY(1,3)*AXY(2,4)
     +           -AXY(1,4)*AXY(2,3))
C
      CALL CKQH2 (AEKQH,AXY,ADB,ADS,TH,RNY,ISHEAR,2,IPRINT,LUNPRF,IERR)
      IF (IERR .NE. 0)                 GO TO 901
C
      CALL FGAH2 (AH,AHH,AXY,TH,RNY,IEL,ISHEAR,IPRINT,LUNPRF,IERR)
      IF (IERR .NE. 0)                 GO TO 902
C
C     TMP(6,12) = AEKQH(6,6)*AHH(6,12)
      CALL DGEMM ('N','N',6,12,6,1.0D0,AEKQH,6,AHH,6,0.0D0,TMP,6)
C     TEM(12,12) = AHH(6,12)^t*TMP(6,12)
      CALL DGEMM ('T','N',12,12,6,1.0D0,AHH,6,TMP,6,0.0D0,TEM,12)
C
      IF (IPRINT .LE. 0)               GO TO  110
C
      WRITE(LUNPRF,6300)
      CALL WRIME (TEM,12,12,LUNPRF)
C
  110 CALL DAXPY (12*12,1.0D0,TEM,1,AEK,1)
C
C-----------------------------------------------------------
C     TRANSFORM THE ELEMENT STIFFNESS MATRIX FROM THE
C     GLOBAL COORDINATE SYSTEM TO THE INCLINED DIRECTION
C-----------------------------------------------------------
C
      CALL DCOPY (12*12,0.0D0,0,TRANM,1)
C
      DO  140 I = 1,4
C
      I3               = 3*I
      COSA             = COS(AALF(I))
      SINA             = SIN(AALF(I))
C
      TRANM(I3-2,I3-2) = 1.0D0
      TRANM(I3-1,I3-1) = COSA
      TRANM(I3-1,I3)   =-SINA
      TRANM(I3,I3-1)   = SINA
      TRANM(I3,I3)     = COSA
C
  140 CONTINUE
C
C-------------------------------------------------------------
C     CALCULATE THE ELEMENT STIFFNESS MATRIX REFERRED TO
C     THE INCLINED DIRECTION
C-------------------------------------------------------------
C
C     TMP(12,12) = AEK(12,12)*TRANM(12,12)
      CALL DGEMM ('N','N',12,12,12,1.0D0,AEK,12,TRANM,12,0.0D0,TMP,12)
C     AEK(12,12) = TRANM(12,12)^t*TMP(12,12)
      CALL DGEMM ('T','N',12,12,12,1.0D0,TRANM,12,TMP,12,0.0D0,AEK,12)
C
C--------------------------------------------------------------
C     PRINT RESULTS
C---------------------------------------------------------------
C
      IF (IPRINT .LE. 0)               GO TO 1000
C
      WRITE(LUNPRF,7000)
      CALL WRIME (TRANM,12,12,LUNPRF)
C
      WRITE(LUNPRF,7100)
      CALL WRIME (AEK,12,12,LUNPRF)
      GO TO 1000
C ------------------
C --- ERROR EXIT ---
C ------------------
  901 IERR = -1
      WRITE(LUNPRF,9000) IERR
      WRITE(LUNPRF,9010)
      GO TO 1000
  902 IERR = -2
      WRITE(LUNPRF,9000) IERR
      WRITE(LUNPRF,9020)
      GO TO 1000
  903 IERR = -3
      WRITE(LUNPRF,9000) IERR
      WRITE(LUNPRF,9030)
      GO TO 1000
C -------------------
C --- R E T U R N ---
C -------------------
 1000 CONTINUE
C
      RETURN
C ---------------
C --- FORMATS ---
C ---------------
 6000 FORMAT (/1H ,'THE ELEMENT NUMBER',I7)
 6100 FORMAT (/1H ,'THE COORDINATES OF NODES IN LOCAL SYSTEM')
 6200 FORMAT (/1H ,'THE STIFFNESS MATRIX -EKC-')
 6300 FORMAT (/1H ,'THE STIFFNESS MATRIX -EKH-')
 6900 FORMAT (/1H ,'THE AREA RATIOS FOR LUMPING UNIFORM LOAD')
C
 7000 FORMAT (/1H ,'TRANSFORMASJONSMATRISEN -TRANM-')
 7100 FORMAT (/1H ,'LOKAL ELEMENTSTIVHETSMATRISE')
C
 9000 FORMAT(/1H ,'*** FEIL VED RETUR FRA RUTINE -QBEK31-',
     +            ', IERR =',I4,' ***'/)
 9010 FORMAT( 1H ,'    IKKE-NULL FEILFLAGG MOTTATT FRA RUTINE',
     +            ' -CKQH2-')
 9020 FORMAT( 1H ,'    IKKE-NULL FEILFLAGG MOTTATT FRA RUTINE',
     +            ' -FGAH2-')
 9030 FORMAT( 1H ,'    NEGATIVT AREAL                        ')
      END
C
C
      SUBROUTINE QBEF2 (AESTRS,AX,AY,AALF,ADB,TH,RNY,
     +                  ISHEAR,IEL,IPRINT,LUNPRF,IERR)
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
C =======================-----------------------------------------------
C === FEDEM                                               MODULE: QBEF2
C =======================-----------------------------------------------
C
C     PURPOSE/
C     METHOD :
C
C     THE SUBROUTINE QBEF2 COMPUTES THE ELEMENT STRESS MATRIX AT
C     FOUR NODES OF AN ELEMENT INCLUDING MOMENTS AND SHEAR FORCES
C     ALL SHAPE FUNCTIONS (INCLUDING RC-MODES AND H-MODES) ARE
C     EXPRESSED IN A CARTESIAN COORDINATE SYSTEM
C
C     THE SHAPE FUNCTIONS EMPLOYED ARE:
C
C     1 , X , Y , X2 , XY , Y2 , X3 , X2Y , XY2 , Y3 , H1, H2.
C
C       H1 = F1C*F2        (F1C = F1*F1*F1)
C       H2 = F2C*F1        (F2C = F2*F2*F2)
C
C     F1 AND F2 ARE FORMED BY THE EDGE EQUATIONS OF THE ELEMENT
C
C     SHAPE DEFORMATION MAY AND MAY NOT BE INCLUDED
C
C     THE EQUATION USED IS:
C
C       ESTRS = D*B*H*T
C
C       D  - ELASTICITY MATRIX
C
C            ! DB  0  !
C       D  = !        !
C            ! 0   DS !
C
C       DB - ELASTICITY MATRIX CORRESPONDING TO BENDING
C       DS - ELASTICITY MATRIX CORRESPONDING TO SHEAR
C
C       B  - STRAIN MATRIX CORRESPONDING TO BENDING DEFORMATION
C
C            !  BB  !
C       B  = !      !
C            !  BS  !
C
C       BB - STRAIN MATRIX CORRESPONDING TO BENDING
C       BS - STRAIN MATRIX CORRESPONDING TO SHEAR
C
C       H  - THE MATRIX H (INVERSION OF THE MATRIX G)
C       T  - THE TRANSFORMATION MATRIX FOR THE DISPLACEMENT VECTOR
C            FROM GLOBAL COORDINATE SYSTEM TO INCLINED DIRECTIONS
C
C     ARGUMENT INPUT:
C
C       AX,AY      - THE ARRAYS CONTAINNING THE COORDINATE
C                    VALUES OF THE NODES REFFERED TO THE
C                    GLOBAL COORDINATE SYSTEM
C       AALF       - THE ARRAY CONTAINNING THE INCLINING
C                    ANGLES OF ELEMENT'S NODES
C       ADB        - ELASTICITY MATRIX CORRESPONDING TO
C                    BENDING
C       LUNPRF     - LOGICAL PRINT UNIT NUMBER
C
C     ARGUMENT OUTPUT:
C
C       AESTRS     - THE STRESS MATRIX AT FOUR NODES OF AN
C                    ELEMENT
C       IERR       - ERROR FLAG
C
C     SWITCHES:
C
C       IPRINT     - PRINT SWITCH
C        0 :          NO PRINT
C        1 :          PRINT
C
C       ISHEAR     - SHEAR DEFORMATION SWITCH
C        0 :          NO SHEAR DEFORMATION INCLUDED
C        1 :          SHEAR DEFORMATION INCLUDED
C
C PROGRAMMED BY: XIUXI WANG
C DATE/VERSION : 82-04-06/1.0
C
C LAST UPDATE  : 86-04-30/1.1 BY: EGIL GIERTSEN - SINTEF DIV.71
C ======================================================================
C
C                                                - EXTERNAL VARIABLES
      DIMENSION AESTRS(20,12),AX(4),AY(4),AALF(4),ADB(3,3)
      DIMENSION AXY(2,4),BB(3,12),BS(2,12),AH(12,12),AHH(6,12),
     +          TRANM(12,12),TEM1(3,12),ANSTRB(3,12),
     +          ANSTRS(2,12),AHMTRA(12,12)
C                                                - COMMON
      COMMON / XPARAX / AREA,S,SI,SI2,SI3,SI4,A1,A2,B1,B2,C1,C2
C
C-------------------------------------------------------------------
C     EXECUTABLE SECTION
C-------------------------------------------------------------------
C
      IERR = 0
      H1X3 = 0.0
      H2X3 = 0.0
      H1Y3 = 0.0
      H2Y3 = 0.0
      H1XY2 = 0.0
      H2XY2 = 0.0
      H1X2Y = 0.0
      H2X2Y = 0.0
      P    = TH*TH/(1.0-RNY)/5.0
C
      CALL CENTPO (AX,AY,AXY)
C
C---------------------------------------------------------------------
C     FORM THE TRANSFORMATION MATRIX TRANM WHICH TRANSFORMS THE
C     DISPLACEMENT FROM THE GLOBAL COORDINATE SYSTEM TO THE INCLINED
C     DIRECTIONS
C---------------------------------------------------------------------
C
      CALL DCOPY (12*12,0.0D0,0,TRANM,1)
C
      DO   40 I = 1,4
C
      COSA               = COS(AALF(I))
      SINA               = SIN(AALF(I))
C
      TRANM(3*I-2,3*I-2) = 1.0
      TRANM(3*I-1,3*I-1) = COSA
      TRANM(3*I-1,3*I)   =-SINA
      TRANM(3*I,3*I-1)   = SINA
      TRANM(3*I,3*I)     = COSA
   40 CONTINUE
C
C--------------------------------------------------------------
C     NORMALIZE THE COORDINATE VALUES OF NODES
C--------------------------------------------------------------
C
      AREA = 0.5*((AXY(1,1)-AXY(1,3))*(AXY(2,2)-AXY(2,4))-
     +            (AXY(1,2)-AXY(1,4))*(AXY(2,1)-AXY(2,3)))
C
      IF (AREA .LE. 0)                 GO TO 903
C
      S    = 0.5*SQRT(AREA)
      SI   = 1.0/S
      SI2  = SI*SI
      SI3  = SI2*SI
      SI4  = SI3*SI
C
      CALL DSCAL (8,SI,AXY,1)
C
      A1 = 0.25*(-AXY(2,1)-AXY(2,2)+AXY(2,3)+AXY(2,4))
      A2 = 0.25*(-AXY(2,1)+AXY(2,2)+AXY(2,3)-AXY(2,4))
      B1 = 0.25*(+AXY(1,1)+AXY(1,2)-AXY(1,3)-AXY(1,4))
      B2 = 0.25*(+AXY(1,1)-AXY(1,2)-AXY(1,3)+AXY(1,4))
      C1 = 0.25*(-AXY(1,1)*AXY(2,4)-AXY(1,2)*AXY(2,3)+AXY(1,3)*AXY(2,2)
     +           +AXY(1,4)*AXY(2,1))
      C2 = 0.25*(-AXY(1,1)*AXY(2,2)+AXY(1,2)*AXY(2,1)+AXY(1,3)*AXY(2,4)
     +           -AXY(1,4)*AXY(2,3))
C
C------------------------------------------------------------------
C     CALCULATE THE PRODUCT OF MATRICES H AND T FOR LATER USE
C------------------------------------------------------------------
C
      CALL FGAH2 (AH,AHH,AXY,TH,RNY,IEL,ISHEAR,IPRINT,LUNPRF,IERR)
      IF (IERR .NE. 0)                 GO TO 901
C
C     AHMTRA(12,12) = AH(12,12)*TRANM(12,12)
      CALL DGEMM ('N','N',12,12,12,1.0D0,AH,12,TRANM,12,0.0D0,AHMTRA,12)
C
      CALL DCOPY (3*12,0.0D0,0,BB,1)
      CALL DCOPY (2*12,0.0D0,0,BS,1)
C
C-----------------------------------------------------------------
C     FORM THE STRAIN MATRIX CORRESPONDING TO BENDING
C-----------------------------------------------------------------
C
      DO  130 I = 1,4
C
      F1    = A1*AXY(1,I)+B1*AXY(2,I)+C1
      F2    = A2*AXY(1,I)+B2*AXY(2,I)+C2
C
      H1X2  = 6.0*A1*F1*(A1*F2+A2*F1)*SI2
      H2X2  = 6.0*A2*F2*(A2*F1+A1*F2)*SI2
      H1Y2  = 6.0*B1*F1*(B1*F2+B2*F1)*SI2
      H2Y2  = 6.0*B2*F2*(B2*F1+B1*F2)*SI2
      H1XY  = 6.0*F1*(2.0*A1*B1*F2+A1*B2*F1+A2*B1*F1)*SI2
      H2XY  = 6.0*F2*(2.0*A2*B2*F1+A2*B1*F2+A1*B2*F2)*SI2
C
      BB(1,4)  = -2.0*SI2
      BB(1,7)  = -6.0*AXY(1,I)*SI2
      BB(1,8)  = -2.0*AXY(2,I)*SI2
      BB(1,11) = -H1X2
      BB(1,12) = -H2X2
      BB(2,6)  = -2.0*SI2
      BB(2,9)  = -2.0*AXY(1,I)*SI2
      BB(2,10) = -6.0*AXY(2,I)*SI2
      BB(2,11) = -H1Y2
      BB(2,12) = -H2Y2
      BB(3,5)  = -2.0*SI2
      BB(3,8)  = -4.0*AXY(1,I)*SI2
      BB(3,9)  = -4.0*AXY(2,I)*SI2
      BB(3,11) = -H1XY
      BB(3,12) = -H2XY
C
      IF (ISHEAR .LE. 0)               GO TO  100
C
      H1X3  = 6.0*A1*A1*(A1*F2+3.0*A2*F1)*SI3
      H2X3  = 6.0*A2*A2*(A2*F1+3.0*A1*F2)*SI3
      H1X2Y = 6.0*A1*(A1*B1*F2+2.0*A2*B1*F1+A1*B2*F1)*SI3
      H2X2Y = 6.0*A2*(A2*B2*F1+2.0*A1*B2*F2+A2*B1*F2)*SI3
      H1XY2 = 6.0*B1*(A1*B1*F2+2.0*A1*B2*F1+A2*B1*F1)*SI3
      H2XY2 = 6.0*B2*(A2*B2*F1+2.0*A2*B1*F2+A1*B2*F2)*SI3
      H1Y3  = 6.0*B1*B1*(B1*F2+3.0*B2*F1)*SI3
      H2Y3  = 6.0*B2*B2*(B2*F1+3.0*B1*F2)*SI3
C
      H1X4  = 24.0*A1*A1*A1*A2*SI4
      H2X4  = 24.0*A2*A2*A2*A1*SI4
      H1X3Y = 6.0*A1*A1*(A1*B2+3.0*A2*B1)*SI4
      H2X3Y = 6.0*A2*A2*(A2*B1+3.0*A1*B2)*SI4
      H1X2Y2= 12.0*A1*B1*(A1*B2+A2*B1)*SI4
      H2X2Y2= 12.0*A2*B2*(A2*B1+A1*B2)*SI4
      H1XY3 = 6.0*B1*B1*(B1*A2+3.0*A1*B2)*SI4
      H2XY3 = 6.0*B2*B2*(B2*A1+3.0*A2*B1)*SI4
      H1Y4  = 24.0*B1*B1*B1*B2*SI4
      H2Y4  = 24.0*B2*B2*B2*B1*SI4
C
      BB(1,11) = BB(1,11) - P*(H1X4+H1X2Y2)
      BB(1,12) = BB(1,12) - P*(H2X4+H2X2Y2)
      BB(2,11) = BB(2,11) - P*(H1Y4+H1X2Y2)
      BB(2,12) = BB(2,12) - P*(H2Y4+H2X2Y2)
      BB(3,11) = BB(3,11) - 2.0*P*(H1X3Y+H1XY3)
      BB(3,12) = BB(3,12) - 2.0*P*(H2X3Y+H2XY3)
C
C------------------------------------------------------------------
C     CALCULATE THE PART OF ELEMENT STRESS MATRIX CORRESPONDING
C     TO MOMENTS
C------------------------------------------------------------------
C
  100 CONTINUE
C     TEM1(3,12) = BB(3,12)*AHMTRA(12,12)
      CALL DGEMM ('N','N',3,12,12,1.0D0,BB,3,AHMTRA,12,0.0D0,TEM1,3)
C     ANSTRB(3,12) = ADB(3,3)*TEM1(3,12)
      CALL DGEMM ('N','N',3,12,3,1.0D0,ADB,3,TEM1,3,0.0D0,ANSTRB,3)
C
      DO  110 L     = 1,3
      DO  110 K     = 1,12
      IP1           = 5*I-5+L
      AESTRS(IP1,K) = ANSTRB(L,K)
  110 CONTINUE
C
C-------------------------------------------------------------
C     CALCULATE THE PART OF THE ELEMENT STRESS MATRIX
C     CORRESPONDING TO SHEAR FORCES
C-------------------------------------------------------------
C
      D         = ADB(1,1)
C
      BS(1,7)   = -6.0*D*SI3
      BS(2,8)   = -2.0*D*SI3
      BS(1,9)   = -2.0*D*SI3
      BS(2,10)  = -6.0*D*SI3
      BS(1,11)  = -D*(H1X3 + H1XY2)
      BS(1,12)  = -D*(H2X3 + H2XY2)
      BS(2,11)  = -D*(H1Y3 + H1X2Y)
      BS(2,12)  = -D*(H2Y3 + H2X2Y)
C
C     ANSTRS(2,12) = BS(2,12)*AHMTRA(12,12)
      CALL DGEMM ('N','N',2,12,12,1.0D0,BS,2,AHMTRA,12,0.0D0,ANSTRS,2)
C
      DO  120 L     = 1,2
      DO  120 K     = 1,12
      IP1           = 5*I-2+L
      AESTRS(IP1,K) = ANSTRS(L,K)
  120 CONTINUE
C
  130 CONTINUE
C
      IF (IPRINT .LE. 0)               GO TO 1000
C
      WRITE(LUNPRF,6000)
      CALL WRIMF2 (AESTRS,20,12,LUNPRF)
C ------------------
C --- ERROR EXIT ---
C ------------------
  901 IERR = -1
      WRITE(LUNPRF,9000) IERR
      WRITE(LUNPRF,9010)
      GO TO 1000
  903 IERR = -3
      WRITE(LUNPRF,9000) IERR
      WRITE(LUNPRF,9030)
      GO TO 1000
C -------------------
C --- R E T U R N ---
C -------------------
 1000 CONTINUE
      RETURN
C
C ---------------
C --- FORMATS ---
C ---------------
 6000 FORMAT (/1H ,'SPENNINGSMATRISEN FOR ELEMENTET')
C
 9000 FORMAT(/1H ,'*** FEIL VED RETUR FRA RUTINE -QBEF2-',
     +            ', IERR =',I4,' ***'/)
 9010 FORMAT( 1H ,'    IKKE-NULL FEILFLAGG MOTTATT FRA RUTINE',
     +            ' -FGAH2-')
 9030 FORMAT( 1H ,'    NEGATIVT AREAK                        ',
     +            ' -FGAH2-')
      END
C
C
      SUBROUTINE FGAH2 (FH,FHH,FXY,TH,RNY,IEL,ISHEAR,IPRINT,
     +                  LUNPRF,IERR)
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
C =======================-----------------------------------------------
C === FEDEM                                              MODULE: FGAH2
C =======================-----------------------------------------------
C
C     PURPOSE/
C     METHOD :
C
C     THE SUBROUTINE FGAH2 FORMS THE MATRIX G, INVERTS IT
C     TO GET THE MATRIX H AND PUTS THE H-PART OF THE MATRIX
C     H INTO THE MATRIX HH
C
C     ALL THE SHAPE FUNCTIONS ARE EXPRESSED IN A CARTESIAN
C     COORDINATE SYSTEM
C
C     THE SHAPE FUNCTIONS USED ARE:
C
C     1 , X , Y , X2 , XY , Y2 , X3 , X2Y , XY2 , Y3 , H1, H2
C
C       H1 = F1C*F2     (F1C = F1*F1*F1)
C       H2 = F2C*F1     (F2C = F2*F2*F2)
C
C     F1 F2 ARE FORMED BY THE EDGE EQUATIONS OF THE ELEMENT.
C
C     SHEAR DEFORMATION MAY AND MAY NOT BE INCLUDED
C
C     ARGUMENT INPUT:
C
C       FXY    - THE ARRAY CONTAINNING THE COORDINATE VALUES
C                OF THE ELEMENT'S NODES , REFERED TO THE LOCAL
C                COORDINATE SYSTEM AND NORMALIZED BY HALF OF
C                SQUARE ROOT OF THE ELEMENT AREA
C       LUNPRF - LOGICAL PRINT UNIT NUMBER
C
C     ARGUMENT OUTPUT:
C
C       FH     - THE ARRAY CONTAINNING THE MATRIX H OF INVERTED
C                MATRIX G
C       FHH    - THE ARRAY CONTAINNING THE H-PART OF INVERTED
C                MATRIX G
C       IERR   - ERROR FLAG
C
C     SWITCHES:
C
C       IPRINT - PRINT SWITCH
C        0 :     NO PRINT
C        1 :     PRINT
C
C       ISHEAR - SHEAR DEFORMATION SWITCH
C        0 :     NO SHEAR DEFORMATION INCLUDED
C        1 :     SHEAR DEFORMATION INCLUDED
C
C PROGRAMMED BY: XIUXI WANG
C DATE/VERSION : 82-04-06/1.0
C
C LAST UPDATE  : 86-04-30/1.1 BY: EGIL GIERTSEN - SINTEF DIV.71
C ======================================================================
C
C                                                - EXTERNAL VARIABLES
      DIMENSION FG(12,12),FH(12,12),FHH(6,12),FXY(2,4)
C                                                - INTERNAL VARIABLES
      DIMENSION IWORK(12)
C                                                - COMMON
      COMMON / XPARAX / AREA,S,SI,SI2,SI3,SI4,A1,A2,B1,B2,C1,C2
C
C----------------------------------------------------------------
C     EXECUTABLE SECTION
C----------------------------------------------------------------
C
      IERR  = 0
C
      IF (ISHEAR .LE. 0)               GO TO  10
C
      P     = TH*TH/(1.0-RNY)/5.0
C
   10 DO 40 I = 1, 4
C
      X     = FXY(1,I)
      Y     = FXY(2,I)
      X2    = X*X
      Y2    = Y*Y
C
      F1    = A1*X+B1*Y+C1
      F2    = A2*X+B2*Y+C2
C
      F1S   = F1*F1
      F2S   = F2*F2
      F1C   = F1S*F1
      F2C   = F2S*F2
C
      H1    = F1C*F2
      H2    = F2C*F1
C
      H1X   = F1S*(3.0*A1*F2+A2*F1)*SI
      H2X   = F2S*(3.0*A2*F1+A1*F2)*SI
      H1Y   = F1S*(3.0*B1*F2+B2*F1)*SI
      H2Y   = F2S*(3.0*B2*F1+B1*F2)*SI
C
      IF (ISHEAR .LE. 0)               GO TO   30
C
      H1X3  = 6.0*A1*A1*(A1*F2+3.0*A2*F1)*SI3
      H2X3  = 6.0*A2*A2*(A2*F1+3.0*A1*F2)*SI3
      H1X2Y = 6.0*A1*(A1*B1*F2+2.0*A2*B1*F1+A1*B2*F1)*SI3
      H2X2Y = 6.0*A2*(A2*B2*F1+2.0*A1*B2*F2+A2*B1*F2)*SI3
      H1XY2 = 6.0*B1*(A1*B1*F2+2.0*A1*B2*F1+A2*B1*F1)*SI3
      H2XY2 = 6.0*B2*(A2*B2*F1+2.0*A2*B1*F2+A1*B2*F2)*SI3
      H1Y3  = 6.0*B1*B1*(B1*F2+3.0*B2*F1)*SI3
      H2Y3  = 6.0*B2*B2*(B2*F1+3.0*B1*F2)*SI3
C
C-----------------------------------------------------------
C     FORM THE PART OF THE MATRIX G CORRESPONDING TO BENDING
C-----------------------------------------------------------
C
   30 I3         = 3*I
C
      FG(I3-2,1) = 1.0
      FG(I3-1,1) = 0.0
      FG(I3,1)   = 0.0
      FG(I3-2,2) = X
      FG(I3-1,2) = 0.0
      FG(I3,2)   =-SI
      FG(I3-2,3) = Y
      FG(I3-1,3) = SI
      FG(I3,3)   = 0.0
      FG(I3-2,4) = X2
      FG(I3-1,4) = 0.0
      FG(I3,4)   =-2.0*X*SI
      FG(I3-2,5) = X*Y
      FG(I3-1,5) = X*SI
      FG(I3,5)   =-Y*SI
      FG(I3-2,6) = Y2
      FG(I3-1,6) = 2.0*Y*SI
      FG(I3,6)   = 0.0
      FG(I3-2,7) = X*X2
      FG(I3-1,7) = 0.0
      FG(I3,7)   =-3.0*X2*SI
      FG(I3-2,8) = X2*Y
      FG(I3-1,8) = X2*SI
      FG(I3,8)   =-2.0*X*Y*SI
      FG(I3-2,9) = X*Y2
      FG(I3-1,9) = 2.0*X*Y*SI
      FG(I3,9)   =-1.0*Y2*SI
      FG(I3-2,10)= Y*Y2
      FG(I3-1,10)= 3.0*Y2*SI
      FG(I3,10)  = 0.0
      FG(I3-2,11)= H1
      FG(I3-1,11)= H1Y
      FG(I3,11)  =-H1X
      FG(I3-2,12)= H2
      FG(I3-1,12)= H2Y
      FG(I3,12)  =-H2X
C
      IF (ISHEAR .LE. 0)               GO TO  20
C
C------------------------------------------------------------
C     FORM THE PART OF THE MATRIX G CORRESPONDING TO SHEAR
C------------------------------------------------------------
C
      FG(I3,7)   = FG(I3,7)   - 6.0*P*SI3
      FG(I3-1,8) = FG(I3-1,8) + 2.0*P*SI3
      FG(I3,9)   = FG(I3,9)   - 2.0*P*SI3
      FG(I3-1,10)= FG(I3-1,10)+ 6.0*P*SI3
      FG(I3-1,11)= FG(I3-1,11)+ P*(H1X2Y + H1Y3 )
      FG(I3,11)  = FG(I3,11)  - P*(H1X3  + H1XY2)
      FG(I3-1,12)= FG(I3-1,12)+ P*(H2X2Y + H2Y3 )
      FG(I3,12)  = FG(I3,12)  - P*(H2X3  + H2XY2)
C
   20 CONTINUE
C
   40 CONTINUE
C
      IF (IPRINT .LE. 0)               GO TO   50
C
      WRITE(LUNPRF,6000)
      CALL WRIMF2 (FG,12,12,LUNPRF)
C
C---------------------------------------------------------
C     INVERTS THE MATRIX G TO GET THE MATRIX H
C---------------------------------------------------------
C
   50 CALL DINV12 (12,FG,FH,IWORK,LUNPRF,IERR)
      IF (IERR .NE. 0)                 GO TO  901
C
      CALL DCOPY (12*12,FG(1,1),1,FH,1)
C
      IF (IPRINT .LE. 0)               GO TO   60
C
      WRITE(LUNPRF,6100)
      CALL WRIMF2 (FH,12,12,LUNPRF)
C
C---------------------------------------------------------
C     PUT THE H-MODE PART OF THE MATRIX H INTO THE MATRIX HH
C---------------------------------------------------------
C
   60 DO  70 I = 1,6
      DO  70 J = 1,12
      FHH(I,J) = FH(I+6,J)
   70 CONTINUE
C
      GO TO 1000
C ------------------
C --- ERROR EXIT ---
C ------------------
  901 IERR = -1
      WRITE(LUNPRF,9000) IERR
      WRITE(LUNPRF,9010) IEL
      GO TO 1000
C -------------------
C --- R E T U R N ---
C -------------------
 1000 CONTINUE
C
      RETURN
C ---------------
C --- FORMATS ---
C ---------------
 6000 FORMAT(/1H ,'MATRISEN -G-')
 6100 FORMAT(/1H ,'MATRISEN -H-')
C
 9000 FORMAT(/1H ,'*** FEIL VED RETUR FRA RUTINE -FGAH2-',
     +            ', IERR =',I4,' ***'/)
 9010 FORMAT( 1H ,'    MATRISEN -G- ER SINGUL[R FOR ELEMENT NUMMER',I7)
      END
C
C
      SUBROUTINE CKQH2 (FEKQH,FXY,FDB,FDS,TH,RNY,ISHEAR,N,IPRINT,
     +                  LUNPRF,IERR)
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
C =======================-----------------------------------------------
C === M O B Y   D I K ===    PART: PLAT2                  MODULE: CKQH2
C =======================-----------------------------------------------
C
C     PURPOSE/
C     METHOD :
C
C     THE SUBROUTINE CKQH2 COMPUTES THE GENERALIZED STIFF-
C     NESS MATRIX CONTRIBUTED FROM H-MODES BY USING NUMERIC
C     INTEGRATION FORMULAS. THE NUMBER OF INTEGRATION POINTS
C     COULD BE 2*2 OR 3*3
C
C     ALL THE SHAPE FUNCTIONS ARE EXPRESSED IN A CARTESIAN
C     COORDINATE SYSTEM
C
C     THE H-MODES EMPLOYED ARE:
C
C       X3 , X2Y , XY2 , Y3 , H1, H2
C
C     SHEAR DEFORMATION MAY AND MAY NOT INCLUDED
C
C     THE FORMULA USED IS:
C
C         EKQH = EKQHB + EKQHS
C              = (INT.(BBHT*DB*BBH))*DA +
C                (INT.(BSHT*DS*BSH))*DA
C
C     ARGUMENT INPUT:
C
C        FXY     - ARRAY CONTAINING THE LOCAL COORDINATE
C                  VALUES OF THE ELEMENT'S NODES REFERRED
C                  TO THE LOCAL COORDINATE SYSTEM AND
C                  NORMALIZED BY HALF OF SQUARE ROOT OF
C                  THE ELEMENT AREA
C        FDB     - THE ELASTICITY MATRIX CORRESPONDING
C                  TO BENDING
C        FDS     - THE ELASTICITY MATRIX CORRESPONDING
C                  TO SHEAR
C        N       - NUMBER OF INTEGRATION POINTS IN ONE
C                  DIRECTION
C        LUNPRF  - LOGICAL PRINT UNIT NUMBER
C
C      ARGUMENT OUTPUT:
C
C        FKQH    - THE GENERALIZED STIFFNESS MATRIX CONTRIBUTED
C                  FORM H-MODES
C        IERR    - ERROR FLAG
C
C      SWITCHES:
C
C        IPRINT  - PRINT SWITCH
C         0 :      NO PRINT
C         1 :      PRINT
C
C        ISHEAR  - SHEAR DEFORMATION SWITCH
C         0 :      NO SHEAR DEFORMATION INCLUDED
C         1 :      SHEAR DEFORMATION INCLUDED
C
C PROGRAMMED BY: XIUXI WANG
C DATE/VERSION : 82-04-25/1.0
C
C LAST UPDATE  : 86-04-30/1.1 BY: EGIL GIERTSEN - SINTEF DIV.71
C ======================================================================
C
C                                                - EXTERNAL VARIABLES
      DIMENSION FXY(2,4),FEKQH(6,6),XIETA(3),WEI(3),
     +          FDB(3,3),FDS(2,2),BH(3,6),BSH(2,6),TMP(18)
C                                                - COMMON
      COMMON / XPARAX / AREA,S,SI,SI2,SI3,SI4,A1,A2,B1,B2,C1,C2
C
C------------------------------------------------------------
C     EXECUTABLE SECTION
C------------------------------------------------------------
C
      IERR     = 0
      H1X3     = 0.0
      H2X3     = 0.0
      H1Y3     = 0.0
      H2Y3     = 0.0
      H1X2Y    = 0.0
      H2X2Y    = 0.0
      H1XY2    = 0.0
      H2XY2    = 0.0
      P        = TH*TH/(1.0-RNY)/5.0
C
      IF (N .LE. 2)                    GO TO 20
C
      XIETA(1) = -0.774597
      XIETA(2) =  0.0
      XIETA(3) =  0.774597
      WEI(1)   =  0.555556
      WEI(2)   =  0.888889
      WEI(3)   =  0.555556
C
      GO TO 30
C
   20 XIETA(1) = -0.577350
      XIETA(2) = +0.577350
      WEI(1)   =  1.0
      WEI(2)   =  1.0
C
   30 CALL DCOPY (6*6,0.0D0,0,FEKQH,1)
      CALL DCOPY (3*6,0.0D0,0,BH,1)
      CALL DCOPY (2*6,0.0D0,0,BSH,1)
C
      DO 150 I = 1,N
      DO 150 J = 1,N
C
C--------------------------------------------------------
C     CALCULATE THE BASIC PARAMETERS 1+XI, 1-XI, 1+ETA
C     AND 1-ETA FOR ALL INTEGRATION POINT
C--------------------------------------------------------
C
      XIP  = 1.0+XIETA(I)
      XIS  = 1.0-XIETA(I)
      ETAP = 1.0+XIETA(J)
      ETAS = 1.0-XIETA(J)
C
C----------------------------------------------------------
C     CALCULATE THE JACOBIAN DETERMINANT J CONNECTED TO
C     NUMERICAL INTEGRATION, MULTIPLIED BY WEIGHT COEFFICIENTS
C----------------------------------------------------------
C
      DXDXI  = (ETAP*FXY(1,1)-ETAP*FXY(1,2)-ETAS*FXY(1,3)+
     +          ETAS*FXY(1,4))/4.0
      DXDETA = (XIP*FXY(1,1)+XIS*FXY(1,2)-XIS*FXY(1,3)-
     +          XIP*FXY(1,4))/4.0
      DYDXI  = (ETAP*FXY(2,1)-ETAP*FXY(2,2)-ETAS*FXY(2,3)+
     +          ETAS*FXY(2,4))/4.0
      DYDETA = (XIP*FXY(2,1)+XIS*FXY(2,2)-XIS*FXY(2,3)-
     +          XIP*FXY(2,4))/4.0
C
      DETJW  = S*S*(DXDXI*DYDETA-DXDETA*DYDXI) * WEI(I)*WEI(J)
C
C----------------------------------------------------------
C     CALCULATE THE COORDINATES X,Y EXPRESSED IN
C     TERMS OF NORMALIZED COORDINATES XI AND ETA
C----------------------------------------------------------
C
      X     = (XIP*ETAP*FXY(1,1)+XIS*ETAP*FXY(1,2)+XIS*ETAS*
     +                  FXY(1,3)+XIP*ETAS*FXY(1,4))/4.0
      Y     = (XIP*ETAP*FXY(2,1)+XIS*ETAP*FXY(2,2)+XIS*ETAS*
     +                  FXY(2,3)+XIP*ETAS*FXY(2,4))/4.0
C
C----------------------------------------------------------
C     FORM THE STRAIN MATRIX BH CONTRIBUTED FROM BENDING
C----------------------------------------------------------
C
      F1    = A1*X+B1*Y+C1
      F2    = A2*X+B2*Y+C2
C
      H1X2  = 6.0*A1*F1*(A1*F2+A2*F1)*SI2
      H2X2  = 6.0*A2*F2*(A2*F1+A1*F2)*SI2
      H1Y2  = 6.0*B1*F1*(B1*F2+B2*F1)*SI2
      H2Y2  = 6.0*B2*F2*(B2*F1+B1*F2)*SI2
      H1XY  = 6.0*F1*(2.0*A1*B1*F2+A1*B2*F1+A2*B1*F1)*SI2
      H2XY  = 6.0*F2*(2.0*A2*B2*F1+A2*B1*F2+A1*B2*F2)*SI2
C
      BH(1,1) = 6.0*X*SI2
      BH(1,2) = 2.0*Y*SI2
      BH(1,5) = H1X2
      BH(1,6) = H2X2
      BH(2,3) = 2.0*X*SI2
      BH(2,4) = 6.0*Y*SI2
      BH(2,5) = H1Y2
      BH(2,6) = H2Y2
      BH(3,2) = 4.0*X*SI2
      BH(3,3) = 4.0*Y*SI2
      BH(3,5) = H1XY
      BH(3,6) = H2XY
C
      IF (ISHEAR .LE. 0)               GO TO 100
C
      H1X3  = 6.0*A1*A1*(A1*F2+3.0*A2*F1)*SI3
      H2X3  = 6.0*A2*A2*(A2*F1+3.0*A1*F2)*SI3
      H1X2Y = 6.0*A1*(A1*B1*F2+2.0*A2*B1*F1+A1*B2*F1)*SI3
      H2X2Y = 6.0*A2*(A2*B2*F1+2.0*A1*B2*F2+A2*B1*F2)*SI3
      H1XY2 = 6.0*B1*(A1*B1*F2+2.0*A1*B2*F1+A2*B1*F1)*SI3
      H2XY2 = 6.0*B2*(A2*B2*F1+2.0*A2*B1*F2+A1*B2*F2)*SI3
      H1Y3  = 6.0*B1*B1*(B1*F2+3.0*B2*F1)*SI3
      H2Y3  = 6.0*B2*B2*(B2*F1+3.0*B1*F2)*SI3
C
      H1X4  = 24.0*A1*A1*A1*A2*SI4
      H2X4  = 24.0*A2*A2*A2*A1*SI4
      H1X3Y = 6.0*A1*A1*(A1*B2+3.0*A2*B1)*SI4
      H2X3Y = 6.0*A2*A2*(A2*B1+3.0*A1*B2)*SI4
      H1X2Y2= 12.0*A1*B1*(A1*B2+A2*B1)*SI4
      H2X2Y2= 12.0*A2*B2*(A2*B1+A1*B2)*SI4
      H1XY3 = 6.0*B1*B1*(B1*A2+3.0*A1*B2)*SI4
      H2XY3 = 6.0*B2*B2*(B2*A1+3.0*A2*B1)*SI4
      H1Y4  = 24.0*B1*B1*B1*B2*SI4
      H2Y4  = 24.0*B2*B2*B2*B1*SI4
C
      BH(1,5) = BH(1,5) + P*(H1X4+H1X2Y2)
      BH(1,6) = BH(1,6) + P*(H2X4+H2X2Y2)
      BH(2,5) = BH(2,5) + P*(H1Y4+H1X2Y2)
      BH(2,6) = BH(2,6) + P*(H2Y4+H2X2Y2)
      BH(3,5) = BH(3,5) + 2.0*P*(H1X3Y+H1XY3)
      BH(3,6) = BH(3,6) + 2.0*P*(H2X3Y+H2XY3)
C
C---------------------------------------------------------
C     CALCULATE THE MATRIX  KQH
C---------------------------------------------------------
C
  100 CONTINUE
C     TMP(3,6) = DETJW*FDB(3,3)*BH(3,6)
      CALL DGEMM ('N','N',3,6,3,DETJW,FDB,3,BH,3,0.0D0,TMP,3)
C     FEKQH(6,6) = FEKQH(6,6) + BH(3,6)^t*TMP(3,6)
      CALL DGEMM ('T','N',6,6,3,1.0D0,BH,3,TMP,3,1.0D0,FEKQH,6)
C
C---------------------------------------------------------
C     FORM THE STRAIN MATRIX BSH CORRESPONDING TO SHEAR
C---------------------------------------------------------
C
      IF (ISHEAR .LE. 0)               GO TO 140
C
      BSH(1,1) = 6.0*P*SI3
      BSH(2,2) = 2.0*P*SI3
      BSH(1,3) = 2.0*P*SI3
      BSH(2,4) = 6.0*P*SI3
      BSH(1,5) = P*(H1X3 + H1XY2)
      BSH(1,6) = P*(H2X3 + H2XY2)
      BSH(2,5) = P*(H1Y3 + H1X2Y)
      BSH(2,6) = P*(H2Y3 + H2X2Y)
C
C     TMP(2,6) = DETJW*FDS(2,2)*BSH(2,6)
      CALL DGEMM ('N','N',2,6,2,DETJW,FDS,2,BSH,2,0.0D0,TMP,2)
C     FEKQH(6,6) = FEKQH(6,6) + BSH(2,6)^t*TMP(2,6)
      CALL DGEMM ('T','N',6,6,2,1.0D0,BSH,2,TMP,2,1.0D0,FEKQH,6)
C
  140 CONTINUE
C
  150 CONTINUE
C
      IF (IPRINT .LE. 0)               GO TO 1000
C
      WRITE(LUNPRF,6000)
      CALL WRIMF2 (FEKQH,6,6,LUNPRF)
C ------------------
C --- ERROR EXIT ---
C ------------------
C -------------------
C --- R E T U R N ---
C -------------------
 1000 CONTINUE
C
      RETURN
C ---------------
C --- FORMATS ---
C ---------------
 6000 FORMAT (/1H ,'THE MATRIX -KQH-')
      END
C
C
      SUBROUTINE FORMLT(FLT,FXY)
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
C =======================-----------------------------------------------
C === M O B Y   D I K ===    PART: PLAT2                  MODULE: FORMLT
C =======================-----------------------------------------------
C
C     PURPOSE/
C     METHOD :
C
C     THE SUBROUTINE FLT FORMS THE MATRIX LT
C
C     THE MATRIX L REPRESENTS THE LUMPING OF TRACTIONS
C     CORRESPONDING TO THE CONSTANT STRESS STATES, TO
C     NODES. IT EXCLUSIVELY DEPENDS ON THE ELEMENT'S
C     GEOMETRY.
C
C     ARGUMENT INPUT:
C      XY     - THE ARRAY CONTAINNING THE LOCAL COORDINATE
C               VALUES OF AN ELEMENT'S NODES
C
C     ARGUMENT OUTPUT:
C      LT     - THE ARRAY
C
C PROGRAMMED BY: XIUXI WANG
C DATE/VERSION : 81-08-17/1.0
C
C LAST UPDATE  : 86-04-30/1.1 BY: EGIL GIERTSEN - SINTEF DIV.71
C ======================================================================
C
C                                                - EXTERNAL VARIABLES
      DIMENSION FLT(3,12), DELTA(2,4), FXY(2,4)
C
C-----------------------------------------------------------
C     CALCULATE THE DELTA MATRIX OF DIFFERENCES FOR THE
C     COORDINATE VALUES
C-----------------------------------------------------------
C
      DELTA(1,1) = (FXY(1,4)-FXY(1,2))/2.0
      DELTA(2,1) = (FXY(2,4)-FXY(2,2))/2.0
C
      DO  30 I   = 2,3
      DELTA(1,I) = (FXY(1,I-1)-FXY(1,I+1))/2.0
      DELTA(2,I) = (FXY(2,I-1)-FXY(2,I+1))/2.0
   30 CONTINUE
C
      DELTA(1,4) = (FXY(1,3)-FXY(1,1))/2.0
      DELTA(2,4) = (FXY(2,3)-FXY(2,1))/2.0
C
C----------------------------------------------------------
C     FORM THE TRANSPOSE OF THE MATRIX L
C----------------------------------------------------------
C
      CALL DCOPY (3*12,0.0D0,0,FLT,1)
C
      DO  50 I     = 1,4
      FLT(2,3*I-1) = DELTA(1,I)
      FLT(3,3*I-1) =-DELTA(2,I)
      FLT(1,3*I)   = DELTA(2,I)
      FLT(3,3*I)   =-DELTA(1,I)
   50 CONTINUE
C ------------------
C --- ERROR EXIT ---
C ------------------
C -------------------
C --- R E T U R N ---
C -------------------
      RETURN
      END
C
C
      SUBROUTINE CENTPO (AX,AY,FXY)
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
C =======================-----------------------------------------------
C === FEDEM                                               MODULE: CENTPO
C =======================-----------------------------------------------
C
C    PURPOSE/
C    METHOD :
C
C    THE SUBROUTINE CENTPO CALCULATES THE CENTROID POSISION
C    OF A QUADRILATERAL ELEMENT AND THE CO-ORDINATE VALUES OF
C    THE NODAL POINTS REFERRED TO THE CENTROID. THE AXES OF THE
C    LOCAL COORDINATE SYSTEM ARE TAKEN PARALLEL TO THE GLOBAL ONES
C
C    ARGUMENT INPUT:
C      AX      - GLOBAL X-COORDINATE VALUES OF THE ELEMENT'S NODES
C      AY      - GLOBAL Y-COORDINATE VALUES OF THE ELEMENT'S NODES
C
C    ARGUMENT OUTPUT:
C      FXY     - ARRAY CONTAINING THE LOCAL COORDINATE
C                VALUES OF THE ELEMENT'S NODES
C
C PROGRAMMED BY: XIUXI WANG
C DATE/VERSION : 82-04-06/1.0
C
C LAST UPDATE  : 86-04-30/1.1 BY: EGIL GIERTSEN - SINTEF DIV.71
C Total rev    : Feb 98  /1.2 Karl Erik Thoresen
C ======================================================================
C
C                                                - EXTERNAL VARIABLES
      DIMENSION AX(4),AY(4),FXY(2,4),XYL(2,4)
C
C------------------------------------------------------------
C     TRANSFORM THE NODAL POINTS FROM GLOBAL COORDINATE
C     SYSTEM TO LOCAL COORDINATE SYSTEM
C-------------------------------------------------------------
C
      DO  10 I = 1,4
      FXY(1,I) = AX(I)-AX(1)
      FXY(2,I) = AY(I)-AY(1)
   10 CONTINUE
C
      DLSQ     = FXY(1,2)**2+FXY(2,2)**2
      DL       = SQRT(DLSQ)
      COSINA   = FXY(1,2)/DL
      SINA     = FXY(2,2)/DL
C
      DO  15 I = 1,4
      XYL(1,I) =  FXY(1,I)*COSINA+FXY(2,I)*SINA
      XYL(2,I) = -FXY(1,I)*SINA+FXY(2,I)*COSINA
   15 CONTINUE
C
      call findCenterofquad(Xyl,Xc,Yc)
C
C     Move the co-ordinate sentre to the centroid
      DO  20 I = 1,4
      FXY(1,I) = FXY(1,I)-XC
      FXY(2,I) = FXY(2,I)-YC
   20 CONTINUE
C
C -------------------
C --- R E T U R N ---
C -------------------
      RETURN
      END
C
      subroutine findCenterofquad(P,Xc,Yc)
C=======================================================================
C
C     Finds the centroid of a quad element by summing
C     the centre of two triangles.
C     The first point P(:,1) must be in origo, the second on the x-axes.
C
C     Programmer: Karl Erik Thoresen Februar 1999
C
C=======================================================================
C
      DOUBLE PRECISION P(2,4), Xc, Yc
      DOUBLE PRECISION A1,A2,A,b,c,h,topX,topY,xc1,yc1,xcl,ycl,xc2,yc2
C
C     triangle 1
      h  = P(2,3)
      b  = P(1,2)
      c  = b - p(1,3)
      A1  = abs(b*h/2)
      Xc1 = (2*b-c)/3
      Yc1 = h/3
C
C     triangle 2
      b    = sqrt(P(1,3)*P(1,3) + P(2,3)*P(2,3))
      if ( b < tiny(b)*10.0D0 ) then
         A2  = 0.0D0
         Xc2 = 0.0D0
         Yc2 = 0.0D0
      else
         topX = (P(1,3)*P(1,4) + P(2,4)*P(2,4))/b
         topY = (P(1,3)*P(2,4) - P(2,3)*P(1,4))/b
         h    = topY
         c    = b - topX
         A2  = abs(b*h/2)
         Xcl = (2*b-c)/3                        ! in local co-ordinates
         Ycl = h/3                              ! in local co-ordinates
         Xc2 = (Xcl * P(1,3) -  Ycl * P(2,3))/b ! in global co-ordinates
         Yc2 = (Xcl * P(2,3) +  Ycl * P(1,3))/b ! in global co-ordinates
      end if
C
C     Sum the two triangles
      A = A1 + A2
      if ( A < tiny(A)*10.0D0 ) then
         Xc = 0.0D0
         Yc = 0.0D0
      else
         Xc = (Xc1 * A1 + Xc2 * A2)/ A
         Yc = (Yc1 * A1 + Yc2 * A2)/ A
      end if
C
      return
      end
C
C
      SUBROUTINE WRIME(A,M,N,LUNPRF)
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
C =======================-----------------------------------------------
C === FEDEM                                                MODUL: WRIME
C =======================-----------------------------------------------
C
C HENSIKT/
C METODE  : SKRIV UT EN DESIMALTALLSMATRISE I E-FORMAT
C
C DATA INN: A(M,N) - MATRISEN SOM SKAL SKRIVES UT
C           M      - ANTALL LINJER I -A(M,N)-
C           N      - ANTALL KOLONNER I -A(M,N)-
C           LUNPRF - ENHETSNUMMER FOR UTSKRIFT
C
C DATA UT : INGEN
C
C KALL    : INGEN
C
C BEGRENSNINGER : INGEN
C
C MASKIN-
C AVHENGIGHETER : INGEN
C
C PROGRAMMERT AV: EGIL GIERTSEN - SINTEF AVD.71
C DATO/VERSJON  : 86-05-08/1.0
C
C SIST OPPDATERT:
C ======================================================================
C
C                                                - EKSTERNE VARIABLE
      INTEGER            M, N, LUNPRF
      DIMENSION          A(M,N)
C                                                - INTERNE VARIABLE
      INTEGER   LTEST1, LTEST2, I, J, K, L
C
C  INITIALISER
C
      LTEST2 = 50
      LTEST1 =  8
C
C  SKRIV MATRISEN -A(M,N)- MED MAKSIMALT 9 TALL
C  P] HVER LINJE
C                                                - L\KKE OVER KOLONNER
      DO 30 I = 1,N,9
C
                    K = I + 8
      IF (K .GT. N) K = N
C                                                - SKRIV OVERSIKT
      WRITE(LUNPRF,6010) (J, J = I,K)
C                                                - L\KKE OVER LINJER
      DO 20 L = 1,M
C
      IF (L .NE. LTEST1)               GO TO 20
      IF (L .EQ. LTEST2)               GO TO 10
C                                                - NY LINJE
      WRITE(LUNPRF,6040)
C
   10 LTEST1 = LTEST1 + 7
C
      IF (L .NE. LTEST2)               GO TO 20
C                                                - SIDESKIFT
      WRITE(LUNPRF,6030)
C
      LTEST2 = LTEST2 + 49
C
C  SKRIV 'DENNE' DELEN AV MATRISEN
C
   20 WRITE(LUNPRF,6020) L,(A(L,J), J = I,K)
C
   30 CONTINUE
C -------------------
C --- R E T U R N ---
C -------------------
      RETURN
C ----------------
C --- FORMATER ---
C ----------------
 6010 FORMAT(/9(10X,I4))
 6020 FORMAT(I5,9E14.7)
 6030 FORMAT(1H1)
 6040 FORMAT(1H )
      END
C
C
      SUBROUTINE WRIMF2(A,M,N,LUNPRF)
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
C =======================-----------------------------------------------
C === FEDEM                                                MODUL: WRIMF2
C =======================-----------------------------------------------
C
C HENSIKT/
C METODE  : SKRIV UT EN MATRISE I F-FORMAT MED 3 DESIMALER
C
C DATA INN: A(M,N) - MATRISE SOM SKAL SKRIVES UT
C           M      - ANTALL LINJER I -A-
C           N      - ANTALL KOLONNER I -A-
C           LUNPRF - ENHETSNUMMER FOR UTSKRIFT
C
C DATA UT : INGEN
C
C KALL    : INGEN
C
C BEGRENSNINGER : INGEN
C
C MASKIN-
C AVHENGIGHETER : INGEN
C
C PROGRAMMERT AV: EGIL GIERTSEN - SINTEF AVD.71
C DATO/VERSJON  : 86-05-08/1.0
C
C SIST OPPDATERT:
C ======================================================================
C
C                                                - EKSTERNE VARIABLE
      INTEGER   M, N, LUNPRF
      DIMENSION A(M,N)
C                                                - INTERNE VARIABLE
      INTEGER   I, J, K, L
C
C  SKRIV UT MATRISEN A(M,N) I F-FORMAT
C
      DO 20 I = 1,N,9
C
                    K = I + 8
      IF (K .GT. N) K = N
C                                                - OVERSKRIFT
      WRITE(LUNPRF,6010) (J, J = I,K)
C                                                - SKRIV DENNE LINJER
      DO 10 L = 1,M
   10 WRITE(LUNPRF,6020) L,(A(L,J), J = I,K)
C
   20 CONTINUE
C ------------------
C --- FEILUTGANG ---
C ------------------
C -------------------
C --- R E T U R N ---
C -------------------
      RETURN
C ----------------
C --- FORMATER ---
C ----------------
 6010 FORMAT(/9(10X,I4))
 6020 FORMAT(I5,9F14.3)
      END
C
C
      SUBROUTINE QBEK35(EM,X,Y,Z,THK,RHO)
C
C***********************************************************************
C
C     QBEK35 DETERMINES A LUMPED MASS MATRIX
C     FOR THE QUADRILATERAL PLATE BENDING ELEMENT FQS.
C
C     PROGRAMMED BY :  KETIL AAMNES
C     DATE/VERSION  :  081085 / 1.0
C
C***********************************************************************
C
      DOUBLE PRECISION  EM(24),EMP(2,9),X(4),Y(4),Z(4),THK,RHO
      DOUBLE PRECISION  AUX1,COSG,SING,SL21,SL31
      DOUBLE PRECISION  X1,X2,X3,Y1,Y2,Y3,Z1,Z2,Z3,XL(4),YL(4)
      DOUBLE PRECISION  AREA,XTP,YTP,H,B,B1,B2,S1,S2
      DOUBLE PRECISION  IXX,IYY1,IYY2,IYY,IZZ1,IZZ2,IZZ
      DOUBLE PRECISION  IX,IY,IZ,M,GAMMA,TR(3,3)
C
C ---------------------------- Local coordinates -------------------
C
      DO 10 I=1,2
C
      X1 = 0.
      Y1 = 0.
      Z1 = 0.
C
      X2 = X(I+1)-X(1)
      Y2 = Y(I+1)-Y(1)
      Z2 = Z(I+1)-Z(1)
C
      X3 = X(I+2)-X(1)
      Y3 = Y(I+2)-Y(1)
      Z3 = Z(I+2)-Z(1)
C
      SL21 = DSQRT(X2*X2 + Y2*Y2 + Z2*Z2)
      SL31 = DSQRT(X3*X3 + Y3*Y3 + Z3*Z3)
C
      COSG = (X3*X2 + Y3*Y2 + Z3*Z2)/(SL31*SL21)
      AUX1 = 1.-COSG*COSG
C
      IF (AUX1) 40,40,50
   40 SING = 0.
      GO TO 60
   50 SING = DSQRT(AUX1)
C
   60 XL(1) = 0.
      XL(2) = SL21
      XL(3) = SL31*COSG
      YL(1) = 0.
      YL(2) = 0.
      YL(3) = SL31*SING
C
C ---------------------------- Element area ------------------------
C
      AREA = (XL(2)*YL(3))/2.
C
      IF (AREA) 80,80,70
C
C ---------------------------- Koordinates for TP ------------------
C
   70 XTP = ((XL(3)+((XL(2)-XL(3))/2.)))*(2./3.)
      YTP = (YL(3))/3.
C
C ------------------------------------------------------------------
C
      H = YL(3)
      B = XL(2)
C
      B1 = XL(3)
      B2 = XL(2)-XL(3)
C
      S1 = DSQRT( (((2./3.)*B1)**2) + (((1./3.)*H)**2) )
      S2 = DSQRT( (((2./3.)*B2)**2) + (((1./3.)*H)**2) )
C
C ---------------------------- Moment of inertia about the X-aksis -
C
      IXX = (RHO*THK*B*H) * (((H**2)/36.) + ((THK**2)/24.))
C
C ---------------------------- Moment of inertia about the Y-aksis -
C
      IYY1 = (RHO*THK*B1*H) * (((B1**2)/36.) + ((THK**2)/24.))
      IYY1 = IYY1 + (((XTP-((2./3.)*B1))**2)*(((B1*H)/2.)*RHO*THK))
C
      IYY2 = (RHO*THK*B2*H) * (((B2**2)/36.) + ((THK**2)/24.))
      IYY2 = IYY2 + (((B1-XTP+((1./3.)*B2))**2)*(((B2*H)/2.)*RHO*THK))
C
      IYY = IYY1 + IYY2
C
C ---------------------------- Moment of inertia about the Z-aksis -
C
      IZZ1 = (RHO*THK*H*B1) * (((B1**2)/4.)+((H**2)/12.))
      IZZ1 = IZZ1 - ((S1**2)*RHO*THK*(B1*H)/2.)
      IZZ1 = IZZ1 + (((XTP-((2./3.)*B1))**2)*((B1*H)/2.)*RHO*THK)
C
      IZZ2 = (RHO*THK*H*B2) * (((B2**2)/4.)+((H**2)/12.))
      IZZ2 = IZZ2 - ((S2**2)*RHO*THK*(B2*H)/2.)
      IZZ2 = IZZ2 + (((B1+((1./3.)*B2)-XTP)**2)*((B2*H)/2.)*RHO*THK)
C
      IZZ = IZZ1 + IZZ2
C
C ---------------------------- The 3 X 3 transformation matrix -----
C
C ---------------------------- V1 ----------------------------------
      TR(1,1) = X2-X1
      TR(2,1) = Y2-Y1
      TR(3,1) = Z2-Z1
      S1 = DSQRT(TR(1,1)**2+TR(2,1)**2+TR(3,1)**2)
      TR(1,1) = TR(1,1)/S1
      TR(2,1) = TR(2,1)/S1
      TR(3,1) = TR(3,1)/S1
C ---------------------------- V3 ----------------------------------
      TR(1,2) = X3-X1
      TR(2,2) = Y3-Y1
      TR(3,2) = Z3-Z1
      TR(1,3) = (TR(2,1)*TR(3,2)-TR(2,2)*TR(3,1))
      TR(2,3) = (TR(1,2)*TR(3,1)-TR(1,1)*TR(3,2))
      TR(3,3) = (TR(1,1)*TR(2,2)-TR(1,2)*TR(2,1))
C
      S1 = DSQRT(TR(1,3)**2+TR(2,3)**2+TR(3,3)**2)
      TR(1,3) = TR(1,3)/S1
      TR(2,3) = TR(2,3)/S1
      TR(3,3) = TR(3,3)/S1
C ---------------------------- V2 ----------------------------------
      TR(1,2) = TR(2,3)*TR(3,1)-TR(2,1)*TR(3,3)
      TR(2,2) = TR(1,1)*TR(3,3)-TR(1,3)*TR(3,1)
      TR(3,2) = TR(1,3)*TR(2,1)-TR(1,1)*TR(2,3)
C
C ---------------------------- Transform to global
C                              moments of inertia ------------------
      IX = ABS(TR(1,1)*IXX+TR(1,2)*IYY+TR(1,3)*IZZ)
      IY = ABS(TR(2,1)*IXX+TR(2,2)*IYY+TR(2,3)*IZZ)
      IZ = ABS(TR(3,1)*IXX+TR(3,2)*IYY+TR(3,3)*IZZ)
C
C ---------------------------- The elements of EM ------------------
C
      M = (RHO*AREA*THK)/3.
      GAMMA = 1./20.
C
      EMP(I,1) = M
      EMP(I,4) = M
      EMP(I,7) = M
C
      EMP(I,2) = GAMMA*IX
      EMP(I,3) = GAMMA*IY
      EMP(I,5) = EMP(I,2)
      EMP(I,6) = EMP(I,3)
      EMP(I,8) = EMP(I,2)
      EMP(I,9) = EMP(I,3)
C
   10 CONTINUE
C
      EM( 3) = EMP(1,1) + EMP(2,1)
      EM( 4) = EMP(1,2) + EMP(2,2)
      EM( 5) = EMP(1,3) + EMP(2,3)
C
      EM( 9) = EMP(1,4)
      EM(10) = EMP(1,5)
      EM(11) = EMP(1,6)
C
      EM(15) = EMP(1,7) + EMP(2,4)
      EM(16) = EMP(1,8) + EMP(2,5)
      EM(17) = EMP(1,9) + EMP(2,6)
C
      EM(21) = EMP(2,7)
      EM(22) = EMP(2,8)
      EM(23) = EMP(2,9)
C
C --- Fordeler masser og treghetsmomenter likt pa de 4 knutepunktene
C
      EM(3) = (EM(3) + EM(9)  + EM(15) + EM(21))/4.0
      EM(4) = (EM(4) + EM(10) + EM(16) + EM(22))/4.0
      EM(5) = (EM(5) + EM(11) + EM(17) + EM(23))/4.0
C
      EM(9)  = EM(3)
      EM(15) = EM(3)
      EM(21) = EM(3)
C
      EM(10) = EM(4)
      EM(16) = EM(4)
      EM(22) = EM(4)
C
      EM(11) = EM(5)
      EM(17) = EM(5)
      EM(23) = EM(5)
C
   80 RETURN
      END
