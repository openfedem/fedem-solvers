C     SPDX-FileCopyrightText: 2023 SAP SE
C
C     SPDX-License-Identifier: Apache-2.0
C
C     This file is part of FEDEM - https://openfedem.org
C
C>    @file beam.f
C>    @brief Subroutines for calculation of beam element matrices.
C
C>    @brief Generates the stiffness matrix for a 12-dof beam element.
      SUBROUTINE BEAM31(EK,X,Y,Z,EP,CA,XS,EFFLEN,PHI,
     +                  IPINA,IPINB,ID,IW,IPSW,IERR)
C
C***********************************************************************
C
C   FINITE ELEMENT LIBRARY :  BEAM31                       SINTEF / NTH
C
C
C     BEAM31 GENERATES THE STIFFNESS MATRIX FOR A 12-DEGREE-OF-FREEDOM
C     BEAM ELEMENT IN A GLOBAL SPACE SYSTEM. THE NODAL POINTS OF THE
C     ELEMENT ARE ECCENTRIC WITH RESPECT TO THE CENTER OF GRAVITY OF THE
C     CROSS SECTION. THE ELEMENT IS STRAIGHT AND THE CROSS SECTION IS
C     CONSTANT. AXIAL DEFORMATIONS, BENDING AND SHEAR ABOUT THE TWO
C     PRINCIPAL AXES AND ST.VENANT TORSION ARE CONSIDERED.
C
C     PROGRAMMED BY :  K.BELL AND F.STENSRUD
C
C     DATE/VERSION  :  23.02.77 / 01
C
C***********************************************************************
C
      IMPLICIT NONE
C
      DOUBLE PRECISION  EK(12,12),X(5),Y(5),Z(5),EP(6),CA(2),XS(2)
      DOUBLE PRECISION  EFFLEN,PHI
      INTEGER           IPINA,IPINB,ID,IW,IPSW,IERR
      DOUBLE PRECISION  BA,BL,BX,BY,BZ,T2(3,3),E1(3),E2(3),W(3)
C
      IERR = 0
C
      IF(IPSW-2) 20,10,10
   10 WRITE(IW,6000) ID
      WRITE(IW,6041) X,Y,Z
      WRITE(IW,6042) EP(1),EP(2),EP(3)
      WRITE(IW,6043) IPINA,IPINB
      WRITE(IW,6044) EFFLEN
      WRITE(IW,6045) PHI
C
C     THE LENGTH OF THE ELEMENT
C
   20 IF(EFFLEN) 22,22,21
   21 BL = EFFLEN
      GO TO 25

   22 BX = X(2) - X(1)
      BY = Y(2) - Y(1)
      BZ = Z(2) - Z(1)
      BL = SQRT(BX*BX + BY*BY + BZ*BZ)
C
C     DEGENERATED ELEMENT ?
C
   25 BA = ABS(X(1))+ABS(X(2))+ABS(Y(1))+ABS(Y(2))+ABS(Z(1))+ABS(Z(2))
      BA = 1.0D-6*BA
      IF(BL-BA) 1010,1010,30
C
C     CHECK THE LOCAL Z-AXIS ALSO
C
   30 BX = X(3) - X(1)
      BY = Y(3) - Y(1)
      BZ = Z(3) - Z(1)
      BZ = SQRT(BX*BX + BY*BY + BZ*BZ)
      IF(BZ-BA) 1020,1020,40
C
C     CHECK THE PROPERTY PARAMETERS
C
   40 if (EP(2) .le. 1.0D-16 .or. EP(3) .le. 1.0D-16) go to 1030
C
C     SET UP GLOBAL-TO-LOCAL TRANSFORMATION MATRIX
C
      CALL DCOS30(T2,X,Y,Z,PHI)
C     Note: This will print the transpose of T2,
C     i.e., the local-to-global transformation matrix
      if (ipsw .gt. 2) write(iw,6046) T2
C
C     STIFFNESS MATRIX IN LOCAL COORDINATES
C
      CALL BELS31(EK,BL,EP,CA,XS,ID,IW,IPSW)
C
C     TRANSFORM TO GLOBAL COORDINATES
C
      if (iPinA .le. 0 .and. iPinB .le. 0) then
C        No pin flags, transform the whole matrix
         CALL MPRO30(EK,T2,12)
      else if (iPinA .le. 0) then
C        Pin flag at end B only, transform the matrix at end A only
         CALL MATTRA(T2,EK,W,12,3,1,IW,IERR)
         CALL MATTRA(T2,EK,W,12,3,4,IW,IERR)
      else if (iPinB .le. 0) then
C        Pin flag at end A only, transform the matrix at end B only
         CALL MATTRA(T2,EK,W,12,3,7,IW,IERR)
         CALL MATTRA(T2,EK,W,12,3,10,IW,IERR)
      end if
C
C     TRANSFORM FROM END POINTS TO NODAL POINTS
C
      E1(1) = X(4)-X(1)
      E1(2) = Y(4)-Y(1)
      E1(3) = Z(4)-Z(1)
      E2(1) = X(5)-X(2)
      E2(2) = Y(5)-Y(2)
      E2(3) = Z(5)-Z(2)
      if (iPinA .gt. 0) CALL VECTRA(T2,E1,W,3,3,1,1,IW,IERR)
      if (iPinB .gt. 0) CALL VECTRA(T2,E2,W,3,3,1,1,IW,IERR)
      if (ipsw  .gt. 2) write(iw,6047) E1,E2
      CALL TRIX30(EK,E1(1),E1(2),E1(3),E2(1),E2(2),E2(3))
C
      GO TO 2000
C
 1010 IERR = -1
      GO TO 1900
 1020 IERR = -2
      GO TO 1900
 1030 IERR = -3
 1900 WRITE(IW,6020)
      WRITE(IW,6030) IERR
      if (IERR .eq. -1) then
         WRITE(IW,6040) ID,'DEGENERATED ELEMENT'
      else if (IERR .eq. -2) then
         WRITE(IW,6040) ID,'ZERO-LENGTH Z-DIRECTION VECTOR'
      else if (IERR .eq. -3) then
         WRITE(IW,6040) ID,'INVALID PROPERTY PARAMETERS'
      end if
      WRITE(IW,6041) X,Y,Z
      WRITE(IW,6042) EP(1),EP(2),EP(3)
C
C     PRINT STIFFNESS MATRIX ?
C
 2000 IF(IPSW-2) 3000,2040,2010
 2010 IF(IERR)   2040,2020,2020
 2020 IF(IPSW-4) 2040,2030,2030
 2030 WRITE(IW,6050)
      CALL MPRT30(EK,12,12,IW)
 2040 WRITE(IW,6060) ID
C
 3000 RETURN
C
C
 6000 FORMAT(///29H ENTERING BEAM31 FOR ELEMENT ,I7,2H :)
 6020 FORMAT(///33H *** ERROR RETURN FROM BEAM31 ***)
 6030 FORMAT(17H     ERROR FLAG =,I3)
 6040 FORMAT(/13H     ELEMENT ,I7,11H   IN ERROR/5X,A)
 6041 FORMAT(9H     XG =,5E15.5/9H     YG =,5E15.5/9H     ZG =,5E15.5)
 6042 FORMAT(9H     E  =, E15.5/9H     G  =, E15.5/9H     A  =, E15.5)
 6043 FORMAT(9H     PA =, I7   /9H     PB =, I7)
 6044 FORMAT(9H     EL =, E15.5)
 6045 FORMAT(10H     PHI =, F10.5)
 6046 FORMAT(43H     LOCAL TO GLOBAL TRANSFORMATION MATRIX:/(4X,3F10.6))
 6047 FORMAT(9H     E1 =, 1P3E13.5 /9H     E2 =, 1P3E13.5)
 6050 FORMAT(//48H     GLOBAL STIFFNESS MATRIX (ECCENTRIC NODES) :)
 6060 FORMAT(/28H LEAVING BEAM31 FOR ELEMENT ,I7)
C
      END
C> @brief Computes stress resultants on a 12-dof beam element.
      SUBROUTINE BEAM32(SR,EK,V,FP,FE,P,X,Y,Z,XS,XP,PHI,ZETA)
C
C***********************************************************************
C
C   FINITE ELEMENT LIBRARY :  BEAM32                       SINTEF / NTH
C
C
C     BEAM32 COMPUTES THE STRESS RESULTANTS AT AN ARBITRARY SECTION
C     OF A 12-DEGREE-OF-FREEDOM BEAM ELEMENT IN A GLOBAL SPACE SYSTEM.
C     THE NODES OF THE ELEMENT ARE ECCENTRIC WITH RESPECT TO THE CENTER
C     OF GRAVITY OF THE CROSS SECTION. THE STIFFNESS MATRIX AND THE
C     NODAL DISPLACEMENTS, THE LOAD VECTOR DUE TO DISTRIBUTED EXTERNAL
C     LOADING AND NODAL FORCES DUE TO INITIAL STRAIN, ARE REQUIRED
C     AS INPUT
C
C     PROGRAMMED BY :  F.STENSRUD AND K.BELL
C
C     DATE/VERSION  :  24.02.77 / 01
C
C***********************************************************************
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION SR(6),EK(12,12),V(12),FP(12),FE(12),P(6),X(5),Y(5),Z(5),
     +          XS(2),XP(2),F(6),FL(6),T2(3,3)
C
      PX1 = P(1)
      PX2 = P(2)
      PY1 = P(3)
      PY2 = P(4)
      PZ1 = P(5)
      PZ2 = P(6)
      YS  = XS(1)
      ZS  = XS(2)
      YP  = XP(1)
      ZP  = XP(2)
      EX1 = X(4)-X(1)
      EY1 = Y(4)-Y(1)
      EZ1 = Z(4)-Z(1)
C
C     CONTROL OF ZETA
C
      IF(ZETA) 1000,10,10
   10 IF(ZETA-1.) 20,20,1000
C
C     THE LENGTH OF THE ELEMENT
C
   20 BX = X(2) - X(1)
      BY = Y(2) - Y(1)
      BZ = Z(2) - Z(1)
      BL = SQRT(BX*BX + BY*BY + BZ*BZ)
C
C     NODAL FORCES IN GLOBAL COORDINATES
C
      DO 40 I=1,6
      FL(I) = 0.
      DO 30 J=1,12
      FL(I) = FL(I) + EK(I,J)*V(J)
   30 CONTINUE
      FL(I) = FL(I) - FP(I) - FE(I)
   40 CONTINUE
C
C     TRANSFORM FROM ECCENTRIC NODES TO ELEMENT ENDS
C
      FL(4) = FL(4) - EZ1*FL(2) + EY1*FL(3)
      FL(5) = FL(5) + EZ1*FL(1) - EX1*FL(3)
      FL(6) = FL(6) - EY1*FL(1) + EX1*FL(2)
C
C     TRANSFORM TO LOCAL COORDINATES
C
      CALL DCOS30(T2,X,Y,Z,PHI)
C
      DO 70 I=1,2
      N = 3*I
      DO 60 J=1,3
      F(J+N-3) = T2(J,1)*FL(N-2) + T2(J,2)*FL(N-1) + T2(J,3)*FL(N)
   60 CONTINUE
   70 CONTINUE
C
C     STRESS RESULTANTS AT NODE 1
C
      SN  = - F(1)
      SQY =   F(2)
      SQZ =   F(3)
      SMT = - F(4) + YS*F(3) - ZS*F(2)
      SMY =   F(5)
      SMZ = - F(6)
C
C     LOAD ECCENTRICITY FROM SHEAR CENTER AND SOME AUXILIARY CONSTANTS
C
      YPS = YP - YS
      ZPS = ZP - ZS
      Z2 = .5*ZETA*ZETA
      Z3 = ZETA*ZETA*ZETA/6.
C
C     STRESS RESULTANTS AT AN ARBITRARY SECTION
C
      SR(1) = SN  - BL*(ZETA*PX1+Z2*(PX2-PX1))
      SR(2) = SQY + BL*(ZETA*PY1+Z2*(PY2-PY1))
      SR(3) = SQZ + BL*(ZETA*PZ1+Z2*(PZ2-PZ1))
      SR(4) = SMT + ZPS*BL*(ZETA*PY1+Z2*(PY2-PY1))
     +            - YPS*BL*(ZETA*PZ1+Z2*(PZ2-PZ1))
      SR(5) = SMY + SQZ*ZETA*BL + BL*BL*(Z2*PZ1+Z3*(PZ2-PZ1))
      SR(6) = SMZ + SQY*ZETA*BL + BL*BL*(Z2*PY1+Z3*(PY2-PY1))
C
      GO TO 2000
C
 1000 DO 1010 I=1,6
 1010 SR(I) = 0.
C
C
 2000 RETURN
C
      END
C> @brief Computes consistent load vector due to line loads on a beam.
      SUBROUTINE BEAM33(F,P,X,Y,Z,EP,CA,XS,XP,PHI)
C
C***********************************************************************
C
C   FINITE ELEMENT LIBRARY :  BEAM33                       SINTEF / NTH
C
C     BEAM33 COMPUTES A KINEMATICALLY CONSISTENT LOAD VECTOR DUE TO A
C     LINEARLY VARYING DISTRIBUTED NORMAL AND/OR TANGENTIAL LOAD ON THE
C     12-DEGREE-OF-FREEDOM  B E A M  ELEMENT
C     THE LOAD COMPONENTS NORMAL TO THE AXIS MAY BE ECCENTRIC
C
C     PROGRAMMED BY :  K.BELL
C     DATE/VERSION  :  78-10-07 / 01
C
C***********************************************************************
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION  F(12),P(6),X(5),Y(5),Z(5),EP(6),CA(2),XS(2),XP(2),
     +           FP(12),T2(3,3)
C                                               ** ELEMENT LENGTH
      BX = X(2) - X(1)
      BY = Y(2) - Y(1)
      BZ = Z(2) - Z(1)
      BL = SQRT(BX*BX + BY*BY + BZ*BZ)
C                                               ** LOCAL LOAD VECTOR
      CALL BELS33(FP,P,BL,EP,CA,XS,XP)
C                                               ** TRANSFORM TO GLOBAL
C                                                  AXES
      CALL DCOS30(T2,X,Y,Z,PHI)
C
      DO 30 I=1,4
      N = 3*I
      DO 20 J=1,3
      F(J+N-3) = T2(1,J)*FP(N-2) + T2(2,J)*FP(N-1) + T2(3,J)*FP(N)
   20 CONTINUE
   30 CONTINUE
C                                               ** TRANSFORM FROM END TO
C                                                  NODAL POINTS
      EX1 = X(4)-X(1)
      EY1 = Y(4)-Y(1)
      EZ1 = Z(4)-Z(1)
      EX2 = X(5)-X(2)
      EY2 = Y(5)-Y(2)
      EZ2 = Z(5)-Z(2)
      F(4) = F(4) + EZ1*F(2) - EY1*F(3)
      F(5) = F(5) - EZ1*F(1) + EX1*F(3)
      F(6) = F(6) + EY1*F(1) - EX1*F(2)
      F(10)= F(10)+ EZ2*F(8) - EY2*F(9)
      F(11)= F(11)- EZ2*F(7) + EX2*F(9)
      F(12)= F(12)+ EY2*F(7) - EX2*F(8)
C
      RETURN
      END
C> @brief Computes consistent load vector due to initial strain on a beam.
      SUBROUTINE BEAM34(F,E0,X,Y,Z,EP)
C
C***********************************************************************
C
C   FINITE ELEMENT LIBRARY :   BEAM34                      SINTEF / NTH
C
C     BEAM34 COMPUTES A VECTOR OF NODAL FORCES DUE TO INITIAL STRAIN FOR
C     THE 12-DEGREE-OF-FREEDOM   B E A M   ELEMENT
C     THE INITIAL STRAIN MAY VARY LINEARLY OVER THE DEPTH AND BREADTH
C     OF THE BEAM, BUT IT MUST BE CONSTANT ALONG THE ELEMENT AXIS
C
C     PROGRAMMED BY :  K.BELL
C     DATE/VERSION  :  78-10-07 / 01
C
C***********************************************************************
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION  F(12),E0(3),X(5),Y(5),Z(5),EP(6),FE(12),T2(3,3)
C
C                                               ** LOCAL FORCE VECTOR
      CALL BELS34(FE,E0,EP)
C                                               ** TRANSFORM TO GLOBAL
C                                                  AXES
      CALL DCOS30(T2,X,Y,Z,0.0D0)
C
      DO 30 I=1,4
      N = 3*I
      DO 20 J=1,3
      F(J+N-3) = T2(1,J)*FE(N-2) + T2(2,J)*FE(N-1) + T2(3,J)*FE(N)
   20 CONTINUE
   30 CONTINUE
C                                               ** TRANSFORM FROM END TO
C                                                  NODAL POINTS
      EX1 = X(4)-X(1)
      EY1 = Y(4)-Y(1)
      EZ1 = Z(4)-Z(1)
      EX2 = X(5)-X(2)
      EY2 = Y(5)-Y(2)
      EZ2 = Z(5)-Z(2)
      F(4) = F(4) + EZ1*F(2) - EY1*F(3)
      F(5) = F(5) - EZ1*F(1) + EX1*F(3)
      F(6) = F(6) + EY1*F(1) - EX1*F(2)
      F(10)= F(10)+ EZ2*F(8) - EY2*F(9)
      F(11)= F(11)- EZ2*F(7) + EX2*F(9)
      F(12)= F(12)+ EY2*F(7) - EX2*F(8)
C
      RETURN
      END
C> @brief Generates a mass matrix for a 12-dof beam element.
      SUBROUTINE BEAM35(EM,X,Y,Z,A,RHO,RIX,IOP,IPINA,IPINB,IW,IPSW,IERR)
C
C***********************************************************************
C
C   FINITE ELEMENT LIBRARY :   BEAM35                      SINTEF / NTH
C
C     BEAM35 GENERATES  A LUMPED (IOP=1), A DIAGONALIZED (IOP=2) OR A
C     CONSISTENT (IOP=3)  MASS MATRIX FOR THE 12-DEGREE-OF-FREEDOM
C     B E A M   ELEMENT
C
C     PROGRAMMED BY :  K.BELL AND I.LANGEN
C     DATE/VERSION  :  78-10-08 / 01
C
C***********************************************************************
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION  EM(12,*),X(5),Y(5),Z(5),T2(3,3)
      DIMENSION  E1(3),E2(3),W(3)
C
      IF(IPSW-2) 12,10,10
   10 WRITE(IW,6000)
      WRITE(IW,6041) X,Y,Z
      WRITE(IW,6042) A,RHO,RIX
      WRITE(IW,6043) IOP,IPINA,IPINB
C
   12 K = 1
      IF (IOP.EQ.3)  K=12
C
      DO 20 I=1,12
      DO 20 J=1,K
      EM(I,J) = 0.
   20 CONTINUE
C                                               ** ELEMENT LENGTH
      BX = X(2)-X(1)
      BY = Y(2)-Y(1)
      BZ = Z(2)-Z(1)
      BL = SQRT(BX*BX + BY*BY + BZ*BZ)
C
C     MASS MATRIX IN LOCAL COORDINATES
C
      AM  = BL*A*RHO
      HM  = AM*0.5D0
      IF (IOP.EQ.1)  GO TO 40
      AM2 = BL*AM
      AM3 = BL*AM2
      RM  = RIX*RHO
      IF (IOP.EQ.3)  GO TO 60
C                                               ** LUMPED MASS MATRIX
   40 EM(1,1) = HM
      EM(2,1) = HM
      EM(3,1) = HM
      EM(7,1) = HM
      EM(8,1) = HM
      EM(9,1) = HM
      IF (IOP.EQ.1)  GO TO 100
C                                               ** DIAGONAL MASS MATRIX
      EM(4,1)  = BL*RM/2.
      EM(5,1)  = AM3/420.
      EM(6,1)  = EM(5,1)
      EM(10,1) = EM(4,1)
      EM(11,1) = EM(5,1)
      EM(12,1) = EM(6,1)
      GO TO 100
C                                               ** CONSIST. MASS MATRIX
   60 EM( 1, 1) = AM/3.
      EM( 1, 7) = AM/6.
      EM( 2, 2) = 13.*AM/35.
      EM( 2, 6) = 11.*AM2/210.
      EM( 2, 8) = 9.*AM/70.
      EM( 2,12) =-13.*AM2/420.
      EM( 3, 3) = EM(2,2)
      EM( 3, 5) =-EM(2,6)
      EM( 3, 9) = EM(2,8)
      EM( 3,11) =-EM(2,12)
      EM( 4, 4) = BL*RM/3.
      EM( 4,10) = BL*RM/6.
      EM( 5, 5) = AM3/105.
      EM( 5, 9) = EM(2,12)
      EM( 5,11) =-AM3/140.
      EM( 6, 6) = EM(5,5)
      EM( 6, 8) = EM(3,11)
      EM( 6,12) = EM(5,11)
      EM( 7, 7) = EM(1,1)
      EM( 8, 8) = EM(2,2)
      EM( 8,12) = EM(3,5)
      EM( 9, 9) = EM(8,8)
      EM( 9,11) = EM(2,6)
      EM(10,10) = EM(4,4)
      EM(11,11) = EM(5,5)
      EM(12,12) = EM(6,6)
C
      DO 80 I=1,12
      DO 80 J=I,12
      EM(J,I) = EM(I,J)
   80 CONTINUE
C
C     TRANSFORM FROM LOCAL TO GLOBAL COORDINATES
C
  100 IF (IOP.EQ.1)  GO TO 200
      CALL DCOS30(T2,X,Y,Z,0.0D0)
      IF (IOP.EQ.3)  GO TO 120
C                                               ** DIAGONAL MASS MATRIX
      X1 = EM(4,1)*T2(1,1)*T2(1,1) +
     +     EM(5,1)*T2(2,1)*T2(2,1) +
     +     EM(6,1)*T2(3,1)*T2(3,1)
      X2 = EM(4,1)*T2(1,2)*T2(1,2) +
     +     EM(5,1)*T2(2,2)*T2(2,2) +
     +     EM(6,1)*T2(3,2)*T2(3,2)
      X3 = EM(4,1)*T2(1,3)*T2(1,3) +
     +     EM(5,1)*T2(2,3)*T2(2,3) +
     +     EM(6,1)*T2(3,3)*T2(3,3)
      if (iPinA .le. 0) then
         EM(4,1)  = X1
         EM(5,1)  = X2
         EM(6,1)  = X3
      end if
      if (iPinB .le. 0) then
         EM(10,1) = X1
         EM(11,1) = X2
         EM(12,1) = X3
      end if
      GO TO 200
C                                               ** CONSIST. MASS MATRIX
  120 if (iPinA .le. 0 .and. iPinB .le. 0) then
C        No pin flags, transform the whole matrix
         CALL MPRO30(EM,T2,12)
      else if (iPinA .le. 0) then
C        Pin flag at end B only, transform the matrix at end A only
         CALL MATTRA(T2,EM,W,12,3,1,IW,IERR)
         CALL MATTRA(T2,EM,W,12,3,4,IW,IERR)
      else if (iPinB .le. 0) then
C        Pin flag at end A only, transform the matrix at end B only
         CALL MATTRA(T2,EM,W,12,3,7,IW,IERR)
         CALL MATTRA(T2,EM,W,12,3,10,IW,IERR)
      end if
C
C     TRANSFORM FROM END POINTS TO NODAL POINTS
C
  200 E1(1) = X(4)-X(1)
      E1(2) = Y(4)-Y(1)
      E1(3) = Z(4)-Z(1)
      E2(1) = X(5)-X(2)
      E2(2) = Y(5)-Y(2)
      E2(3) = Z(5)-Z(2)
      if (iPinA .gt. 0) CALL VECTRA(T2,E1,W,3,3,1,1,IW,IERR)
      if (iPinB .gt. 0) CALL VECTRA(T2,E2,W,3,3,1,1,IW,IERR)
      IF (IOP.EQ.3)  GO TO 220
C                                               ** LUMPED OR DIAGONAL
C                                                  MASS MATRIX
      EM(4,1)  = EM(4,1)  + HM*(E1(2)*E1(2) + E1(3)*E1(3))
      EM(5,1)  = EM(5,1)  + HM*(E1(1)*E1(1) + E1(3)*E1(3))
      EM(6,1)  = EM(6,1)  + HM*(E1(1)*E1(1) + E1(2)*E1(2))
      EM(10,1) = EM(10,1) + HM*(E2(2)*E2(2) + E2(3)*E2(3))
      EM(11,1) = EM(11,1) + HM*(E2(1)*E2(1) + E2(3)*E2(3))
      EM(12,1) = EM(12,1) + HM*(E2(1)*E2(1) + E2(2)*E2(2))
      GO TO 300
C                                               ** CONSIST. MASS MATRIX
  220 CALL TRIX30(EM,E1(1),E1(2),E1(3),E2(1),E2(2),E2(3))
C
C     PRINT MASS MATRIX ?
C
  300 IF(IPSW-2) 390,340,320
  320 IF(IPSW-4) 340,330,330
  330 WRITE(IW,6050)
      IF(IOP-2) 331,331,332
  331 CALL MPRT30(EM,12,1,IW)
      GO TO 340
  332 CALL MPRT30(EM,12,12,IW)
  340 WRITE(IW,6060)
C
  390 RETURN
C
 6000 FORMAT(///18H ENTERING BEAM35 :)
 6041 FORMAT(9H     XG =,5E15.5/9H     YG =,5E15.5/9H     ZG =,5E15.5)
 6042 FORMAT(9H     A  =, E15.5/9H     RO =, E15.5/9H     IX =, E15.5)
 6043 FORMAT(9H     OP =, I7   /9H     PA =, I7   /9H     PB =, I7)
 6050 FORMAT(//43H     GLOBAL MASS MATRIX (ECCENTRIC NODES) :)
 6060 FORMAT(/15H LEAVING BEAM35)
C
      END
C> @brief Generates a geometric stiffness matrix for a 12-dof beam element.
      SUBROUTINE BEAM36(EKG,X,Y,Z,EP,CA,XS,PHI,S1,S2)
C
C***********************************************************************
C
C   FINITE ELEMENT LIBRARY :   BEAM36                      SINTEF / NTH
C
C     BEAM36 GENERATES A GEOMETRICAL STIFFNESS MATRIX FOR THE 12-DEGREE-
C     OF-FREEDOM   B E A M   ELEMENT ON THE BASIS OF A LINEARLY VARYING
C     AXIAL FORCE (N=S1-X*(S1+S2)/L)
C
C     PROGRAMMED BY :  K.BELL
C     DATE/VERSION  :  78-10-08 / 01
C
C***********************************************************************
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION EKG(12,12),X(5),Y(5),Z(5),EP(6),CA(2),XS(2),T2(3,3)
C
C                                               ** ELEMENT LENGTH
      BX = X(2)-X(1)
      BY = Y(2)-Y(1)
      BZ = Z(2)-Z(1)
      BL = SQRT(BX*BX + BY*BY + BZ*BZ)
C                                               ** GEOMETRIC STIFFNESS
C                                                  IN LOCAL COORDINATES
      CALL BELS36(EKG,S1,S2,BL,EP,CA,XS)
C                                               ** TRANSFORM TO GLOBAL
C                                                  AXES
      CALL DCOS30(T2,X,Y,Z,PHI)
      CALL MPRO30(EKG,T2,12)
C                                               ** TRANSFORM FROM END TO
C                                                  NODAL POINTS
      EX1 = X(4)-X(1)
      EY1 = Y(4)-Y(1)
      EZ1 = Z(4)-Z(1)
      EX2 = X(5)-X(2)
      EY2 = Y(5)-Y(2)
      EZ2 = Z(5)-Z(2)
      CALL TRIX30(EKG,EX1,EY1,EZ1,EX2,EY2,EZ2)
C
      RETURN
      END
      SUBROUTINE BELS31(EK,BL,EP,CA,XS,ID,IW,IPSW)
C
C***********************************************************************
C
C   FINITE ELEMENT LIBRARY :  BELS31                       SINTEF / NTH
C
C
C     BELS31 GENERATES THE STIFFNESS MATRIX FOR A 12-DEGREE-OF-FREEDOM
C     BEAM ELEMENT IN A LOCAL SPACE SYSTEM. THE ELEMENT IS STRAIGHT
C     AND HAS A CONSTANT CROSS SECTION. AXIAL DEFORMATIONS, BENDING
C     AND SHEAR ABOUT THE TWO PRINCIPAL AXES AND ST.VENANT TORSION
C     ARE CONSIDERED.
C
C     PROGRAMMED BY : F.STENSRUD
C
C     DATE/VERSION  : 25.05.73 / 01
C
C***********************************************************************
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION EK(12,12),EP(6),CA(2),XS(2)
C
      E   = EP(1)
      G   = EP(2)
      A   = EP(3)
      RIY = EP(4)
      RIZ = EP(5)
      RIT = EP(6)
      CAY = CA(1)
      CAZ = CA(2)
      YS  = XS(1)
      ZS  = XS(2)
C
C     SOME AUXILIARY CONSTANTS
C
      if (E .lt. -1.0d-16) then
C        A negative Young's modulus (E) is used to flag
C        alternative specification of the structural
C        cross section properties, as follows:
C        Axial stiffness E*A is given in A
         EA  = A
C        Bending stiffnesses E*I are given in RIY and RIZ
         EIY = RIY
         EIZ = RIZ
C        Shear stiffnesses G*As are given in CAY and CAZ
         if (CAY .gt. 1.0d-16) then
            ALY = 12.0d0*EIY/(CAY*BL*BL)
         else
C           Infinite shear stiffness (no shear deformation)
            ALY = 0.0d0
         end if
         if (CAZ .gt. 1.0d-16) then
            ALZ = 12.0d0*EIZ/(CAZ*BL*BL)
         else
C           Infinite shear stiffness (no shear deformation)
            ALZ = 0.0d0
         end if
C        Torsional stiffness G*It is given as RIT
         GIT = RIT
      else
         EA  = E*A
         EIY = E*RIY
         EIZ = E*RIZ
         ALY = 12.0d0*CAY*EIY/(A*G*BL*BL)
         ALZ = 12.0d0*CAZ*EIZ/(A*G*BL*BL)
         GIT = G*RIT
      end if
C
C     PRINT ?
C
      IF(IPSW-2) 20,10,10
   10 WRITE(IW,6000) ID
      WRITE(IW,6010) BL
      WRITE(IW,6020) E
      WRITE(IW,6030) G
      WRITE(IW,6040) A
      WRITE(IW,6050) RIY
      WRITE(IW,6060) RIZ
      WRITE(IW,6070) RIT
      WRITE(IW,6080) CAY
      WRITE(IW,6090) CAZ
      WRITE(IW,6100) YS
      WRITE(IW,6110) ZS
      WRITE(IW,6120) EA
      WRITE(IW,6130) EIY
      WRITE(IW,6140) EIZ
      WRITE(IW,6150) GIT
      WRITE(IW,6160) ALY
      WRITE(IW,6170) ALZ
C
   20 DO 30 I=1,12
      DO 30 J=1,12
   30 EK(I,J) = 0.
C
C     UPPER TRIANGLE OF THE STIFFNESS MATRIX IN CASE OF DOUBLE SYMMETRY
C
      EK(1,1)   =   EA/BL
      EK(2,2)   =   12.*EIY/(BL*BL*BL*(1.+ALY))
      EK(3,3)   =   12.*EIZ/(BL*BL*BL*(1.+ALZ))
      EK(4,4)   =   GIT/BL
      EK(3,5)   = - .5*BL*EK(3,3)
      EK(2,6)   =   .5*BL*EK(2,2)
      EK(5,5)   =   EIZ*(4.+ALZ)/(BL*(1.+ALZ))
      EK(6,6)   =   EIY*(4.+ALY)/(BL*(1.+ALY))
      EK(1,7)   = - EK(1,1)
      EK(7,7)   =   EK(1,1)
      EK(2,8)   = - EK(2,2)
      EK(6,8)   = - EK(2,6)
      EK(8,8)   =   EK(2,2)
      EK(3,9)   = - EK(3,3)
      EK(5,9)   = - EK(3,5)
      EK(9,9)   =   EK(3,3)
      EK(4,10)  = - EK(4,4)
      EK(10,10) =   EK(4,4)
      EK(3,11)  =   EK(3,5)
      EK(5,11)  =   EIZ*(2.-ALZ)/(BL*(1.+ALZ))
      EK(9,11)  = - EK(3,5)
      EK(11,11) =   EK(5,5)
      EK(2,12)  =   EK(2,6)
      EK(6,12)  =   EIY*(2.-ALY)/(BL*(1.+ALY))
      EK(8,12)  = - EK(2,6)
      EK(12,12) =   EK(6,6)
C
C     LOWER TRIANGLE OG EK FROM SYMMETRY
C
      DO 40 I=1,12
      DO 40 J=1,I
   40 EK(I,J) = EK(J,I)
C
C     DOUBLE SYMMETRICAL CROSS SECTION ?
C
      BE = ABS(YS) + ABS(ZS)
      BA = 0.00001D0*SQRT(A)
      IF(BA-BE) 50,1000,1000
C
C     ADJUSTMENT OF EK DUE TO SINGLE OR NO SYMMETRY
C
   50 CONTINUE
      DO 60 I=1,12
      EK( 4,I) = EK( 4,I) - ZS*EK(2,I) + YS*EK(3,I)
   60 EK(10,I) = EK(10,I) - ZS*EK(8,I) + YS*EK(9,I)
C
      DO 70 I=1,12
      EK(I, 4) = EK(I, 4) - EK(I,2)*ZS + EK(I,3)*YS
   70 EK(I,10) = EK(I,10) - EK(I,8)*ZS + EK(I,9)*YS
C
C     PRINT ?
C
 1000 CONTINUE
      IF(IPSW-2) 3000,2010,2000
 2000 WRITE(IW,6180)
      CALL MPRT30(EK,12,12,IW)
 2010 WRITE(IW,6190) ID
C
 3000 RETURN
C
C
 6000 FORMAT(///29H ENTERING BELS31 FOR ELEMENT ,I7,2H :)
 6010 FORMAT(14H     L      = ,1PE12.5)
 6020 FORMAT(14H     E      = ,1PE12.5)
 6030 FORMAT(14H     G      = ,1PE12.5)
 6040 FORMAT(14H     A      = ,1PE12.5)
 6050 FORMAT(14H     IY     = ,1PE12.5)
 6060 FORMAT(14H     IZ     = ,1PE12.5)
 6070 FORMAT(14H     IT     = ,1PE12.5)
 6080 FORMAT(14H     KAPPAY = ,1PE12.5)
 6090 FORMAT(14H     KAPPAZ = ,1PE12.5)
 6100 FORMAT(14H     YSHEAR = ,1PE12.5)
 6110 FORMAT(14H     ZSHEAR = ,1PE12.5)
 6120 FORMAT(14H     EA     = ,1PE12.5)
 6130 FORMAT(14H     EIY    = ,1PE12.5)
 6140 FORMAT(14H     EIZ    = ,1PE12.5)
 6150 FORMAT(14H     GIT    = ,1PE12.5)
 6160 FORMAT(14H     ALFAY  = ,1PE12.5)
 6170 FORMAT(14H     ALFAZ  = ,1PE12.5)
 6180 FORMAT(//37H     LOCAL ELEMENT STIFFNESS MATRIX :)
 6190 FORMAT(/28H LEAVING BELS31 FOR ELEMENT ,I7)
      WRITE(IW,6120) EA
      WRITE(IW,6130) EIY
      WRITE(IW,6140) EIZ
      WRITE(IW,6150) GIT
      WRITE(IW,6160) ALY
      WRITE(IW,6170) ALZ
C
      END
      SUBROUTINE BELS32(SR,EK,V,FP,FE,P,BL,XS,XP,ZETA)
C
C***********************************************************************
C
C   FINITE ELEMENT LIBRARY :  BELS32                       SINTEF / NTH
C
C
C     BELS32 COMPUTES THE STRESS RESULTANTS AT AN ARBRITRARY SECTION OF
C     A 12-DEGREE-OF-FREEDOM BEAM ELEMENT IN A LOCAL SPACE SYSTEM. THE
C     STIFFNESS MATRIX AND THE NODAL DISPLACEMENTS, THE LOAD VECTOR DUE
C     TO DISTRIBUTED EXTERNAL LOADING AND NODAL FORCES DUE TO INITIAL
C     STRAIN ARE REQUIRED AS INPUT
C
C     PROGRAMMED BY : F.STENSRUD
C
C     DATE/VERSION  : 28.05.73 / 01
C
C***********************************************************************
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION SR(6),EK(12,12),V(12),FP(12),FE(12),P(6),
     +          XS(2),XP(2),F(6)
C
      PX1 = P(1)
      PX2 = P(2)
      PY1 = P(3)
      PY2 = P(4)
      PZ1 = P(5)
      PZ2 = P(6)
      YS  = XS(1)
      ZS  = XS(2)
      YP  = XP(1)
      ZP  = XP(2)
C
C     CONTROL OF ZETA
C
      IF(ZETA) 50,10,10
   10 IF(ZETA-1.) 20,20,50
   20 CONTINUE
C
C     NODAL FORCES AT NODE 1 IN CASE OF DOUBLE SYMMETRY
C
      DO 40 I=1,6
      F(I) = 0.
      DO 30 J=1,12
   30 F(I) = F(I) + EK(I,J)*V(J)
   40 F(I) = F(I) - FP(I) - FE(I)
C
C     STRESS RESULTANTS AT NODE 1
C
      SN  = - F(1)
      SQY =   F(2)
      SQZ =   F(3)
      SMT = - F(4) + YS*F(3) - ZS*F(2)
      SMY =   F(5)
      SMZ = - F(6)
C
C     LOAD ECCENTRICITY FROM SHEAR CENTER AND SOME AUXILIARY CONSTANTS
C
      YPS = YP - YS
      ZPS = ZP - ZS
      Z2 = .5*ZETA*ZETA
      Z3 = ZETA*ZETA*ZETA/6.
C
C     STRESS RESULTANTS AT AN ARBITRARY SECTION
C
      SR(1) = SN  - BL*(ZETA*PX1+Z2*(PX2-PX1))
      SR(2) = SQY + BL*(ZETA*PY1+Z2*(PY2-PY1))
      SR(3) = SQZ + BL*(ZETA*PZ1+Z2*(PZ2-PZ1))
      SR(4) = SMT + ZPS*BL*(ZETA*PY1+Z2*(PY2-PY1))
     +            - YPS*BL*(ZETA*PZ1+Z2*(PZ2-PZ1))
      SR(5) = SMY + SQZ*ZETA*BL + BL*BL*(Z2*PZ1+Z3*(PZ2-PZ1))
      SR(6) = SMZ + SQY*ZETA*BL + BL*BL*(Z2*PY1+Z3*(PY2-PY1))
C
      GO TO 100
C
   50 DO 60 I=1,6
   60 SR(I) = 0.
C
C
  100 RETURN
C
      END
      SUBROUTINE BELS33(F,P,BL,EP,CA,XS,XP)
C
C***********************************************************************
C
C   FINITE ELEMENT LIBRARY :  BELS33                       SINTEF / NTH
C
C
C     BELS33 COMPUTES A KINEMATICALLY CONSISTENT LOAD VECTOR DUE TO
C     DISTRIBUTED EXTERNAL LOADING FOR A 12-DEGREE-OF-FREEDOM BEAM
C     ELEMENT IN A LOCAL SPACE SYSTEM. THE DISTRIBUTED LOADING MAY
C     VARY LINEARLY ALONG THE ELEMENT AXIS. THE COMPONENTS NORMAL TO
C     THE BEAM AXIS MAY BE ECCENTRIC.
C
C     PROGRAMMED BY : F.STENSRUD
C
C     DATE/VERSION  : 28.05.73 / 01
C
C***********************************************************************
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION F(12),P(6),EP(6),CA(2),XS(2),XP(2)
C
      PX1 = P(1)
      PX2 = P(2)
      PY1 = P(3)
      PY2 = P(4)
      PZ1 = P(5)
      PZ2 = P(6)
      E   = EP(1)
      G   = EP(2)
      A   = EP(3)
      RIY = EP(4)
      RIZ = EP(5)
      CAY = CA(1)
      CAZ = CA(2)
      YS  = XS(1)
      ZS  = XS(2)
      YP  = XP(1)
      ZP  = XP(2)
C
C     THE SHEAR CONSTANTS ALY AND ALZ
C
      B = 12.*E/(A*G*BL*BL)
      ALY = B*CAY*RIY
      ALZ = B*CAZ*RIZ
C
C     LOAD ECCENTRICITY FROM SHEAR CENTER
C
      YPS = YP - YS
      ZPS = ZP - ZS
C
C     SOME AUXILIARY CONSTANTS
C
      CX = BL/6.
      CY = BL/(120.*(1.+ALY))
      CZ = BL/(120.*(1.+ALZ))
C
C     THE LOAD VECTOR IN CASE OF DOUBLE SYMMETRY
C
      F(1)  =   CX*(2.*PX1+PX2)
      F(2)  =   CY*((42.+40.*ALY)*PY1 + (18.+20.*ALY)*PY2)
      F(3)  =   CZ*((42.+40.*ALZ)*PZ1 + (18.+20.*ALZ)*PZ2)
      F(4)  = - CX*(ZPS*(2.*PY1+PY2) - YPS*(2.*PZ1+PZ2))
      F(5)  = - CZ*BL*((6.+5.*ALZ)*PZ1 + (4.+5.*ALZ)*PZ2)
      F(6)  =   CY*BL*((6.+5.*ALY)*PY1 + (4.+5.*ALY)*PY2)
      F(7)  =   CX*(PX1+2.*PX2)
      F(8)  =   CY*((18.+20.*ALY)*PY1 + (42.+40.*ALY)*PY2)
      F(9)  =   CZ*((18.+20.*ALZ)*PZ1 + (42.+40.*ALZ)*PZ2)
      F(10) = - CX*(ZPS*(PY1+2.*PY2) - YPS*(PZ1+2.*PZ2))
      F(11) =   CZ*BL*((4.+5.*ALZ)*PZ1 + (6.+5.*ALZ)*PZ2)
      F(12) = - CY*BL*((4.+5.*ALY)*PY1 + (6.+5.*ALY)*PY2)
C
C     ADJUSTMENT FOR SINGLE OR NO SYMMETRY
C
      F(4)  = F(4)  - ZS*F(2) + YS*F(3)
      F(10) = F(10) - ZS*F(8) + YS*F(9)
C
      RETURN
C
      END
      SUBROUTINE BELS34(F,E0,EP)
C
C***********************************************************************
C
C   FINITE ELEMENT LIBRARY :  BELS34                       SINTEF / NTH
C
C
C     BELS34 COMPUTES A VECTOR OF NODAL FORCES DUE TO INITIAL STRAIN
C     FOR A 12-DEGREE-OF-FREEDOM BEAM ELEMENT IN A LOCAL SPACE SYSTEM.
C     THE INITIAL STRAIN MAY VARY LINEARLY OVER THE CROSS SECTION OF
C     THE BEAM, BUT IT MUST BE CONSTANT ALONG THE ELEMENT AXIS.
C
C     PROGRAMMED BY : F.STENSRUD
C
C     DATE/VERSION  : 28.05.73 / 01
C
C***********************************************************************
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION F(12),E0(3),EP(6)
C
      EX  = E0(1)
      EY  = E0(2)
      EZ  = E0(3)
      E   = EP(1)
      A   = EP(3)
      RIY = EP(4)
      RIZ = EP(5)
C
      DO 10 I=1,12
   10 F(I) = 0.
C
C     NODAL FORCE VECTOR
C
      F(1)  = - E*A*EX
      F(5)  = - E*RIZ*EZ
      F(6)  =   E*RIY*EY
      F(7)  = - F(1)
      F(11) = - F(5)
      F(12) = - F(6)
C
C
      RETURN
C
      END
      SUBROUTINE BELS36(EKG,S1,S2,BL,EP,CA,XS)
C
C***********************************************************************
C
C   FINITE ELEMENT LIBRARY :  BELS36                       SINTEF / NTH
C
C
C     BELS36 GENERATES A GEOMETRICAL STIFFNESS MATRIX, BASED ON A
C     LINEARLY VARYING INITIAL FORCE, FOR A 12-DEGREE-OF-FREEDOM BEAM
C     ELEMENT IN A LOCAL SPACE SYSTEM
C
C     PROGRAMMED BY : F.STENSRUD
C
C     DATE/VERSION  : 21.06.73 / 01
C
C***********************************************************************
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION EKG(12,12),EP(6),CA(2),XS(2)
C
      E   = EP(1)
      G   = EP(2)
      A   = EP(3)
      RIY = EP(4)
      RIZ = EP(5)
      CAY = CA(1)
      CAZ = CA(2)
      YS  = XS(1)
      ZS  = XS(2)
C
C     SOME AUXILIARY CONSTANTS
C
      B = 12.*E/(A*G*BL*BL)
      ALY = B*CAY*RIY
      ALZ = B*CAZ*RIZ
      CY = 1./(120.*BL*(1.+ALY)**2)
      CZ = 1./(120.*BL*(1.+ALZ)**2)
C
      DO 10 I=1,12
      DO 10 J=1,12
   10 EKG(I,J) = 0.
C
C     UPPER TRIANGLE OF THE GEOMETRICAL STIFFNESS MATRIX
C
      EKG(2,2)   = CY*(72.+120.*ALY+60.*ALY*ALY)*(S1-S2)
      EKG(2,6)   = - CY*BL*((16.*ALY+10.*ALY*ALY)*(S1+S2)+12.*S2)
      EKG(2,8)   = - EKG(2,2)
      EKG(2,12)  = CY*BL*(12.*S1+(16.*ALY+10.*ALY*ALY)*(S1+S2))
      EKG(3,3)   = CZ*(72.+120.*ALZ+60.*ALZ*ALZ)*(S1-S2)
      EKG(3,5)   = CZ*BL*((16.*ALZ+10.*ALZ*ALZ)*(S1+S2)+12.*S2)
      EKG(3,9)   = - EKG(3,3)
      EKG(3,11)  = - CZ*BL*(12.*S1+(16.*ALZ+10.*ALZ*ALZ)*(S1+S2))
      EKG(5,5)   = CZ*BL*BL*((12.+14.*ALZ+5.*ALZ*ALZ)*S1-
     +                       (4.+6.*ALZ+5.*ALZ*ALZ)*S2)
      EKG(5,9)   = - EKG(3,5)
      EKG(5,11)  = CZ*BL*BL*(2.+10.*ALZ+5.*ALZ*ALZ)*(-S1+S2)
      EKG(6,6)   = CY*BL*BL*((12.+14.*ALY+5.*ALY*ALY)*S1-
     +                       (4.+6.*ALY+5.*ALY*ALY)*S2)
      EKG(6,8)   = - EKG(2,6)
      EKG(6,12)  = CY*BL*BL*(2.+10.*ALY+5.*ALY*ALY)*(-S1+S2)
      EKG(8,8)   =   EKG(2,2)
      EKG(8,12)  = - EKG(2,12)
      EKG(9,9)   =   EKG(3,3)
      EKG(9,11)  = - EKG(3,11)
      EKG(11,11) = CZ*BL*BL*((4.+6.*ALZ+5.*ALZ*ALZ)*S1-
     +                       (12.+14.*ALZ+5.*ALZ*ALZ)*S2)
      EKG(12,12) = CY*BL*BL*((4.+6.*ALY+5.*ALY*ALY)*S1-
     +                       (12.+14.*ALY+5.*ALY*ALY)*S2)
C
C     LOWER TRIANGLE OF EKG FROM SYMMETRY
C
      DO 20 I=1,12
      DO 20 J=1,I
   20 EKG(I,J) = EKG(J,I)
C
C     DOUBLE SYMMETRIC CROSS SECTION ?
C
      BE = ABS(YS) + ABS(ZS)
      BA = 0.00001*SQRT(A)
      IF(BA-BE) 30,1000,1000
C
C     ADJUSTMENT OF EKG DUE TO SINGLE OR NO SYMMETRY
C
   30 CONTINUE
      DO 40 I=1,12
      EKG( 4,I) = EKG( 4,I) - ZS*EKG(2,I) + YS*EKG(3,I)
   40 EKG(10,I) = EKG(10,I) - ZS*EKG(8,I) + YS*EKG(9,I)
C
      DO 50 I=1,12
      EKG(I, 4) = EKG(I, 4) - EKG(I,2)*ZS + EKG(I,3)*YS
   50 EKG(I,10) = EKG(I,10) - EKG(I,8)*ZS + EKG(I,9)*YS
C
C
 1000 RETURN
C
      END
