      SUBROUTINE SDVINI(STATEV,COORDS,NSTATV,NCRDS,NOEL,NPT,
     1 LAYER,KSPT)
C
      INCLUDE 'ABA_PARAM.INC'
C
      DIMENSION STATEV(NSTATV),COORDS(NCRDS)
      REAL(8) PDF, X, B_TOTAL
      PARAMETER PI=3.1416,MEAN =2.2D0,N=493,BARWIDTH=0.007,X0=1.05D0

      STATEV(1) = -1.D0 !Total Time
      B_TOTAL=0.D0
      
      DO I=1,N
          X=X0 + (I-1)*BARWIDTH
          PDF = 32.D0*(X-X0)**2.D0/(PI**2)/(MEAN-X0)**3.D0
     1    *EXP(-4.D0/PI*(X-X0)**2/(MEAN-X0)**2)
          STATEV(I+1) = PDF * BARWIDTH
          B_TOTAL = B_TOTAL + STATEV(I+1)
      END DO
      STATEV(N+2) = B_TOTAL
      
      RETURN
      END
      
      
      SUBROUTINE UMAT(STRESS,STATEV,DDSDDE,SSE,SPD,SCD,
     1 RPL,DDSDDT,DRPLDE,DRPLDT,
     2 STRAN,DSTRAN,TIME,DTIME,TEMP,DTEMP,PREDEF,DPRED,CMNAME,
     3 NDI,NSHR,NTENS,NSTATV,PROPS,NPROPS,COORDS,DROT,PNEWDT,
     4 CELENT,DFGRD0,DFGRD1,NOEL,NPT,LAYER,KSPT,JSTEP,KINC)
C
      INCLUDE 'ABA_PARAM.INC'
C
      CHARACTER*8 CMNAME
      DIMENSION STRESS(NTENS),STATEV(NSTATV),
     1 DDSDDE(NTENS,NTENS),DDSDDT(NTENS),DRPLDE(NTENS),
     2 STRAN(NTENS),DSTRAN(NTENS),TIME(2),PREDEF(1),DPRED(1),
     3 PROPS(NPROPS),COORDS(3),DROT(3,3),DFGRD0(3,3),DFGRD1(3,3)
      
CCC   Change Here!!! (BI(N))
      DIMENSION BBAR(6),DISTGR(3,3),BI(494),CHAINLENGTH(493),BETA_I(493)
      PARAMETER(ZERO=0.D0, ONE=1.D0, TWO=2.D0, THREE=3.D0, SIX=6.D0, 
     1          FOUR=4.D0, NINE=9.D0, N=493, BARWIDTH=0.007, X0=1.05D0)
      REAL(8) I1BAR, SUM
      
C=====================================================================
C     This file calculates and returns Cauchy Stress tensor (STRESS) 
C      and Jacobian(DDSDDE)
C     DISTGR - Distorted tensor
C     BBAR - Deviatoric Left Cauchy-Green Deformation Tensor
C     DET  - Determinant of deformation gradient
C=====================================================================
C=====================================================================
C     Define Parameters
      DTILDE = 1.0D-6
C=====================================================================
C     Initialize chain length
      DO I=1,N
          CHAINLENGTH(I)=X0 + (I-1)*BARWIDTH
      END DO
C=====================================================================
C     Calculate determinant of Deformation Gradient, DISTGR and BBAR
      
      DET = DFGRD1(1,1) * DFGRD1(2,2) * DFGRD1(3,3)
     1    - DFGRD1(1,2) * DFGRD1(2,1) * DFGRD1(3,3)
     2    + DFGRD1(1,2) * DFGRD1(2,3) * DFGRD1(3,1)
     3    + DFGRD1(1,3) * DFGRD1(3,2) * DFGRD1(2,1)
     4    - DFGRD1(1,3) * DFGRD1(3,1) * DFGRD1(2,2)
     5    - DFGRD1(2,3) * DFGRD1(3,2) * DFGRD1(1,1)

      SCALE = DET**(-ONE/THREE)
      DO K1 = 1, 3
          DO K2 = 1, 3
              DISTGR(K2,K1) = SCALE * DFGRD1(K2,K1)
          END DO
      END DO

C=================================================================
C     CALCULATE DEVIATORIC LEFT CAUCHY-GREEN DEFORMATION TENSOR
      
      BBAR(1) = DISTGR(1,1)**2.+DISTGR(1,2)**2.+DISTGR(1,3) ** 2.
      BBAR(2) = DISTGR(2,1)**2.+DISTGR(2,2)**2.+DISTGR(2,3) ** 2.
      BBAR(3) = DISTGR(3,3)**2.+DISTGR(3,1)**2.+DISTGR(3,2) ** 2.
      BBAR(4) = DISTGR(1,1)*DISTGR(2,1) + DISTGR(1,2) * DISTGR(2,2)
     1        + DISTGR(1,3)*DISTGR(2,3)
      BBAR(5) =  DISTGR(1,1) * DISTGR(3,1) + DISTGR(1,2) * DISTGR(3,2)
     1        + DISTGR(1,3) * DISTGR(3,3)
      BBAR(6) =  DISTGR(2,1) * DISTGR(3,1) + DISTGR(2,2) * DISTGR(3,2)
     1        + DISTGR(2,3) * DISTGR(3,3)
      
      I1BAR = BBAR(1) + BBAR(2) + BBAR(3)
C====================================================================
CCC   Change Here!!! (DO I=1,N)
      SUM = 0.D0
      DO I = 1,N
          BI(I) = STATEV(I+1)
          SUM = SUM + BI(I)
      END DO
      BI(N+1) = SUM
      
C===================================================================
C     Set Coefficients
      C1 = 0.0D0
      C2 = 0.0D0
C     Solve Inverse Langevin
      
      
C===================================================================
      
CCC   Change BETA
CCC   Change DO I=1,1
      DO I=1,N
          X = SQRT(I1BAR/(THREE*CHAINLENGTH(I)))
          IF (X .LT. 1.D0) THEN
              CALL INVLANGEVIN(X, OUTPUT)
              BETA_I(I) = OUTPUT
              C1 = C1 + BI(I)/(TWO * I1BAR)*(ONE/((ONE/BETA_I(I))**TWO -
     1            (TWO /(EXP(BETA_I(I))-EXP(-BETA_I(I))))**TWO))-
     2            ONE/(SIX*CHAINLENGTH(I))
     3            *BI(I)*BETA_I(I)*(I1BAR/(THREE*CHAINLENGTH(I)))
     4            **(-THREE/TWO)
          
              C2 = C2 + BI(I)*((I1BAR/(THREE * CHAINLENGTH(I)))
     1        **(-0.5D0))* BETA_I(I)
          END IF
      END DO
      C1 = TWO/(THREE * DET)*C1
      C2 = ONE/(THREE * DET)*C2


C=====================================================================
C     Update Cauchy Stress
      STRESS(1) = C2 * BBAR(1) - I1BAR/THREE + ONE/DTILDE * 
     1            SUM * (DET - ONE/DET)
      STRESS(2) = C2 * BBAR(2) - I1BAR/THREE + ONE/DTILDE *
     1            SUM * (DET - ONE/DET)
      STRESS(3) = C2 * BBAR(3) - I1BAR/THREE + ONE/DTILDE * 
     1            SUM * (DET - ONE/DET)
      STRESS(4) = C2 * BBAR(4)
      STRESS(5) = C2 * BBAR(5)
      STRESS(6) = C2 * BBAR(6)


C================================================================
C     Update Stiffness
      DDSDDE(1,1) = C1 * (BBAR(1) - ONE/THREE * I1BAR)**2.
     1            + C2 * (TWO/THREE*BBAR(1) + TWO/NINE * I1BAR) 
     2            + TWO*DET*SUM/DTILDE
      DDSDDE(1,2) = C1 * (BBAR(1) - ONE/THREE * I1BAR)*(BBAR(2) 
     1            - ONE/THREE*I1BAR)+C2*(-TWO/THREE*BBAR(2)
     2            - TWO/THREE*BBAR(1) + TWO/NINE*I1BAR) 
     3            + TWO*DET*SUM/DTILDE
      DDSDDE(1,3) = C1 * (BBAR(1) - ONE/THREE * I1BAR)*(BBAR(3) 
     1            - ONE/THREE*I1BAR)+C2*(-TWO/THREE*BBAR(3)
     2            - TWO/THREE*BBAR(1) + TWO/NINE*I1BAR) 
     3            + TWO*DET*SUM/DTILDE
      DDSDDE(1,4) = C1 * (BBAR(1) - ONE/THREE * I1BAR) * BBAR(4)
     1            + C2 * (ONE/THREE*BBAR(4))
      DDSDDE(1,5) = C1 * (BBAR(1) - ONE/THREE * I1BAR) * BBAR(5)
     1            + C2 * (ONE/THREE * BBAR(5))
      DDSDDE(1,6) = C1 * (BBAR(1) - ONE/THREE * I1BAR) * BBAR(6)
     1            - C2 * (TWO/THREE * BBAR(6))
      
      DDSDDE(2,2) = C1 * (BBAR(2) - ONE/THREE * I1BAR)**2.
     1            + C2 * (TWO/THREE*BBAR(2) + TWO/NINE * I1BAR) 
     2            + TWO*DET*SUM/DTILDE
      DDSDDE(2,3) = C1 * (BBAR(2) - ONE/THREE * I1BAR)*(BBAR(3) 
     1            - ONE/THREE*I1BAR)+C2*(-TWO/THREE*BBAR(2)
     2            - TWO/THREE*BBAR(3) + TWO/NINE*I1BAR) 
     3            + TWO*DET*SUM/DTILDE
      DDSDDE(2,4) = C1 * (BBAR(2) - ONE/THREE * I1BAR) * BBAR(4)
     1            + C2 * (ONE/THREE*BBAR(4))
      DDSDDE(2,5) = C1 * (BBAR(2) - ONE/THREE * I1BAR) * BBAR(5)
     1            - C2 * (TWO/THREE * BBAR(5))
      DDSDDE(2,6) = C1 * (BBAR(2) - ONE/THREE * I1BAR) * BBAR(6)
     1            + C2 * (ONE/THREE*BBAR(6))
      
      DDSDDE(3,3) = C1 * (BBAR(3) - ONE/THREE * I1BAR)**2.
     1            + C2 * (TWO/THREE*BBAR(3) + TWO/NINE * I1BAR) 
     2            + TWO*DET*SUM/DTILDE
      
      DDSDDE(3,4) = C1 * (BBAR(3) - ONE/THREE * I1BAR) * BBAR(4)
     1            - C2 * (TWO/THREE * BBAR(4))
      DDSDDE(3,5) = C1 * (BBAR(3) - ONE/THREE * I1BAR) * BBAR(5)
     1            + C2 * (ONE/THREE*BBAR(5))
      DDSDDE(3,6) = C1 * (BBAR(3) - ONE/THREE * I1BAR) * BBAR(6)
     1            + C2 * (ONE/THREE*BBAR(6))
      
      DDSDDE(4,4) = C1 * BBAR(4)**2. + C2 * (ONE/TWO*(BBAR(1)+BBAR(2)))
      DDSDDE(4,5) = C1 * BBAR(4)*BBAR(5) + C2 * BBAR(6)/TWO
      DDSDDE(4,6) = C1 * BBAR(4)*BBAR(6) + C2 * BBAR(5)/TWO
      
      DDSDDE(5,5) = C1 * BBAR(5)**2. + C2 * (ONE/TWO*(BBAR(1)+BBAR(3)))
      DDSDDE(5,6) = C1 * BBAR(5)*BBAR(6) + C2 * BBAR(4)/TWO
      
      DDSDDE(6,6) = C1 * BBAR(6)**2. + C2 * (ONE/TWO*(BBAR(2)+BBAR(3)))
      
      DO I=1, NTENS
          DO J=1,I-1
              DDSDDE(I,J) = DDSDDE(J,I)
          END DO
      END DO

C====================================================================
C     Update BI
CCC   Change Here!!!
      NR = 1.D0
      RS = 5.D4
      GAMMA_L = ONE/THREE
      SUM=0.D0
      IF (TIME(2) .NE. STATEV(1)) THEN
          STATEV(1) = TIME(2)
          DO I=1,N
              IF (BI(I) .NE. 0.D0) THEN
                  BI(I) = BI(I)-DTIME*NR*CHAINLENGTH(I)
     1                    *BI(I)/RS*EXP(GAMMA_L*BETA_I(I))
                  IF (BI(I) .LT. 0.D0) THEN
                      BI(I) = 0.D0
                  END IF
              END IF
              STATEV(I+1) = BI(I)
              SUM=SUM + BI(I)
          END DO
          STATEV(N+2)=SUM
      END IF
      RETURN
      END
C=====================================================================

C      FUNCTION INVLANGEVIN(X)
C      IMPLICIT NONE
C      REAL(8) INVLANGEVIN
C      REAL(8), INTENT(IN):: X
C      REAL(8) :: Y, Y0, ERROR, J
C      IF (X .EQ. 0.0D0) THEN
C         Y = 0.0D0
C         INVLANGEVIN = Y
C      ELSE IF (DABS(X) .LT. 0.999D0) THEN
C         ERROR = 1.0D0
C         Y0 = 1.0D-1
C         DO WHILE(ERROR .GT. 1.0D-12)
C            J = -(2.0D0/(EXP(Y0)-EXP(-Y0)))**2.0D0 + 1.0D0/(Y0**2.0D0)
C            Y = -1.0D0/J * ((EXP(Y0) + EXP(-Y0))/(EXP(Y0) - EXP(-Y0)) -
C     1           1.0D0/Y0 - X) + Y0
C            ERROR = DABS(Y - Y0)
C            Y0 = Y
C         END DO        
C         INVLANGEVIN = Y
C      ELSE
C         INVLANGEVIN = SIGN(1.0D0/(1 - DABS(X)),X)
C      END IF
C      RETURN
C      END FUNCTION
C======================================================================
      SUBROUTINE INVLANGEVIN(INPUT, OUTPUT)
      
      INCLUDE 'ABA_PARAM.INC'
      REAL(8) INPUT, OUTPUT, DISCREPANCY, TEMP, JACOBI   
      
      IF (INPUT .EQ. 0.D0) THEN
         OUTPUT = 0.D0
      ELSE IF (ABS(INPUT) .LT. 0.99D0) THEN
         DISCREPANCY = 1.D0
         TEMP = 1.D-1
         DO WHILE(DISCREPANCY .GT. 1.D-7)
            JACOBI = -(2.D0/(EXP(TEMP)-EXP(-TEMP)))**2.D0 + 
     1           1.D0/(TEMP ** 2.D0)
            OUTPUT = -1.D0/JACOBI*((EXP(TEMP)+EXP(-TEMP))/(EXP(TEMP) 
     1      - EXP(-TEMP)) - 1.D0/TEMP - INPUT) + TEMP
            DISCREPANCY = ABS(OUTPUT - TEMP)
            TEMP = OUTPUT
         END DO
      ELSE
         OUTPUT = SIGN(1.0D0/(1 - ABS(INPUT)),INPUT)
      END IF
      END
      
