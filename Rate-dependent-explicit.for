      subroutine vumat(
C Read only (unmodifiable)variables -
     1  nblock, ndir, nshr, nstatev, nfieldv, nprops, lanneal,
     2  stepTime, totalTime, dt, cmname, coordMp, charLength,
     3  props, density, strainInc, relSpinInc,
     4  tempOld, stretchOld, defgradOld, fieldOld,
     5  stressOld, stateOld, enerInternOld, enerInelasOld,
     6  tempNew, stretchNew, defgradNew, fieldNew,
C Write only (modifiable) variables -
     7  stressNew, stateNew, enerInternNew, enerInelasNew )
C
      include 'vaba_param.inc'
C
      dimension props(nprops), density(nblock), coordMp(nblock,*),
     1  charLength(nblock), strainInc(nblock,ndir+nshr),
     2  relSpinInc(nblock,nshr), tempOld(nblock),
     3  stretchOld(nblock,ndir+nshr),
     4  defgradOld(nblock,ndir+nshr+nshr),
     5  fieldOld(nblock,nfieldv), stressOld(nblock,ndir+nshr),
     6  stateOld(nblock,nstatev), enerInternOld(nblock),
     7  enerInelasOld(nblock), tempNew(nblock),
     8  stretchNew(nblock,ndir+nshr),
     8  defgradNew(nblock,ndir+nshr+nshr),
     9  fieldNew(nblock,nfieldv),
     1  stressNew(nblock,ndir+nshr), stateNew(nblock,nstatev),
     2  enerInternNew(nblock), enerInelasNew(nblock)
C
      character*80 cmname
      dimension F(3,3),eye(3,3),xU(3,3),B(3,3),xUbar(3,3),
     1          Bst(3,3),DevBst(3,3),coSig(3,3),alak(3,3),
     2          Bi(494),chainLength(493),beta(493)
      
      parameter (zero=0.d0, one=1.d0, two=2.d0, three=3.d0,
     1            third=one/three, n=493,barWidth=0.007)
      
C=============material parameters=================================
      dtilde=1.d-6

C=============identity matrix=====================================
      eye=zero
      do i=1,3
          eye(i,i)=one
      end do
C=================initialize chain length=========================      
      do i=1,n
          chainLength(i)=1.05d0 + (i-1)*barWidth
      end do
C=============main program starts here============================
      do km = 1,nblock
C       deformation gradient
        F(1,1)=defgradNew(km,1)
        F(2,2)=defgradNew(km,2)
        F(3,3)=defgradNew(km,3)
        F(1,2)=defgradNew(km,4)
        F(2,3)=defgradNew(km,5)
        F(3,1)=defgradNew(km,6)
        F(2,1)=defgradNew(km,7)
        F(3,2)=defgradNew(km,8)
        F(1,3)=defgradNew(km,9)
C       strerch tensor
        xU(1,1)=stretchNew(km,1)
        xU(2,2)=stretchNew(km,2)
        xU(3,3)=stretchNew(km,3)
        xU(1,2)=stretchNew(km,4)
        xU(2,3)=stretchNew(km,5)
        xU(3,1)=stretchNew(km,6)
        xU(2,1)=xU(1,2)
        xU(3,2)=xU(2,3)
        xU(1,3)=xU(3,1)
C       Jacobian
        det=F(1,1)*F(2,2)*F(3,3)
     1     +F(1,2)*F(2,3)*F(3,1)
     2     +F(1,3)*F(2,1)*F(3,2)
     3     -F(1,1)*F(2,3)*F(3,2)
     4     -F(1,2)*F(2,1)*F(3,3)
     5     -F(1,3)*F(2,2)*F(3,1)
        
C       State Variable Bi(i)
        sum = 0.d0
        do i = 1,n
          Bi(i) = stateOld(km, i)
          sum = sum + Bi(i)
        end do
        Bi(n+1) = sum
        
C       Left Cauchy-Green strain
        B=MATMUL(F,TRANSPOSE(F))
        Bst=B*(det**(-two/three))
        xI1bar=Bst(1,1)+Bst(2,2)+Bst(3,3)
        xUbar=xU*(det**(-third))
C=====================================
        coeff=zero
        do i=1,n
          x = sqrt(xI1bar/(three*chainLength(i)))
          if (x .lt. 1.d0) then
              CALL invlangevin(x, output)
              beta(i) = output
              
              coeff = coeff + Bi(i)*((xI1Bar/(three*chainLength(i)))
     1        **(-0.5d0))* beta(i)
              
          end if
        end do
        
        coeff=100.d0/(three * det)*coeff
C=====================================
        
C       Corotational Cauchy stress
        coSig=coeff*(MATMUL(xUbar,xUbar)-third*xI1bar*eye)
     1   + (det-one/det)*sum*eye/dtilde
        
C       update stress
        stressNew(km,1)=coSig(1,1)
        stressNew(km,2)=cosig(2,2)
        stressNew(km,3)=cosig(3,3)
        stressNew(km,4)=cosig(1,2)
        stressNew(km,5)=cosig(2,3)
        stressNew(km,6)=cosig(3,1)
        
C       update bt
        nr=one
        rs=5.d4
        gamma=third
        sum=zero
        do i=1,n
            if (Bi(i) .ne. zero) then
                Bi(i)=Bi(i)
     1                -dt*nr*chainLength(i)*Bi(i)/rs*exp(gamma*beta(i))
                if (Bi(i) .lt. 1.d-7) then
                    Bi(i)=zero
                end if
            end if
            stateNew(km,i)=Bi(i)
            sum=sum+Bi(i)
        end do
        stateNew(km,n+1)=sum
          
      end do
        
      return
      end
C==========================================================
      subroutine invlangevin(input, output)
      
      include 'vaba_param.inc'
      real(8) input, output, discrepancy, temp, jacobi   
      
      if (input .eq. 0.d0) THEN
         output = 0.d0
      else if (abs(input) .lt. 0.99d0) then
         discrepancy = 1.d0
         temp = 1.d-1
         do while(discrepancy .gt. 1.d-7)
            jacobi = -(2.d0/(exp(temp)-exp(-temp)))**2.d0 + 
     1           1.d0/(temp ** 2.d0)
            output = -1.d0/jacobi*((exp(temp)+exp(-temp))/(exp(temp) 
     1      - exp(-temp)) - 1.d0/temp - input) + temp
            discrepancy = abs(output - temp)
            temp = output
         end do
      else
         output = sign(1.0d0/(1 - abs(input)),input)
      end if
      end
      