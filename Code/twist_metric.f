	!!! this code calculate the quantum metric 
	module lanczos1
       integer, parameter :: nx0=32,ny0=32
       integer, parameter :: kx0=nx0,ky0=ny0
        integer, parameter::nq1=21,nq2=21,nl=2,nk1=kx0,nk2=ky0
        integer,parameter :: nsp=nq1*nq2*nl, nsk=nk1*nk2
	integer, parameter :: nnp=nsp, np1=nq1, np2=nq2
        integer, parameter:: nnp1=nsp
        logical chernband
       real*8 tw_angle
	end module lanczos1

        module parameter
        use lanczos1
       integer, parameter :: ns0=ny0*nx0,ns=ns0*2
	integer mnt_d(ns0, ns0), kp1, kp2,kp11
        integer  mq0(2,2)
        integer  mg0(2,2)
	integer mnt(kx0*ky0*2,2), ind_mnt(0:kx0,0:ky0)
        complex*16 com, com1,com2,com3,com4,com11,com12
        real*8  valp(nnp),valf(nsk),tw1
        double complex ci, glist(7),k_plus,k_minus,g1,g2,kg1,kg2
        double complex q1i(6), q2i(6)
        real*8 vq, ee1, ee2, area, dgate
        integer ix0, iy0, ng1, nn11
	parameter(pi=3.14159265358979)
        integer q11, q12, q1, q2, q3, q4, q5, q6
        double precision phi,V1,V2,W1,W2, VD1
        integer, parameter :: neibt=12, kb=3
        integer  neib(0:nsp,30), gr0(100)
        complex*16 hamp(nnp,nnp),vecp(nnp,nnp)
        complex*16 vecf0(nnp,nsk), vecf1(nnp*10, 4*nk1,4*nk2)
        double precision mass,a0,aM,factor
        double precision angle,x,y, g, prefac
        double complex tensor(0:30,2,2), ck, ck0, ckshift
            real*8 me, m00, hbar, rJ, prefactor
	end module parameter

        subroutine band_spectrum_para()
        use parameter
        integer m1(2), m2(2)

               me= 9.1093837e-31
        m00 = 0.620d0*me
        hbar=1.054571817e-34
        rJ=6.582119569509066e-13/hbar
        prefactor = hbar**2/(2.0d0*m00)*rJ*1.0e20
        ci=dcmplx(0.0d0, 1.0d0)

        a0=3.52d0
        factor=prefactor
        tw1=tw_angle*180.0d0/pi
        aM=a0/(2.0d0*sin(tw_angle/2.0d0))


        neib=0
        do ix=1, np1
        do jy=1, np2
        ii=jy+(ix-1)*np2
        ix1=ix+1
        if(ix==np1)ix1=1
        jy1=jy+1
        if(jy==np2)jy1=1
        ix2=ix-1
        if(ix==1)ix2=np1
        if(ix1.ne.0)neib(ii,1)=jy+(ix1-1)*np2
        if(jy1.ne.0)neib(ii,2)=jy1+(ix-1)*np2
        if(jy1*ix2.ne.0)neib(ii,3)=jy1+(ix2-1)*np2
        j11=jy+(ix1-1)*np2
        j12=jy1+(ix-1)*np2
        j13=jy1+(ix2-1)*np2
        neib(j11,neibt/2+1)=ii
        neib(j12,neibt/2+2)=ii
        neib(j13,neibt/2+3)=ii
        enddo
        enddo 

        do i=1, np1*np2
                j1=neib(i,1)
                if(j1.ne.0)then
                j2=neib(j1,2)
                neib(i,4)=j2
                neib(j2, 4+neibt/2)=i
                endif

                j1=neib(i,2)
                j2=neib(j1,3)
                        if(j1*j2.ne.0)then
                neib(i,5)=j2
                neib(j2, 5+neibt/2)=i
                        endif
                
                j1=neib(i,3)
                j2=neib(j1,1+neibt/2)
                if(j1*j2.ne.0)then
                neib(i,6)=j2
                neib(j2, 6+neibt/2)=i
                        endif
        enddo

33      format(8i7)


        do j=1, 7
        angle=(j-1)*pi/3.0d0
        prefac=4*pi/(dsqrt(3.0d0)*aM)
        glist(j)=prefac*cdexp(ci*angle)
        enddo

        g1=glist(1)
        g2=glist(3)  
        m1(1)=1
        m1(2)=0   
        m2(1)=-1
        m2(2)=1

        mg0=0
        mq0=0
        mq0(1,1)=m1(1)   
        mq0(1,2)=m1(2)  
        mq0(2,1)=m2(1) 
        mq0(2,2)=m2(2)  
        mg0(1,1)=dfloat(m2(2))/(m2(2)*m1(1)-m2(1)*m1(2))
        mg0(1,2)=dfloat(-m1(2))/(m2(2)*m1(1)-m2(1)*m1(2))
        mg0(2,1)=dfloat(m2(1))/(m2(1)*m1(2)-m1(1)*m2(2))
        mg0(2,2)=dfloat(-m1(1))/(m2(1)*m1(2)-m1(1)*m2(2))

        if(nx0==2)then
        g1=glist(1)
        g2=2*glist(1)+glist(2)
        mq0(1,1)=1
        mq0(1,2)=0
        mq0(2,1)=2
        mq0(2,2)=1

        mg0(1,1)=1
        mg0(1,2)=0
        mg0(2,1)=-2
        mg0(2,2)=1
        m1(1)=1
        m1(2)=0    
        m2(1)=2
        m2(2)=1
        mq0(1,1)=m1(1)   
        mq0(1,2)=m1(2)   
        mq0(2,1)=m2(1)   
        mq0(2,2)=m2(2)  
        mg0(1,1)=dfloat(m2(2))/(m2(2)*m1(1)-m2(1)*m1(2))
        mg0(1,2)=dfloat(-m1(2))/(m2(2)*m1(1)-m2(1)*m1(2))
        mg0(2,1)=dfloat(m2(1))/(m2(1)*m1(2)-m1(1)*m2(2))
        mg0(2,2)=dfloat(-m1(1))/(m2(1)*m1(2)-m1(1)*m2(2))
        endif

        k_minus= (glist(1)+glist(6))/3.0d0
        k_plus=(glist(1)+glist(2))/3.0d0
        kg1=(glist(2)+glist(3))/3.0d0
        return
        end


        subroutine band_spectrum()
        use parameter 
        integer dlayer, kspin, ks1
        dlayer=1

        tensor=0
        tensor(0,1,2)=W1 
        tensor(0,2,1)=W1 
        do j1=1,3
        tensor(j1,1,1)=tensor(j1,1,1)-V1*cdexp(ci*(-1)**(j1-1)*phi)
        tensor(j1,2,2)=tensor(j1,2,2)-V1*cdexp(ci*(-1)**(j1-1)*(-phi))
        enddo

        do j1=4,6
        tensor(j1,1,1)=tensor(j1,1,1)-V2 
        tensor(j1,2,2)=tensor(j1,2,2)-V2  
        enddo

        do j1=2,3
        tensor(j1,2,1)=tensor(j1,2,1)+W1
        enddo
        j1=1
        tensor(j1,2,1)=tensor(j1,2,1)+W2   
        tensor(j1,1,2)=tensor(j1,1,2)+W2  
        do j1=5, 5
        tensor(j1,2,1)=tensor(j1,2,1)+W2
        enddo

        kspin=1
        ks1=1
        if(kspin==2)ks1=-1
        do k1=1,kx0+1
        do k2=1,ky0+1 
        do kp11=0,1

        dx=36.0d0*2.0d0
        ckshift=kp11*g1/nx0  
        ckshift=ckshift/dx
        ck0=(k1-1)/dfloat(kx0)*g1+(k2-1)/dfloat(ky0)*g2
        ki=k2+(k1-1)*ky0
         do i=1,nnp
         do j=1,nnp
           hamp(i,j)=cmplx(0.d0,0.d0)        
         enddo
         enddo

        do ix=1, np1
        do iy=1, np2
        ii=iy+(ix-1)*np2     
        ix0=np1/2+1
        iy0=np2/2+1
        ck=ck0+glist(1)*(ix-ix0)+glist(2)*(iy-iy0)+ckshift
        hamp(2*ii-1,2*ii-1)=hamp(2*ii-1,2*ii-1)+
     &    factor*cdabs(ck-ks1*k_plus)**2-VD1/2.0d0
        hamp(2*ii,2*ii)=hamp(2*ii,2*ii)+factor*cdabs(ck-ks1*k_minus)**2
     &     +VD1/2.0d0
        hamp(2*ii-1,2*ii)=hamp(2*ii-1,2*ii)+tensor(0,1,2)
        hamp(2*ii,2*ii-1)=hamp(2*ii,2*ii-1)+tensor(0,2,1)
        do ik=1,6
        jj=neib(ii,ik)
        if(jj.ne.0)then
        hamp(2*ii-1,2*jj-1)=hamp(2*ii-1,2*jj-1)+tensor(ik,1,1)
        hamp(2*ii,2*jj)=hamp(2*ii,2*jj)+tensor(ik,2,2)
        hamp(2*ii-1,2*jj)=hamp(2*ii-1,2*jj)+tensor(ik,1,2)
        hamp(2*ii,2*jj-1)=hamp(2*ii,2*jj-1)+tensor(ik,2,1)
        
        hamp(2*jj-1,2*ii-1)=hamp(2*jj-1,2*ii-1)+
     &    conjg(tensor(ik,1,1))
        hamp(2*jj,2*ii)=hamp(2*jj,2*ii)+conjg(tensor(ik,2,2))
        hamp(2*jj,2*ii-1)=hamp(2*jj,2*ii-1)+conjg(tensor(ik,1,2))
        hamp(2*jj-1,2*ii)=hamp(2*jj-1,2*ii)+conjg(tensor(ik,2,1))
        endif
        enddo
        end do
        end do
        do i=1, nnp-1
        do j=i+1, nnp
        if(hamp(i,j).ne.conjg(hamp(j,i)))then
        write(*,*)hamp(i,j),hamp(j,i), 'wrong', i,j
        endif
        enddo
        enddo
        call eigen1(nnp, hamp,valp, vecp)
           do ii=1, kb
        vecf1(1+(ii-1)*nnp:ii*nnp,k1+(nk1+1)*kp11,k2)
     &    =vecp(1:nnp,ii)
        enddo
        ck=ck0
        enddo                 !!!for maxnum loop
        enddo
        enddo

        open(91,file='quantum_metric.dat', access='append')
        open(81,file='T_metric.dat', access='append')

        do ii=1, kb
        gr1=0.0d0
        do k1=1,ny0
        do k2=1,nx0
        nj=k1+(k2-1)*ny0
        com=0.0
        com1=0.0
             do j=1+(ii-1)*nnp,nnp*ii
        com11=(vecf1(j,k1,k2)-vecf1(j,k1+nk1+1,k2))*dx
        com=com+dconjg(com11)*com11
        com1=com1+dconjg(com11)*vecf1(j,k1,k2)
        enddo
        gr1=gr1+(com-dconjg(com1)*com1)*dsqrt(3.0d0)/2.0d0
        write(91,24)ii,tw_angle*180/pi, gr1
24      format(i6, 8f18.9)
        enddo
        enddo
        gr0(ii)=gr1
        enddo
25      format(8f18.9)
        write(81,25)tw_angle*180/pi,V2,gr0(1:kb)/pi
        close(91)
        close(81)
        return
        end


	subroutine eigen1(n,a,val,vec)
	use lanczos1   
      parameter (n3=nnp1)
      PARAMETER (NMAX=n3,LDA=NMAX)
      PARAMETER (LWORK=NMAX*NMAX+2*NMAX,LIWORK=3+5*NMAX,
     +LRWORK=4*NMAX*NMAX)
      CHARACTER JOB,UPLO
	complex*16 vec(n,n)
      real*8 val(n)
      complex*16 A(LDA,NMAX),WORK(LWORK)
      real*8     RWORK(LRWORK)
      INTEGER IWORK(LIWORK)
      EXTERNAL zheevd!!!!cheevd

      UPLO='U'
      JOB='V'
      IF (UPLO.EQ.'U') THEN
      END IF
      CALL zheevd (JOB,UPLO,N,A,LDA,val,WORK,LWORK,RWORK,LRWORK,IWORK,
     +LIWORK,INFO)

      IF (INFO.GT.0)THEN
      WRITE (*,*) 'Failure to converge'
      ENDIF


	do i=1,n
	do j=1,n
	vec(i,j)=a(i,j)
	enddo
	enddo
	return
	end
c---------------------------------------------------------
	!!! paragram  main
	use parameter
        chernband=.true.
        ckshift=0.0d0

        VD1=0.0d0
        V1=17.50d0
        W1=-6.50d0
        W2=11.750d0
        phi=-1.020843d0
        V20=-10.50
        do 1005  kt=0, 8
        do 1005   kv2=-11, 5
        V2=V20+kv2*0.5
        ckshift=(kp1-1)*g1/nx0+(kp2-1)*g2/ny0
        tw_angle=(1.400d0+kt*0.1)/180.0d0*pi

         call band_spectrum_para()
        call band_spectrum()
1005    continue
	end
