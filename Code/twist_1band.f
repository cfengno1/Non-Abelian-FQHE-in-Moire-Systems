	!!! this code  works for the case we project holes into one band
	module lanczos1
        integer, parameter :: nx0=2, nx=nx0
       integer, parameter :: ny0=14,ny=ny0, nphase=1
!! change the system size nx0, ny0, nphase=1 no Chern, and nphase=13
!with Chern number in a mesh of 12x12 in phase space.

        integer, parameter::n16=1000, n24=114000,n28=1480000
        integer*8, parameter::n32=19000000
        integer,parameter::k16=200000,k24=20000000
        integer*8,parameter::k28=250000000,k32=4600000000
       integer, parameter :: ns0=ny0*nx0,ns=ns0*2
        integer*8,parameter::nlanc=n32*(ns0/32)+n28*(ns0/28)+n24
        integer*8,parameter::ksave1=k32*(ns0/32)+k28*(ns0/28)+k24
	integer, parameter::  ll1=1
       integer, parameter::  jchonx=1, jchony=1
	parameter(ne=ns0/2,nt=1,nt1=1,nt11=1)
        parameter(ne1d=(ne+1)/2,ne2d=ne-ne1d)!nbasis,depends nd ns only
	parameter(pi=3.14159265358979,ma=2,numd=ne)
        integer, parameter::nq1=21,nq2=21,nl=2
        integer, parameter:: nnp1=nq1*nq2*nl   
        integer,parameter :: nsp=nq1*nq2*nl
	integer, parameter :: nnp=nsp, np1=nq1, np2=nq2
	end module lanczos1

	module param1
	use lanczos1
	integer icomb
        integer model
        logical chernband
       integer, parameter :: kx0=nx0,ky0=ny0,nst=ns0
        integer, parameter::nk1=kx0,nk2=ky0, nsk=nk1*nk2
        integer, parameter :: nx1=(nx0+1)/2, ns00=(ns0+1)/2 !!!!ny0*nx1
       integer, parameter :: ns01=ns0-ns00, nbasis01=2**ns01
          integer, parameter :: ms0=nq1*6, nsq=ms0*2
	integer jcho(2),jcho1(2),kmn(10),kin, nbasisf
        double precision phi,V1,V2,W1,W2,VD1,tw_angle
        double precision epsilon, dgate
	integer mnt_d(ns0, ns0), kp1, kp2
        integer  mq0(2,2)
        integer  mg0(2,2)
	integer mnt(kx0*ky0*2,2), ind_mnt(0:kx0,0:ky0)
        complex*16 :: fq(ns0,ns0,-nq1-10:nq1+10,-nq2-10:nq2+10)
        complex*16 :: fq0(ns0,ns0,-nq1-10:nq1+10,-nq2-10:nq2+10)
        complex*16 :: fq00(ns0,ns0,-nq1-10:nq1+10,-nq2-10:nq2+10)
        complex*16 ::fq10(ns0,ns0,-nq1-10:nq1+10,-nq2-10:nq2+10)
         real*8 density(ns, 6)
        complex*16 vecb(nnp,ns0, nphase, nphase)
        real*8  valp(nnp),valf(nsk),tw1,val_hf(nsk),val_hf1(nsk)
        double complex ci, glist(7),k_plus,k_minus,g1,g2,kg1,kg2
        double complex q1i(6), q2i(6)
        double precision gk1(-ms0:ms0, -ms0:ms0, 2)
        real*8 vq, ee1, ee2, area
        integer ix0, iy0, ng1, nn11

       integer, parameter::ndimen=(nst+1)*(nst+1)*(nst+1)
	real aa1
	integer i0s(10*ns)
	integer*8 jc1
        parameter(nbasis=nlanc)
	integer*8 kjjt,kjjt1,kjj,kjj1
        integer q11, q12, q1, q2, q3, q4, q5, q6
	integer ja,jb,krr,nseed
	integer, parameter::ntot3=ns*ns*ns
        complex*16 v(1,1:ndimen),vin(-ndimen:ndimen)
          integer kzone(-ms0:ms0,-ms0:ms0)
        complex*16 vsq(ndimen,-ms0:ms0, -ms0:ms0)
        complex*16 vsq1(ndimen,-ms0:ms0, -ms0:ms0)
           complex*16 spin(-ms0:ms0, -ms0:ms0,ma)
           complex*16 spin0(-ms0:ms0, -ms0:ms0,ma)
               complex*16 temps(ndimen,6)
        integer ind_jj(1:ns0,1:ns0,1:ns0,1:ns0,1:ns0)
        parameter (nbasis0=2**ns00) !!! nbasis0=2**ns0
	integer*8, parameter::nn=nlanc! nn for nt nt1 for chern number
	real*8 ham(nn)
        real*8 hop
	real*8 val(nn), val1(nn)
	end module param1

	
	module chern0
	use param1
	real cxy(ma),bt1(ma),bt(ma),bt0(ma)
     & ,btm(ma),b0(nt1,ma), dbt1(ma),dbt2(ma),dbt11(ma),dbt12(ma)
     & ,btm1(ma),btm2(ma)
	integer icxy(nn)
	complex vec0(nn,ma,nt),vec00(nn,ma)
	real bt00(nn)
	end module chern0


	module interact
	use param1
	integer, parameter:: n11=1140
        integer, allocatable, dimension(:) ::  nvin,nvin1
        integer*2, allocatable, dimension(:):: nvi 
        integer, allocatable, dimension(:,:):: label,label1
	integer*8 kjjs(0:nn),kjjss(0:nn),kjjs1(0:nn+1)
	integer kjjs2(0:nn+1) 
        integer nns1(1:n11),nns2(1:n11),nne(1:n11)
	integer nns(n11)
	integer ki(ns0*3),ki1(ns0*3),numb(nbasis0,5)
	integer nnsb(nbasis0,ns0),nneb(nbasis0,ns0)
	integer newbb(nbasis,2)!!!!!!,newb1(nbasisd)
	integer njm(nbasis0),njm1(nbasis0),newj0(nbasis0)
	integer numr(nbasis0),numr1(nbasis0)
        real*8 vec(nn*2*ma),vecph(nn*2*ma,nphase,nphase)
        real*8 vec1(nn*2*ma) !, vecph(nn*8, nphase, nphase)
	end module interact

        module tmd
        use param1
	integer,parameter::lwork=2*nnp+nnp**2
        integer,parameter:: lrwork=1+5*nnp+4*nnp**2
        integer,parameter::liwork=3+5*nnp
	real*8 :: rwork(1+5*nnp+4*nnp**2)
        integer, parameter :: neibt=12
        integer  neib(0:nsp,30)
        integer  iwork(3+5*nnp)
        complex*16 hamp(nnp,nnp),vecp(nnp,nnp),vecf(nnp,nsk)
        complex*16 vecf0(nnp,nsk)
        complex*16 work(2*nnp+nnp**2)
        double precision mass,a0,aM,factor
        double precision angle,x,y, g, prefac
        double complex tensor(0:30,2,2), ck, ck0, ckshift
            real*8 me, m00, hbar, rJ, prefactor
        end module tmd
           module constant
        double precision d11(10)
        end module constant

        subroutine band_spectrum_para()
        use tmd
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

!!!! newer
        g1=glist(1)
        g2=glist(3)  ! g(2)-g(1)
        m1(1)=1
        m1(2)=0    ! for g1
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
        write(222,*)real(k_minus), imag(k_minus)
        write(222,*)real(k_plus), imag(k_plus)
        write(222,*)real(kg1), imag(kg1)
        write(222,*)-real(k_minus), -imag(k_minus)
        write(222,*)-real(k_plus), -imag(k_plus)
        write(222,*)-real(kg1), -imag(kg1)
        write(222,*)real(k_minus), imag(k_minus)

         kzone=0
        do iy=-ms0, ms0
        do ix=-ms0, ms0
           kg1=g1/nx0*ix+g2/ny0*iy
        gk1(ix, iy,1)=real(kg1)
        gk1(ix, iy,2)=imag(kg1)
        if(cdabs(k_plus**2)+0.0001.ge.(gk1(ix,iy,1)**2+
     &     gk1(ix,iy,2)**2))then
        endif
        if(cdabs(glist(1)**2/4)-0.000001.gt.(gk1(ix,iy,1)**2+
     &     gk1(ix,iy,2)**2))then
        kzone(ix,iy)=1
                else
        if(cdabs(glist(1)**2/4.+0.00028).gt.(gk1(ix,iy,1)**2+
     &     gk1(ix,iy,2)**2))then
        if(gk1(ix,iy,1).gt.-0.032)then
        if(gk1(ix,iy,2).lt.-0.031.or.
     &  gk1(ix,iy,2).gt.-0.028)then
        kzone(ix,iy)=1
        endif
        endif
        endif
        endif
        enddo
        enddo


        return
        end

        subroutine band_spectrum()
        use tmd
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
        do kspin=1, dlayer
        ks1=1
        if(kspin==2)ks1=-1
        if(kspin.eq.1)then
        do k1=1,kx0
        do k2=1,ky0 
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
        do i=1,1 
        write(22,*)k1-1,k2-1, valp(1:2)
        enddo
 
        vecf(1:nnp, ki)=vecp(1:nnp,1+ll1)
        vecf0(1:nnp, ki)=vecp(1:nnp,1)
        vecb(1:nnp,ki, kp1,kp2)=vecp(1:nnp,1+ll1)
        valf(ki)=valp(1+ll1)
        if(ki==1)then
        do i=1,np1
        do j=1,np2
        ii=j+(i-1)*np2
        if(cdabs(vecf(2*ii-1,1))+cdabs(vecf(2*ii,1)).gt.0.00001)then
        write(102,16)i,j,ii,real(vecf(2*ii-1,1)),imag(vecf(2*ii,1))
     &     ,real(vecf(2*ii,1)),imag(vecf(2*ii,1))
        endif
        enddo
        enddo
        endif
16      format(3i6, 4f18.12)
        ck=ck0
        enddo                 !!!for maxnum loop
        enddo

        fq=0.0d0
        fq0=0.0d0
        fq00=0.0d0
        fq10=0.0d0
        do k11=1, ns0
        do k12=1, ns0
        do q1=1, nq1  
        do q2=1, nq2  
        do q11=1,nq1
        do q12=1,nq2
        iq1=q2+(q1-1)*nq2
        iq2=q12+(q11-1)*nq2

        ia1=q1-q11  
        ib1=q2-q12
        fq(k11,k12,ia1,ib1)=fq(k11,k12,ia1,ib1)+
     & dconjg(vecf(2*iq1-1,k11))*vecf(2*iq2-1,k12)
     & +dconjg(vecf(2*iq1,k11))*vecf(2*iq2,k12)

        fq0(k11,k12,ia1,ib1)=fq0(k11,k12,ia1,ib1)+
     & dconjg(vecf0(2*iq1-1,k11))*vecf0(2*iq2-1,k12)
     & +dconjg(vecf0(2*iq1,k11))*vecf0(2*iq2,k12)

        fq10(k11,k12,ia1,ib1)=fq10(k11,k12,ia1,ib1)+
     & dconjg(vecf(2*iq1-1,k11))*vecf0(2*iq2-1,k12)
     & +dconjg(vecf(2*iq1,k11))*vecf0(2*iq2,k12)
        enddo
        enddo
        enddo
        enddo
        enddo
        enddo
        endif

        if(kspin.eq.2)then
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
        if(kspin==1)tensor(j1,2,1)=tensor(j1,2,1)+W1
        if(kspin==2)tensor(j1,1,2)=tensor(j1,1,2)+W1
        enddo
        j1=1
        tensor(j1,2,1)=tensor(j1,2,1)+W2 
        tensor(j1,1,2)=tensor(j1,1,2)+W2 
        do j1=5, 5
        if(kspin==1)tensor(j1,2,1)=tensor(j1,2,1)+W2
        if(kspin==2)tensor(j1,1,2)=tensor(j1,1,2)+W2
        enddo
        do k1=1,kx0
        do k2=1,ky0 
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
        vecf0(1:nnp, ki)=vecp(1:nnp,1)
        if(ki==1)then
        do i=1,np1
        do j=1,np2
        ii=j+(i-1)*np2
        if(cdabs(vecf(2*ii-1,1))+cdabs(vecf(2*ii,1)).gt.0.00001)then
        write(102,16)i,j,ii,real(vecf(2*ii-1,1)),imag(vecf(2*ii,1))
     &     ,real(vecf(2*ii,1)),imag(vecf(2*ii,1))
        endif
        enddo
        enddo
        endif
        ck=ck0
        enddo                 
        enddo

        fq00=0.0d0
        do k11=1, ns0
        do k12=1, ns0
        do q1=1, nq1  
        do q2=1, nq2  
        do q11=1,nq1
        do q12=1,nq2
        iq1=q2+(q1-1)*nq2
        iq2=q12+(q11-1)*nq2

        ia1=q1-q11  !! momentum transfer-x
        ib1=q2-q12

       fq00(k11,k12,ia1,ib1)=fq00(k11,k12,ia1,ib1)+
     & dconjg(vecf0(2*iq1-1,k11))*vecf0(2*iq2-1,k12)
     & +dconjg(vecf0(2*iq1,k11))*vecf0(2*iq2,k12)
        enddo
        enddo
        enddo
        enddo
        enddo
        enddo
        endif
        enddo




14      format(4i6, 2f18.11)
        return
        end


	subroutine coulomb()                !  ns ~ ab/2*pi*ell^2
	use param1
        use tmd
	complex*16 com1,com11,com12,com3,com4
        real*8 vq0
        complex*16 u, cqx, cqx1 
        integer qd1, qd2

        ee1=1.60217663*8.9875517923d0*1e+03  

        area=1.0d0
        ee1=ee1*2.0d0*pi/area/epsilon
        area=aM*aM*dsqrt(3.0d0)/2.0d0*kx0*ky0
        ee1=ee1/area
	mnt=0
	ind_mnt=0
	do kx=0,kx0-1
	do ky=0,ky0-1
	j=ky+1+kx*ky0
        mnt(j,1)=kx
	mnt(j,2)=ky
	mnt(j+ns0,1)=kx
	mnt(j+ns0,2)=ky
	ind_mnt(kx,ky)=j
	enddo
	enddo
        vc=0.0d0

	do kx=0,kx0-1
	do ky=0,ky0-1

	j1=ind_mnt(kx,ky)

	do kx1=0,kx0-1
	do ky1=0,ky0-1
	j2=ind_mnt(kx1,ky1)

	jx1=mod(kx-kx1+kx0,kx0)
	jy1=mod(ky-ky1+ky0,ky0)

	mnt_d(j1,j2)=ind_mnt(jx1,jy1)
	enddo
	enddo
	enddo
	enddo

        vsq=0
        vsq1=0.0d0
	v=dcmplx(0.d0,0.d0)
        u=dcmplx(0.0d0, 1.0d0)
	
	do j1=1,ns0
	do j2=1,ns0
	j12=j1+(j2-1)*ns0
	do j3=1,ns0
	do j41=1,1!!!ns0
	j4x=mod(mnt(j1,1)+mnt(j2,1)-mnt(j3,1)+kx0,kx0)
	j4y=mod(mnt(j1,2)+mnt(j2,2)-mnt(j3,2)+ky0,ky0)
	j4=ind_mnt(j4x, j4y)
	jj=j1+(j2-1)*ns0+(j3-1)*ns0**2
	n=jj

        vq0=0.0d0
        do q1=-nq1,nq1
        do q2=-nq2, nq2
        qd1=mnt(j1,1)+mnt(j2,1)-mnt(j3,1)-mnt(j4,1)
        qd2=mnt(j1,2)+mnt(j2,2)-mnt(j3,2)-mnt(j4,2)
        qd1=qd1/kx0
        qd2=qd2/ky0
        q3=-q1-qd1*mq0(1,1)-qd2*mq0(2,1)
        q4=-q2-qd1*mq0(1,2)-qd2*mq0(2,2)

        cqx=(mnt(j1,1)-mnt(j4,1))/dfloat(kx0)*g1+q1*glist(1)
        cqx=cqx+(mnt(j1,2)-mnt(j4,2))/dfloat(ky0)*g2
        cqx=cqx+q2*glist(2)

        vq0=cdabs(cqx)  
        if(vq0.ge.1e-15)then
        if(dgate.gt.0.1)then
         vq0=tanh(dgate*vq0)/vq0
        else
         vq0=1.0d0/vq0
        endif
        v(1,n)=v(1,n)+vq0*fq(j1,j4,q1,q2)*fq(j2,j3,q3,q4)
        iaa=mnt(j1,1)-mnt(j4,1) !!!!+ns0+(q1+nq1)*ns0*2
        ibb=mnt(j1,2)-mnt(j4,2) !!!!+ns0+(q1+nq1)*ns0*2
        iaa=iaa+(q1*mg0(1,1)+q2*mg0(2,1))*kx0
        ibb=ibb+(q1*mg0(1,2)+q2*mg0(2,2))*ky0
        cqx1=iaa/dfloat(kx0)*g1+ibb/dfloat(ky0)*g2

        if(abs(iaa).le.ms0.and.abs(ibb).le.ms0)then
        vsq(n,iaa,ibb)=vsq(n,iaa,ibb)+fq(j1,j4,q1,q2)*fq(j2,j3,q3,q4)
        vsq1(n,iaa,ibb)=vsq1(n,iaa,ibb)+vq0*fq(j1,j4,q1,q2)*
     &     fq(j2,j3,q3,q4)
        endif
        endif
	enddo
	enddo
	enddo
	enddo
	enddo
        end do
        v=v*ee1/2.0d0
        vsq1=vsq1*ee1/2.0d0

111	continue
1110	continue

        val_hf=0.0d0
        val_hf1=0.0d0
        do j1=1,ns0
        do j2=1,ns0
        do ij=1,2
        if(ij==1)then
        j3=j2
        j4=j1
        else
        j3=j1
        j4=j2
        endif
        vq0=0.0d0
        do q1=-nq1,nq1
        do q2=-nq2, nq2
        q3=-q1
        q4=-q2
        cqx=(mnt(j1,1)-mnt(j4,1))/dfloat(kx0)*g1+q1*glist(1)
        cqx=cqx+(mnt(j1,2)-mnt(j4,2))/dfloat(ky0)*g2
        cqx=cqx+q2*glist(2)
        vq0=cdabs(cqx)  
        if(vq0.ge.1e-15)then
        if(dgate.gt.0.1)then
         vq0=tanh(dgate*vq0)/vq0
        else
         vq0=1.0d0/vq0
        endif
        if(ij==1)then
        val_hf(j2)=val_hf(j2)+2.0d0*vq0*fq0(j1,j4,q1,q2)*
     &     fq(j2,j3,q3,q4)
        val_hf(j2)=val_hf(j2)+2.0d0*vq0*fq00(j1,j4,q1,q2)*
     &     fq(j2,j3,q3,q4)
        else
        val_hf(j1)=val_hf(j1)-2.0d0*vq0*abs(fq10(j1,j4,q1,q2)**2)
        endif
        endif
        enddo
        enddo
        enddo 
        enddo
        enddo
        val_hf=val_hf*ee1/2.0d0

        do j1=1,ns0
        write(23,*)j1, val_hf(j1),val_hf1(j1)
        enddo


	kk=0
	vin=dcmplx(0.d0,0.d0)
	do j1=1,nst
	do j2=1,nst
	j12=mod(j1-1,ns0)+1+mod(j2-1,ns0)*ns0
	do j3=1,nst
	do j41=1,1!!!nst
	j4x=mod(mnt(j1,1)+mnt(j2,1)-mnt(j3,1)+kx0,kx0)
	j4y=mod(mnt(j1,2)+mnt(j2,2)-mnt(j3,2)+ky0,ky0)
	j4=ind_mnt(j4x, j4y)
	j34=mod(j3-1,ns0)+1+mod(j4-1,ns0)*ns0
	ji=j1+(j2-1)*nst+(j3-1)*nst*nst!!+(j4-1)*ns*ns*ns
	jk=j2+(j1-1)*nst+(j3-1)*nst*nst!!+(j4-1)*ns*ns*ns
	jk1=j1+(j2-1)*nst+(j4-1)*nst*nst!!+(j3-1)*ns*ns*ns
	jk0=j2+(j1-1)*nst+(j4-1)*nst*nst!!+(j3-1)*ns*ns*ns

	if(icomb.eq.1)vin(ji)=v(1,ji)!-v(1,jk)*(icomb-1)
	if(icomb.eq.2)vin(ji)=v(1,ji)-v(1,jk)
	if(icomb.eq.4)vin(ji)=v(1,ji)+v(1,jk0)-v(1,jk)-v(1,jk1)
	vin(-ji)=-vin(ji)
	enddo
	enddo
	enddo
	enddo
        return
        end 

        subroutine den_sq()!
	use interact
        complex*16 temp0,temp1
        real*8 eng(10)
        complex*16 c11, c12

        spin= dcmplx(0.d0, 0.d0)
        temps=0.0d0

        do 600 n2a=1,nn11
	n2=n2a
	nb1=newbb(n2a,1)
	nb2=newbb(n2a,2)
	ne1=numb(nb1,1)
	ne2=numb(nb2,1)
	if(numb(nb1,1)+numb(nb2,1).ne.ne)write(*,*)'numb wrong',ne,jcho
	do j=1,ns00
	nns1(j)=nnsb(nb1,j)
	nns1(j+ns00)=nnsb(nb2,j)
	enddo

	do j=1,ne
	ia=numb(nb1,1)
	ia1=numb(nb2,1)
	if(j.le.ia)nne(j)=nneb(nb1,j)
	if(j.gt.ia)nne(j)=nneb(nb2,j-ia)+ns00
	enddo

        do it=1,ns
        nns2(it)=nns1(it)
        end do

        do 55 i1=1,ne!!!!-1
        do 54 i2=1,ne
        j1=nne(i1)
        j2=nne(i2)  !!! this chosen an order for j1, j2
	if(j1.eq.j2)go to  54
	j12=mod(j1-1,ns0)+1+mod(j2-1,ns0)*ns0!!!ind_mnt(j1,j2) !!!j1+(j2-1)*ns
        nns2(j1)=0
        nns2(j2)=0
        do 53 j3=1,nst
        if(nns2(j3).eq.1) go to 53 !!! j3 can be the same as j1, or j2
	do 531 j41=1,1!!!nst
	j4x=mod(mnt(j1,1)+mnt(j2,1)-mnt(j3,1)+kx0,kx0)
	j4y=mod(mnt(j1,2)+mnt(j2,2)-mnt(j3,2)+ky0,ky0)
	j4=ind_mnt(j4x, j4y)
        if(nns2(j4).eq.1.or.j4.eq.j3) go to 532

        nns2(j3)=1
        nns2(j4)=1
	nb11=numr1(nb1)
	nb22=numr1(nb2)
	ne1=numb(nb1,1)
	ne2=numb(nb2,1)
	ne11=ne1!
	ne12=ne2
	
	if(j1.le.ns00)then
	nb11=nb11-ki1(j1)	
	ne11=ne11-1
	else
	nb22=nb22-ki1(j1-ns00)
	ne12=ne12-1
	endif
	if(j2.le.ns00)then
	nb11=nb11-ki1(j2)
	ne11=ne11-1
	else
	nb22=nb22-ki1(j2-ns00)
	ne12=ne12-1
	endif
	if(j3.le.ns00)then
	nb11=nb11+ki1(j3)	 
	ne11=ne11+1
	else
	nb22=nb22+ki1(j3-ns00)
	ne12=ne12+1
	endif
	if(j4.le.ns00)then
	nb11=nb11+ki1(j4)
	ne11=ne11+1
	else
	nb22=nb22+ki1(j4-ns00)
	ne12=ne12+1
	endif

	if(abs(ne11-ne1d).gt.numd.and.abs(ne11-ne2d).gt.numd)go to 6006
	nb11=numr(nb11)
	nb22=numr(nb22)
	ne10=numb(nb11,1)
	je1=numb(nb11,1)+numb(nb22,1)
	jm1=numb(nb11,2)+numb(nb22,4)

	if(je1.ne.ne.or.ne10.ne.ne11)then
	write(*,*)n2,'new wrong1',ne11,je1,nb11,nb22,je1,jm1
	write(*,*)j1,j2,j3,j4
	endif

	nb1a=newj0(nb11)
        jcf=njm(nb1a)+njm1(nb22)
	nb1x=newbb(jcf,1)
	nb2x=newbb(jcf,2)
	if(nb11.ne.nb1x.or.nb22.ne.nb2x)write(*,*)'new wrong',nb11,nb1x
     &   ,nb22,nb2x
	
	n1=jcf
	n1a=n1
	nd=j1+(j2-1)*nst+(j3-1)*nst*nst
        np11=0
	if(j1.gt.j2)np11=np11+1
	if(j4.gt.j3)np11=np11+1
	
        do it=min(j3,j4)+1,max(j3,j4)-1
        if(nns2(it).eq.1) np11=np11+1
        end do
        do it=min(j1,j2)+1,max(j1,j2)-1 
        if(nns1(it).eq.1) np11=np11+1
        end do
	nd11=nd
          if(n1.ne.0)then
        temp0=dcmplx(vec(n1),-vec(n1+nn11))*dcmplx(vec(n2),vec(n2+nn11))
        if(mod(np11,2).ne.0)temp0=-temp0
        temps(nd,1)=temps(nd,1)+temp0!!!!*(-1)**np1a
        temp0=dcmplx(vec(n1+2*nn11),-vec(n1+3*nn11))*dcmplx
     &   (vec(n2+2*nn11),vec(n2+3*nn11))
        if(mod(np11,2).ne.0)temp0=-temp0
        temps(nd,2)=temps(nd,2)+temp0!!!!*(-1)**np1a
        temp0=dcmplx(vec(n1+4*nn11),-vec(n1+5*nn11))*dcmplx
     &   (vec(n2+4*nn11), vec(n2+5*nn11))
        if(mod(np11,2).ne.0)temp0=-temp0
        temps(nd,3)=temps(nd,3)+temp0!!!!*(-1)**np1a
        endif

6006	continue
        nns2(j3)=0
        nns2(j4)=0
532	continue
531     continue	
53      continue
        nns2(j1)=1
        nns2(j2)=1
54      continue
55      continue

600	continue



        eng=0.0d0
        spin=0.0d0
	do j1=1,ns0
	do j2=1,ns0
	j12=j1+(j2-1)*ns0
	do j3=1,ns0
	do j41=1,1!!!ns0
	j4x=mod(mnt(j1,1)+mnt(j2,1)-mnt(j3,1)+kx0,kx0)
	j4y=mod(mnt(j1,2)+mnt(j2,2)-mnt(j3,2)+ky0,ky0)
	j4=ind_mnt(j4x, j4y)
	jj=j1+(j2-1)*ns0+(j3-1)*ns0**2
	n=jj
        do q1=-nq1,nq1
        do q2=-nq2, nq2
        qd1=mnt(j1,1)+mnt(j2,1)-mnt(j3,1)-mnt(j4,1)
        qd2=mnt(j1,2)+mnt(j2,2)-mnt(j3,2)-mnt(j4,2)
        qd1=qd1/kx0
        qd2=qd2/ky0
        cqx=(mnt(j1,1)-mnt(j4,1))/dfloat(kx0)*g1+q1*glist(1)
        cqx=cqx+(mnt(j1,2)-mnt(j4,2))/dfloat(ky0)*g2
        cqx=cqx+q2*glist(2)
        nq1a=0
        nq1b=0
        iaa=mnt(j1,1)-mnt(j4,1) !!!!+ns0+(q1+nq1)*ns0*2
        ibb=mnt(j1,2)-mnt(j4,2) !!!!+ns0+(q1+nq1)*ns0*2
        iaa=iaa+(q1*mg0(1,1)+q2*mg0(2,1))*kx0
        ibb=ibb+(q1*mg0(1,2)+q2*mg0(2,2))*ky0
        
        do ij=1, ma
        if(abs(iaa).le.ms0.and.abs(ibb).le.ms0)then
        spin(iaa,ibb,ij)=spin(iaa,ibb,ij)+vsq(n,iaa,ibb)*temps(n,ij)
        eng(ij)=eng(ij)+vsq1(n,iaa,ibb)*temps(n, ij)
        if(j1==j3)then
        spin(iaa, ibb,ij)=spin(iaa,ibb,ij)+vsq(n,iaa,ibb)*density(j1,ij)     !!!!!temps(n,iaa,ibb)
        endif
        endif
	enddo
	enddo
	enddo

	enddo
	enddo
	enddo
        end do


                spin=spin/(nx0*ny0)
        temp1=dcmplx(0.d0,0.d0)
        spin0=spin0+spin/2.0d0
        open(25,file='Sq.dat', access='append')
        do iy=-ms0, ms0
        do ix=-ms0, ms0
17       format(3i8, 9f18.6)
        if(kzone(ix,iy)==1)then
        write(25,17)ix,iy,ns0,gk1(ix,iy,1:2),real(spin(ix,iy,1:3))
     &     ,imag(spin(ix,iy,1:3))
        endif
        enddo
        enddo
        close(25)
        close(28)
        return
        end


	subroutine get_particle_ent(Na)	
	use interact 
	integer jj,k1,k2,lwork,kdim,info,i1,i2, jth, i3, i, j,js,ik,n_h
        integer is, na, nh1(nbasis0), nh, nh0
	complex*16  wave(nbasis0*4)
        REAL*8 aa, a00, a10, a20, a30
        integer*8, parameter :: nlanc1=nlanc*100
        integer  :: indx_c(nlanc1)
        integer pick(nlanc1), pick1(nlanc),nw1(nlanc), nw2
        integer index(nlanc,400), index1(nlanc, 400)
        integer index2(nlanc,400),index3(nlanc, 400)
        complex*16,allocatable::rho(:,:),rho1(:,:),rho2(:,:)
        real*8,allocatable, dimension(:)::den_val
        complex*16,allocatable, dimension(:)::den_work
         real*8, allocatable :: rwork(:)
        integer, allocatable ::iwork(:)
        integer nm1(2), nm2(2), kn(2), jm(nlanc)
        complex*16 c11, c22, c12
        external zheevd

        if(Na==3)then
        open(100+Na,file="PES_Na3.dat", access="append")
        endif
        if(Na==4)then
        open(100+Na,file="PES_Na4.dat", access="append")
        endif
        
        a00=dfloat(4*3*2)/ne/(ne-1)/(ne-2)/(ne-3)
        if(na==3)a00=dfloat(3*2)/ne/(ne-1)/(ne-2)
        a10=0.0d0
        a20=0.0d0
        a30=0.0d0
        aa=0.0d0
	entropy=0.d0
	entropy1=0.d0

        nw1=0
          kn=0
           jj=0
            jm=-1
        if(Na==4)then
        do j1=1,ns0-3
                do j2=j1+1, ns0-2
                do j3=j2+1, ns0-1
                do j4=j3+1, ns0
                       do ia=1,2
	kn(ia)=mnt(j1,ia)+mnt(j2,ia)+mnt(j3,ia)+mnt(j4,ia)
                        enddo
          kn(1)=mod(kn(1),kx0)
          kn(2)=mod(kn(2),ky0)
	ja=ind_mnt(kn(1), kn(2))
                                jj=jj+1
                                jm(jj)=ja
                                        enddo
                                        enddo
                                        enddo
                                        enddo
                        else
                        if(Na==3)then
        do j1=1,ns0-2
                do j2=j1+1, ns0-1
                do j3=j2+1, ns0
                       do ia=1,2
	kn(ia)=mnt(j1,ia)+mnt(j2,ia)+mnt(j3,ia)
                        enddo
          kn(1)=mod(kn(1),kx0)
          kn(2)=mod(kn(2),ky0)
	ja=ind_mnt(kn(1), kn(2))
                                jj=jj+1
                                jm(jj)=ja
                                        enddo
                                        enddo
                                        enddo
                        endif
                        endif


        do ik1=1, ns0
                jj1=0
                jjt=0
                pick=0
        pick1=0
                        if(Na==4)then
                  do j1=1, ns0-3
                  do j2=j1+1, ns0-2
                  do j3=j2+1, ns0-1
                  do j4=j3+1, ns0
                        jjt=jjt+1
         if(jm(jjt)==ik1)then
                        if(pick(jjt)==0)then
                 jj1=jj1+1
                              pick(jjt)=jj1  
                        endif
        endif
        enddo
        enddo
        enddo
        enddo

                else
                if(Na==3)then
                  do j1=1, ns0-2
                  do j2=j1+1, ns0-1
                  do j3=j2+1, ns0
                        jjt=jjt+1
         if(jm(jjt)==ik1)then
                        if(pick(jjt)==0)then
                 jj1=jj1+1
                              pick(jjt)=jj1  
                        endif
        endif
        enddo
        enddo
        enddo
                endif
                endif

                k10=mnt(ik1,1)
                k20=mnt(ik1,2)
        nh=0
        nh0=0
                index=0
                index1=0
                index2=0
                        kk=0
                                jjt=0
                                jj=0
        ikt=1
                do ika=1,ikt
              nns=0
        nns1=0
        
        if(Na==4)then
                do i1=1,ns0-9
                do i2=i1+1, ns0-8
                do i3=i2+1, ns0-7
                do i4=i3+1, ns0-6
                do i5=i4+1, ns0-5
                do i6=i5+1, ns0-4
                do i7=i6+1, ns0-3
                do i8=i7+1, ns0-2
                do i9=i8+1, ns0-1
                do i10=i9+1, ns0
                        do ia=1,2
	kn(ia)=mnt(i1,ia)+mnt(i2,ia)+mnt(i3,ia)+mnt(i4,ia)
        kn(ia)=kn(ia)+mnt(i5,ia)+mnt(i6,ia)+mnt(i7,ia)
        kn(ia)=kn(ia)+mnt(i8,ia)+mnt(i9,ia)+mnt(i10,ia)
                        enddo
                        nns=0
               nns(i1)=1
               nns(i2)=1
               nns(i3)=1
               nns(i4)=1
               nns(i5)=1
               nns(i6)=1
               nns(i7)=1
               nns(i8)=1
               nns(i9)=1
               nns(i10)=1
        nns1=nns
                nb10=1
                        nb20=1
                  do ia=1,ns00
         if(nns(ia)==1)then
        nb10=nb10+ki1(ia)         !   jr(3)=jr(1)+ki1(j3)+ki1(j4)
                endif
         if(nns(ia+ns00)==1)then
        nb20=nb20+ki1(ia)
        endif
                enddo

            kn(1)=mod(kn(1)+k10+ns0, nx0)
            kn(2)=mod(kn(2)+k20+ns0, ny0)
        if(ind_mnt(kn(1),kn(2))==kmn(ika))then
                        kk=kk+1

                jj=0
                jjt=0
                  do j1=1, ns0-3
                  do j2=j1+1, ns0-2
                  do j3=j2+1, ns0-1
                  do j4=j3+1, ns0
                        jjt=jjt+1
                        nns=nns1
         if(jm(jjt)==ik1)then
               if(nns(j1)+nns(j2)+nns(j3)+nns(j4)==0)then

                                nns(j1)=1
                                nns(j2)=1
                                nns(j3)=1
                                nns(j4)=1
                     nb11=nb10
                     nb22=nb20
           if(j1.le.ns00)then 
        nb11=nb11+ki1(j1)         
               else
        nb22=nb22+ki1(j1-ns00)
        endif
           if(j2.le.ns00)then
        nb11=nb11+ki1(j2)         
                else
        nb22=nb22+ki1(j2-ns00)
        endif
           if(j3.le.ns00)then 
        nb11=nb11+ki1(j3)         
                else
        nb22=nb22+ki1(j3-ns00)
        endif
           if(j4.le.ns00)then 
        nb11=nb11+ki1(j4)         
                else
        nb22=nb22+ki1(j4-ns00)
        endif
        nb11=numr(nb11)
        nb22=numr(nb22)
          ja1=newj0(nb11)
        nj=njm(ja1)+njm1(nb22)
        if(nj.ne.0)then
        nb1=newbb(nj,1)
        nb2=newbb(nj,2)
        if(nb1.ne.nb11.or.nb2.ne.nb22)write(*,*)'wrong2',
     &     nb1,nb11,nb2,nb22
        pick1(nj)=pick1(nj)+1
        endif
                        jj=jj+1
                        index(kk,jj)=nj+(ika-1)*nlanc
                index1(kk,jj)=pick(jjt)
         nsin=1
        if(j2.ne.j1+1)then
        do it=j1+1,j2-1
        if(nns(it).eq.1) nsin=-nsin
        end do
        endif
        if(j4.ne.j3+1)then
        do it=j3+1, j4-1
        if(nns(it).eq.1) nsin=-nsin
        end do
        endif
                index2(kk,jj)=nsin

                        
                                endif
                                endif
                                enddo
                                enddo
                                enddo
                                enddo
                        endif
                                enddo
                                nw1(kk)=jj
                        enddo
                        enddo
                        enddo
                        enddo
                        enddo
                        enddo
                        enddo
                        enddo
                        enddo
                        endif


	if(Na==3)then
                do i1=1,ns0-10
                do i2=i1+1, ns0-9
                do i3=i2+1, ns0-8
                do i4=i3+1, ns0-7
                do i5=i4+1, ns0-6
                do i6=i5+1, ns0-5
                do i7=i6+1, ns0-4
                do i8=i7+1, ns0-3
                do i9=i8+1, ns0-2
                do i10=i9+1, ns0-1
                do i11=i10+1, ns0
                        do ia=1,2
	kn(ia)=mnt(i1,ia)+mnt(i2,ia)+mnt(i3,ia)+mnt(i4,ia)
        kn(ia)=kn(ia)+mnt(i5,ia)+mnt(i6,ia)+mnt(i7,ia)
        kn(ia)=kn(ia)+mnt(i8,ia)+mnt(i9,ia)+mnt(i10,ia)
        kn(ia)=kn(ia)+mnt(i11,ia)
                        enddo
                        nns=0
               nns(i1)=1
               nns(i2)=1
               nns(i3)=1
               nns(i4)=1
               nns(i5)=1
               nns(i6)=1
               nns(i7)=1
               nns(i8)=1
               nns(i9)=1
               nns(i10)=1
               nns(i11)=1
        nns1=nns
                nb10=1
                        nb20=1
                  do ia=1,ns00
         if(nns(ia)==1)then
        nb10=nb10+ki1(ia)         !   jr(3)=jr(1)+ki1(j3)+ki1(j4)
                endif
         if(nns(ia+ns00)==1)then
        nb20=nb20+ki1(ia)
        endif
                enddo

            kn(1)=mod(kn(1)+k10+ns0, nx0)
            kn(2)=mod(kn(2)+k20+ns0, ny0)
        if(ind_mnt(kn(1),kn(2))==kmn(ika))then
                        kk=kk+1

                jj=0
                jjt=0
                  do j1=1, ns0-2
                  do j2=j1+1, ns0-1
                  do j3=j2+1, ns0
                        jjt=jjt+1
                        nns=nns1
         if(jm(jjt)==ik1)then
               if(nns(j1)+nns(j2)+nns(j3)==0)then
                                nns(j1)=1
                                nns(j2)=1
                                nns(j3)=1
                     nb11=nb10
                     nb22=nb20
           if(j1.le.ns00)then 
        nb11=nb11+ki1(j1)         
               else
        nb22=nb22+ki1(j1-ns00)
        endif
           if(j2.le.ns00)then
        nb11=nb11+ki1(j2)         
                else
        nb22=nb22+ki1(j2-ns00)
        endif
           if(j3.le.ns00)then 
        nb11=nb11+ki1(j3)         
                else
        nb22=nb22+ki1(j3-ns00)
        endif
        nb11=numr(nb11)
        nb22=numr(nb22)
          ja1=newj0(nb11)
        nj=njm(ja1)+njm1(nb22)
        if(nj.ne.0)then
        nb1=newbb(nj,1)
        nb2=newbb(nj,2)
        if(nb1.ne.nb11.or.nb2.ne.nb22)write(*,*)'wrong2',
     &     nb1,nb11,nb2,nb22
        pick1(nj)=pick1(nj)+1
        endif
                        jj=jj+1
                        index(kk,jj)=nj+(ika-1)*nlanc
                index1(kk,jj)=pick(jjt)
         nsin=1
        if(j2.ne.j1+1)then
        do it=j1+1,j2-1
        if(nns(it).eq.1) nsin=-nsin
        end do
        endif
        do it=1, j3-1
        if(nns(it).eq.1) nsin=-nsin
        end do
                index2(kk,jj)=nsin
                                endif
                                endif
                                enddo
                                enddo
                                enddo
                        endif
                                enddo
                                nw1(kk)=jj
                        enddo
                        enddo
                        enddo
                        enddo
                        enddo
                        enddo
                        enddo
                        enddo
                        enddo
                        enddo
                        endif

                        enddo

                                nw2=jj1
                                     nh=nw2
                          
        kj=nh
        if(kj.ne.0)then
	kdim=kj
        lrwork=4*kdim*kdim
	lwork=lrwork

	if(allocated(den_val))deallocate(den_val)
	if(allocated(den_work))deallocate(den_work)
	if(allocated(rho))deallocate(rho,rho1,rho2)
	if(allocated(rwork))deallocate(rwork, iwork)
        allocate(den_val(kdim),den_work(lrwork),rho(nh,nh))
        allocate(rho1(nh,nh),rho2(nh,nh))
        allocate(rwork(lwork),iwork(lrwork))

	rho=0.d0
        rho1=0.0d0
        rho2=0.0d0
        den_val=0

                do ka=1, kk
                                do ia=1,nw1(ka)
                                k1=index(ka,ia)
        im=(k1-1)/nlanc+1
        k1=mod(k1-1,nlanc)+1
        ma1=ma
        do ij=0,ma1-1
         c11=dcmplx(vec(k1+ij*2*nn11),-vec(k1+(ij*2+1)*nn11))
                                do ja=1, nw1(ka)
                                k2=index(ka,ja)
        k2=mod(k2-1,nlanc)+1
         c22=dcmplx(vec(k2+ij*2*nn11),vec(k2+(ij*2+1)*nn11))
                                ia1=index1(ka, ia)
                                ia2=index1(ka, ja)
        
        if(index2(ka,ja).ne.index2(ka,ia))c22=-c22
              rho(ia1,ia2)=rho(ia1,ia2)+c11*c22
                                  enddo
                                        enddo
                                                enddo
        enddo


	!!! now we have reduced density matrix!!
	!! get the eigen values for density matrix
	rho=-rho/ikt/ma
        do i=1,nh
        a10=a10+rho(i,i)
        enddo

        call ZHEEV('V','U',nh,rho,nh,den_val,den_work,lwork,rwork,
     &    lrwork, iwork, lrwork,info)


        do i=1,kdim
                if(dabs(den_val(i))>=1.0e-15) then
              entropy=entropy+dlog(dabs(den_val(i)))*dabs(den_val(i))
              aa=aa+den_val(i)
                endif
                write(100+Na,29)i,-den_val(i)*a00,
     &     -dlog(dabs(den_val(i)*a00)),ik1,ns0,Na
29      format(i6, 2f18.9, 3i6)
        enddo
112	format(2i6,6f18.6)

        endif
        enddo
        
        close(100+Na)
	return
	end
	


	subroutine eigen1(n,a,val,vec)
	use lanczos1    !!paramsizew
      parameter (n3=nnp1)!!!*2)!!!!!nn)                 !!!!!!  n3==n=nx*nx
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

c	write(*,*)'before'
      UPLO='U'
      JOB='V'
      IF (UPLO.EQ.'U') THEN
      END IF
c      CALL cheevd (JOB,UPLO,N,A,LDA,val,WORK,LWORK,RWORK,LRWORK,IWORK,

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

	function ran2(idum)
	integer idum, im1,im2,imm1,ia1,ia2,iq1,iq2,ir1,ir2,ntab,ndiv
	real ran2,am,eps,rnmx
	parameter (im1=2147483563,im2=2147483399,am=1./im1,imm1=im1-1,
     &  ia1=40014,ia2=40692,iq1=53668,iq2=52774,ir1=12211,
     &  ir2=3791,ntab=32,ndiv=1+imm1/ntab,eps=1.2e-7,rnmx=1.-eps)

	integer idum2,j,k,iv(ntab),iy
	save iv,iy,idum2
	data idum2/123456789/,iv/ntab*0/,iy/0/

	if(idum.le.0)then
	idum=max(-idum,1)
	idum2=idum
	do 11 j=ntab+8,1,-1
	k=idum/iq1
	idum=ia1*(idum-k*iq1)-k*ir1
	if(idum.lt.0)idum=idum+im1
	if(j.le.ntab)iv(j)=idum
11	continue

	iy=iv(1)
	endif

	k=idum/iq1
	idum=ia1*(idum-k*iq1)-k*ir1
	if(idum.lt.0)idum=idum+im1
	k=idum2/iq2
	idum2=ia2*(idum2-k*iq2)-k*ir2
	if(idum2.lt.0)idum2=idum2+im2
	j=1+iy/ndiv
	iy=iv(j)-idum2
	iv(j)=idum
	if(iy.lt.1)iy=iy+imm1
	ran2=min(am*iy,rnmx)
	return
	end

	subroutine getvec()
	use interact
        integer*8 n
        real*8 ph, ph1, ph00(120), a1
	complex*16  v1a, v2a, v3a, v4a, vcom(4)
	complex*16 com1,com2,com3,com,com4, com11, com12,com13,com10
        hop=1.0d0

	krr=0
	km=1
        call coulomb()
        if(krr.eq.0)call CMATV1()

	ng1=ma
	iseed=99982221
	nty=2
        n=nbasisf
	call eivalnew(n,ng1,val,vec,iseed,nty)
        
        aa1=0.0d0
        density=0.0d0
             open(14,file="density.dat", access='append')
	do n1a=1,n
	nb1=newbb(n1a,1)
	nb2=newbb(n1a,2)

	do j=1,ns00
	j11=nnsb(nb1,j)
	j12=nnsb(nb2,j)
        aa1=vec(n1a)**2+vec(n1a+n)**2
        if(j11.eq.1)density(j,1)=density(j,1)+aa1
        if(j12.eq.1)density(j+ns00,1)=density(j+ns00,1)+aa1
        aa1=vec(n1a+2*n)**2+vec(n1a+n+2*n)**2
        if(j11.eq.1)density(j,2)=density(j,2)+aa1
        if(j12.eq.1)density(j+ns00,2)=density(j+ns00,2)+aa1
        aa1=vec(n1a+4*n)**2+vec(n1a+n+4*n)**2
        if(j11.eq.1)density(j,3)=density(j,3)+aa1
        if(j12.eq.1)density(j+ns00,3)=density(j+ns00,3)+aa1
	enddo
        enddo
        a1=tw_angle*180.0d0/pi
        aa1=0.0d0
        do j=1, ns0
        write(14,51)Ns0, Ne, j,a1,density(j,1:3)
        aa1=aa1+density(j,1)+density(j,2)+density(j,3)
        enddo
        write(14,*)aa1, jcho(1:2),'kx, ky'
51      format(3i6, 8f20.12)
        close(14)


        if(jcho(1)==nx0/2.and.jcho(2)==0)then
        if(nphase.eq.1)call den_sq()
              kmn(1)=ind_mnt(jcho(1),jcho(2))
        if(nphase.eq.1)call get_particle_ent(3)
        if(nphase.eq.1)call get_particle_ent(4)
        endif
        
        vecph(1:n*2*ng1,kp1, kp2)=vec(1:2*n*ng1) 

	do i=1,ng1
11      format(2i6, 18f19.8)
	enddo

	aa=0.0
        if(kp1+kp2.eq.2*nphase.and.nphase.ne.1)then
        open(92,file='Chern_num.dat', access='append')
        open(91,file='Chern_curve.dat', access='append')
        do ii=0, (ng1-1)*2,2
        ph00=0.0
        ph1=0.0
        do k1=1,nphase-1
        do k2=1,nphase-1
        nj=k2+(k1-1)*nphase
        com=0.0
        com1=0.0
        com2=0.0
        com3=0.0 
        com4=0.0d0
        
             do j=1,n  
            v1a=dcmplx(vecph(j+n*ii,k1,k2),vecph(j+n*ii+n,k1,k2))
            v2a=dcmplx(vecph(j+n*ii,k1+1,k2),vecph(j+n*ii+n,k1+1,k2))
         v3a=dcmplx(vecph(j+n*ii,k1+1,k2+1),vecph(j+n*ii+n,k1+1,k2+1))
            v4a=dcmplx(vecph(j+n*ii,k1,k2+1),vecph(j+n*ii+n,k1,k2+1))
                n1a=j
        nb1=newbb(n1a,1)
        nb2=newbb(n1a,2)
        vcom=1.0d0
        do ja=1,ns00
        j11=nnsb(nb1,ja)
        j12=nnsb(nb2,ja)
        if(j11.eq.1)then
        com10=0.0d0
        com11=0.0d0
        com12=0.0d0
        com13=0.0d0
        do jk=1, nnp
        com10=com10+dconjg(vecb(jk,ja,k1,k2))*vecb(jk,ja,k1+1,k2)
        com11=com11+dconjg(vecb(jk,ja,k1+1,k2))*vecb(jk,ja,k1+1,k2+1)
        com12=com12+dconjg(vecb(jk,ja,k1+1,k2+1))*vecb(jk,ja,k1,k2+1)
        com13=com13+dconjg(vecb(jk,ja,k1,k2+1))*vecb(jk,ja,k1,k2)
        enddo
        vcom(1)=vcom(1)*com10
        vcom(2)=vcom(2)*com11
        vcom(3)=vcom(3)*com12
        vcom(4)=vcom(4)*com13
        endif
        if(j12.eq.1)then
        com10=0.0d0
        com11=0.0d0
        com12=0.0d0
        com13=0.0d0
        do jk=1, nnp
        jb=ja+ns00
        com10=com10+dconjg(vecb(jk,jb,k1,k2))*vecb(jk,jb,k1+1,k2)
        com11=com11+dconjg(vecb(jk,jb,k1+1,k2))*vecb(jk,jb,k1+1,k2+1)
        com12=com12+dconjg(vecb(jk,jb,k1+1,k2+1))*vecb(jk,jb,k1,k2+1)
        com13=com13+dconjg(vecb(jk,jb,k1,k2+1))*vecb(jk,jb,k1,k2)
        enddo
        vcom(1)=vcom(1)*com10
        vcom(2)=vcom(2)*com11
        vcom(3)=vcom(3)*com12
        vcom(4)=vcom(4)*com13
        endif
        enddo
        vcom=dconjg(vcom)

        com=dconjg(v1a)*v2a*vcom(1)+com
        com1=dconjg(v2a)*v3a*vcom(2)+com1
        com2=dconjg(v3a)*v4a*vcom(3)+com2
        com3=dconjg(v4a)*v1a*vcom(4)+com3
        com4=dconjg(v4a)*v4a+com4
        enddo
        com=com*com1*com2*com3
21      format(2i6, 2f18.9)
        ph00(nj)=imag(cdlog(com))
        if(ph00(nj).gt.pi)then
        ph00(nj)=ph00(nj)-2*pi
        endif
        if(ph00(nj).le.-pi)then
        ph00(nj)=ph00(nj)+2*pi
        endif
        ph1=ph1+ph00(nj)
        write(91,23)ii/2,k1,k2,nj,ph1, ph00(nj), tw1
23      format(4i6, 5f18.9)
        enddo
        enddo
        write(92,24)tw1,ii/2,jcho(1:2),ph1
24      format(f18.6,3i6, 2f18.9)
        enddo
        close(91)
        close(92)


        endif

	return
	end

	!!! paragram  main
	use chern0
        use tmd
        
	complex*16 ccc1
	icomb=4
        chernband=.true.
        ckshift=0.0d0
        
        tw_angle=1.800d0/180.0d0*pi
        epsilon=2.0d0
        dgate=0.0d0
        VD1=0.0d0
        V1=17.50d0
        W1=-6.50d0
        V2=-10.50
        W2=11.750d0
        phi=-1.020843d0
        call band_spectrum_para()

        open(520,file= 'Energy.dat',access='append')
	do 1004 iter2= 0, ny0-1
	do 1004 iter1= 0, nx0-1

        
        do 1004 kp1=1,  nphase
        do 1004 kp2=1,  nphase
        ckshift=(kp1-1)*g1/nx0+(kp2-1)*g2/ny0
        if(nphase.ne.1)then
        ckshift=ckshift/(nphase-1)
        endif
        call band_spectrum()
        call coulomb()

	jcho(1)=iter1
	jcho(2)=iter2

	kin=1
        mmax=1
	do iter=1,mmax

        krr=0
	call getvec()

	write(520,122)ne,jcho(1:2),kp1,kp2,val(1:ng1),V1,W1,V2,W2
122     format(5i5, 15f21.7)
	enddo


1009	continue
1000    continue
1004 	continue
        close(520)
	end




c-------------------------------------------------------------------

	subroutine CMATV1()!
	use interact
        complex*16 temp0,temp1
        integer jr(4), jr1(4)
        if(krr.eq.0)then  !!!! no reading or using the memory
	ham(1:nn)=0.d0!-1000.d0
        kjjs=0
        kjjs1=0
        hdia=0.0d0

	jc=0
!!!! for 16 electron problem on 32 sites. Separate sites to two
!!!! groups each has ns1=16 states.
	nnsb=0
	nneb=0

      if(allocated(nvin))then  !!,work(lwork),rwork(lwork))
      deallocate(nvin,label)
     	endif 
	allocate(nvin(ksave1),label(1,ksave1))
	jc=0
	do i=1,ns00
	ki(i)=0
	ki1(i)=2**(i-1)
	enddo

	i0=0
	i0s=0
	i0s(1:ns00)=1

	do i24=0,i0s(24)
	ki(24)=0
	if(i24.eq.1)ki(24)=ki1(24)
	do i23=0,i0s(23)
	ki(23)=0
	if(i23.eq.1)ki(23)=ki1(23)
	do i22=0,i0s(22)
	ki(22)=0
	if(i22.eq.1)ki(22)=ki1(22)
	do i21=0,i0s(21)
	ki(21)=0
	if(i21.eq.1)ki(21)=ki1(21)
	do i20=0,i0s(20)
	ki(20)=0
	if(i20.eq.1)ki(20)=ki1(20)
	do i19=0,i0s(19)
	ki(19)=0
	if(i19.eq.1)ki(19)=ki1(19)
	do i18=0,i0s(18)
	ki(18)=0
	if(i18.eq.1)ki(18)=ki1(18)
	do i17=0,i0s(17)
	ki(17)=0
	if(i17.eq.1)ki(17)=ki1(17)

	do i16=0,i0s(16)
	ki(16)=0
	if(i16.eq.1)ki(16)=ki1(16)
	do i15=0,i0s(15)
	ki(15)=0
	if(i15.eq.1)ki(15)=ki1(15)

	do i14=0,i0s(14)
	ki(14)=0
	if(i14.eq.1)ki(14)=ki1(14)
	do i13=0,i0s(13)
	ki(13)=0
	if(i13.eq.1)ki(13)=ki1(13)
	do i12=0,i0s(12)
	ki(12)=0
	if(i12.eq.1)ki(12)=ki1(12)
	do i11=0,i0s(11)
	ki(11)=0
	if(i11.eq.1)ki(11)=ki1(11)
	do i10=0,i0s(10)
	ki(10)=0
	if(i10.eq.1)ki(10)=ki1(10)
	do i9=0,i0s(9)
	ki(9)=0
	if(i9.eq.1)ki(9)=ki1(9)
	do i8=0,i0s(8)
	ki(8)=0
	if(i8.eq.1)ki(8)=ki1(8)
	do i7=0,i0s(7)
	ki(7)=0
	if(i7.eq.1)ki(7)=ki1(7)
	do i6=0,i0s(6)
	ki(6)=0
	if(i6.eq.1)ki(6)=ki1(6)
	do i5=0,i0s(5)
	ki(5)=0
	if(i5.eq.1)ki(5)=ki1(5)
	do i4=0, i0s(4)
	ki(4)=0
	if(i4.eq.1)ki(4)=ki1(4)
	do i3=0, i0s(3)
	ki(3)=0
	if(i3.eq.1)ki(3)=ki1(3)
	do i2=0,i0s(2)
	ki(2)=0
	if(i2.eq.1)ki(2)=ki1(2)
	do i1=0,i0s(1)
	ki(1)=0
	if(i1.eq.1)ki(1)=ki1(1)
	ni1=i1+i2+i3+i4+i5+i6+i7+i8+i9+i10+i11+i12+i13+i14+i15+i16
     &   +i17+i18+i19+i20+i21+i22+i23+i24

	jm1=0
	jm2=0
	jc1=1
	jm3=0
	jm4=0
	do i=1,ns00
	jc1=jc1+ki(i)
	if(ki(i).ne.0)then
        if(chernband)then
	jm1=jm1+mnt(i,1) !! momentum between 0 and kx0-1
	jm2=jm2+mnt(i,2)
	jm3=jm3+mnt(i+ns00,1) !! momentum between 0 and kx0-1
	jm4=jm4+mnt(i+ns00,2)
        else
	jm1=jm1+i !! momentum between 0 and kx0-1
	jm2=jm1
	jm3=jm3+i+ns00 !! momentum between 0 and kx0-1
	jm4=jm3
        endif 
        endif
	enddo
        if(.not.chernband)then
        jm1=mod(jm1,ns0)
        jm2=mod(jm2,ns0)
        jm3=mod(jm3,ns0)
        jm4=mod(jm4,ns0)
        endif
        

	if(abs(ni1-ne1d).le.numd.or.abs(ni1-ne2d).le.numd)then
	if(ni1.le.ne)then
	jc=jc+1
	numb(jc,1)=ni1   ! number of electrons
        if(chernband)then
	numb(jc,2)=mod(jm1,kx0)    ! total momentum-x
	numb(jc,3)=mod(jm2,ky0)   ! total momentum for up states
	numb(jc,4)=mod(jm3,kx0)   ! total momentum for up states
	numb(jc,5)=mod(jm4,ky0)   ! total momentum for up states
        else
	numb(jc,2)=jm1    ! total momentum-x
	numb(jc,3)=jm2   ! total momentum for up states
	numb(jc,4)=jm3   ! total momentum for up states
	numb(jc,5)=jm4   ! total momentum for up states
        endif
	numr(jc1)=jc
	numr1(jc)=jc1

	nns(1)=i1
	nns(2)=i2
	nns(3)=i3
	nns(4)=i4
	nns(5)=i5
	nns(6)=i6
	nns(7)=i7
	nns(8)=i8
	nns(9)=i9
	nns(10)=i10
	nns(11)=i11
	nns(12)=i12
	nns(13)=i13
	nns(14)=i14

	nns(15)=i15
	nns(16)=i16
	nns(17)=i17
	nns(18)=i18
	nns(19)=i19
	nns(20)=i20
	nns(21)=i21
	nns(22)=i22
	nns(23)=i23
	nns(24)=i24

	do j=1,ns00
	nne(j)=0
	enddo
	ia=0
	do ie=1,ns00
	if(nns(ie).ne.0)then
	ia=ia+1
	nne(ia)=ie
	endif
	enddo
	

	nb10=jc
	jt=0

	ne11=ni1
	do j=1,ns00
	nnsb(nb10,j)=nns(j)
	if(j.le.ne11)nneb(nb10,j)=nne(j)
	enddo

	
	endif
	endif

           if(jc1==nbasis01)then
        nbasis0n1=jc
        endif


	enddo
	enddo
	enddo
	enddo
	enddo
	enddo
	enddo
	enddo

	enddo
	enddo
	enddo
	enddo
	enddo
	enddo
	enddo
	enddo

	enddo
	enddo
	enddo
	enddo

	enddo
	enddo
	enddo
	enddo
! figure out the representation for translation operator
	nbasis0n=jc
	do j=1,nbasis0n
	newj0(j)=0
	enddo

	jc=0
	njm(1)=0
	njm(2)=0
	do nb10=1,nbasis0n
	njm(nb10)=0
	enddo

	ja1=0
	do  nb10=1,nbasis0n   !nbasis0
	jc1=0
	ne1=numb(nb10,1)
	if(ne1.le.ns0.and.ne1.le.ne)then!!!!!eq.ne1)then
	do nb20=1,nbasis0n1        !    nb10-1
	ne2=numb(nb20,1)
	je=ne1+ne2
	if(je.ne.ne)go to 6002
	jm=numb(nb10,2)+numb(nb20,4)
	jm=mod(jm,kx0)
	if(jm.ne.jcho(1).and.jchonx.eq.1)go to 6002
	jm=numb(nb10,3)+numb(nb20,5)
	jm=mod(jm,ky0)

	if(jm.ne.jcho(2).and.jchony.eq.1)go to 6002
	!if(ne-ne1.gt.numd.and.jcut.eq.1)go to 6002
	if(abs(ne1-ne1d).gt.numd.or.abs(ne2-ne2d).gt.numd)go to 6002
	
     	jc=jc+1
	jc1=jc1+1
	newbb(jc,1)=nb10
	newbb(jc,2)=nb20
6002	continue
	enddo    ! end of loop for nb20
	if(jc1.ne.0)then
	ja1=ja1+1
	newj0(nb10)=ja1
	if(nb10.ne.nbasis0)njm(ja1+1)=jc
	endif    
	endif   !!!!! gt.  ns0-nq0
	enddo
!!!!! including they are the same since there is no symmetry there
	nbasisf=jc
         write(*,*)jc,'new ham  basis Hamilt dimens'


!!!!! this part tells njm1, given nb10,nb20, one could find back to jc
	njm1(1)=1
	do nb2=2,nbasis0n1
	njm1(nb2)=1
	je1=numb(nb2,1)
	jm1=numb(nb2,4)
	km1=numb(nb2,5)
	do nb21=1,nb2-1
	je=numb(nb21,1)
	jm=numb(nb21,4)
	km=numb(nb21,5)
	if(je1.eq.je)then!!!!.and.jm1.eq.jm)then
	if(jm1.eq.jm.or.jchonx.eq.0)then!!!!.and.jm1.eq.jm)then
	if(km1.eq.km.or.jchony.eq.0)then!!!!.and.jm1.eq.jm)then
	njm1(nb2)=njm1(nb2)+1
	endif
	endif
	endif
	enddo
	enddo

	do jc=1,nbasisf    !nbasis
	nb10=newbb(jc,1)
	nb20=newbb(jc,2)
	ja1=newj0(nb10)
	nj=njm(ja1)+njm1(nb20)
	if(nj.ne.jc)write(*,*)'wrong nj matches jc?',jc,nj
	enddo

	nn11=nbasisf

	kjj=0
        kjj1=0
	kjjt=0
        kjjt1=0
        do 600 n2a=1,nn11
	kjjs1(n2a)=kjjs1(n2a-1)+kjjs(n2a-1)
	kjjs2(n2a)=kjjs2(n2a-1)+kjjss(n2a-1)
	n2=n2a
	nb1=newbb(n2a,1)
	nb2=newbb(n2a,2)
	ne1=numb(nb1,1)
	ne2=numb(nb2,1)
	if(numb(nb1,1)+numb(nb2,1).ne.ne)write(*,*)'numb wrong',ne,jcho
	do j=1,ns00
	nns1(j)=nnsb(nb1,j)
	nns1(j+ns00)=nnsb(nb2,j)
	enddo


	do j=1,ne
	ia=numb(nb1,1)
	ia1=numb(nb2,1)
	if(j.le.ia)nne(j)=nneb(nb1,j)
	if(j.gt.ia)nne(j)=nneb(nb2,j-ia)+ns00
	enddo


        do it=1,ns
        nns2(it)=nns1(it)
        end do

	do i1=1,ne !!! looking at qn-k of each electrons
	j1=nne(i1)
	ham(n2)=ham(n2)+hop*(valf(j1)+val_hf(j1))  !!!!!!!eng_dia(j1)*0

	enddo

           
	nb11=numr1(nb1)
	nb22=numr1(nb2)

        do 55 i1=1,ne!!!!-1
        do 54 i2=1,ne
        j1=nne(i1)
        j2=nne(i2)  !!! this chosen an order for j1, j2
	if(j1.eq.j2)go to  54
	j12=mod(j1-1,ns0)+1+mod(j2-1,ns0)*ns0!!!ind_mnt(j1,j2) !!!j1+(j2-1)*ns
        nns2(j1)=0
        nns2(j2)=0
        do 53 j3=1,nst
        if(nns2(j3).eq.1) go to 53 !!! j3 can be the same as j1, or j2
        !j4x=mod(mnt(j1,1)+mnt(j2,1)-mnt(j3,1)+kx0,kx0)
        !j4y=mod(mnt(j1,1)+mnt(j2,1)-mnt(j3,1)+kx0,kx0)
	!!! j4=ind_mnt(j4x, j4y)
	do 531 j41=1,1!!!nst
	j4x=mod(mnt(j1,1)+mnt(j2,1)-mnt(j3,1)+kx0,kx0)
	j4y=mod(mnt(j1,2)+mnt(j2,2)-mnt(j3,2)+ky0,ky0)
	j4=ind_mnt(j4x, j4y)
        if(.not.chernband)j4=mod(j1+j2-j3+ns0-1, ns0)+1
	j4r=(j4-1)/ns0
	j34=mod(j3-1,ns0)+1+mod(j4-1,ns0)*ns0!!!ind_mnt(j1,j2) !!!j1+(j2-1)*ns
        if(nns2(j4).eq.1.or.j4.eq.j3) go to 532
        if(j4.ge.j3.and.icomb.eq.4) go to 532
        if(j1.ge.j2.and.icomb.ge.2) go to 532
        nns2(j3)=1
        nns2(j4)=1
	nb11=numr1(nb1)
	nb22=numr1(nb2)
	ne1=numb(nb1,1)
	ne2=numb(nb2,1)
	ne11=ne1!
	ne12=ne2
	
	if(j1.le.ns00)then
	nb11=nb11-ki1(j1)	!     jr(3)=jr(1)+ki1(j3)+ki1(j4)
	ne11=ne11-1
	else
	nb22=nb22-ki1(j1-ns00)
	ne12=ne12-1
	endif
	if(j2.le.ns00)then
	nb11=nb11-ki1(j2)	 !    jr(3)=jr(1)+ki1(j3)+ki1(j4)
	ne11=ne11-1
	else
	nb22=nb22-ki1(j2-ns00)
	ne12=ne12-1
	endif
	if(j3.le.ns00)then
	nb11=nb11+ki1(j3)	  !   jr(3)=jr(1)+ki1(j3)+ki1(j4)
	ne11=ne11+1
	else
	nb22=nb22+ki1(j3-ns00)
	ne12=ne12+1
	endif
	if(j4.le.ns00)then
	nb11=nb11+ki1(j4)	  !   jr(3)=jr(1)+ki1(j3)+ki1(j4)
	ne11=ne11+1
	else
	nb22=nb22+ki1(j4-ns00)
	ne12=ne12+1
	endif

	!if(ne-ne11.gt.numd.or.go to 6006
	if(abs(ne11-ne1d).gt.numd.and.abs(ne11-ne2d).gt.numd)go to 6006
	nb11=numr(nb11)
	nb22=numr(nb22)
	ne10=numb(nb11,1)
	je1=numb(nb11,1)+numb(nb22,1)
	jm1=numb(nb11,2)+numb(nb22,4)

	if(je1.ne.ne.or.ne10.ne.ne11)then!!!!.or.mod(jm1,ns0).ne.0)then
	write(*,*)'new wrong1',ne11,je1,nb11,nb22,je1,jm1
	endif

	nb1a=newj0(nb11)
        jcf=njm(nb1a)+njm1(nb22)
	nb1x=newbb(jcf,1)
	nb2x=newbb(jcf,2)
	if(nb11.ne.nb1x.or.nb22.ne.nb2x)write(*,*)'new wrong',nb11,nb1x
     &   ,nb22,nb2x
	
	!nb11=nb1
	!nb22=nb2
	!!!!n1a=newb1(n1)
        !!jcf1=njm(nb11)+njm1(nb22)
	n1=jcf
	n1a=n1
	nd=j1+(j2-1)*nst+(j3-1)*nst*nst!!+(j4-1)*ns*ns*ns
        !!temp0=temp0-v(j4r+1,nd)
        np11=0
	if(j1.gt.j2)np11=np11+1
	if(j4.gt.j3)np11=np11+1
	
!!!!!!!!!!	if(j4.lt.j3)np1=1
        do it=min(j3,j4)+1,max(j3,j4)-1
        if(nns2(it).eq.1) np11=np11+1
        end do
        do it=min(j1,j2)+1,max(j1,j2)-1 !!!!!!!!!!!j2+1,j1  !!!!+1,j2
        if(nns1(it).eq.1) np11=np11+1
        end do
        !!!temp0=2.0d0*dble((-1)**np)*temp0
	!!!!!!!!!!if(j4r.eq.0)np0=np1
	nd=nd*(-1)**np11!j1+(j2-1)*ns+(j3-1)*ns*ns
	nd11=nd+j4r*2*ndimen1*0
	temp0=vin(nd11)
ccccccccccccccccccccccccccccccccccccccccc        b(n2)=b(n2)+temp*a(n1)
ccc        hr(n1,n2)=hr(n1,n2)+dreal(temp)
ccc        hi(n1,n2)=hi(n1,n2)+dimag(temp)

	if(n1.eq.n2)then
	ham(n1)=ham(n1)+temp0
	endif

        if(n1.gt.n2)then
	if(cdabs(temp0).gt.1e-9)then
	kjj=kjj+1
        kjjt=kjjt+1
	if(kjj.le.ksave1)then
	kjjs(n2a)=kjjs(n2a)+1
	label(1,kjj)=n1a
	!!!!label(2,kjj)=0!n1a
	nvin(kjj)=nd11  !!!!(-1)**np*(nd1+(nd-1)*ndim1)
	endif
        end if
	endif

	ks=0
	if(kjj.eq.-ksave1)then
	ks=ks+1
	kjj=0
	endif

6006	continue
        nns2(j3)=0
        nns2(j4)=0
532	continue
531     continue	
	!!!!enddo !this loop is for multiple j4--reduced conservation law
53      continue
        nns2(j1)=1
        nns2(j2)=1
54      continue
55      continue

        
        do it=1,ns
        nns2(it)=nns1(it)
        end do

600	continue
	kjjs1(nn11+1)=kjjs1(nn11)+kjjs(nn11)
	if(kjj.ne.0.and.ks.ne.0)then
	ks=ks+1
	endif


        nni=min(nn,10)
	write(*,*)'ham-count',jcho,nn11,kjjt,kjjt1,ham(1:nni)
	write(*,*)kjjt, ksave1, 'kjjt-ksave1'
                if(ksave1.lt.kjjt)stop
                
	kjj=0
	ks=0



	endif

c-------------------------------------------
c Calculate Vim*A
c-------------------------------------------
c-------------------------------------------
902 	continue
        return
        end


c-------------------------------------------------------------------
	subroutine MATVEC(a,b,n, ntype, add, acc)
	use interact
        integer*8 n, n2, n1, im, i
        double precision a(n,ntype),b(n,ntype),sum,a1ij,a2ij,zero
        double complex acc
        logical add
        parameter (zero=0.0d0)
         complex*16 aa
	complex*16 temp0
        if(.not.add)then
        do i=1,n
        b(i,1:ntype)=zero
        enddo
        endif
        do n2=1,n
        b(n2,1)=b(n2,1)+ham(n2)*a(n2,1)
        if(ntype.eq.2)b(n2,2)=b(n2,2)+ham(n2)*a(n2,2)
        end do
	kjj0=0
	do n1=1,n
	do kjj=kjjs1(n1)+1,kjjs1(n1+1)
         im=kjj
        i=n1
        j1=label(1,im)
	ji=nvin(kjj)
	!ji=nvin(kjj)
	if(j1.ne.0)then
	temp0=vin(ji)
	a1ij=dreal(temp0)
	a2ij=dimag(temp0)
	jcff=j1
	jca=i
	b(jcff,1)=b(jcff,1)+a1ij*a(jca,1)
        if(ntype.eq.2)then
        b(jcff,1)=b(jcff,1)-a2ij*a(jca,2)
        b(jcff,2)=b(jcff,2)+a1ij*a(jca,2)!!+a2ij*a(1,jca)
        b(jcff,2)=b(jcff,2)+a2ij*a(jca,1)
        endif
        if(jcff.ne.jca)then
        a2ij=-a2ij
        b(jca,1)=b(jca,1)+a1ij*a(jcff,1)!-a2ij*a(2,jcff)
        if(ntype.eq.2)then
        b(jca,1)=b(jca,1)-a2ij*a(jcff,2)
        b(jca,2)=b(jca,2)+a1ij*a(jcff,2)!+a2ij*a(1,jcff)
        b(jca,2)=b(jca,2)+a2ij*a(jcff,1)
        endif
	endif
	endif
	enddo
	enddo
      continue
      return
      end

c------------------------------------------------------------------
	subroutine eivalnew(nbasis,nev,val,vec,iseed,ctype)
      integer*8 NBASIS
        integer  i, ii
        integer CTYPE,abtype,nev,iseed
      logical ADD
     	double precision val(nev),vec(2*nbasis*nev) 
      
	abtype=0
	i=0
	ii=0
	
	call eivalha(nbasis,nev,ctype,abtype,i,ii,iseed,val,vec)
	return
	end

      subroutine MATVEC11(VIN,VOUT,NBASIS,ctype,ADD,A)
      integer*8 NBASIS
        integer CTYPE
      logical ADD
      double precision VIN(NBASIS,CTYPE), VOUT(NBASIS,CTYPE)
      integer i,j
      double precision  ar, ai, zero
      parameter (zero = 0.0d0)
      
      double complex A, ac
      external A
      
         do i = 1,NBASIS
            if(.not.ADD) then
               VOUT(i,1) = zero
               VOUT(i,2) = zero
            endif         

            do j = i,i   !1,NBASIS
               ac = A(i,j)
               ar = dreal(ac)
               ai = dimag(ac)
               VOUT(i,1) = VOUT(i,1) + ar*VIN(j,1) - ai*VIN(j,2) 
               VOUT(i,2) = VOUT(i,2) + ar*VIN(j,2) + ai*VIN(j,1)
            enddo
         enddo
      
      return
      end
      
      
	subroutine eivalha(nbasis,nev,ctype,abtype,i,ii,iseed,eval,evec)
      implicit none
      integer*8 nbasis, i1
      integer ctype,abtype,i,ii,nev,maxvec,iseed,j,maxstep,nstep
      integer nevec,id,laninfo, ierr,npolish
        
      parameter (maxvec=6 , npolish = 5)
      integer*8, parameter ::maxdim =20000000
! nbasis == < maxdim    nev == < maxvec
      logical internal, limit, relative

        integer*8, parameter :: maxdim1=2*maxdim
        integer*8, parameter :: maxdim2=2*4*maxdim
      double precision zero, a1, a2, shift, dd, one,start(maxdim1)
      parameter (zero = 0.0d0, one = 1.0d0, maxstep=1000)
      double precision eval(nev), fluct(maxvec),evec(2*nev*nbasis)
      external matva,matvb
      integer iwork(maxstep) 
        integer*8 nev0, imax
      double precision gap, range, ev(maxdim)
      common/cd/dd(maxdim)
      dimension id(maxdim)
      integer*8 maxwork
      parameter (maxwork = 9*maxstep + 10*maxdim)
      double precision rwork(maxwork), testeval
      external testeval
      double complex amatrix

        nev0=nev
      limit = .false.
      internal = .false.
      relative = .false.

 201  format(' Test of lanczab (on a diagonal  matrix A, with B=I!)')

      write(6,202) maxdim
 202  format(' give matrix dimension ( <',i5,')')

      if(nbasis.gt.maxdim) then
         write(6,101) nbasis, maxdim
 101     format(' nbasis = ',i10,' must not exceed maxdim=',i10)
         stop
      endif

      if(nbasis.le.0) then
         write(6,99) nbasis
 99      format(' invalid nbasis=',i5)
         stop
      endif

      write(6,203)
 203  format(' give CTYPE = 1 for real or 2 for complex:')


      if(ctype.ne.1.and.ctype.ne.2) then
         write(6,102) ctype
 102     format(' valid CTYPE are 1 (real) or 2(complex)')
         stop
      endif

      write(6,*)'give ABTYPE = 0  for standard or 1 generalized'

      if(abtype.ne.0.and.abtype.ne.1) then
         write(6,*) abtype
         stop
      endif


      write(6,*)'(INTERNAL = .false.) for lowesteig,1'

      internal = .false.
      if (i.eq.1) then
         internal = .true.
      endif
      
      limit = .false.
      
      if(.not.internal) then
         write(6,1040)
 1040    format(' enter 0 for no eigenvalue limits,'
     &        /' enter 1 for to set an absolute upper bound',
     &        /' enter 2 to set an upper bound relative to',
     &        ' the lowest eigenvalue:')
         if(ii.eq.0) then
            limit = .false.
            write(6,1060) 
 1060       format(' give the desired number of eigenvalues NEV')
         else 
            write(6,1050)
 1050       format(' give the desired number of eigenvalues NEV',
     &           /' and the upper bound ')
            limit = .true.
            if(ii.eq.2) then
               relative = .true.
            endif
         endif
      endif
      
      
      if(internal) then
         write(6,2)
 2       format(' give desired number of eigenvalues nev',
     &        ' and range [a1,a2]')
      endif


      write(6,4) 
 4    format(' to see lanczos progress every N steps, give N > 0:')
	laninfo=10
      do i1 = 1,nbasis
         dd(i1) =  testeval(i1,nbasis)
      enddo
c-------------------------------------

        a1=0.0d0
      if(.not.internal) then 
         call indexd(nbasis,dd,id)
         a2 = a1
         if(relative) a2 = a2 + dd(id(1))
         imax= nev
        
         if(nev0.gt.nbasis) imax= nbasis
         do i = 1,imax
            if(.not.limit.or.dd(id(i)).le.a2) then
               ev(i) = dd(id(i))
               nev0 = i
            endif
         enddo
         do i = 1,nev0
            id(i) = i
         enddo
      else
         do i = 1,nbasis
            ev(i) = (dd(i)-a1)*(dd(i)-a2)
         enddo
         call indexd(nbasis,ev,id)
         do i = 1,nev
            if(ev(id(i)).le.zero) nev0 = i
         enddo
         do i = 1,nev0
            ev(i) = dd(id(i))
         enddo
         call indexd(nev0,ev,id)
      endif
c-------------------------------------

      do i = 1,nev
         call vstart(nbasis,ctype,evec(1+(i-1)*ctype*nbasis),iseed)
      enddo


      
      eval(1) = a1
      eval(2) = a2
      if(limit) then
         eval(2) = a1
         if(relative)  eval(2) = eval(2) - one
      endif
      ierr = laninfo
      if(internal.or.limit) then
         nev = -nev
      endif
      if(abtype.eq.0) then
      call lanczos('H',ctype, matva,matva,nbasis,nev,eval,
     &        evec,fluct,rwork,maxwork,iwork,maxstep,ierr)
      else
         call lanczos('G',ctype, matva,matvb,nbasis,nev,eval,
     &        evec,fluct,rwork,maxwork,iwork,maxstep,ierr)
      endif

 
      if(internal) then
 113     format(' desired spectral range was ['
     &        ,f12.4,',',f12.4,']')
      endif

      if(limit) then
         if(relative) then
 131        format(' eigenvalues requested in range eval(i) <  eval(1)',
     &           ' +',f15.8)
         else
 132        format(' eigenvalues requested in range eval(i) < ',f15.8)
         endif
      endif

 112  format(' eigenvalues found: nev = ',i5)
      write(6,120) (eval(i),i=1,nev)
 120  format(3f25.15)

      write(6,135) nev0
 135  format(' expected answer: nev = ',i5)

	return
      end

      
      

      subroutine MATVA(VIN,VOUT,NBASIS,CTYPE,ADD)
      integer*8 NBASIS,maxdim
        integer ctype
      logical ADD
      double precision VIN(NBASIS,CTYPE), VOUT(NBASIS,CTYPE)
      double complex amatrix
      external amatrix

      call MATVEC(VIN,VOUT,NBASIS,CTYPE,ADD,amatrix)

      return
      end


      double complex function amatrix(i,j)
      integer i,j
      double precision d
      double complex zero,one
      parameter (zero = (0.0d0,0.0d0), one = (1.0d0,0.0d0))
      parameter (maxdim=500)
      common/cd/d(maxdim)

      amatrix = zero
      if(i.eq.j) then
         amatrix = one*d(i)
      endif
      return
      end


      subroutine MATVB(VIN,VOUT,NBASIS,CTYPE,ADD)
      integer*8 NBASIS
        integer CTYPE
      logical ADD
      double precision VIN(NBASIS,CTYPE), VOUT(NBASIS,CTYPE)
      double complex bmatrix
      external bmatrix

      call MATVEC(VIN,VOUT,NBASIS,CTYPE,ADD,bmatrix)

      return
      end


      double complex function bmatrix(i,j)
      integer i,j
      double complex zero,one
      parameter (zero = (0.0d0,0.0d0), one = (1.0d0,0.0d0))

      bmatrix = zero
      if(i.eq.j) bmatrix = one

      return
      end


      subroutine lanczab(ctype,abtype,nbasis,matva,matvb,
     &     polish,start,shift,maxstep,internal,a1,a2,
     &     nevec,eval,evec,fluct,
     &     work,rwork,iwork, d,e,e2, lanvec,
     &     gap,range,nstep,ierr)
        use constant
      implicit none
      logical polish,internal
      external matva, matvb
      integer maxstep, nbasis, nevec, nstep, ierr,ctype,abtype
      double precision a1, a2, shift, range, gap
      double precision eval(nevec),fluct(nevec)
      double precision d(maxstep), e(maxstep), e2(maxstep),
     &   rwork(5*maxstep)
      double precision lanvec(maxstep)
      integer iwork(maxstep)
      double precision evec(nbasis,ctype*nevec), start(nbasis,ctype)
      double precision work(nbasis,ctype*(2+abtype))


      double precision tiny, small
      parameter (tiny = 1.0d-14, small = 1.0d-6)

      double precision zero, one
      parameter (zero = 0.0d0, one = 1.0d0)

      double precision a, ac(2), oldev0, anorm, ev0, evshift  
      double precision  proj0, eval0, ev0lanc,err
      integer ix, iy, iz, is, it, i, j, k, info, id, idi 
      logical done, test, add, general
      integer ninfo



      evshift = zero
      id = 1 + ctype*(nevec-1)
      info = 0
      if(ierr.gt.0) then
         info = 1
         ninfo = ierr
      end if
      ierr = 0

      if (ctype.eq.1) then
         ix = 1
         iz = 2
         if(abtype.eq.0) then 
            it = 3
         else
            iy = 3
            it = 4
            is = 5
         endif
      else if (ctype.eq.2) then
         ix = 1 
         iz = 3
         if(abtype.eq.0) then
            it = 5
         else
            iy = 5
            it = 7
            is = 9
         endif
      else
         ierr = -6
         return
      end if

      if (abtype.eq.0) then
         general = .false.
      else if (abtype.eq.1) then
         general = .true.
      else
         ierr = -6
         return
      end if


      if(polish) then
         call vcopy(nbasis,ctype,evec(1,id),work(1,ix))
      else
         call vcopy(nbasis,ctype,start(1,1),work(1,ix))
         if (general) then
            call matvb(work(1,ix),work(1,iy),nbasis,ctype,.false.)
         else
            iy = ix
         end if
         do i = 1,nevec-1
            idi = 1 + ctype*(i-1)
            call vdot(nbasis,ctype,evec(1,idi),work(1,iy),ac)
            do j = 1,ctype
               ac(j) = -ac(j)
            end do
            call vaxpy(nbasis,ctype,ac,evec(1,idi),work(1,iy))
         end do
         
         if (general) then
            call matvb(work(1,ix),work(1,iy),nbasis,ctype,.false.)
            call vdot(nbasis,ctype,work(1,ix),work(1,iy),ac)
            a = ac(1)
            if(a.lt.-tiny) then
               write(6,10) 
 10            format('invalid MATVB: not positive definite')   
               ierr = -3
               return
            end if
            a = dsqrt(a)
         else
            call vdnrm2(nbasis,ctype,work(1,ix),a)
         end if

         if(a.lt.tiny) then
            write(6,20)
 20         format('nul entry into lanczab ')
            ierr = -5
            return
         end if

         if (general.and.ctype.eq.2) then
            if(dabs(ac(2)).gt.10*tiny) then
               write(6,11) ac(2),tiny
 11            format('invalid MATVB: not Hermitian')
               ierr = -3
               return
            end if
         end if

         a = 1.0d0/a
         call vdscal(nbasis,ctype,a,work(1,ix))
         call vcopy(nbasis,ctype,work(1,ix),evec(1,id))
      end if


      nstep = 0
      if (general) then
         call matvb(work(1,ix),work(1,iy),nbasis,ctype,.false.)
      else
         iy = ix
      end if
      e(1) = zero
      e2(1) = zero

      done = .false.
      if(ctype.eq.2) then
         test = .true.
      else
         test = .false.
      end if
      add = .false.

      do while (.not.done)
         nstep = nstep + 1

         if(.not.internal) then
            call matva(work(1,iy),work(1,iz),nbasis,ctype,add)
         else
            if(abtype.eq.0) then
               call matva(work(1,iy),work(1,it),nbasis,ctype,.false.)
            else
               call matva(work(1,iy),work(1,is),nbasis,ctype,.false.)
               call matvb(work(1,is),work(1,it),nbasis,ctype,.false.)
            endif
            ac(1) = -(a1+a2)
            ac(2) = zero
            call vaxpy(nbasis,ctype,ac,work(1,iy),work(1,it))
            call matva(work(1,it),work(1,iz),nbasis,ctype,add)
         endif
         
         add = .true.
         if(test) then
            if(.not.internal) then
               call vdot(nbasis,ctype,work(1,iz),work(1,iy),ac)
            else if(abtype.eq.0) then
               call vdot(nbasis,ctype,work(1,it),work(1,iy),ac)
            else
               call vdot(nbasis,ctype,work(1,is),work(1,iy),ac)
            endif
            if(dabs(ac(2)).gt.100*tiny) then
               ierr = -3
               write(6,12) dabs(ac(2))
 12            format('invalid MATVA: not Hermitian', f25.15)
               return
            endif
            test = .false.
         end if
      
         do i = 1,nevec-1
            if(.not.internal) then
               a = eval(nevec-1) - eval(i)
            else
               a = -(eval(i)-a1)*(eval(i) -a2)
            endif

            if (.not.internal) then
               a = a+shift
            else
               if(shift.gt.((a1-a2)/2)**2) then
                  a = a + shift
               else
                  a = a + ((a1-a2)/2)**2
               endif
            endif

            if(a.gt.zero) then
               idi = 1 + ctype*(i-1)
               call vdot(nbasis,ctype,evec(1,idi),work(1,iy),ac)
               do j = 1,ctype
                  ac(j) =  a*ac(j)
               end do
               call vaxpy(nbasis,ctype,ac,evec(1,idi),work(1,iz))
            endif
         end do
         call vdot(nbasis,ctype,work(1,iy),work(1,iz),ac)
         d(nstep) = ac(1)


         if(polish) then
            if(nstep.eq.1) then
               evshift = one - d(1)
            end if
            d(nstep) = d(nstep) + evshift
         end if

         if(nstep.gt.1) then
            call teval(nstep,d,e,e2,ev0,gap,range,rwork)
          d11(1:6)=rwork(1:6)
         else
            ev0 = d(1)
          d11(1:6)=d(1:6)
            gap = zero
            range = zero
         end if

         ev0lanc = ev0-evshift
          d11(1:6)=d11(1:6)-evshift
         if(internal) ev0lanc =  ev0lanc + a1*a2
         if(internal) d11(1:6) =d11(1:6)+ a1*a2
         

         if(info.eq.1.and.(mod(nstep-1,ninfo).eq.0))then
             write(6,30) nstep,ev0lanc,d(nstep),e(nstep)
             write(6,*) 'all eiv', d11(1:6)
                endif
 30      format(i7,f30.15,2(1pd12.3))
         
         if(nstep.ge.maxstep) then
            write(6,35) maxstep,eval0
 35         format(' lanczos1: ',/,' limit MAXSTEP =',i6,
     &           '  on number of lanczos steps is too small',
     &           ' for convergence',/' upper bound on energy :',f20.15)
            ierr = -1
            return
         end if
         
         
         
         if(nstep.gt.1.and.ev0.ge.oldev0) then
            if(.not.polish.or.info.eq.1)
     &           write(6,40) nevec,nstep
 40         format(' convergence : nevec =',i5,' t-matrix length:',i10)
            nstep = nstep - 1
            done = .true.
         end if


         if (.not.done) then
               ac(1) = evshift - d(nstep)
               ac(2) = zero
            call vaxpy(nbasis,ctype,ac,work(1,ix),work(1,iz))
            if (general) then
               call matvb(work(1,iz),work(1,iy),nbasis,ctype,.false.)
               call vdot(nbasis,ctype,work(1,iz),work(1,iy),ac)
               a = ac(1)
               if(a.lt.-tiny) then
                  write(6,10)
                  ierr = -3
                  return
               end if
               if (a.lt.zero) then
                  write(6,10)
                  ierr = -4
                  return
               end if
               a = dsqrt(a)
            else
               call vdnrm2(nbasis,ctype,work(1,iz),a)
            end if
         
         
            e2(nstep+1) = a**2
            e(nstep+1) = a 
            
            i = ix
            ix = iz
            iz = i
            oldev0 = ev0
            
       
            a = dsqrt(d(nstep)**2 + e2(nstep) + e2(nstep+1))
            if(e(nstep+1)/a.lt.tiny) then
               write(6,50)
 50            format(' T-matrix terminates ')
               done = .true.
            else
               a = 1/e(nstep+1)
               call vdscal(nbasis,ctype,a,work(1,ix))
               if (general) then
                  call vdscal(nbasis,ctype,a,work(1,iy))
               else
                  iy = ix
               end if
               a = -e(nstep+1)
               call vdscal(nbasis,ctype,a,work(1,iz))
            end if
         end if
      end do
      
      
      if (nstep.gt.1) then
         call tevec(nstep,d,e,e2,oldev0,lanvec,rwork,iwork)
      else
         lanvec(1) = one
      end if
      

      call vcopy(nbasis,ctype,evec(1,id),work(1,ix))
      if (general) then
         call matvb(work(1,ix),work(1,iy),nbasis,ctype,.false.)
      else
         iy = ix
      end if

      call vdscal(nbasis,ctype,zero,evec(1,id))
      add = .false.
      do k = 1,nstep
         if(info.eq.1.and.(mod(k-1,ninfo).eq.0))
     &     write(6,60) k,nstep
 60      format(' lanczos: reconstructing vector ',i6,' out of ',i6)
         ac(1) = lanvec(k)
         ac(2) = zero
         call vaxpy(nbasis,ctype,ac,work(1,ix),evec(1,id))

         if(k.lt.nstep) then
            if(.not.internal) then
               call matva(work(1,iy),work(1,iz),nbasis,ctype,add)
            else
               if(abtype.eq.0) then
                  call matva(work(1,iy),work(1,it),nbasis,ctype,.false.)
               else
                  call matva(work(1,iy),work(1,is),nbasis,ctype,.false.)
                  call matvb(work(1,is),work(1,it),nbasis,ctype,.false.)
               endif
               ac(1) = -(a1+a2)
               ac(2) = zero
               call vaxpy(nbasis,ctype,ac,work(1,iy),work(1,it))
               call matva(work(1,it),work(1,iz),nbasis,ctype,add)
            endif
            add = .true.
            
            do  i = 1,nevec-1
               if(.not.internal) then
                  a = eval(nevec-1) - eval(i)
               else
                  a =  -(eval(i)-a1)*(eval(i) -a2)
               endif

               if(.not.internal) then
                  a = a+shift
               else
                  if(shift.gt.((a1-a2)/2)**2) then
                     a = a + shift
                  else
                     a = a + ((a1-a2)/2)**2
                  endif
               endif

               if(a.gt.zero) then
                  idi = 1 + ctype*(i-1)
                  call vdot(nbasis,ctype,evec(1,idi),work(1,iy),ac)
                  do j = 1,ctype
                     ac(j) = a*ac(j)
                  end do
                  call vaxpy(nbasis,ctype,ac,evec(1,idi),work(1,iz))
               endif
            end do
            
            call vdot(nbasis,ctype,work(1,iy),work(1,iz),ac)
            a = ac(1)
            ac(1) = -a
            ac(2) = zero
            call vaxpy(nbasis,ctype,ac,work(1,ix),work(1,iz))
            
            if (general) then
               call matvb(work(1,iz),work(1,iy),nbasis,ctype,.false.)
               call vdot(nbasis,ctype,work(1,iz),work(1,iy),ac)
               a = ac(1)
               if(a.lt.-tiny) then
                  write(6,10)
                  ierr = -3
                  return
               end if
               if (a.le.zero) then
                  ierr = -4
                  return
               end if
               a = dsqrt(a)
            else
               call vdnrm2(nbasis,ctype,work(1,iz),a)
            end if

            a = -a
            i = ix
            ix = iz
            iz = i
            call vdscal(nbasis,ctype,a,work(1,iz))
            a = -1/a
            call vdscal(nbasis,ctype,a,work(1,ix))
            if (general) then
               call vdscal(nbasis,ctype,a,work(1,iy))
            else
               iy = ix
            end if
         end if
      end do
      


      if (general) then
         call matvb(evec(1,id),work(1,1),nbasis,ctype,.false.)
         call vdot(nbasis,ctype,evec(1,id),work(1,1),ac)
         a = ac(1)
         if (a.lt.-tiny) then
            write(6,10)
            ierr = -3
            return
         end if
         if(a.lt.zero) then
            ierr = -4
            return
         end if
         a = dsqrt(a)
      else
         call vdnrm2(nbasis,ctype,evec(1,id),a)
      end if

      anorm = 1/a
      call vdscal(nbasis,ctype,anorm,evec(1,id))

      if(nevec.gt.1) then
         i = 1
         idi = 1 + ctype*(i-1) 
         if(general) then
            call vdot(nbasis,ctype,evec(1,idi),work(1,1),ac)
            do j = 1,ctype
               ac(j) = ac(j)*anorm
            end do
         else
            call vdot(nbasis,ctype,evec(1,idi),evec(1,id),ac)
         end if
         
         proj0 = zero
         do j = 1,ctype
            proj0 = proj0 + ac(j)**2
         end do
         
         if(proj0.ge.small) then
            call vdscal(nbasis,ctype,zero,evec(1,id))
            write(6,70) proj0,i, shift
 70         format(' next eigenvector has finite projection ',f10.6,
     *           ' on eigenvector',i5,
     *           /' LANCZAB parameter SHIFT = ',f20.10,' is too small:'
     *           /' to find more eigenvectors,',  
     *           /' call LANCZAB with larger value of  SHIFT') 
            ierr = -2
            return
         end if
      end if

      if (general) then
         call matvb(evec(1,id),work(1,1),nbasis,ctype,.false.)
      else
         call vcopy(nbasis,ctype,evec(1,id),work(1,1))
      end if

      call matva(work(1,1),work(1,1+ctype),nbasis,ctype,.false.)
      call vdot(nbasis,ctype,work(1,1),work(1,1+ctype),ac)
      eval0  = ac(1)

      if (general) then
         call matvb(work(1,1+ctype),work(1,1+2*ctype),nbasis,ctype,
     &        .false.)
         call vdot(nbasis,ctype,work(1,1+ctype),work(1,1+2*ctype),ac)
      else
         call vdot(nbasis,ctype,work(1,1+ctype),work(1,1+ctype),ac)
      end if

      err = ac(1)- eval0**2
      fluct(nevec) = err
      eval(nevec) = eval0

      if(internal) then
         if(polish) ev0 = ev0 -evshift
         if(dabs(ev0-eval0*(eval0-a1-a2)).gt.small) then
            eval(nevec) = (a1+a2)/2 - dsqrt(((a1+a2)/2)**2 + ev0)
            ierr = 3
            if((eval(nevec)-a1)*(eval(nevec)-a2).gt.zero) 
     &           ierr = 1
            return
         endif
      endif


      if(internal.and.(eval0-a1)*(eval0-a2).gt.zero) ierr = 1

      if(internal.and.nevec.gt.1.and.ierr.ne.1.and.
     &     eval0.lt.(eval(1)+eval(nevec-1))/2) then
         call vcopy(nbasis,ctype,evec(1,id),work(1,1))
         do i = 1,nevec-1
            idi = 1 +  ctype*(nevec - i -1)
            call vcopy(nbasis,ctype,evec(1,idi),evec(1,idi+ctype))
            eval(nevec-i+1) = eval(nevec-i)
            fluct(nevec-i+1) = fluct(nevec-i)
         enddo
         call vcopy(nbasis,ctype,work(1,1),evec(1,1))
         eval(1) = eval0
         fluct(1) = err
      endif

      return
      end

      subroutine lanczos(eigtype,ctype,matveca,matvecb,n,nev,eval,
     &     evec,fluct,rwork,maxwork,iwork,maxstep,ierr)
      implicit none
      integer*8 n,  maxwork 
      integer  nev, ierr,  maxstep, ctype
      double precision eval(nev), evec(n,ctype,nev), fluct(nev)
      double precision rwork(maxwork)
      integer iwork(maxstep)
      external matveca,matvecb
      character eigtype

        integer*8 imem
      integer  npolish, laninf,nshift, nevec,i
      parameter ( nshift = 10, npolish = 5)
      double precision a,b,scale, small,shift, zero, one, range, gap
      integer*8 d1, e1, e21, lanvec1,rwork1,start1,work1,imax,icount 
      integer      nstep, abtype
      logical done, limit, internal, relative
      parameter (zero=0.0d0, one = 1.0d0, small = 1.0d-4)
      double precision evalmin,evalmax,deltamax,a1,a2,c(2)
     

      if(eigtype == 'H') then
         abtype = 0
      else if (eigtype == 'G') then
         abtype = 1
      else
         write(6,1) eigtype
 1       format('LANCZOS: allowed values of EIGTYPE are H or G',
     &        /' input value was: ',a1)
         stop
      endif


      if(n.lt.0.or.maxstep.le.0.or.maxwork.le.0.or.
     &     (ctype.ne.1.and.ctype.ne.2)) then
         write(6,2) n,ctype,maxstep,maxwork
 2       format('LANCZOS: called with invalid arguments:',
     &        'N=',i10,' CTYPE=',i3,' MAXSTEP=',i10,' MAXWORK=',i10) 
         stop
      endif
    



      if(n.eq.0)  nev = 0
      
      if(nev.eq.0) then
         ierr = 0
         return
      endif


      internal = .false.
      limit = .false.
      relative = .false.

      if(nev.lt.0) then
         nev = -nev
         evalmax = eval(1)
         if(eval(1).eq.eval(2)) then
            limit = .true.
         else if(eval(1).gt.eval(2)) then
            limit = .true.
            relative = .true.
         else
            internal = .true.
            evalmin = eval(1)
            evalmax = eval(2)
         endif
      endif


      if(nev.gt.n) nev = n
      if (n.eq.1) then
         evec(1,1,1) = one
         if(ctype.eq.2) evec(1,2,1) = zero
         call matveca(evec,rwork(start1),1,ctype,.false.)
         a = rwork(start1)
         if(abtype.eq.0) then
            eval(1) = a
         else 
            call matvecb(evec,rwork(start1),1,ctype,.false.)
            b = rwork(start1)
            eval(1) = a*b
            if(eigtype == 'A') then
               evec(1,1,1) = one/dsqrt(b)
            else
               evec(1,1,1) = dsqrt(b)
            endif
         endif
         ierr = 0
         if(.not.internal) then
            if(limit.and..not.relative.and.eval(1).ge.evalmax) then
               nev = 0
               evec(1,1,1) = zero
            endif
            return
         else 
            if(eval(1).ge.evalmin.and.eval(1).le.evalmax) return
            nev = 0
            evec(1,1,1) = zero
            return
         endif
      endif

      imem = 1
      d1 = imem
      imem = imem + maxstep
      lanvec1 = imem
      imem = imem +  maxstep
      e1 = imem
      imem = imem + maxstep
      e21 = imem
      imem = imem + maxstep
      rwork1 = imem
      imem = imem + 5*maxstep
      start1 = imem
      imem = imem + n*ctype
      work1 = imem
      imem = imem + (2+abtype)*ctype*n
      if(internal) then
         imem = imem + (1+abtype)*ctype*n
      endif
         
      if(imem.gt.maxwork) then
         write(6,10) imem, maxwork
 10      format('LANCZOS: workspace array RWORK must have size ',i12,
     &        ' but reported size (MAXWORK) =',i12)
         ierr = -100
         return
      endif

      

      if(abtype.eq.0) then
         call matveca(evec,rwork(start1),n,ctype,.false.)
         call vdot(n,ctype,evec,evec,c)
         b = c(1)
         call vdot(n,ctype,evec,rwork(start1),c)
         a = -c(1)/b
         call vaxpy(n,ctype,a,evec,rwork(start1))
         call vdot(n,ctype,rwork(start1),rwork(start1),c)
         a = c(1)
      else
         call matvecb(evec,rwork(start1),n,ctype,.false.)
         call matveca(rwork(start1),rwork(work1),n,ctype,.false.)
         call vdot(n,ctype,evec,rwork(start1),c)
         b = c(1)
         call vdot(n,ctype,rwork(start1),rwork(work1),c)
         a = -c(1)/b
         call vaxpy(n,ctype,a,rwork(start1),rwork(work1))
         call matvecb(rwork(work1),rwork(start1),n,ctype,.false.)
         call vdot(n,ctype,rwork(work1),rwork(start1),c)
         a = c(1)
      endif
      
      if(a.ge.zero.and.b.gt.zero) then
         scale = dsqrt(a/b)
      else
         write(6,6) a, b
 6       format('LANCZOS: incorrect MATVECB, not positive ',2d12.3)
      endif
      

c--------------------------------------------
      laninf = ierr
      imax = nev

      nevec = 1
        evalmin=0.0
        evalmax=0.0
      do while (nevec.le.imax)
         call vcopy(n,ctype,evec(1,1,nevec),rwork(start1))
         a1 = evalmin 
         a2 = evalmax
         shift = scale
         
         done = .false.
         icount = 0
         do while (icount.lt.nshift.and.(.not.done))            
            icount = icount +1
            if(icount.ne.1) then
               write(6,28) icount, shift
 28            format(' repeat lanczab:',i5,' shift=',d12.3)
            endif
            ierr = laninf 
            call lanczab(ctype,abtype,n,matveca,matvecb,.false.,
     *           rwork(start1),
     *           shift,maxstep,internal,a1,a2,
     *           nevec,eval,evec,fluct,rwork(work1),
     *           rwork(rwork1),iwork,
     *           rwork(d1),rwork(e1),rwork(e21),rwork(lanvec1),
     *           gap,range,nstep,ierr)
  

            if(ierr.eq.1) then
               nev = nevec-1
               write(6,555) nev
 555           format(' ierr=1 return, nev = ',i6)
               return
            endif

            if(ierr.eq.-2) then
               shift = 10*shift
            else
               done = .true.
            endif
         enddo
         
         if (ierr.lt.0) then
            write(6,30) ierr
 30         format(' lanczab failed, ierr=',i4)
            return
         else  if(laninf.gt.0) then
            if(ierr.eq.0) then
               write(6,20) nstep, eval(nevec),fluct(nevec)
 20            format('lanczab: convergence required',i10,
     &              ' lanczos steps:',/' eigenvalue:',f25.15,
     &              ' fluctuation ',d12.2)
            else
               write(6,21) nstep,ierr
 21            format(' lanczab: nstep=',i8,' ierr=',i8)
            endif
               
         endif
         

         if(internal) then
            if(eval(nevec).lt.(a1+a2)/2) then
               a2 = 2*eval(nevec) -a1
            else
               a1 = 2*eval(nevec) - a2
            endif
         endif

         done = .false.
         icount = 0
         do while (icount.lt.npolish.and.(.not.done))
            ierr=laninf
            icount = icount + 1
            call lanczab(ctype,abtype,n,matveca,matvecb,.true.,
     *           rwork(start1),
     *           shift,maxstep,internal,a1,a2,
     *           nevec,eval,evec,fluct,rwork(work1),
     *           rwork(rwork1),iwork,
     *           rwork(d1),rwork(e1),rwork(e21),rwork(lanvec1),
     *           gap,range,nstep,ierr)

            if(laninf.gt.0) then
               if(ierr.eq.0) then
                  write(6,36) nstep,eval(nevec)
 36               format(' polishing eigenvalue: nstep=',i8,f25.16)
               else
                  write(6,37) nstep,ierr
 37               format(' polishing eigenvalue: nstep=',i8,' ierr=',i8)
               endif
            endif

            if(nstep.lt.2) done = .true.
         enddo

         if(nevec.eq.1.and.ierr.eq.0.and.limit.and.relative)
     &        evalmax = evalmax + eval(1)

         if(laninf.gt.0) then
            write(6,40) eval(nevec),fluct(nevec)
 40         format('lanczab: polished eigenvalue:',f25.15,
     &           ' fluctuation ',d12.2)
         endif

         if(ierr.eq.1.and.
     &        (eval(nevec)-evalmin)*(eval(nevec)-evalmax).gt.zero)
     &        then
            nev = nevec-1
            return
         endif

         if(ierr.ge.0) then
            if(.not.limit.or.eval(nevec).lt.evalmax) then
               nev = nevec
            else
               return
            endif

         else if (ierr.lt.0) then
            write(6,50) ierr
 50         format('lanczab: polishing step failed, ierr=',i5)
            return
         endif
         nevec = nevec+1
      enddo

      return
      end

      DOUBLE PRECISION FUNCTION EPSLON (X) 
      DOUBLE PRECISION X 
      DOUBLE PRECISION A,B,C,EPS 
      A = 4.0D0/3.0D0 
   10 B = A - 1.0D0 
      C = B + B + B 
      EPS = DABS(C-1.0D0) 
      IF (EPS .EQ. 0.0D0) GO TO 10 
      EPSLON = EPS*DABS(X) 
      RETURN 
      END 
      DOUBLE PRECISION FUNCTION PYTHAG(A,B) 
      DOUBLE PRECISION A,B 
      DOUBLE PRECISION P,R,S,T,U 
      P = DMAX1(DABS(A),DABS(B)) 
      IF (P .EQ. 0.0D0) GO TO 20 
      R = (DMIN1(DABS(A),DABS(B))/P)**2 
   10 CONTINUE 
         T = 4.0D0 + R 
         IF (T .EQ. 4.0D0) GO TO 20 
         S = R/T 
         U = 1.0D0 + 2.0D0*S 
         P = U*P 
         R = (S/U)**2 * R 
      GO TO 10 
   20 PYTHAG = P 
      RETURN 
      END 

      subroutine teval(n,d,e,e2,ev0,gap,range,rwork)
      implicit none
      integer n
      double precision d(n),e(n),e2(n),rwork(n,2),ev0,gap,range
      integer i,ierr
      character *6 name
      if(n.le.1) then
         write(6,10) n
 10      format(" invalid call to teval with n = ",i5)
         stop
      end if
      

      do i = 1,n
         rwork(i,1) = d(i)
         rwork(i,2) = e2(i)
      end do
      write(name,16)
 16   format('TQLRAT')
      call tqlrat(n,rwork(1,1),rwork(1,2),ierr)
      if(ierr.ne.0) then
         write(6,20) name,n,ierr
 20     format(' failure of ',a6,' in LANCZAB; NSTEP=',i5,
     &        ' ierr = ',i5)
         stop
      end if

      ev0 = rwork(1,1)
      range = rwork(n,1) - rwork(1,1)
      gap   = rwork(2,1) - rwork(1,1)

      return
      end

      subroutine tevec(n,d,e,e2,eval,evec,rwork,iwork)
      implicit none
      integer n
      integer iwork(n)
      double precision d(n),e(n),e2(n),eval,evec(n),rwork(n,5)
      character*6 name
      double precision en(1)
      integer ind(1),isplit(1),ifail(1),ierr
      en(1) = eval
      ind(1) = 1
      isplit(1) = n

      if(n.le.1) then
         write(6,10) n
 10      format(" invalid call to tevec with n = ",i5)
         stop
      end if

      write(name,15)
 15   format('TINVIT')
      call tinvit(n,n,d,e,e2,1,en,ind,evec,ierr,
     &     rwork(1,1),rwork(1,2),rwork(1,3),
     &     rwork(1,4),rwork(1,5))

      if(ierr.ne.0) then
         write(6,20) name,n,ierr
 20     format(' failure of ',a6,' in TEVEC: N=',i5,
     &        ' ierr=',i5)
         stop
      end if 
      return
      end

      SUBROUTINE TINVIT(NM,N,D,E,E2,M,W,IND,Z, 
     X                  IERR,RV1,RV2,RV3,RV4,RV6) 
      INTEGER I,J,M,N,P,Q,R,S,II,IP,JJ,NM,ITS,TAG,IERR,GROUP 
      DOUBLE PRECISION D(N),E(N),E2(N),W(M),Z(NM,M), 
     X       RV1(N),RV2(N),RV3(N),RV4(N),RV6(N) 
      DOUBLE PRECISION U,V,UK,XU,X0,X1,EPS2,EPS3,EPS4,NORM,ORDER,EPSLON, 
     X       PYTHAG 
      INTEGER IND(M) 
      IERR = 0 
      IF (M .EQ. 0) GO TO 1001 
      TAG = 0 
      ORDER = 1.0D0 - E2(1) 
      Q = 0 
  100 P = Q + 1 
      DO 120 Q = P, N 
         IF (Q .EQ. N) GO TO 140 
         IF (E2(Q+1) .EQ. 0.0D0) GO TO 140 
  120 CONTINUE 
  140 TAG = TAG + 1 
      S = 0 
      DO 920 R = 1, M 
         IF (IND(R) .NE. TAG) GO TO 920 
         ITS = 1 
         X1 = W(R) 
         IF (S .NE. 0) GO TO 510 
         XU = 1.0D0 
         IF (P .NE. Q) GO TO 490 
         RV6(P) = 1.0D0 
         GO TO 870 
  490    NORM = DABS(D(P)) 
         IP = P + 1 
         DO 500 I = IP, Q 
  500    NORM = DMAX1(NORM, DABS(D(I))+DABS(E(I))) 
         EPS2 = 1.0D-3 * NORM 
         EPS3 = EPSLON(NORM) 
         UK = Q - P + 1 
         EPS4 = UK * EPS3 
         UK = EPS4 / DSQRT(UK) 
         S = P 
  505    GROUP = 0 
         GO TO 520 
  510    IF (DABS(X1-X0) .GE. EPS2) GO TO 505 
         GROUP = GROUP + 1 
         IF (ORDER * (X1 - X0) .LE. 0.0D0) X1 = X0 + ORDER * EPS3 
  520    V = 0.0D0 
         DO 580 I = P, Q 
            RV6(I) = UK 
            IF (I .EQ. P) GO TO 560 
            IF (DABS(E(I)) .LT. DABS(U)) GO TO 540 
            XU = U / E(I) 
            RV4(I) = XU 
            RV1(I-1) = E(I) 
            RV2(I-1) = D(I) - X1 
            RV3(I-1) = 0.0D0 
            IF (I .NE. Q) RV3(I-1) = E(I+1) 
            U = V - XU * RV2(I-1) 
            V = -XU * RV3(I-1) 
            GO TO 580 
  540       XU = E(I) / U 
            RV4(I) = XU 
            RV1(I-1) = U 
            RV2(I-1) = V 
            RV3(I-1) = 0.0D0 
  560       U = D(I) - X1 - XU * V 
            IF (I .NE. Q) V = E(I+1) 
  580    CONTINUE 
         IF (U .EQ. 0.0D0) U = EPS3 
         RV1(Q) = U 
         RV2(Q) = 0.0D0 
         RV3(Q) = 0.0D0 
  600    DO 620 II = P, Q 
            I = P + Q - II 
            RV6(I) = (RV6(I) - U * RV2(I) - V * RV3(I)) / RV1(I) 
            V = U 
            U = RV6(I) 
  620    CONTINUE 
         IF (GROUP .EQ. 0) GO TO 700 
         J = R 
         DO 680 JJ = 1, GROUP 
  630       J = J - 1 
            IF (IND(J) .NE. TAG) GO TO 630 
            XU = 0.0D0 
            DO 640 I = P, Q 
  640       XU = XU + RV6(I) * Z(I,J) 
            DO 660 I = P, Q 
  660       RV6(I) = RV6(I) - XU * Z(I,J) 
  680    CONTINUE 
  700    NORM = 0.0D0 
         DO 720 I = P, Q 
  720    NORM = NORM + DABS(RV6(I)) 
         IF (NORM .GE. 1.0D0) GO TO 840 
         IF (ITS .EQ. 5) GO TO 830 
         IF (NORM .NE. 0.0D0) GO TO 740 
         RV6(S) = EPS4 
         S = S + 1 
         IF (S .GT. Q) S = P 
         GO TO 780 
  740    XU = EPS4 / NORM 
         DO 760 I = P, Q 
  760    RV6(I) = RV6(I) * XU 
  780    DO 820 I = IP, Q 
            U = RV6(I) 
            IF (RV1(I-1) .NE. E(I)) GO TO 800 
            U = RV6(I-1) 
            RV6(I-1) = RV6(I) 
  800       RV6(I) = U - RV4(I) * RV6(I-1) 
  820    CONTINUE 
         ITS = ITS + 1 
         GO TO 600 
  830    IERR = -R 
         XU = 0.0D0 
         GO TO 870 
  840    U = 0.0D0 
         DO 860 I = P, Q 
  860    U = PYTHAG(U,RV6(I)) 
         XU = 1.0D0 / U 
  870    DO 880 I = 1, N 
  880    Z(I,R) = 0.0D0 
         DO 900 I = P, Q 
  900    Z(I,R) = RV6(I) * XU 
         X0 = X1 
  920 CONTINUE 
      IF (Q .LT. N) GO TO 100 
 1001 RETURN 
      END 
      SUBROUTINE TQLRAT(N,D,E2,IERR) 
      INTEGER I,J,L,M,N,II,L1,MML,IERR 
      DOUBLE PRECISION D(N),E2(N) 
      DOUBLE PRECISION B,C,F,G,H,P,R,S,T,EPSLON,PYTHAG 
      IERR = 0 
      IF (N .EQ. 1) GO TO 1001 
      DO 100 I = 2, N 
  100 E2(I-1) = E2(I) 
      F = 0.0D0 
      T = 0.0D0 
      E2(N) = 0.0D0 
      DO 290 L = 1, N 
         J = 0 
         H = DABS(D(L)) + DSQRT(E2(L)) 
         IF (T .GT. H) GO TO 105 
         T = H 
         B = EPSLON(T) 
         C = B * B 
  105    DO 110 M = L, N 
            IF (E2(M) .LE. C) GO TO 120 
  110    CONTINUE 
  120    IF (M .EQ. L) GO TO 210 
  130    IF (J .EQ. 30) GO TO 1000 
         J = J + 1 
         L1 = L + 1 
         S = DSQRT(E2(L)) 
         G = D(L) 
         P = (D(L1) - G) / (2.0D0 * S) 
         R = PYTHAG(P,1.0D0) 
         D(L) = S / (P + DSIGN(R,P)) 
         H = G - D(L) 
         DO 140 I = L1, N 
  140    D(I) = D(I) - H 
         F = F + H 
         G = D(M) 
         IF (G .EQ. 0.0D0) G = B 
         H = G 
         S = 0.0D0 
         MML = M - L 
         DO 200 II = 1, MML 
            I = M - II 
            P = G * H 
            R = P + E2(I) 
            E2(I+1) = S * R 
            S = E2(I) / R 
            D(I+1) = H + S * (H + D(I)) 
            G = D(I) - E2(I) / G 
            IF (G .EQ. 0.0D0) G = B 
            H = G * P / R 
  200    CONTINUE 
         E2(L) = S * G 
         D(L) = H 
         IF (H .EQ. 0.0D0) GO TO 210 
         IF (DABS(E2(L)) .LE. DABS(C/H)) GO TO 210 
         E2(L) = H * E2(L) 
         IF (E2(L) .NE. 0.0D0) GO TO 130 
  210    P = D(L) + F 
         IF (L .EQ. 1) GO TO 250 
         DO 230 II = 2, L 
            I = L + 2 - II 
            IF (P .GE. D(I-1)) GO TO 270 
            D(I) = D(I-1) 
  230    CONTINUE 
  250    I = 1 
  270    D(I) = P 
  290 CONTINUE 
      GO TO 1001 
 1000 IERR = L 
 1001 RETURN 
      END 
      subroutine  vaxpy(n,ctype,a,x,y)
      implicit none
      integer*8 n
        integer  ctype
      double precision a(ctype),x(n,ctype),y(n,ctype)
      integer*8 i
      
      if(n.le.0.or.(ctype.ne.1.and.ctype.ne.2)) then
         write(6,1) n,ctype
 1       format('VAXPY: invalid N, CTYPE = ',2i12)
         stop
      endif
      if(ctype.eq.1) then
         do i = 1,n
            y(i,1) = y(i,1) + a(1)*x(i,1)
         end do
      else
         do i = 1,n
            y(i,1) = y(i,1) +  a(1)*x(i,1) - a(2)*x(i,2)
            y(i,2) = y(i,2) +  a(2)*x(i,1) + a(1)*x(i,2)
         end do
      end if

      return
      end
      
      subroutine  vcopy(n,ctype,x,y)
      implicit none
      integer*8 n
        integer  ctype
      double precision x(n,ctype),y(n,ctype)
      integer*8 i
      
      if(n.le.0.or.(ctype.ne.1.and.ctype.ne.2)) then
         write(6,1) n,ctype
 1       format('VCOPY: invalid N, CTYPE = ',2i12)
         stop
      endif
      if(ctype.eq.1) then
         do i = 1,n
            y(i,1) = x(i,1)
         end do
      else
         do i = 1,n
            y(i,1) = x(i,1)
            y(i,2) = x(i,2)
         end do
      end if

      return
      end


      subroutine  vdnrm2(n,ctype,x,a)
      implicit none
      integer n, ctype
      double precision a,x(n,ctype)
      integer i
      double precision zero
      parameter (zero = 0.0d0)

      if(n.le.0.or.(ctype.ne.1.and.ctype.ne.2)) then
         write(6,1) n,ctype
 1       format('VDNRM2: invalid N, CTYPE = ',2i12)
         stop
      endif
      a = zero
      if(ctype.eq.1) then
         do i = 1,n
            a = a + x(i,1)**2
         enddo
      else
         do i = 1,n
            a = a + x(i,1)**2 + x(i,2)**2
         end do
      end if
      a = dsqrt(a)
      return
      end


      subroutine vdot(n,ctype,x,y,a)
      implicit none
      integer*8 n
        integer  ctype
      double precision x(n,ctype),y(n,ctype),a(ctype)
      
      integer*8 i
      double precision zero
      parameter (zero=0.0d0)
      
      if(n.le.0.or.(ctype.ne.1.and.ctype.ne.2)) then
         write(6,1) n,ctype
 1       format('VDOT: invalid n, ctype = ',2i12)
         stop
      endif
      if(ctype.eq.1) then
         a(1) = zero
         do i = 1,n
            a(1) = a(1) + x(i,1)*y(i,1)
         end do
      else
         a(1) = zero
         a(2) = zero
         do i = 1,n
            a(1) = a(1) + x(i,1)*y(i,1) + x(i,2)*y(i,2)
            a(2) = a(2) + x(i,1)*y(i,2) - x(i,2)*y(i,1)
         enddo
      endif

      return
      end
      subroutine vdscal(n,ctype,a,x)
      implicit none
      integer*8 n
        integer  ctype
      double precision a, x(n,ctype)
      integer*8 i

       if(n.le.0.or.(ctype.ne.1.and.ctype.ne.2)) then
         write(6,1) n,ctype
 1       format('VDSCAL: invalid N, CTYPE = ',2i12)
         stop
      endif
      if(ctype.eq.1) then
         do i = 1,n
            x(i,1) = a*x(i,1)
         end do
      else
         do i = 1,n
            x(i,1) = a*x(i,1)
            x(i,2) = a*x(i,2)
         end do
      end if

      return
      end

      SUBROUTINE INDEXD(N,ARRIN,INDX)
      implicit none
      integer*8 n,i,j,l,ir,indxt
      double precision arrin(n),q
      integer indx(n)



      if(n.le.0) return
      indx(1) = 1
      if(n.eq.1) return

      DO 11 J=1,N
        INDX(J)=J
11    CONTINUE
      L=N/2+1
      IR=N
10    CONTINUE
        IF(L.GT.1)THEN
          L=L-1
          INDXT=INDX(L)
          Q=ARRIN(INDXT)
        ELSE
          INDXT=INDX(IR)
          Q=ARRIN(INDXT)
          INDX(IR)=INDX(1)
          IR=IR-1
          IF(IR.EQ.1)THEN
            INDX(1)=INDXT
            RETURN
          ENDIF
        ENDIF
        I=L
        J=L+L
20      IF(J.LE.IR)THEN
          IF(J.LT.IR)THEN
            IF(ARRIN(INDX(J)).LT.ARRIN(INDX(J+1)))J=J+1
          ENDIF
          IF(Q.LT.ARRIN(INDX(J)))THEN
            INDX(I)=INDX(J)
            I=J
            J=J+J
          ELSE
            J=IR+1
          ENDIF
        GO TO 20
        ENDIF
        INDX(I)=INDXT
      GO TO 10
      END

      double precision function testeval(i,nbasis)
      integer*8 i, nbasis,maxbasis,minbasis
      double precision one
      parameter (one = 1.0d0, minbasis=0,maxbasis = 0)

      if(maxbasis.gt.0.and.nbasis.gt.manbasis) then
         write(6,1) nbasis, maxbasis
 1       format('TESTEVAL:  nbasis =',i10,
     &        ' exceeds maxbasis=',i10)
         stop
      else if (minbasis.gt.0.and.nbasis.lt.minbasis) then
         write(6,2) nbasis, maxbasis
 2       format('TESTEVAL:  nbasis =',i10,
     &        ' is less than minbasis=',i10)
         stop
      endif


      testeval =  2*one*((i+1)/2)

      return
      end
      subroutine vstart(n,ctype,x,iseed)
      integer*8 n, i
        integer  iseed, ctype
      double precision x(n,ctype)
      real ran2
      external ran2
      integer j
 
      if(n.le.0.or.(ctype.ne.1.and.ctype.ne.2)) then
         write(6,1) n, ctype
 1       format('VSTART: call with invalid n, ctype =',2i10)
         stop
      end if

      do j = 1,ctype
         do i = 1,n
            x(i,j) = dble(ran2(iseed))
         enddo
      end do
      return
      end



