           
           
           subroutine search(key,labelbasis,countbasis,flag4,location)
            use global
            implicit none

             integer :: countbasis,flag4,location
             integer(kind=8) :: labelbasis(1:MAXbasis)


             integer(kind=8) :: key     ! 所要寻找的值
             integer :: up      ! 记录每一个小组的类型起始位置
             integer :: down    ! 记录每一个小组的类型结束位置
             integer :: middle  ! 记录每一个小组的类型中间位置

                                ! 一开始的小组范围就是整个类型
             down=1
             up=countbasis
             middle=int((down+up)/2)

                                ! 如果key值超出范围, 铁定不存在类型中
             if((key>labelbasis(down)).or.(key<labelbasis(up))) then
               flag4=0
               return
             end if


             do while(down<=up)
              if(key>labelbasis(middle)) then
                                ! 如果key>中间值，那数据就落在上半部
                up=middle-1
                middle=int((down+up)/2)
              else if(key<labelbasis(middle)) then
                                ! 如果 key < 中间值，那数据就落在下半部
                down=middle+1
                middle=int((down+up)/2)
              else if(key==labelbasis(middle)) then
                location=middle
                flag4=1
                return
              end if
             end do

             flag4=0
            return
           end subroutine search











        subroutine make_band()
        use global
        implicit none


        !integer  kp1, kp2

        integer,parameter :: nsp=nq1*nq2*nl, nsk=nk1*nk2
        integer, parameter :: nnp=nsp, np1=nq1, np2=nq2


        integer, parameter :: neibt=12
        integer :: neib(1:nq1*nq2,12)

        double precision :: VD1
        double precision :: angle, x,y, g, prefac
        double complex :: tensor(0:30,2,2), ck, ck0, ckshift

        double precision phi,V1,V2,W1,W2

        real*8  valp(nnp)
        complex*16 hamp(nnp,nnp),vecp(nnp,nnp) !,vecf(nnp,nsk)
        !complex*16 vecf0(nnp,nsk),vecf3(nnp,nsk)

        integer  kspin, ks1

        double complex ci, k_plus,k_minus, kg1,kg2
        integer m1(2), m2(2)

        double precision :: mom_line(0:100,2)


        integer :: i,j,j1,j2,k1,k2,ki,ik

        integer :: ix,iy,ix0,iy0,ix1,iy1,ix2,iy2

        integer :: jj,j11,j12,j13,jy,jy1

        integer :: k11,k12,q1,q2,q11,q12,iq1,iq2,ia1,ib1

        character *2  chx,chy

        double complex :: curvature,tmp1,tmp2,tmp3,tmp4
        double precision :: C1,C2,C3
        double precision, external :: Converttheta

        double complex :: metric(1:nsk,1:nsk),Dxu(1:nnp,1:nsk),Dyu(1:nnp,1:nsk)

        character *3 cV


        !-----------------------------------------
        me= 9.1093837e-31
        m00 = 0.620d0*me
        hbar=1.054571817e-34
        rJ=6.582119569509066e-13/hbar

        prefactor = hbar**2/(2.0d0*m00)*rJ*1.0e20
        ci=cmplx(0.0d0, 1.0d0)
        VD1= 0.0d0

        a0=3.52d0
        factor=prefactor
        !theta=pi/3.0d0
        !antip=0
        tw_angle=1.95d0/180.0d0*pi !1.800d0/180.0d0*pi
        tw1=tw_angle*180.0d0/pi
        aM=a0/(2.0d0*sin(tw_angle/2.0d0))

        do j=1, 7
        angle=(j-1)*pi/3.0d0
        prefac=4*pi/(dsqrt(3.0d0)*aM)
        glist(j)=prefac*exp(ci*angle)
        enddo



        g1=glist(1)
        g2=glist(3)  !! g(2)-g(1)
        m1(1)=1
        m1(2)=0    !!! for g1
        m2(1)=-1
        m2(2)=1

        !mg0=0
        mq0=0
        mq0(1,1)=m1(1)   !1
        mq0(1,2)=m1(2)   !0
        mq0(2,1)=m2(1)   !-1
        mq0(2,2)=m2(2)   ! 1
        !mg0(1,1)=dfloat(m2(2))/(m2(2)*m1(1)-m2(1)*m1(2))
        !mg0(1,2)=dfloat(-m1(2))/(m2(2)*m1(1)-m2(1)*m1(2))
        !mg0(2,1)=dfloat(m2(1))/(m2(1)*m1(2)-m1(1)*m2(2))
        !mg0(2,2)=dfloat(-m1(1))/(m2(1)*m1(2)-m1(1)*m2(2))

        k_minus= (glist(1)+glist(6))/3.0d0
        k_plus=(glist(1)+glist(2))/3.0d0
        kg1=(glist(2)+glist(3))/3.0d0



        neib=0
        do ix=1, np1
        do jy=1, np2
        i= jy+(ix-1)*np2

        ix1=ix+1
        if(ix==np1)ix1=1
        jy1=jy+1
        if(jy==np2)jy1=1
        ix2=ix-1
        if(ix==1)ix2=np1

        if(ix1.ne.0)neib(i,1)=jy+(ix1-1)*np2
        if(jy1.ne.0)neib(i,2)=jy1+(ix-1)*np2
        if(jy1*ix2.ne.0)neib(i,3)=jy1+(ix2-1)*np2
        j11=jy+(ix1-1)*np2
        j12=jy1+(ix-1)*np2
        j13=jy1+(ix2-1)*np2
        neib(j11,neibt/2+1)=i
        neib(j12,neibt/2+2)=i
        neib(j13,neibt/2+3)=i
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

        !do i=1,np1*np2
        !do ik=1,6
        !if(abs(neib(i,ik))>1000.or.neib(i,ik)==0) print*,i,ik,neib(i,ik)
        !end do
        !end do
        !print*,neib(6,1)
        !stop

        V1=17.50d0
        W1=-6.50d0
        V2=-10.500d0
        W2=11.750d0
        phi= -0.97 !-1.02084 !-0.9859365d0  !-56.49d0/180.0d0*pi !-91.d0/180.0d0*pi !-0.9859365d0

        !V1=11.2d0
        !W1=-13.3d0
        !V2=0.0d0
        !W2=0.0d0
        !phi= -91.d0/180.0d0*pi 

        !V2=-10.0
        !do while(V2>=-24.0d0)


        tensor=0.0d0
        tensor(0,1,2)=W1 !+W2
        tensor(0,2,1)=W1 !+W2
        do j1=1,3
        tensor(j1,1,1)=tensor(j1,1,1)-V1*cdexp(ci*(-1)**(j1-1)*phi)
        tensor(j1,2,2)=tensor(j1,2,2)-V1*cdexp(ci*(-1)**(j1-1)*(-phi))
        enddo

        do j1=4,6
        tensor(j1,1,1)=tensor(j1,1,1)-V2 !*cdexp(ci*(-1)**(j1-1)*phi)
        tensor(j1,2,2)=tensor(j1,2,2)-V2  !*cdexp(ci*(-1)**(j1-1)*(-phi))
        enddo

        do j1=2,3
        tensor(j1,2,1)=tensor(j1,2,1)+W1
        enddo
        j1=1
        tensor(j1,2,1)=tensor(j1,2,1)+W2   !! *2
        tensor(j1,1,2)=tensor(j1,1,2)+W2   !! *2
        do j1=5, 5
        tensor(j1,2,1)=tensor(j1,2,1)+W2
        enddo


        do k1=1,40
          mom_line(k1,1)=0.0
          mom_line(k1,2)=k1*0.5/40
        end do
        do k1=41,65
          mom_line(k1,1)=(k1-40)*0.333333/25
          mom_line(k1,2)=0.5+(k1-40)*(0.666666-0.5)/25
        end do
        do k1=66,100
          mom_line(k1,1)=0.333333-(k1-65)*0.333333/35
          mom_line(k1,2)=0.666666-(k1-65)*0.666666/35
        end do



        kspin=1
        ks1=1
        !!if(kspin==2)ks1=-1
        !!if(kspin.eq.1)then
        k1=1
        k2=1
        !do k1=1,kx0
        !do k2=1,ky0 !!nk2
        !ki=k2+(k1-1)*ky0

        ! ck0=(k1-1)/dble(kx0)*g1+(k2-1)/dble(ky0)*g2
        ! !!ck0=(k1-1)/dble(kx0)*glist(1)+(k2-1)/dble(ky0)*glist(2)

        do ki=100,1,-1

         ck0=mom_line(ki,1)*g1 + mom_line(ki,2)*g2


         do i=1,nnp
         do j=1,nnp
           hamp(i,j)=cmplx(0.d0,0.d0)
         enddo
         enddo


        ix0=(np1-1)/2+1
        iy0=(np2-1)/2+1 

        ix=ix0
        iy=iy0
        do ix=1, np1
        do iy=1, np2

        i=iy+(ix-1)*np2          !!iy+(ix-1)*ny
        !ix0=(np1-1)/2+1
        !iy0=(np2-1)/2+1

        ck=ck0+glist(1)*(ix-ix0)+glist(2)*(iy-iy0)  !+ckshift

        hamp(2*i-1,2*i-1)=hamp(2*i-1,2*i-1)+factor*abs(ck-ks1*k_plus)**2-VD1/2.0d0
        hamp(2*i,2*i)    =hamp(2*i,2*i)    +factor*abs(ck-ks1*k_minus)**2+VD1/2.0d0
        hamp(2*i-1,2*i)  =hamp(2*i-1,2*i)  +W1  !tensor(0,1,2)
        hamp(2*i,2*i-1)  =hamp(2*i,2*i-1)  +W1  !tensor(0,2,1)

        go to 14
        do ik=1,6
        jj=neib(i,ik)
        !if(abs(jj)>1000) print*,i,ik,jj
        !!if(j.ne.0)then
        !print*, i,ik,jj
        !if(ik==5) print*, ix,iy, int(jj/np2)+1,mod(jj,np2)

        hamp(2*i-1,2*jj-1)=hamp(2*i-1,2*jj-1) +tensor(ik,1,1)
        hamp(2*i,2*jj)=hamp(2*i,2*jj) +tensor(ik,2,2)
        hamp(2*i-1,2*jj)=hamp(2*i-1,2*jj) +tensor(ik,1,2)
        hamp(2*i,2*jj-1)=hamp(2*i,2*jj-1) +tensor(ik,2,1)

        hamp(2*jj-1,2*i-1)=hamp(2*jj-1,2*i-1) +conjg(tensor(ik,1,1))
        hamp(2*jj,2*i)=hamp(2*jj,2*i) +conjg(tensor(ik,2,2))
        hamp(2*jj,2*i-1)=hamp(2*jj,2*i-1) +conjg(tensor(ik,1,2))
        hamp(2*jj-1,2*i)=hamp(2*jj-1,2*i) +conjg(tensor(ik,2,1))
        !!endif
        enddo!ik
        14 continue

        !go to 16
        !-------------------------------------------------------
        !neighbor 1
        ix1=ix+1
        iy1=iy
        if(ix1<=np1) then
        jj=iy1+(ix1-1)*np2

        !print*, ix1,iy1
        hamp(2*i-1,2*jj-1)=hamp(2*i-1,2*jj-1) -V1*cdexp(ii*phi) !
        hamp(2*i,2*jj)=hamp(2*i,2*jj) -V1*cdexp(-ii*phi)

        hamp(2*jj-1,2*i-1)=hamp(2*jj-1,2*i-1) -V1*cdexp(-ii*phi) 
        hamp(2*jj,2*i)=hamp(2*jj,2*i) -V1*cdexp(ii*phi)

        hamp(2*i,2*jj-1)=hamp(2*i,2*jj-1) + W2
        hamp(2*jj-1,2*i)=hamp(2*jj-1,2*i) + W2

        hamp(2*i-1,2*jj)=hamp(2*i-1,2*jj) + W2
        hamp(2*jj,2*i-1)=hamp(2*jj,2*i-1) + W2
        end if

        !neighbor 2
        ix1=ix
        iy1=iy+1
        if(iy1<=np2) then
        jj=iy1+(ix1-1)*np2

        !print*, ix1,iy1
        hamp(2*i-1,2*jj-1)=hamp(2*i-1,2*jj-1) -V1*cdexp(-ii*phi) !-V1*cdexp(ii*phi)
        hamp(2*i,2*jj)=hamp(2*i,2*jj) -V1*cdexp(ii*phi)

        hamp(2*jj-1,2*i-1)=hamp(2*jj-1,2*i-1) -V1*cdexp(ii*phi) !-V1*cdexp(ii*phi)
        hamp(2*jj,2*i)=hamp(2*jj,2*i) -V1*cdexp(-ii*phi)

        hamp(2*i,2*jj-1)=hamp(2*i,2*jj-1) + W1
        hamp(2*jj-1,2*i)=hamp(2*jj-1,2*i) + W1
        end if

        !neighbor 3
        ix1=ix-1
        iy1=iy+1
        if(ix1>=1.and.iy1<=np2) then
        jj=iy1+(ix1-1)*np2

        !print*, ix1,iy1
        hamp(2*i-1,2*jj-1)=hamp(2*i-1,2*jj-1) -V1*cdexp(ii*phi) !-V1*cdexp(-ii*phi)
        hamp(2*i,2*jj)=hamp(2*i,2*jj) -V1*cdexp(-ii*phi)

        hamp(2*jj-1,2*i-1)=hamp(2*jj-1,2*i-1) -V1*cdexp(-ii*phi) !-V1*cdexp(-ii*phi)
        hamp(2*jj,2*i)=hamp(2*jj,2*i) -V1*cdexp(ii*phi)

        hamp(2*i,2*jj-1)=hamp(2*i,2*jj-1) + W1
        hamp(2*jj-1,2*i)=hamp(2*jj-1,2*i) + W1
        end if


        !go to 16
        !next-neighbor 1
        ix1=ix+1
        iy1=iy+1
        if(ix1<=np1.and.iy1<=np2) then
        jj=iy1+(ix1-1)*np2

        !print*, ix1,iy1
        hamp(2*i-1,2*jj-1)=hamp(2*i-1,2*jj-1) -V2  !-V1*cdexp(ii*phi)
        hamp(2*i,2*jj)=hamp(2*i,2*jj) -V2 

        hamp(2*jj-1,2*i-1)=hamp(2*jj-1,2*i-1) -V2  !-V1*cdexp(ii*phi)
        hamp(2*jj,2*i)=hamp(2*jj,2*i) -V2 
        end if

        !next-neighbor 2
        ix1=ix-1
        iy1=iy+2
        if(ix1>=1.and.iy1<=np2) then
        jj=iy1+(ix1-1)*np2

        !print*, ix1,iy1
        hamp(2*i-1,2*jj-1)=hamp(2*i-1,2*jj-1) -V2  !-V1*cdexp(ii*phi)
        hamp(2*i,2*jj)=hamp(2*i,2*jj) -V2 

        hamp(2*jj-1,2*i-1)=hamp(2*jj-1,2*i-1) -V2  !-V1*cdexp(ii*phi)
        hamp(2*jj,2*i)=hamp(2*jj,2*i) -V2 

        hamp(2*i,2*jj-1)=hamp(2*i,2*jj-1) + W2
        hamp(2*jj-1,2*i)=hamp(2*jj-1,2*i) + W2
        end if

        !next-neighbor 3
        ix1=ix-2
        iy1=iy+1
        if(ix1>=1.and.iy1<=np2) then
        jj=iy1+(ix1-1)*np2

        hamp(2*i-1,2*jj-1)=hamp(2*i-1,2*jj-1) -V2  !-V1*cdexp(-ii*phi)
        hamp(2*i,2*jj)=hamp(2*i,2*jj) -V2 

        hamp(2*jj-1,2*i-1)=hamp(2*jj-1,2*i-1) -V2 !-V1*cdexp(-ii*phi)
        hamp(2*jj,2*i)=hamp(2*jj,2*i) -V2         
        end if
        16 continue

        end do!iy
        end do!ix
        !-------------------------------------
     


        do i=1, nnp-1
        do j=i+1, nnp
        if(hamp(i,j).ne.conjg(hamp(j,i)))then
        write(*,*)hamp(i,j),hamp(j,i), 'non-hermitian', i,j
        stop
        endif
        ! if(abs(hamp(i,j))>0.1d0) write(150,*) i,j,hamp(i,j)
        enddo
        enddo
        !stop
        !call eigen1(nnp, hamp, valp, vecp)
        call mkl_zheev(nnp, hamp, valp)


             ck=ck0
             open(11,file="val.dat", access='append')
             do i=1,10
               write(11,'(I4,I6,I6,4X,I4, F24.16)') ki, k1,k2,i,valp(i) !real(ck), aimag(ck), i,valp(i)
             end do
             write(11,*) ' '
             close(11)

         
         vecf(1:nnp, ki)= hamp(1:nnp,2)
         vecf0(1:nnp, ki)= hamp(1:nnp,1)
         vecf3(1:nnp, ki)= hamp(1:nnp,3)
         valf(ki,1:10) = valp(1:10)

         !write(800,'(I4,2F8.4,2X,10F12.4)') k1,mom_line(k1,1),mom_line(k1,2),valp(1:10)

        !enddo                 !!!for maxnum loop
        !enddo
        enddo !ki

        return

        !-------------------------------------------

        !----------------------------------------------
        !!form factor

        fq22=0.0d0
        fq11=0.0d0
        fq00=0.0d0
        fq01=0.0d0
        fq10=0.0d0
        do k11=1, ns0 !momentum in mBZ 
        do k12=1, ns0

        do q1=1, nq1  !larger momentum transfer V(q)
        do q2=1, nq2  !LARGER MOMENTUM TRANSFER V(a)
         do q11=1,nq1
         do q12=1,nq2

         iq1=q2+(q1-1)*nq2
         iq2=q12+(q11-1)*nq2

         ia1=q1-q11  !! momentum transfer-x
         ib1=q2-q12

         fq11(k11,k12,ia1,ib1)=fq11(k11,k12,ia1,ib1) &
              +conjg(vecf(2*iq1-1,k11))*vecf(2*iq2-1,k12)&
              +conjg(vecf(2*iq1,k11))*vecf(2*iq2,k12)

         fq00(k11,k12,ia1,ib1)=fq00(k11,k12,ia1,ib1) &
              +conjg(vecf0(2*iq1-1,k11))*vecf0(2*iq2-1,k12)&
              +conjg(vecf0(2*iq1,k11))*vecf0(2*iq2,k12)

         fq22(k11,k12,ia1,ib1)=fq22(k11,k12,ia1,ib1) &
              +conjg(vecf3(2*iq1-1,k11))*vecf3(2*iq2-1,k12)&
              +conjg(vecf3(2*iq1,k11))*vecf3(2*iq2,k12)

         fq10(k11,k12,ia1,ib1)=fq10(k11,k12,ia1,ib1) &
              +conjg(vecf(2*iq1-1,k11))*vecf0(2*iq2-1,k12) &
              +conjg(vecf(2*iq1,k11))*vecf0(2*iq2,k12)

         fq01(k11,k12,ia1,ib1)=fq01(k11,k12,ia1,ib1) &
              +conjg(vecf0(2*iq1-1,k11))*vecf(2*iq2-1,k12) &
              +conjg(vecf0(2*iq1,k11))*vecf(2*iq2,k12)      

         fq20(k11,k12,ia1,ib1)=fq20(k11,k12,ia1,ib1) &
              +conjg(vecf3(2*iq1-1,k11))*vecf0(2*iq2-1,k12) &
              +conjg(vecf3(2*iq1,k11))*vecf0(2*iq2,k12)

         fq02(k11,k12,ia1,ib1)=fq02(k11,k12,ia1,ib1) &
              +conjg(vecf0(2*iq1-1,k11))*vecf3(2*iq2-1,k12) &
              +conjg(vecf0(2*iq1,k11))*vecf3(2*iq2,k12) 

         fq12(k11,k12,ia1,ib1)=fq12(k11,k12,ia1,ib1) &
              +conjg(vecf(2*iq1-1,k11))*vecf3(2*iq2-1,k12) &
              +conjg(vecf(2*iq1,k11))*vecf3(2*iq2,k12)

         fq21(k11,k12,ia1,ib1)=fq21(k11,k12,ia1,ib1) &
              +conjg(vecf3(2*iq1-1,k11))*vecf(2*iq2-1,k12) &
              +conjg(vecf3(2*iq1,k11))*vecf(2*iq2,k12) 

         enddo
         enddo
        enddo
        enddo

        enddo
        enddo
        print*, 'form factor.'


        return
        end subroutine make_band




        subroutine make_harteefock()                !  ns ~ ab/2*pi*ell^2
        use global
        implicit none

        !double complex :: Vn(0:Nx*Ny-1,0:Nx*Ny-1,0:Nx*Ny-1,0:Nx*Ny-1)


        double precision :: vq0

        integer :: j1,j2,j3,j4,jj,jx1,jx4,jy1,jy4

        integer :: kx,ky,kx1,ky1

        integer :: ni,nj,n,i,j,j41,j12,iaa,ibb,ij,int_m,int_n,int_l

        double complex :: cqx ,cqx1

        integer :: q1,q2,q3,q4, qd1,qd2

        double complex :: Ma(Nb,Nb)
        double precision :: wa(Nb)

        double complex :: VH(1:Nb,1:Nb,1:ns0,1:ns0),VF(1:Nb,1:Nb,1:ns0,1:ns0)
        double complex :: err0,err1,tmp
        real :: rand

        integer :: ia1,ib1,iq1,iq2,k11,k12,q11,q12 

        double complex :: curvature,tmp1,tmp2,tmp3,tmp4
        double precision :: C1,C2,C3
        double precision, external :: Converttheta

        integer,parameter :: nsp=nq1*nq2*nl, nsk=nk1*nk2
        integer, parameter :: nnp=nsp, np1=nq1, np2=nq2

        integer :: ix,iy,ix0,iy0,ix1,iy1,ix2,iy2
        integer :: k1,k2,j11,j13

        !double precision :: val_hf(1:ns0)

        ee1=1.60217663*8.9875517923d0*1e+03
        epsilon=5.0d0


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


        print* ,'begin hartee-fock ... '

        go to 110
        !Hartee-Fock  
        val_hf=0.0d0
        do j1=1,ns0
        do j2=1,ns0

        do ij=1,2

        if(ij==1)then !Hartee
        j3=j2
        j4=j1
        else   !Fock
        j3=j1
        j4=j2
        endif

        vq0=0.0d0
        do q1=-nq1,nq1
        do q2=-nq2,nq2
        q3=-q1
        q4=-q2

        cqx= (mnt(j1,1)-mnt(j4,1))/dble(kx0)*g1+q1*glist(1)+&
             (mnt(j1,2)-mnt(j4,2))/dble(ky0)*g2+q2*glist(2)
        vq0=cdabs(cqx)


        if(vq0.ge.1e-15)then
        !if(dgate.gt.0.1)then
        ! vq0=tanh(dgate*vq0)/vq0
        !else
         vq0=1.0d0/vq0
        !endif
        if(ij==1)then
!!!!!        val_hf(j2)=val_hf(j2)+2.0d0*vq0*fq0(j1,j4,q1,q2)*fq(j2,j3,q3,q4)
        !val_hf(j2)=val_hf(j2)+2.0d0*vq0*fq00(j1,j4,q1,q2)*fq(j2,j3,q3,q4)
        else
!!!!!        val_hf(j2)=val_hf(j2)-2.0d0*vq0*abs(fq10(j1,j4,q1,q2)**2)
        endif
        endif
        enddo !q1
        enddo !q2
        enddo !ij 
!!!switching
        enddo !j2
        enddo !j1
        val_hf=val_hf*ee1/2.0d0

        !do j1=1,ns0
        !write(23,*)j1, val_hf(j1)
        !enddo
        110 continue




               !random initialization 
               rho=0.0d0
               do int_m=1,Nb
               do int_n=1,Nb
               do j1=1,Nx*Ny
                   if(int_m==1.and.int_n==1) then
                    rho(int_m,int_n,j1)=1.0d0
                   !else
                   ! !call random_number(rand)
                   ! rho(int_m,int_n,j1)=1.0d0 !rand !/(Nx*Ny)
                   end if
               end do
               end do
               end do




        !do int_l=1,10
        val_hf=0.0d0

        VH=0.0d0
        do j1=1,ns0
        do j2=1,ns0

        !Hartee term
        j3=j2
        j4=j1

        vq0=0.0d0
        do q1=-nq1,nq1 !difference between k1 and k4 
        do q2=-nq2,nq2

        q3=-q1 !difference between k2 and k3
        q4=-q2

        cqx= (mnt(j1,1)-mnt(j4,1))/dble(kx0)*g1+q1*glist(1)+&
             (mnt(j1,2)-mnt(j4,2))/dble(ky0)*g2+q2*glist(2)
        vq0=cdabs(cqx)


        if(vq0.ge.1e-15)then
         vq0=1.0d0/vq0

         VH(1,2,j1,j2)=VH(1,2,j1,j2)+1.0d0*vq0*fq00(j1,j4,q1,q2)*fq11(j2,j3,q3,q4)
         VH(1,3,j1,j2)=VH(1,3,j1,j2)+1.0d0*vq0*fq00(j1,j4,q1,q2)*fq22(j2,j3,q3,q4)
         VH(2,3,j1,j2)=VH(2,3,j1,j2)+1.0d0*vq0*fq11(j1,j4,q1,q2)*fq22(j2,j3,q3,q4)

         VH(2,1,j1,j2)=VH(2,1,j1,j2)+1.0d0*vq0*fq11(j1,j4,q1,q2)*fq00(j2,j3,q3,q4)
         VH(3,1,j1,j2)=VH(3,1,j1,j2)+1.0d0*vq0*fq22(j1,j4,q1,q2)*fq00(j2,j3,q3,q4)
         VH(3,2,j1,j2)=VH(3,2,j1,j2)+1.0d0*vq0*fq22(j1,j4,q1,q2)*fq11(j2,j3,q3,q4)
         !print*, j1,j2, VH(j1,j2)
     !    val_hf(j2)=val_hf(j2)+2.0d0*vq0*fq00(j1,j4,q1,q2)*fq11(j2,j3,q3,q4)
         endif

        enddo !q1
        enddo !q2

        enddo !j2
        enddo !j1
        VH=VH*ee1/2.0d0
        !print*, 'VH term'


        VF=0.0d0
        do j1=1,ns0
        do j2=1,ns0

        !Fock term
        j3=j1
        j4=j2

        vq0=0.0d0
        do q1=-nq1,nq1 !difference between k1 and k4 
        do q2=-nq2,nq2

        q3=-q1 !difference between k2 and k3
        q4=-q2

        cqx= (mnt(j1,1)-mnt(j4,1))/dble(kx0)*g1+q1*glist(1)+&
             (mnt(j1,2)-mnt(j4,2))/dble(ky0)*g2+q2*glist(2)
        vq0=cdabs(cqx)


        if(vq0.ge.1e-15)then
         vq0=1.0d0/vq0

         VF(1,2,j1,j2)=VF(1,2,j1,j2)+1.0d0*vq0*fq01(j1,j4,q1,q2)*fq10(j2,j3,q3,q4)
         !VF(3,2,j1,j2)=VF(3,2,j1,j2)+1.0d0*vq0*fq12(j1,j4,q1,q2)*fq21(j2,j3,q3,q4)
         !VF(3,1,j1,j2)=VF(3,1,j1,j2)+1.0d0*vq0*fq02(j1,j4,q1,q2)*fq20(j2,j3,q3,q4)

         VF(2,1,j1,j2)=VF(2,1,j1,j2)+1.0d0*vq0*fq10(j1,j4,q1,q2)*fq01(j2,j3,q3,q4)
         !VF(2,3,j1,j2)=VF(2,3,j1,j2)+1.0d0*vq0*fq21(j1,j4,q1,q2)*fq12(j2,j3,q3,q4)
         !VF(1,3,j1,j2)=VF(1,3,j1,j2)+1.0d0*vq0*fq20(j1,j4,q1,q2)*fq02(j2,j3,q3,q4)

         !VF(2,1,j1,j2)=VF(2,1,j1,j2)+ vq0*abs(fq10(j1,j4,q1,q2)**2)
         !VF(1,2,j1,j2)=VF(1,2,j1,j2)+ vq0*abs(fq10(j1,j4,q1,q2)**2)

         !val_hf(j1)=val_hf(j1)-2.0d0*vq0*abs(fq10(j1,j4,q1,q2)**2)
         val_hf(j1)=val_hf(j1)-2.0d0*vq0*fq10(j1,j4,q1,q2)*fq01(j2,j3,q3,q4)
         !print*, j1,j2, VH(j1,j2)
        endif

        enddo !q1
        enddo !q2

        enddo !j2
        enddo !j1
        VF=VF*ee1/2.0d0
        !print*, 'VF term'

        !VF=0.0d0
        val_hf=val_hf*ee1/2.0d0

        !vq0=0.0d0
        !do j1=1,ns0
        !  vq0=vq0+ val_hf(j1)/ns0
        !end do
        !val_hf=val_hf-vq0
        !vq0=0.0d0
        !do j1=1,ns0
        !  vq0=vq0+ val_hf_wrong(j1)/ns0
        !end do
        !val_hf_wrong=val_hf_wrong-vq0


        do kx=1,kx0
        do ky=1,ky0
        j2=ky+(kx-1)*ky0

        tmp=0.0d0
        do j1=1,ns0
           tmp=tmp-2*VF(2,1,j2,j1)
        end do

        write(23,'(I4,I4,I4,2X,F12.6,F12.6)')j2, kx,ky, val_hf(j2), dble(tmp)
        write(24,'(I4,I4,I4,2X,F12.6,F12.6)')j2, kx,ky, val_hf(j2)+valf(j2,2)
        enddo
        write(23,*) ' '
        write(24,*) ' '
        end do







             do j2=1,ns0

                  Ma(:,:)=0.0d0

                  !!Hartee
                  !tmp=0.0d0
                  !do j1=1,ns0 !Nx*Ny 
                  !  tmp = tmp + VH(j2,j1)*rho11(j1) !*rho00(j2)
                  !end do
                  !Ma(1,1)=tmp + valf(j2,1) !*rho00(j2)

                  !tmp=0.0d0
                  !do j1=1,ns0 !Nx*Ny
                  !  tmp = tmp + VH(j1,j2)*rho00(j1) !*rho00(j2)
                  !end do
                  !Ma(2,2)=tmp + valf(j2,2) !*rho11(j2)

                  do i=1,Nb
                   Ma(i,i)=valf(j2,i)
                  end do

                  do i=1,Nb

                    tmp=0.0d0
                    do j=1,Nb
                    if(j/=i) then
                     tmp=0.0d0
                     do j1=1,ns0 ! 
                      tmp = tmp + 2.0*VH(i,j,j1,j2)*rho(i,i,j1) 
                     end do

                     Ma(j,j)=Ma(j,j)+tmp
                    end if
                    end do

                  end do !i


                  !!Fock
                  !tmp=0.0d0
                  !do j1=1,ns0
                  !  tmp = tmp + VF(j1,j2)*rho01(j1)
                  !end do
                  !Ma(2,1)= tmp

                  !tmp=0.0d0
                  !do j1=1,ns0
                  !  tmp = tmp + VF(j2,j1)*rho10(j1)
                  !end do
                  !Ma(1,2)= tmp


                  i=1; j=2
                  do i=1,Nb

                    tmp=0.0d0
                    do j=1,Nb
                    if(j/=i) then
                     tmp=0.0d0
                     do j1=1,ns0 !
                      tmp = tmp - 2.0*VF(j,i,j2,j1) *rho(i,i,j1) 
                      !tmp=tmp - 2.0d0* VF(i,j,j2,j1)*rho(i,i,j1)
                     end do
                     !print*, j2, tmp
                     Ma(j,j)=Ma(j,j)+tmp
                    end if
                    end do

                  end do !i

                  !print*, Ma(2,2)

                  !call mkl_zheev(Nb,Ma,wa)

                  !bandk(1:Nb,j2)=wa(1:Nb)
                  !eigveck(1:Nb,1:Nb,j2)=Ma(1:Nb,1:Nb)

                  bandk(2,j2)=Ma(2,2)
                  !print*, bandk(2,j2)
             end do


             err0=0.0d0;err1=0.0d0
             do j1=1,Nx*Ny

                do int_m=1,Nb
                do int_n=1,Nb
                  tmp=0.0d0
                  tmp=conjg(eigveck(int_m,1,j1))*eigveck(int_n,1,j1)+&
                      conjg(eigveck(int_m,2,j1))*eigveck(int_n,2,j1)
                  err0=err0+abs(rho(int_m,int_n,j1)-tmp)
                  rho(int_m,int_n,j1)=tmp
                end do
                end do

             end do
             !write(*,'(I4,F12.6,F12.6)') int_l, abs(err0) !,abs(err1) 

           do j1=1,Nx*Ny
             write(490,'(I4,F8.4,F8.4,F8.4,F8.4,F8.4,F8.4,F8.4)') j1,abs(rho(1,1,j1)),abs(rho(2,2,j1)),abs(rho(3,3,j1)), &
                abs(rho(1,2,j1)),abs(rho(2,1,j1)),abs(rho(1,3,j1)),abs(rho(2,3,j1))
           end do

           !  if(abs(err0)<0.001.and.abs(err1)<0.001) go to 119
           !end do !int_l
           119 continue

           
           do k1=1,kx0
           do k2=1,ky0
            j1=k2+(k1-1)*ky0
             !write(490,'(I4,F8.4,F8.4,F8.4,F8.4)') j1,abs(rho(1,1,j1)),abs(rho(2,2,j1)),abs(rho(1,2,j1)),abs(rho(2,1,j1))
write(390,'(I4,I4,I4,2X,F12.6,F12.6,F12.6,4X,F12.6,F12.6,F12.6)') j1,k1,k2,bandk(1,j1), bandk(2,j1),bandk(3,j1), valf(j1,1), valf(j1,2),valf(j1,3)
           end do
           end do

           return

        do j1=1,ns0
        vecf(:,j1)=eigveck(1,2,j1)*vecf0(:,j1) +eigveck(2,2,j1)*vecf(:,j1)
        vecf0(:,j1)=eigveck(1,1,j1)*vecf0(:,j1) +eigveck(2,1,j1)*vecf(:,j1)
        end do


        !go to 40
        !-------------------------------------------
        !Berry phase
        C1=0.0d0; C2=0.0d0; C3=0.0d0
        do k1=1,kx0
        do k2=1,ky0

          !ki=k2+(k1-1)*ky0
          ix0=k1
          iy0=k2
          jj=iy0+(ix0-1)*ky0

          ix=ix0+1
          if(ix0==kx0) ix=1
          iy=iy0
          j11=iy+(ix-1)*ky0

          ix1=ix0+1
          if(ix0==kx0) ix1=1
          iy1=iy0+1
          if(iy0==ky0) iy1=1
          j12=iy1+(ix1-1)*ky0

          ix2=ix0
          iy2=iy0+1
          if(iy0==ky0) iy2=1
          j13=iy2+(ix2-1)*ky0

          !---------------------------------- 
          tmp1=0.0d0 !     
          do i=1,nnp
            tmp1=tmp1+ conjg(vecf0(i,jj))*vecf0(i,j11)
          end do
          tmp2=0.0d0 !
          do i=1,nnp
            tmp2=tmp2+ conjg(vecf0(i,j11))*vecf0(i,j12)
          end do
          tmp3=0.0d0 !
          do i=1,nnp
            tmp3=tmp3+ conjg(vecf0(i,j12))*vecf0(i,j13)
          end do
          tmp4=0.0d0 !
          do i=1,nnp
            tmp4=tmp4+ conjg(vecf0(i,j13))*vecf0(i,jj)
          end do

          curvature=tmp1*tmp2*tmp3*tmp4

          C1=C1+aimag(log(curvature)) !Converttheta(curvature)
          !---------------------------------- 
          tmp1=0.0d0 !     
          do i=1,nnp
            tmp1=tmp1+ conjg(vecf(i,jj))*vecf(i,j11)
          end do
          tmp2=0.0d0 !
          do i=1,nnp
            tmp2=tmp2+ conjg(vecf(i,j11))*vecf(i,j12)
          end do
          tmp3=0.0d0 !
          do i=1,nnp
            tmp3=tmp3+ conjg(vecf(i,j12))*vecf(i,j13)
          end do
          tmp4=0.0d0 !
          do i=1,nnp
            tmp4=tmp4+ conjg(vecf(i,j13))*vecf(i,jj)
          end do

          curvature=tmp1*tmp2*tmp3*tmp4

          C2=C2+aimag(log(curvature))  !Converttheta(curvature)

        end do
        end do

        print*, C1/2/pi, C2/2/pi !, C3/2/pi
        40 continue


        fq11=0.0d0
        !fq00=0.0d0
        !fq01=0.0d0
        !fq10=0.0d0
        do k11=1, ns0 !momentum in mBZ 
        do k12=1, ns0

        do q1=1, nq1  !larger momentum transfer V(q)
        do q2=1, nq2  !LARGER MOMENTUM TRANSFER V(a)
         do q11=1,nq1
         do q12=1,nq2

         iq1=q2+(q1-1)*nq2
         iq2=q12+(q11-1)*nq2

         ia1=q1-q11  !! momentum transfer-x
         ib1=q2-q12

         fq11(k11,k12,ia1,ib1)=fq11(k11,k12,ia1,ib1) &
              +conjg(vecf(2*iq1-1,k11))*vecf(2*iq2-1,k12)&
              +conjg(vecf(2*iq1,k11))*vecf(2*iq2,k12)

         !fq00(k11,k12,ia1,ib1)=fq00(k11,k12,ia1,ib1) &
         !     +conjg(vecf0(2*iq1-1,k11))*vecf0(2*iq2-1,k12)&
         !     +conjg(vecf0(2*iq1,k11))*vecf0(2*iq2,k12)

         !fq10(k11,k12,ia1,ib1)=fq10(k11,k12,ia1,ib1) &
         !     +conjg(vecf(2*iq1-1,k11))*vecf0(2*iq2-1,k12) &
         !     +conjg(vecf(2*iq1,k11))*vecf0(2*iq2,k12)

         !fq01(k11,k12,ia1,ib1)=fq01(k11,k12,ia1,ib1) &
         !     +conjg(vecf0(2*iq1-1,k11))*vecf(2*iq2-1,k12) &
         !     +conjg(vecf0(2*iq1,k11))*vecf(2*iq2,k12)
         enddo
         enddo
        enddo
        enddo

        enddo
        enddo


        return
        end subroutine make_harteefock








