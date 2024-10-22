


         subroutine constructbasis(labelbasis,countbasis,two)
          use global
          implicit none
          
           integer(kind=8) :: two(0:Nx*Ny-1)

           integer(kind=8) :: labelbasis(1:MAXbasis)


           integer :: int_i,int_l
           integer :: i1,i2,i3,i4,i5,i6,i7,i8,i9,i10,i11,i12

           integer :: countbasis

           integer :: orb(0:Nx*Ny-1,2),kx,ky

           character *2 ci


        do kx=0,Nx-1
        do ky=0,Ny-1

         orb(ky+kx*Ny ,1)=kx
         orb(ky+kx*Ny ,2)=ky

        enddo
        enddo


           labelbasis(:)=0

           two(Nx*Ny-1)=1
           do int_i=Nx*Ny-2,0,-1
            two(int_i)=two(int_i+1)*2
           end do

           countbasis=0
           do i1=0,Nx*Ny-1
            do i2=i1+1,Nx*Ny-1
             do i3=i2+1,Nx*Ny-1
              do i4=i3+1,Nx*Ny-1
               do i5=i4+1,Nx*Ny-1
                do i6=i5+1,Nx*Ny-1
                 do i7=i6+1,Nx*Ny-1
                  do i8=i7+1,Nx*Ny-1
              !     do i9=i8+1,Nx*Ny-1
              !      do i10=i9+1,Nx*Ny-1
              !       do i11=i10+1,Nx*Ny-1
              !        do i12=i11+1,Nx*Ny-1


                   kx= orb(i1,1)+orb(i2,1)+orb(i3,1)+orb(i4,1)+&
                       orb(i5,1)+orb(i6,1)+orb(i7,1)+orb(i8,1) !+&
                       !orb(i9,1)+orb(i10,1)+orb(i11,1)+orb(i12,1)     
                   ky= orb(i1,2)+orb(i2,2)+orb(i3,2)+orb(i4,2)+&
                       orb(i5,2)+orb(i6,2)+orb(i7,2)+orb(i8,2) !+&
                       !orb(i9,2)+orb(i10,2)+orb(i11,2)+orb(i12,2) 

                   if(mod(kx,Nx)==Total_Kx .and. mod(ky,Ny)==Total_Ky) then
                     countbasis=countbasis+1
                if(countbasis>MAXbasis) print *, 'countbasis overflow.'
                     !basis(countbasis,0)=int_i
                     !basis(countbasis,1)=int_j
                     !basis(countbasis,2)=int_k
                     !basis(countbasis,3)=int_p
                     !basis(countbasis,4)=int_q
                     !basis(countbasis,5)=int_r
                     !basis(countbasis,6)=int_s
                     !basis(countbasis,7)=int_t
                     !basis(countbasis,8)=int_o
                     !basis(countbasis,9)=int_u

                 labelbasis(countbasis)=two(i1)+two(i2)+two(i3)+two(i4)+&
                                        two(i5)+two(i6)+two(i7)+two(i8) !+&
                                        !two(i9)+two(i10)+two(i11)+two(i12) 
                   end if

               !       end do
               !      end do
               !     end do
               !    end do
                  end do
                 end do
                end do
               end do
              end do
             end do
            end do
           end do

           Write(*,'(A18,A2,I4,A2,I4,A2,I8)') 'Total basis for','(',Total_Kx,',',Total_Ky,')',countbasis


           !open(10,file='basis.dat')
           !do int_l=0,countbasis
           ! write(10,'(I8,A2,I4,I4,I4,I4,I4,I4,A2)') int_l,&
           ! '(',basis(int_l,0),basis(int_l,1),basis(int_l,2),basis(int_l,3),&
           !     basis(int_l,4),basis(int_l,5),&
           !     !basis(int_l,6),basis(int_l,7),&
           !     !basis(int_l,8),basis(int_l,9),&
           !  ')'
           !end do
           !close(10)  !write basis on disk

           !**********************************************************************
            !if(dimen/=countbasis) stop
            dimen=countbasis

            !two(Nx*Ny-1)=1
            !do int_i=Nx*Ny-2,0,-1
            ! two(int_i)=two(int_i+1)*2
            !end do

            !labelbasis(:)=0
            !do int_i=0,countbasis
            ! do int_j=0,Ne-1
            !  labelbasis(int_i)=labelbasis(int_i)+two(basis(int_i,int_j))
            ! end do
            ! if(int_i>0.and.labelbasis(int_i)>=labelbasis(int_i-1)) print *, int_i,labelbasis(int_i)
            !end do


           !**********************************************************************

            !write(ci,'(I2.2)') Total_Kx*Ny+Total_Ky
            !open(11,file='labelbasis'//ci//'.dat')
            ! do int_l=1,countbasis
            !  write(11,'(I8,I12)') int_l,labelbasis(int_l)
            ! end do
            !close(11)  !write basis on disk

         return
         end subroutine constructbasis


