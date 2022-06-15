      program generate
      implicit none

      integer i,j,k,l,indice
      double precision a1(2),a2(2),a1r(2),a2r(2),r(2)
      double precision tau1(2),tau2(2),tau1r(2),tau2r(2),theta
c rectangular unit cell:
      double precision ra1(2),ra2(2),rt1(2),rt2(2),rt3(2),rt4(2)
      double precision ra1r(2),ra2r(2),
     .                 rt1r(2),rt2r(2),rt3r(2),rt4r(2)
      double precision Lx,Ly

c AA' stacking from
c http://research.physics.berkeley.edu/zettl/pdf/
c 544.Alternative%20Stacking%20Sequences%20in%20hbn.pdf

c swap A and B atoms on top layer

      a1(1)= 2.504545468d0*dsqrt(3.d0)/2.d0
      a1(2)= 2.504545468d0/2.d0
      a2(1)= a1(1)
      a2(2)=-a1(2)
      tau1(1)=0.d0
      tau1(2)=0.d0
      tau2(1)= 2.504545468d0/dsqrt(3.d0)
      tau2(2)= 0.d0
      theta=5.0d0*(4.d0*datan(1.d0)/180.d0)

      a1r(1)= dcos(theta)*a1(1)+dsin(theta)*a1(2)
      a1r(2)=-dsin(theta)*a1(1)+dcos(theta)*a1(2)
      a2r(1)= dcos(theta)*a2(1)+dsin(theta)*a2(2)
      a2r(2)=-dsin(theta)*a2(1)+dcos(theta)*a2(2)

      tau1r(1)=tau1(1)
      tau1r(2)=tau1(2)
      tau2r(1)= dcos(theta)*(tau2(1))+dsin(theta)*(tau2(2))
      tau2r(2)=-dsin(theta)*(tau2(1))+dcos(theta)*(tau2(2))

      Lx=2.504545468d0*dsqrt(3.d0)*19.d0
      Ly=2.504545468d0*19.d0

      write(*,*) 'Lx',Lx,'Ly',Ly

      ra1(1)=2.504545468d0*dsqrt(3.d0)
      ra1(2)=0.d0
      ra2(1)=0.d0
      ra2(2)=2.504545468d0

      ra1r(1)= dcos(theta)*ra1(1)+dsin(theta)*ra1(2)
      ra1r(2)=-dsin(theta)*ra1(1)+dcos(theta)*ra1(2)
      ra2r(1)= dcos(theta)*ra2(1)+dsin(theta)*ra2(2)
      ra2r(2)=-dsin(theta)*ra2(1)+dcos(theta)*ra2(2)


      rt1(1)=tau1(1)
      rt1(2)=tau1(2)
      rt2(1)=tau2(1)
      rt2(2)=tau2(2)
      rt3(1)=tau1(1)+a1(1)
      rt3(2)=tau1(2)+a1(2)
      rt4(1)=tau2(1)+a1(1)
      rt4(2)=tau2(2)+a1(2)

      rt1r(1)= dcos(theta)*rt1(1)+dsin(theta)*rt1(2)
      rt1r(2)=-dsin(theta)*rt1(1)+dcos(theta)*rt1(2)
      rt2r(1)= dcos(theta)*rt2(1)+dsin(theta)*rt2(2)
      rt2r(2)=-dsin(theta)*rt2(1)+dcos(theta)*rt2(2)
      rt3r(1)= dcos(theta)*rt3(1)+dsin(theta)*rt3(2)
      rt3r(2)=-dsin(theta)*rt3(1)+dcos(theta)*rt3(2)
      rt4r(1)= dcos(theta)*rt4(1)+dsin(theta)*rt4(2)
      rt4r(2)=-dsin(theta)*rt4(1)+dcos(theta)*rt4(2)



      write(*,*) a1(1),a1r(1)
      write(*,*) a1(2),a1r(2)

      indice=1

      open(10,file="coords.data",status="replace")
      open(20,file="1p0.xyz",status="replace")

c bottom layer:
      DO i=1,19
       Do j=1,19
        do k=1,2
         r(k)=dble(i-1)*ra1(k)+dble(j-1)*ra2(k)+rt1(k)
        enddo
c bring atom back to box:
c        if(r(1).lt.0.d0)then
c          r(1)=r(1)+Lx
c        elseif(r(1).gt.Lx)then
c          r(1)=r(1)-Lx
c        endif
c        if(r(2).lt.0.d0)then
c          r(2)=r(2)+Ly
c        elseif(r(2).gt.Ly)then
c          r(2)=r(2)-Ly
c        endif
        write(10,'(3I10, 1F8.3, 3F17.8)') indice, 1, 1,-0.420, 
     .                       r(1), r(2), 10.d0
        write(20,'(A5, 3F12.5)') 'N', r(1),r(2),10.d0
        indice=indice+1

        do k=1,2
         r(k)=dble(i-1)*ra1(k)+dble(j-1)*ra2(k)+rt2(k)
        enddo
c bring atom back to box:
c        if(r(1).le.0.d0)then
c          r(1)=r(1)+Lx
c        elseif(r(1).gt.Lx)then
c          r(1)=r(1)-Lx
c        endif
c        if(r(2).le.0.d0)then
c          r(2)=r(2)+Ly
c        elseif(r(2).gt.Ly)then
c          r(2)=r(2)-Ly
c        endif
        write(10,'(3I10, 1F8.3, 3F17.8)') indice, 1, 2,0.420, 
     .                       r(1),r(2), 10.d0
        write(20,'(A5, 3F12.5)') 'B', r(1),r(2),10.d0
        indice=indice+1

        do k=1,2
         r(k)=dble(i-1)*ra1(k)+dble(j-1)*ra2(k)+rt3(k)
        enddo
        write(10,'(3I10, 1F8.3, 3F17.8)') indice, 1, 1,-0.420,
     .                       r(1),r(2), 10.d0
        write(20,'(A5, 3F12.5)') 'N', r(1),r(2),10.d0
        indice=indice+1

        do k=1,2
         r(k)=dble(i-1)*ra1(k)+dble(j-1)*ra2(k)+rt4(k)
        enddo
        write(10,'(3I10, 1F8.3, 3F17.8)') indice, 1, 2,0.420,
     .                       r(1),r(2), 10.d0
        write(20,'(A5, 3F12.5)') 'B', r(1),r(2),10.d0
        indice=indice+1

       EndDo
      ENDDO

c bottom layer+1:
      DO i=1,19
       Do j=1,19
        do k=1,2
         r(k)=dble(i-1)*ra1(k)+dble(j-1)*ra2(k)+rt1(k)
        enddo
c bring atom back to box:
c        if(r(1).lt.0.d0)then
c          r(1)=r(1)+Lx
c        elseif(r(1).gt.Lx)then
c          r(1)=r(1)-Lx
c        endif
c        if(r(2).lt.0.d0)then
c          r(2)=r(2)+Ly
c        elseif(r(2).gt.Ly)then
c          r(2)=r(2)-Ly
c        endif
        write(10,'(3I10, 1F8.3, 3F17.8)') indice, 2, 3, 0.420, 
     .                       r(1), r(2), 13.33d0
        write(20,'(A5, 3F12.5)') 'B', r(1),r(2),13.33d0
        indice=indice+1

        do k=1,2
         r(k)=dble(i-1)*ra1(k)+dble(j-1)*ra2(k)+rt2(k)
        enddo
c bring atom back to box:
c        if(r(1).le.0.d0)then
c          r(1)=r(1)+Lx
c        elseif(r(1).gt.Lx)then
c          r(1)=r(1)-Lx
c        endif
c        if(r(2).le.0.d0)then
c          r(2)=r(2)+Ly
c        elseif(r(2).gt.Ly)then
c          r(2)=r(2)-Ly
c        endif
        write(10,'(3I10, 1F8.3, 3F17.8)') indice, 2, 4,-0.420, 
     .                       r(1),r(2), 13.33d0
        write(20,'(A5, 3F12.5)') 'N', r(1),r(2),13.33d0
        indice=indice+1

        do k=1,2
         r(k)=dble(i-1)*ra1(k)+dble(j-1)*ra2(k)+rt3(k)
        enddo
        write(10,'(3I10, 1F8.3, 3F17.8)') indice, 2, 3, 0.420,
     .                       r(1),r(2), 13.33d0
        write(20,'(A5, 3F12.5)') 'B', r(1),r(2),13.33d0
        indice=indice+1

        do k=1,2
         r(k)=dble(i-1)*ra1(k)+dble(j-1)*ra2(k)+rt4(k)
        enddo
        write(10,'(3I10, 1F8.3, 3F17.8)') indice, 2, 4,-0.420,
     .                       r(1),r(2), 13.33d0
        write(20,'(A5, 3F12.5)') 'N', r(1),r(2),13.33d0
        indice=indice+1

       EndDo
      ENDDO

c bottom layer+3:
      DO i=1,19
       Do j=1,19
        do k=1,2
         r(k)=dble(i-1)*ra1r(k)+dble(j-1)*ra2r(k)+rt1(k)
        enddo
c bring atom back to box:
        if(r(1).lt.0.d0)then
          r(1)=r(1)+Lx
        elseif(r(1).gt.Lx)then
          r(1)=r(1)-Lx
        endif
        if(r(2).lt.0.d0)then
          r(2)=r(2)+Ly
        elseif(r(2).gt.Ly)then
          r(2)=r(2)-Ly
        endif
        write(10,'(3I10, 1F8.3, 3F17.8)') indice, 4, 7,-0.420, 
     .                       r(1), r(2), 19.99d0
        write(20,'(A5, 3F12.5)') 'N', r(1),r(2),19.99d0
        indice=indice+1

        do k=1,2
         r(k)=dble(i-1)*ra1r(k)+dble(j-1)*ra2r(k)+rt2(k)
        enddo
c bring atom back to box:
        if(r(1).le.0.d0)then
          r(1)=r(1)+Lx
        elseif(r(1).gt.Lx)then
          r(1)=r(1)-Lx
        endif
        if(r(2).le.0.d0)then
          r(2)=r(2)+Ly
        elseif(r(2).gt.Ly)then
          r(2)=r(2)-Ly
        endif
       write(10,'(3I10, 1F8.3, 3F17.8)') indice, 4, 8,0.420, 
     .                       r(1),r(2), 19.99d0
        write(20,'(A5, 3F12.5)') 'B', r(1),r(2),19.99d0
        indice=indice+1

        do k=1,2
         r(k)=dble(i-1)*ra1r(k)+dble(j-1)*ra2r(k)+rt3(k)
        enddo
c bring atom back to box:
        if(r(1).le.0.d0)then
          r(1)=r(1)+Lx
        elseif(r(1).gt.Lx)then
          r(1)=r(1)-Lx
        endif
        if(r(2).le.0.d0)then
          r(2)=r(2)+Ly
        elseif(r(2).gt.Ly)then
          r(2)=r(2)-Ly
        endif
        write(10,'(3I10, 1F8.3, 3F17.8)') indice, 4, 7, -0.420,
     .                       r(1),r(2), 19.99d0
        write(20,'(A5, 3F12.5)') 'N', r(1),r(2),19.99d0
        indice=indice+1

        do k=1,2
         r(k)=dble(i-1)*ra1r(k)+dble(j-1)*ra2r(k)+rt4(k)
        enddo
c bring atom back to box:
        if(r(1).le.0.d0)then
          r(1)=r(1)+Lx
        elseif(r(1).gt.Lx)then
          r(1)=r(1)-Lx
        endif
        if(r(2).le.0.d0)then
          r(2)=r(2)+Ly
        elseif(r(2).gt.Ly)then
          r(2)=r(2)-Ly
        endif
        write(10,'(3I10, 1F8.3, 3F17.8)') indice, 4, 8, 0.420,
     .                       r(1),r(2), 19.99d0
        write(20,'(A5, 3F12.5)') 'B', r(1),r(2),19.99d0
        indice=indice+1

       EndDo
      ENDDO

! Adding the flake with 0.4 degree angle

c flake layer:
      DO i=1,13
       Do j=1,13
        do k=1,2
         r(k)=dble(i-1)*ra1r(k)+dble(j-1)*ra2r(k)+rt1r(k)+3.0d0*ra1(k)
     .                                                   +4.0d0*ra2(k)
        enddo
c bring atom back to box:
c        if(r(1).lt.0.d0)then
c          r(1)=r(1)+Lx
c        elseif(r(1).gt.Lx)then
c          r(1)=r(1)-Lx
c        endif
c        if(r(2).lt.0.d0)then
c          r(2)=r(2)+Ly
c        elseif(r(2).gt.Ly)then
c          r(2)=r(2)-Ly
c        endif
        write(10,'(3I10, 1F8.3, 3F17.8)') indice, 3, 5, 0.420,
     .                       r(1),r(2), 16.66d0
        write(20,'(A5, 3F12.5)') 'B', r(1),r(2),16.66d0
        indice=indice+1

        do k=1,2
         r(k)=dble(i-1)*ra1r(k)+dble(j-1)*ra2r(k)+rt2r(k)+3.0d0*ra1(k)+
     .                                                    4.0d0*ra2(k)
        enddo
c bring atom back to box:
c        if(r(1).le.0.d0)then
c          r(1)=r(1)+Lx
c        elseif(r(1).gt.Lx)then
c          r(1)=r(1)-Lx
c        endif
c        if(r(2).le.0.d0)then
c          r(2)=r(2)+Ly
c        elseif(r(2).gt.Ly)then
c          r(2)=r(2)-Ly
c        endif
        write(10,'(3I10, 1F8.3, 3F17.8)') indice, 3, 6,-0.420,
     .                       r(1),r(2), 16.66d0
       write(20,'(A5, 3F12.5)') 'N',r(1),r(2),16.66d0
        indice=indice+1

        do k=1,2
         r(k)=dble(i-1)*ra1r(k)+dble(j-1)*ra2r(k)+rt3r(k)+3.0d0*ra1(k)+
     .                                                    4.0d0*ra2(k)
        enddo
c bring atom back to box:
c        if(r(1).lt.0.d0)then
c          r(1)=r(1)+Lx
c        elseif(r(1).gt.Lx)then
c          r(1)=r(1)-Lx
c        endif
c        if(r(2).lt.0.d0)then
c          r(2)=r(2)+Ly
c        elseif(r(2).gt.Ly)then
c          r(2)=r(2)-Ly
c        endif
        write(10,'(3I10, 1F8.3, 3F17.8)') indice, 3, 5, 0.420,
     .                      r(1),r(2), 16.66d0
        write(20,'(A5, 3F12.5)') 'B', r(1),r(2),16.66d0
        indice=indice+1

        do k=1,2
         r(k)=dble(i-1)*ra1r(k)+dble(j-1)*ra2r(k)+rt4r(k)+3.0d0*ra1(k)+
     .                                                    4.0d0*ra2(k)
        enddo
c bring atom back to box:
c        if(r(1).le.0.d0)then
c          r(1)=r(1)+Lx
c        elseif(r(1).gt.Lx)then
c          r(1)=r(1)-Lx
c        endif
c        if(r(2).le.0.d0)then
c          r(2)=r(2)+Ly
c        elseif(r(2).gt.Ly)then
c          r(2)=r(2)-Ly
c        endif
        write(10,'(3I10, 1F8.3, 3F17.8)') indice, 3, 6, -0.420,
     .                       r(1),r(2), 16.66d0
       write(20,'(A5, 3F12.5)') 'N',r(1),r(2),16.66d0
        indice=indice+1

       EndDo
      ENDDO



      close(10)
      close(20)

      end program
