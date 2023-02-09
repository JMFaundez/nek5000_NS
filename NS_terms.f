c-----------------------------------------------------------------------
c-----------------------------------------------------------------------

      subroutine comp_res_main
      include 'SIZE'
      include 'SOLN'
      include 'MASS'
      include 'INPUT'
      include 'MYNS'

      real*8 dudt(lx1,ly1,lz1,lelv), dvdt(lx1,ly1,lz1,lelv),
     &       dwdt(lx1,ly1,lz1,lelv)
      real*8 fux(lx1,ly1,lz1,lelv), fuy(lx1,ly1,lz1,lelv),
     &       fuz(lx1,ly1,lz1,lelv)
      real*8 dpdx(lx1,ly1,lz1,lelv), dpdy(lx1,ly1,lz1,lelv),
     &       dpdz(lx1,ly1,lz1,lelv)
      real*8 convu(lx1,ly1,lz1,lelv), convv(lx1,ly1,lz1,lelv),
     &       convw(lx1,ly1,lz1,lelv)
      real*8 lapu(lx1,ly1,lz1,lelv), lapv(lx1,ly1,lz1,lelv),
     &       lapw(lx1,ly1,lz1,lelv)
      real*8 resu(lx1,ly1,lz1,lelv), rhsu(lx1,ly1,lz1,lelv),
     &       resv(lx1,ly1,lz1,lelv), rhsv(lx1,ly1,lz1,lelv),
     &       resw(lx1,ly1,lz1,lelv), rhsw(lx1,ly1,lz1,lelv)
      real*8 normres(lx1,ly1,lz1,lelv), absdudt(lx1,ly1,lz1,lelv)
      real*8 allvar(lx1*ly1*lz1*lelv,10)
      real*8 wn1
      ntot1 = lx1*ly1*lz1*nelv


      if(nid.eq.0.and.verb1)  write(*,*) 'Computing dudt'
      call dveldt(dudt,dvdt,dwdt)
      if(nid.eq.0.and.verb1)  write(*,*) 'Computing user forcing'
      call ext_forc(fux,fuy,fuz)
      if(nid.eq.0.and.verb1)  write(*,*) 'Computing dpdx'
      call gradp(dpdx,dpdy,dpdz)
      if(nid.eq.0.and.verb1)  write(*,*) 'Computing convective term'
      call conv_term(convu,convv,convw)
      if(nid.eq.0.and.verb1)  write(*,*) 'Computing Laplacian'
      call laplu(lapu,lapv,lapw)

      call compute_rhs(rhsu,dpdx,lapu,convu,fux)
      call compute_rhs(rhsv,dpdy,lapv,convv,fuy)
      if (if3d) then
        call compute_rhs(rhsw,dpdz,lapw,convw,fuz)
      endif

      if(nid.eq.0.and.verb1)  write(*,*) 'Computing residuals'
      call sub3(resu,dudt,rhsu,ntot1)
      call sub3(resv,dvdt,rhsv,ntot1)
      if (if3d) then
        call sub3(resw,dwdt,rhsw,ntot1)
        call vsq3d(normres,resu,resv,resw,ntot1)
        call vsq3d(absdudt,dudt,dvdt,dwdt,ntot1)
      else
        call rzero(resw,ntot1)
        call vsq2d(normres,resu,resv,ntot1)
        call vsq2d(absdudt,dudt,dvdt,ntot1)
      endif
      call vsqrt(normres,ntot1) 
      wn1 = gl2norm(absdudt,ntot1)
      call cmult(normres,1./wn1,ntot1) 

      if(nid.eq.0.and.verb1)  write(*,*) 'Saving values'
      call copy(allvar(1,1),normres,ntot1)  !temp
      call copy(allvar(1,2),rhsu,ntot1)  !s1
      call copy(allvar(1,3),rhsv,ntot1)  !s2
      call copy(allvar(1,4),rhsw,ntot1)  !s3
      call copy(allvar(1,5),dudt,ntot1)  !s4
      call copy(allvar(1,6),dvdt,ntot1)  !s5
      call copy(allvar(1,7),dwdt,ntot1)  !s6
      call copy(allvar(1,8),convu,ntot1) !s7 
      call copy(allvar(1,9),convv,ntot1) !s8 
      call copy(allvar(1,10),convw,ntot1) !s9 
      if(nid.eq.0.and.verb1)  write(*,*) 'Values copied'
      call outpost2(vx,vy,vz,pr,allvar,10,'res')
      if(nid.eq.0.and.verb1)  write(*,*) 'File saved'

      !vmax = glmax(resu,ntot1)
      !vmin = glmin(resu,ntot1)
      !tnorm = gl2norm(resu,ntot1)

      !wn1 = gl2norm(dudt,ntot1)
      !wn2 = gl2norm(fux,ntot1)
      !wn3 = gl2norm(dpdx,ntot1)
      !wn4 = gl2norm(convu,ntot1)
      !wn5 = gl2norm(lapu,ntot1)

      !if(nid.eq.0.and.verb1) write(*,*) 'Weak Max res:', vmax
      !if(nid.eq.0.and.verb1) write(*,*) 'Weak Min res:', vmin
      !if(nid.eq.0.and.verb1) write(*,*) 'Weak norm res:', tnorm
      !if(nid.eq.0.and.verb1) write(*,*) 'Weak norm dudt:', wn1
      !if(nid.eq.0.and.verb1) write(*,*) 'Weak norm fux:', wn2
      !if(nid.eq.0.and.verb1) write(*,*) 'Weak norm dpdx:', wn3
      !if(nid.eq.0.and.verb1) write(*,*) 'Weak norm convu:', wn4
      !if(nid.eq.0.and.verb1) write(*,*) 'Weak norm lapu:', wn5

      return
      end subroutine

c-----------------------------------------------------------------------
c-----------------------------------------------------------------------

      subroutine compute_rhs(rhs,v1,v2,v3,v4)
      include 'SIZE'
      real*8 rhs(lx1,ly1,lz1,lelv)
      real*8 v1(lx1,ly1,lz1,lelv),v2(lx1,ly1,lz1,lelv),
     &       v3(lx1,ly1,lz1,lelv),v4(lx1,ly1,lz1,lelv)
      
      ntot1 = lx1*ly1*lz1*nelv
      !           (1)       (2)       (3)     (4)
      ! dudt = M-1.DT.p + M-1.A.u - conv(u) + fux
      call add3(rhs,v1,v2,ntot1) 
      call sub2(rhs,v3,ntot1)
      call add2(rhs,v4,ntot1)

      end subroutine

c-----------------------------------------------------------------------
c-----------------------------------------------------------------------

      subroutine ext_forc(fux,fuy,fuz)
      include 'SIZE'
      include 'SOLN' ! BFX, BFY, BFZ
      real*8 fux(lx1,ly1,lz1,lelv), fuy(lx1,ly1,lz1,lelv),
     &       fuz(lx1,ly1,lz1,lelv)
      ntot1 = lx1*ly1*lz1*nelv
      call makeuf
      call opcopy(fux,fuy,fuz,bfx,bfy,bfz)
      
      end subroutine

c-----------------------------------------------------------------------
c-----------------------------------------------------------------------

      subroutine gradp(dpdx,dpdy,dpdz)
      include 'SIZE'
      include 'SOLN'
      include 'MASS'
      include 'INPUT'
      include 'MYNS'
      
      common /scrns/ w1(lx1,ly1,lz1,lelv)
     $ ,             w2(lx1,ly1,lz1,lelv)
     $ ,             w3(lx1,ly1,lz1,lelv)
      common /scrvh/ h2(lx1,ly1,lz1,lelv)
      common /scrhi/ h2inv(lx1,ly1,lz1,lelv)

      real*8 dpdx(lx1,ly1,lz1,lelv), dpdy(lx1,ly1,lz1,lelv),
     &       dpdz(lx1,ly1,lz1,lelv)
      
      ! Copied from subroutine incomprn in NEK's core
      call copy    (h2,vtrans(1,1,1,1,1),ntot1)
      call invers2 (h2inv,h2,ntot1)

      call opgradt(w1,w2,w3,pr) !Gradient of p from pressure mesh to velocity mesh
      call opbinv (dpdx,dpdy,dpdz,w1,w2,w3,h2inv) !Inverted mass matrix

      end subroutine

c-----------------------------------------------------------------------
c-----------------------------------------------------------------------

      subroutine laplu(lapu,lapv,lapw)
      include 'SIZE'
      include 'SOLN'
      include 'INPUT'
      include 'MASS'
      include 'MYNS'
      real*8 lapu(lx1,ly1,lz1,lelv), lapv(lx1,ly1,lz1,lelv),
     &       lapw(lx1,ly1,lz1,lelv)
      real*8 helm2(lx1,ly1,lz1,lelv)

      common /scrns/ w1(lx1,ly1,lz1,lelv), w2(lx1,ly1,lz1,lelv),
     $               w3(lx1,ly1,lz1,lelv)
      common /scrvh/ h2(lx1,ly1,lz1,lelv)
      common /scrhi/ h2inv(lx1,ly1,lz1,lelv)

      ntot1 = lx1*ly1*lz1*nelv
      ntot = nx1*ny1*nz1*nelv
      
      imesh=1
      call rzero(h2,ntot1)
      if (nid.eq.0) write(*,*) 'VDIFF:', vdiff(1,1,1,1,1)
      call wlaplacian(w1,vx,vdiff,1)
      call wlaplacian(w2,vy,vdiff,1)
      if (if3d) call wlaplacian(w3,vz,vdiff,1)

      call opbinv (lapu,lapv,lapw,w1,w2,w3,h2inv) !Inverted mass matrix
      end subroutine


c-----------------------------------------------------------------------
c-----------------------------------------------------------------------

      subroutine dveldt(dudt,dvdt,dwdt)
      include 'SIZE'
      include 'TSTEP'
      include 'SOLN'
      include 'INPUT'
      ! using these as working arrays
      common /SCRNS/ h2(lx1,ly1,lz1,lelv),
     &              TA1(lx1,ly1,lz1,lelv),
     &              TA2(lx1,ly1,lz1,lelv),
     &              TA3(lx1,ly1,lz1,lelv)
      integer ilag
      real*8 mybd(3)

      ntot1 = lx1*ly1*lz1*nelv
      const = 1./dt
!      if(nid.eq.0) write(*,*) 'const', const
      
      mybd(1) = 3./2
      mybd(2) = -4./2
      mybd(3) = 1./2
      ! h2=rho/dt
      call cmult2(h2,vtrans,const,ntot1)

      call copy(dudt,vx,ntot1)
      call copy(dvdt,vy,ntot1)
      call cmult(dudt,mybd(1),ntot1)
      call cmult(dvdt,mybd(1),ntot1)
      if (if3d) then
        call copy(dwdt,vz,ntot1)
        call cmult(dwdt,mybd(1),ntot1)
      endif
       
      do ilag=1,2
        call copy(TA1,vxlag(1,1,1,1,ilag),ntot1)
        call copy(TA2,vylag(1,1,1,1,ilag),ntot1)
        call cmult(TA1,mybd(ilag+1),ntot1)
        call cmult(TA2,mybd(ilag+1),ntot1)
        if (if3d) then
          call copy(TA3,vzlag(1,1,1,1,ilag),ntot1)
          call cmult(TA3,mybd(ilag+1),ntot1)
        else
          call rzero(TA3,ntot1)
        endif
        call opadd2(dudt,dvdt,dwdt,TA1,TA2,TA3)
      enddo
      call opcolv(dudt,dvdt,dwdt,h2)

      end subroutine
      
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------

      subroutine conv_term(convu,convv,convw)
      include 'SIZE'
      include 'SOLN'
      include 'INPUT'
      real*8 convu(lx1,ly1,lz1,lelv), convv(lx1,ly1,lz1,lelv),
     &       convw(lx1,ly1,lz1,lelv)
      ntot1 = lx1*ly1*lz1*nelv

      call convop(convu,vx)
      call convop(convv,vy)
      if (if3d) call convop(convw,vz)

      end subroutine                                     
                                                                       
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------

      subroutine vsq2d(a,b,c,n)
      REAL A(1),B(1),C(1),D(1)
    
      include 'OPCTR'
    
      DO 100 i=1,n
          a(i)=b(i)**2+c(i)**2
  100 CONTINUE
      return
      END

      subroutine vsq3d(a,b,c,d,n)
      REAL A(1),B(1),C(1),D(1)

      include 'OPCTR'

      DO 100 i=1,n
          a(i)=b(i)**2+c(i)**2+d(i)**2
  100 CONTINUE
      return
      END