c-----------------------------------------------------------------------
c-----------------------------------------------------------------------

      subroutine strong_ns
      include 'SIZE'
      include 'SOLN'

      real*8 dudt(lx1,ly1,lz1,lelv), dvdt(lx1,ly1,lz1,lelv),
     &         dwdt(lx1,ly1,lz1,lelv)
      real*8 fux(lx1,ly1,lz1,lelv), fuy(lx1,ly1,lz1,lelv),
     &         fuz(lx1,ly1,lz1,lelv)
      real*8 dpdx(lx1,ly1,lz1,lelv), dpdy(lx1,ly1,lz1,lelv),
     &         dpdz(lx1,ly1,lz1,lelv)
      real*8 convu(lx1,ly1,lz1,lelv), convv(lx1,ly1,lz1,lelv),
     &         convw(lx1,ly1,lz1,lelv)
      real*8 lapu(lx1,ly1,lz1,lelv), lapv(lx1,ly1,lz1,lelv),
     &        lapw(lx1,ly1,lz1,lelv)
      real*8 resu(lx1,ly1,lz1,lelv), rhsu(lx1,ly1,lz1,lelv),
     &       resv(lx1,ly1,lz1,lelv), rhsv(lx1,ly1,lz1,lelv),
     &       resw(lx1,ly1,lz1,lelv), rhsw(lx1,ly1,lz1,lelv)
      real*8 allvar(lx1*ly1*lz1*lelv,7)
      logical verb1

      
      !write(*,*) nbd
      !write(*,*) nfield
      ntot1 = lx1*ly1*lz1*nelv
      verb1 = .true.

      if(nid.eq.0.and.verb1)  write(*,*) 'Computing dudt'
      call sdveldt(dudt,dvdt,dwdt)
      if(nid.eq.0.and.verb1)  write(*,*) 'Computing user forcing'
      call sext_forc(fux,fuy,fuz)
      if(nid.eq.0.and.verb1)  write(*,*) 'Computing dpdx'
      call sdpdx(dpdx,dpdy,dpdz)
      if(nid.eq.0.and.verb1)  write(*,*) 'Computing convective term'
      call sconvu(convu,convv,convw)
      if(nid.eq.0.and.verb1)  write(*,*) 'Computing Laplacian'
      call slapu(lapu,lapv,lapw)

      ! dudt=-dpdx+lapu+fux-con(u)
      call compute_srhs(rhsu,dpdx,lapu,fux,convu)
      call sub3(resu,dudt,rhsu,ntot1)
      if(nid.eq.0.and.verb1)  write(*,*) 'strong-NS computed'
      
      call copy(allvar(1,1),resu,ntot1)  !temp
      call copy(allvar(1,2),rhsu,ntot1)  !s1
      call copy(allvar(1,3),dudt,ntot1)  !s2
      call copy(allvar(1,4),fux,ntot1)   !s3 
      call copy(allvar(1,5),dpdx,ntot1)  !s4
      call copy(allvar(1,6),convu,ntot1) !s5 
      call copy(allvar(1,7),lapu,ntot1)  !s6 
      if(nid.eq.0.and.verb1)  write(*,*) 'Files copied'
      call outpost2(vx,vy,vz,pr,allvar,7,'res')
      !call outpost(lapu,convu,dpdx,pr,t,'res')
      if(nid.eq.0.and.verb1)  write(*,*) 'file saved'

      vmax = glmax(resu,ntot1)
      vmin = glmin(resu,ntot1)
      tnorm = gl2norm(resu,ntot1)

      if(nid.eq.0.and.verb1) write(*,*) 'Strong Max res:', vmax
      if(nid.eq.0.and.verb1) write(*,*) 'Strong Min res:', vmin
      if(nid.eq.0.and.verb1) write(*,*) 'Strong norm res:', tnorm

      return
      end subroutine

c-----------------------------------------------------------------------
c-----------------------------------------------------------------------

c-----------------------------------------------------------------------
c-----------------------------------------------------------------------

      subroutine weak_ns
      include 'SIZE'
      include 'SOLN'

      real*8 dudt(lx1,ly1,lz1,lelv), dvdt(lx1,ly1,lz1,lelv),
     &         dwdt(lx1,ly1,lz1,lelv)
      real*8 fux(lx1,ly1,lz1,lelv), fuy(lx1,ly1,lz1,lelv),
     &         fuz(lx1,ly1,lz1,lelv)
      real*8 dpdx(lx1,ly1,lz1,lelv), dpdy(lx1,ly1,lz1,lelv),
     &         dpdz(lx1,ly1,lz1,lelv)
      real*8 convu(lx1,ly1,lz1,lelv), convv(lx1,ly1,lz1,lelv),
     &         convw(lx1,ly1,lz1,lelv)
      real*8 lapu(lx1,ly1,lz1,lelv), lapv(lx1,ly1,lz1,lelv),
     &        lapw(lx1,ly1,lz1,lelv)
      real*8 resu(lx1,ly1,lz1,lelv), rhsu(lx1,ly1,lz1,lelv),
     &       resv(lx1,ly1,lz1,lelv), rhsv(lx1,ly1,lz1,lelv),
     &       resw(lx1,ly1,lz1,lelv), rhsw(lx1,ly1,lz1,lelv)
      real*8 allvar(lx1*ly1*lz1*lelv,7)
      logical verb1
      !write(*,*) nbd
      !write(*,*) nfield
      ntot1 = lx1*ly1*lz1*nelv
      verb1 = .true.


      if(nid.eq.0.and.verb1)  write(*,*) 'Computing dudt'
      call wdveldt(dudt,dvdt,dwdt)
      if(nid.eq.0.and.verb1)  write(*,*) 'Computing user forcing'
      call wext_forc(fux,fuy,fuz)
      if(nid.eq.0.and.verb1)  write(*,*) 'Computing dpdx'
      call wdpdx(dpdx,dpdy,dpdz)
      if(nid.eq.0.and.verb1)  write(*,*) 'Computing convective term'
      call wconvu(convu,convv,convw)
      if(nid.eq.0.and.verb1)  write(*,*) 'Computing Laplacian'
      call wlapu(lapu,lapv,lapw)

      ! Mdudt=Dp-Ku+fux-Cu
      call compute_wrhs(rhsu,dpdx,lapu,fux,convu)
      call sub3(resu,dudt,rhsu,ntot1)
      if(nid.eq.0.and.verb1)  write(*,*) 'weak-NS computed'

      call copy(allvar(1,1),resu,ntot1)  !temp
      call copy(allvar(1,2),rhsu,ntot1)  !s1
      call copy(allvar(1,3),dudt,ntot1)  !s2
      call copy(allvar(1,4),fux,ntot1)   !s3 
      call copy(allvar(1,5),dpdx,ntot1)  !s4
      call copy(allvar(1,6),convu,ntot1) !s5 
      call copy(allvar(1,7),lapu,ntot1)  !s6 
      call outpost2(vx,vy,vz,pr,allvar,7,'res')
      if(nid.eq.0.and.verb1)  write(*,*) 'file saved'

      vmax = glmax(resu,ntot1)
      vmin = glmin(resu,ntot1)
      tnorm = gl2norm(resu,ntot1)

      if(nid.eq.0.and.verb1) write(*,*) 'Weak Max res:', vmax
      if(nid.eq.0.and.verb1) write(*,*) 'Weak Min res:', vmin
      if(nid.eq.0.and.verb1) write(*,*) 'Weak norm res:', tnorm

      return
      end subroutine

c-----------------------------------------------------------------------
c-----------------------------------------------------------------------

      subroutine compute_wrhs(rhs,v1,v2,v3,v4)
      include 'SIZE'
      real*8 rhs(lx1,ly1,lz1,lelv)
      real*8 v1(lx1,ly1,lz1,lelv),v2(lx1,ly1,lz1,lelv),
     &       v3(lx1,ly1,lz1,lelv),v4(lx1,ly1,lz1,lelv)
      
      ntot1 = lx1*ly1*lz1*nelv
      ! Mdudt=Dp-Ku+fux-Cu

      !call chsign(v1,ntot1) 
      call chsign(v2,ntot1) 
      call chsign(v4,ntot1) 
      call add3(rhs,v1,v3,ntot1) 
      call add2(rhs,v2,ntot1)
      call add2(rhs,v4,ntot1)


      end subroutine

c-----------------------------------------------------------------------
c-----------------------------------------------------------------------

      subroutine compute_srhs(rhs,v1,v2,v3,v4)
      include 'SIZE'
      real*8 rhs(lx1,ly1,lz1,lelv)
      real*8 v1(lx1,ly1,lz1,lelv),v2(lx1,ly1,lz1,lelv),
     &       v3(lx1,ly1,lz1,lelv),v4(lx1,ly1,lz1,lelv)
      
      ntot1 = lx1*ly1*lz1*nelv
      ! dudt=-dpdx+lapu+fux-con(u)

      !call chsign(v1,ntot1) 
      call add3(rhs,v3,v2,ntot1) 
      call sub2(rhs,v1,ntot1)
      call sub2(rhs,v4,ntot1)


      end subroutine

c-----------------------------------------------------------------------
c-----------------------------------------------------------------------

      subroutine wdveldt(dudt,dvdt,dwdt)
      include 'SIZE'
      include 'TSTEP'
      include 'SOLN'
      include 'MASS'
      ! using these as working arrays
!      common /SCRNS/ h2(lx1,ly1,lz1,lelv),
!     &              TA1(lx1,ly2,lz1,lelv),
!     &              TA2(lx1,ly2,lz1,lelv),
!     &              TA3(lx1,ly2,lz1,lelv)
      integer ilag
      real*8 dudt(lx1,ly1,lz1,lelv), dvdt(lx1,ly1,lz1,lelv),
     &         dwdt(lx1,ly1,lz1,lelv)
      real*8 TA1(lx1,ly1,lz1,lelv), TA2(lx1,ly1,lz1,lelv),
     &         TA3(lx1,ly1,lz1,lelv), h2(lx1,ly1,lz1,lelv)
      real*8 dumm1(lx1,ly1,lz1,lelv)
      real*8 mybd(3)



      ntot1 = lx1*ly1*lz1*nelv
      const = 1./dt
!      ! check what's going with vlag
!      call sub3(dumm1,vx,vx,ntot1)
!      velnorm = gl2norm(dumm1,ntot1)
!      if(nid.eq.0) write(*,*) 'norm vel dif:', velnorm
!      
!      call sub3(dumm1,vx,vxlag(1,1,1,1,1),ntot1)
!      velnorm = gl2norm(dumm1,ntot1)
!      if(nid.eq.0) write(*,*) 'norm vel dif:', velnorm
!
!      call sub3(dumm1,vxlag(1,1,1,1,1),vxlag(1,1,1,1,2),ntot1)
!      velnorm = gl2norm(dumm1,ntot1)
!      if(nid.eq.0) write(*,*) 'norm vel dif:', velnorm

      mybd(1) = 3./2
      mybd(2) = -4./2
      mybd(3) = 1./2
      ! h2=rho/dt
      call cmult2(h2,vtrans(1,1,1,1,1),const,ntot1)

      call opcolv3c(dudt,dvdt,dwdt,vx,
     &                             vy,
     &                             vz,
     &                             bm1, mybd(1))
       
      !if(nid.eq.0) write(*,*) bd(1)
      do ilag=1,2
      !   if(nid.eq.0) write(*,*) bd(ilag+1)
         call opcolv3c(TA1,TA2,TA3,vxlag(1,1,1,1,ilag),
     &                                 vylag(1,1,1,1,ilag),
     &                                 vzlag(1,1,1,1,ilag),
     &                                 bm1, mybd(ilag+1))
         call opadd2(dudt,dvdt,dwdt,TA1,TA2,TA3)
      enddo
      call opcolv(dudt,dvdt,dwdt,h2)

      end subroutine

c-----------------------------------------------------------------------
c-----------------------------------------------------------------------

      subroutine wext_forc(fux,fuy,fuz)
      include 'SIZE'
      include 'SOLN' ! BFX, BFY, BFZ
      real*8 fux(lx1,ly1,lz1,lelv), fuy(lx1,ly1,lz1,lelv),
     &         fuz(lx1,ly1,lz1,lelv)
      ntot1 = lx1*ly1*lz1*nelv
      call makeuf
      call opcopy(fux,fuy,fuz,bfx,bfy,bfz)
      
      end subroutine

c-----------------------------------------------------------------------
c-----------------------------------------------------------------------

      subroutine wdpdx(dpdx,dpdy,dpdz)
      include 'SIZE'
      include 'SOLN'
      real*8 dpdx(lx1,ly1,lz1,lelv), dpdy(lx1,ly1,lz1,lelv),
     &         dpdz(lx1,ly1,lz1,lelv)
      ntot1 = lx1*ly1*lz1*nelv

      call opgradt(dpdx,dpdy,dpdz,pr)

      end subroutine
      
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------

      subroutine wconvu(convu,convv,convw)
      include 'SIZE'
      include 'SOLN'
      include 'MASS'
      include 'INPUT'
      real*8 convu(lx1,ly1,lz1,lelv), convv(lx1,ly1,lz1,lelv),
     &         convw(lx1,ly1,lz1,lelv)
      ntot1 = lx1*ly1*lz1*nelv

      call convop(convu,vx)
      call col2(convu,bm1,ntot1)
      call convop(convv,vy)
      call col2(convv,bm1,ntot1)
      if (if3d) then
      call convop(convw,vz)
      call col2(convw,bm1,ntot1)
      endif

      end subroutine

c-----------------------------------------------------------------------
c-----------------------------------------------------------------------

      subroutine wlapu(lapu,lapv,lapw)
      include 'SIZE'
      include 'SOLN'
      include 'INPUT'
      real*8 lapu(lx1,ly1,lz1,lelv), lapv(lx1,ly1,lz1,lelv),
     &        lapw(lx1,ly1,lz1,lelv)
      ntot1 = lx1*ly1*lz1*nelv

      call wlaplacian(lapu,vx,vdiff,1)
      call wlaplacian(lapv,vy,vdiff,1)
      if (if3d) call wlaplacian(lapw,vz,vdiff,1)

      end subroutine


c-----------------------------------------------------------------------
c-----------------------------------------------------------------------

      subroutine sdveldt(dudt,dvdt,dwdt)
      include 'SIZE'
      include 'TSTEP'
      include 'SOLN'
      include 'INPUT'
      ! using these as working arrays
!      common /SCRNS/ h2(lx1,ly1,lz1,lelv),
!     &              TA1(lx1,ly2,lz1,lelv),
!     &              TA2(lx1,ly2,lz1,lelv),
!     &              TA3(lx1,ly2,lz1,lelv)
      integer ilag
      real*8 dudt(lx1,ly1,lz1,lelv), dvdt(lx1,ly1,lz1,lelv),
     &         dwdt(lx1,ly1,lz1,lelv)
      real*8 TA1(lx1,ly1,lz1,lelv), TA2(lx1,ly1,lz1,lelv),
     &         TA3(lx1,ly1,lz1,lelv), h2(lx1,ly1,lz1,lelv)
      real*8 mybd(3)

      ntot1 = lx1*ly1*lz1*nelv
      const = 1./dt
!      if(nid.eq.0) write(*,*) 'const', const
      
      mybd(1) = 3./2
      mybd(2) = -4./2
      mybd(3) = 1./2
      ! h2=rho/dt
      call cmult2(h2,vtrans(1,1,1,1,1),const,ntot1)

      call copy(dudt,vx,ntot1)
      call copy(dvdt,vy,ntot1)
      call cmult(dudt,mybd(1),ntot1)
      call cmult(dvdt,mybd(1),ntot1)
      if (if3d) then
      call copy(dwdt,vz,ntot1)
      call cmult(dwdt,mybd(1),ntot1)
      endif
       
      !if(nid.eq.0) write(*,*) bd(1)
      do ilag=1,2
!      !   if(nid.eq.0) write(*,*) bd(ilag+1)
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
!      call cmult(dudt,const,ntot1)
!      call cmult(dvdt,const,ntot1)
      !call cmult(dwdt,const,ntot1)

      end subroutine

c-----------------------------------------------------------------------
c-----------------------------------------------------------------------

      subroutine sext_forc(fux,fuy,fuz)
      include 'SIZE'
      include 'SOLN' ! BFX, BFY, BFZ
      real*8 fux(lx1,ly1,lz1,lelv), fuy(lx1,ly1,lz1,lelv),
     &         fuz(lx1,ly1,lz1,lelv)
      ntot1 = lx1*ly1*lz1*nelv

      call nekuf(fux,fuy,fuz)
      
      end subroutine

c-----------------------------------------------------------------------
c-----------------------------------------------------------------------

      subroutine sdpdx(dpdx,dpdy,dpdz)
      include 'SIZE'
      include 'SOLN'
      real*8 dpdx(lx1,ly1,lz1,lelv), dpdy(lx1,ly1,lz1,lelv),
     &         dpdz(lx1,ly1,lz1,lelv)
      real*8 TA1(lx1,ly1,lz1,lelv), TA2(lx1,ly1,lz1,lelv),
     &         TA3(lx1,ly1,lz1,lelv), h2(lx1,ly1,lz1,lelv)
      ! using these as working arrays
!      common /SCRNS/ TA1(lx1,ly2,lz1,lelv),
!     &               TA2(lx1,ly2,lz1,lelv),
!     &               TA3(lx1,ly2,lz1,lelv)
      ntot1 = lx1*ly1*lz1*nelv
      
      call mappr(TA1,pr,TA2,TA3)
      call gradm1(dpdx,dpdy,dpdz,TA1)

      end subroutine
      
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------

      subroutine sconvu(convu,convv,convw)
      include 'SIZE'
      include 'SOLN'
      include 'INPUT'
      real*8 convu(lx1,ly1,lz1,lelv), convv(lx1,ly1,lz1,lelv),
     &         convw(lx1,ly1,lz1,lelv)
      ntot1 = lx1*ly1*lz1*nelv

      call convop(convu,vx)
      call convop(convv,vy)
      if (if3d) call convop(convw,vz)

      end subroutine

c-----------------------------------------------------------------------
c-----------------------------------------------------------------------

      subroutine slapu(lapu,lapv,lapw)
      include 'SIZE'
      include 'SOLN'
      include 'INPUT'
      real*8 lapu(lx1,ly1,lz1,lelv), lapv(lx1,ly1,lz1,lelv),
     &        lapw(lx1,ly1,lz1,lelv)
      ntot1 = lx1*ly1*lz1*nelv

      call laplacian(lapu,vx)
      call laplacian(lapv,vy)
			call col2(lapu,vdiff,ntot1)
			call col2(lapv,vdiff,ntot1)
			if (if3d) then
      call laplacian(lapw,vz)
			call col2(lapw,vdiff,ntot1)
			endif

      end subroutine

c-----------------------------------------------------------------------
c-----------------------------------------------------------------------

      subroutine laplacian(lapu,up)                                      
      implicit none                                                      
!                                                                        
      include 'SIZE'                                                     
      include 'INPUT'           ! if3d                                   
      include 'DXYZ'            ! dxm1,d[xy]tm1                          
      include 'GEOM'            ! r[xy]m1,s[xy]m1,t[xy]m1,jacmi          
!                                                                        
      ! argument list                                                    
      real up (lx1*ly1*lz1*lelv)       ! perturbation velocity componen
!                                                                        
      ! output                                                           
      real lapu (lx1*ly1*lz1,lelv)                                       
!                                                                        
      ! local variables                                                  
      real ux  (lx1*ly1*lz1,lelv)                                        
     $ ,   uy  (lx1*ly1*lz1,lelv)                                        
     $ ,   uz  (lx1*ly1*lz1,lelv)                                        
     $ ,   ur  (lx1*ly1*lz1)                                             
     $ ,   us  (lx1*ly1*lz1)                                             
     $ ,   ut  (lx1*ly1*lz1)                                             
!                                                                        
      common /ctmp1/ ur,us,ut                                            
!                                                                        
      integer e,i,lxyz,nel,N                                             
!----------------------------------------------------------------------- 
      lxyz = lx1*ly1*lz1                                                  
      nel = lx1-1                                                         
      call gradm1(ux,uy,uz,up)                                            
      do e=1,lelt                                                         
        if (if3d) then                                                    
          call local_grad3(ur,us,ut,ux,nel,e,dxm1,dxtm1)                  
          do i=1,lxyz                                                     
            lapu(i,e) = jacmi(i,e)*(  ur(i)*rxm1(i,1,1,e)                 
     $                              + us(i)*sxm1(i,1,1,e)                 
     $                              + ut(i)*txm1(i,1,1,e) )               
          enddo                                                           
          call local_grad3(ur,us,ut,uy,nel,e,dxm1,dxtm1)                  
          do i=1,lxyz                                                     
            lapu(i,e) = lapu(i,e) + jacmi(i,e)*(  ur(i)*rym1(i,1,1,e)     
     $                                          + us(i)*sym1(i,1,1,e)     
     $                                          + ut(i)*tym1(i,1,1,e) )   
          enddo                                                           
          call local_grad3(ur,us,ut,uz,nel,e,dxm1,dxtm1)                  
          do i=1,lxyz                                                     
            lapu(i,e) = lapu(i,e) + jacmi(i,e)*(  ur(i)*rzm1(i,1,1,e)     
     $                                          + us(i)*szm1(i,1,1,e)     
     $                                          + ut(i)*tzm1(i,1,1,e) )   
          enddo                                                           
        else ! 2D                                                         
          call local_grad2(ur,us,ux,nel,e,dxm1,dytm1)                     
          do i=1,lxyz                                                     
            lapu(i,e) = jacmi(i,e)*(ur(i)*rxm1(i,1,1,e)                   
     $                            + us(i)*sxm1(i,1,1,e) )                 
          enddo                                                           
          call local_grad2(ur,us,uy,nel,e,dxm1,dytm1)                     
          do i=1,lxyz                                                     
            lapu(i,e) = lapu(i,e)                                         
     $                  + jacmi(i,e)*(ur(i)*rym1(i,1,1,e)                 
     $                              + us(i)*sym1(i,1,1,e) )               
          enddo                                                           
        endif ! if3d                                                      
      enddo                                                               
!                                                                         
      return                                                              
      end subroutine laplacian                                            
!                                                                         







