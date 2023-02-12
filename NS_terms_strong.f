c-----------------------------------------------------------------------
c-----------------------------------------------------------------------

      subroutine strong_ns
      include 'SIZE'
      include 'SOLN'
      include 'MYNS'

      real*8 dudt(lx1,ly1,lz1,lelv), dvdt(lx1,ly1,lz1,lelv),
     &         dwdt(lx1,ly1,lz1,lelv)
      real*8 fux(lx1,ly1,lz1,lelv), fuy(lx1,ly1,lz1,lelv),
     &         fuz(lx1,ly1,lz1,lelv)
      real*8 bcx(lx1,ly1,lz1,lelv), bcy(lx1,ly1,lz1,lelv),
     &         bcz(lx1,ly1,lz1,lelv)
      real*8 dpdx(lx1,ly1,lz1,lelv), dpdy(lx1,ly1,lz1,lelv),
     &         dpdz(lx1,ly1,lz1,lelv)
      real*8 convu(lx1,ly1,lz1,lelv), convv(lx1,ly1,lz1,lelv),
     &         convw(lx1,ly1,lz1,lelv)
      real*8 lapu(lx1,ly1,lz1,lelv), lapv(lx1,ly1,lz1,lelv),
     &        lapw(lx1,ly1,lz1,lelv)
      real*8 resu(lx1,ly1,lz1,lelv), rhsu(lx1,ly1,lz1,lelv),
     &       resv(lx1,ly1,lz1,lelv), rhsv(lx1,ly1,lz1,lelv),
     &       resw(lx1,ly1,lz1,lelv), rhsw(lx1,ly1,lz1,lelv)
      real*8 allvar(lx1*ly1*lz1*lelv,8)
      real*8 sn1,sn2,sn3,sn4,sn5
      real*8 L1,L2

      
      !write(*,*) nbd
      !write(*,*) nfield
      ntot1 = lx1*ly1*lz1*nelv

      if(nid.eq.0.and.verb1)  write(*,*) 'Computing dudt'
      call sdveldt(dudt,dvdt,dwdt)
      if(nid.eq.0.and.verb1)  write(*,*) 'Computing user forcing'
      call sext_forc(fux,fuy,fuz)
      if(nid.eq.0.and.verb1)  write(*,*) 'Computing BC forcing'
      call sBC_term(bcx,bcy,bcz)
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

      vmax = glmax(resu,ntot1)
      vmin = glmin(resu,ntot1)
      tnorm = gl2norm(resu,ntot1)
      sn1 = gl2norm(dudt,ntot1)
      sn2 = gl2norm(fux,ntot1)
      sn3 = gl2norm(dpdx,ntot1)
      sn4 = gl2norm(convu,ntot1)
      sn5 = gl2norm(lapu,ntot1)

      const = 1./sn1
      call cmult(resu,const,ntot1)
      
      call copy(allvar(1,1),resu,ntot1)  !temp
      call copy(allvar(1,2),rhsu,ntot1)  !s1
      call copy(allvar(1,3),dudt,ntot1)  !s2
      call copy(allvar(1,4),fux,ntot1)   !s3 
      call copy(allvar(1,5),dpdx,ntot1)  !s4
      call copy(allvar(1,6),convu,ntot1) !s5 
      call copy(allvar(1,7),lapu,ntot1)  !s6 
      call copy(allvar(1,8),bcx,ntot1)   !s7 
      if(nid.eq.0.and.verb1)  write(*,*) 'Files copied'
      call outpost2(vx,vy,vz,pr,allvar,8,'res')
      !call outpost(lapu,convu,dpdx,pr,t,'res')
      if(nid.eq.0.and.verb1)  write(*,*) 'file saved'

      call my_norm(L1,L2,rhsu)

      if(nid.eq.0) write(*,*)'L1,L2:', L1,L2


      if(nid.eq.0.and.verb1) write(*,*) 'Strong Max res:', vmax
      if(nid.eq.0.and.verb1) write(*,*) 'Strong Min res:', vmin
      if(nid.eq.0.and.verb1) write(*,*) 'Strong norm res:', tnorm
      if(nid.eq.0.and.verb1) write(*,*) 'Strong norm dudt:', sn1
      if(nid.eq.0.and.verb1) write(*,*) 'Strong norm fux:', sn2
      if(nid.eq.0.and.verb1) write(*,*) 'Strong norm dpdx:', sn3
      if(nid.eq.0.and.verb1) write(*,*) 'Strong norm convu:', sn4
      if(nid.eq.0.and.verb1) write(*,*) 'Strong norm lapu:', sn5

      return
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

      subroutine sBC_term(fux,fuy,fuz)
      include 'SIZE'
      include 'SOLN' 
      include 'INPUT'
      real*8 fux(lx1,ly1,lz1,lelv), fuy(lx1,ly1,lz1,lelv),
     &         fuz(lx1,ly1,lz1,lelv)
      real*8 v1(lx1,ly1,lz1,lelv), v2(lx1,ly1,lz1,lelv),
     &         v3(lx1,ly1,lz1,lelv)
      real*8 w1(lx1,ly1,lz1,lelv), w2(lx1,ly1,lz1,lelv),
     &         w3(lx1,ly1,lz1,lelv)
      real*8 h1(lx1,ly1,lz1,lelv), h2(lx1,ly1,lz1,lelv)
      real*8 h2inv(lx1,ly1,lz1,lelv)
      ntot1 = lx1*ly1*lz1*nelv

      call rzero(v1,ntot1)
      call rzero(v2,ntot1)
      if(if3d) call rzero(v3,ntot1)
      !intype=0 ! I'm not sure if one has to include time derivatives
      !call sethlm(h1,h2,intype)
      call bcdirvc(v1,v2,v3,v1mask,v2mask,v3mask)
      call opcopy(fux,fuy,fuz,v1,v2,v3)
      !call ophx (w1,w2,w3,v1,v2,v3,h1,h2)
      !call rone(h2inv,ntot1)
      !scl = 1
      !call opbinv (fux,fuy,fuz,w1,w2,w3,h2inv) !Inverted mass matrix
      !call opbinv1 (fux,fuy,fuz,w1,w2,w3,scl) !Inverted mass matrix
      
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

c-----------------------------------------------------------------------
c-----------------------------------------------------------------------

      subroutine make_cont(u)
      include 'SIZE'
      include 'MASS'
      real*8 u(lx1,ly1,lz1,lelv)

      ntot = nx1*ny1*nz1*lelv
      call col2(u,bm1,ntot)
      call dssum(u,ntot)
      call col2(u,binvm1,ntot)
      end subroutine

c-----------------------------------------------------------------------
c-----------------------------------------------------------------------


      subroutine my_norm(L1,L2,u)
      include 'SIZE'
      include 'MASS'
      real*8 u(lx1,ly1,lz1,lelv)
      real*8 ubm(lx1,ly1,lz1,lelv)
      real*8 ubm2(lx1,ly1,lz1,lelv)
      real*8 L1,L2

      ntot = nx1*ny1*nz1*lelv
      call col3(ubm,u,bm1,ntot)
      L1 = glsum(ubm,ntot)/VOLVM1
      call col3(ubm2,u,u,ntot)
      call col2(ubm2,bm1,ntot)
      L2 = glsum(ubm2,ntot)/VOLVM1
      end subroutine




