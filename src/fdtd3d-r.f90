!----------------------------------------------------------------------
!
!  module: fdtd3d(-r) / max3d
!
!  fdtd core algorithm.
!
!----------------------------------------------------------------------


module fdtd3d
 
  use grid
  use constant
  implicit none
  save
  
contains

  subroutine StepH

    implicit none
    
    real(8) :: dtx, dty, dtz
    real(8) :: Exh,Eyh,Ezh
    integer :: i, j, k
    
    dtx = DT/Sx
    dty = DT/Sy
    dtz = DT/Sz


    ! H in GT+1/2DT

    do k=KSIG, KEIG-1
       do j=JSIG, JEIG-1
          do i=ISIG, IEIG-1

             Exh=Ex(i,j,k)
             Eyh=Ey(i,j,k)
             Ezh=Ez(i,j,k)

             Hx(i,j,k) =  Hx(i,j,k) &
                  - dty*( Ez(i,j+1,k) - Ezh ) &
                  + dtz*( Ey(i,j,k+1) - Eyh )
             Hy(i,j,k) = Hy(i,j,k) &
                  - dtz*( Ex(i,j,k+1) - Exh ) &
                  + dtx*( Ez(i+1,j,k) - Ezh )
             Hz(i,j,k) = Hz(i,j,k) &
                  - dtx*( Ey(i+1,j,k) - Eyh ) &
                  + dty*( Ex(i,j+1,k) - Exh )

          enddo
       enddo
    enddo
        
  end subroutine StepH


  subroutine StepE

    implicit none
    
    real(8) :: dtx, dty, dtz
    real(8) :: Hxh, Hyh, Hzh
    real(8) :: epsinvx, epsinvy, epsinvz
    integer :: i, j, k
    
    dtx = DT/Sx
    dty = DT/Sy
    dtz = DT/Sz
    
    ! E in GT + DT
    do k=KSIG, KEIG-1
       do j=JSIG, JEIG-1
          do i=ISIG, IEIG-1

             Hxh=Hx(i,j,k)
             Hyh=Hy(i,j,k)
             Hzh=Hz(i,j,k) 

             epsinvx = 0.5*(EPSINV(i,j,k) +  EPSINV(i+1,j,k))
             epsinvy = 0.5*(EPSINV(i,j,k) +  EPSINV(i,j+1,k))
             epsinvz = 0.5*(EPSINV(i,j,k) +  EPSINV(i,j,k+1))

             Ex(i,j,k) =  Ex(i,j,k) +  epsinvx* &
                  ( dty*( Hzh - Hz(i,j-1,k) ) &
                  - dtz*( Hyh - Hy(i,j,k-1) ))
             Ey(i,j,k) =  Ey(i,j,k) +  epsinvy* &
                  ( dtz*( Hxh - Hx(i,j,k-1) ) &
                  - dtx*( Hzh - Hz(i-1,j,k) ))
             Ez(i,j,k) =  Ez(i,j,k) +  epsinvz* &
                  ( dtx*( Hyh - Hy(i-1,j,k) )  &
                  - dty*( Hxh - Hx(i,j-1,k) ))

          enddo
       enddo
    enddo
    
  end subroutine StepE

end module fdtd3d
