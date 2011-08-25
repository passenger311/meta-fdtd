!-*- F90 -*------------------------------------------------------------
!
!  module: diagebalmode_outgpl / meta
!
!  this module handles GPL output of data related to the diagebalmode
!  module.
!
!----------------------------------------------------------------------


module diagebalmode_outgpl

  use constant
  use strings
  use reglist
  use outlist
  use buflist
  use mpiworld
  use grid 
  use fdtd
  use fdtd_calc
  use diagebalmode
  use out_calc
  use matlorentz

  implicit none
  private
  save

  ! --- Module Identifier

  character(len=STRLNG), parameter :: modname = 'DIAGEBALMODE_OUTGPL'

 ! --- Public Methods

  public :: InitializeDiagebalModeOutgplObj
  public :: FinalizeDiagebalModeOutgplObj
  public :: WriteDataDiagebalModeOutgplObj

contains

!----------------------------------------------------------------------

  subroutine InitializeDiagebalModeOutgplObj(out)

    type (T_OUT) :: out
    
    out%numnodes = 1    

  end subroutine InitializeDiagebalModeOutgplObj


!----------------------------------------------------------------------

  subroutine FinalizeDiagebalModeOutgplObj(out)

    type (T_OUT) :: out

  end subroutine FinalizeDiagebalModeOutgplObj


!----------------------------------------------------------------------


  subroutine WriteDataDiagebalModeOutgplObj(out, mode)

    type (T_OUT) :: out
    logical :: mode

    if ( .not. mode ) return

    select case (out%fn)
    case('En') ! differential
       call WriteValuesDiagebalMode(out, 1 )
    case('EnI') ! time integrated
       call WriteValuesDiagebalMode(out, 2 )
    case('DS') ! differential
       call WriteValuesDiagebalMode(out, 3 )
    case('DSI') ! time integrated
       call WriteValuesDiagebalMode(out, 4 )
    case('Um')
       call WriteMatrixDiagebalMode(out, 1 )   
    case('Sm')
       call WriteMatrixDiagebalMode(out, 2 )   
    case('Jm')
       call WriteMatrixDiagebalMode(out, 3 )   
    case default
       write(out%funit,*) "OUTPUT FUNCTION NOT IMPLEMENTED" 
    end select
    
  end subroutine WriteDataDiagebalModeOutgplObj

  subroutine WriteMatrixDiagebalMode(out, mode)

    type (T_OUT) :: out
    integer :: mode
    type (T_DIAGEBALMODE) :: diag
    type (T_REG) :: reg
    integer :: i,j,k, m_P, n_P, m4_m
    real(kind=8) :: v
    integer :: err

    diag = diagebalmodeobj(out%objidx)
    
    M4_WRITE_DBG({"WriteValues!"})
    M4_IFELSE_DBG({call EchoRegObj(regobj(out%regidx))})

    M4_MODOBJ_GETREG(out,reg)
    k = ((reg%ke-reg%ks)/2)+reg%ks
    
    do j=reg%js, reg%je, reg%dj
       do i=reg%is, reg%ie, reg%di
          select case ( mode )
             case (1)
!                v = abs(diag%buf_E(diag%h_pos_E,i,j,k,0))
                v = M4_VOLEX(i,j,k) / epsinvx(i,j,k) * &
                     (REALPART(diag%buf_E(diag%h_pos_E,i,j,k,0))**2+&
                     IMAGPART(diag%buf_E(diag%h_pos_E,i,j,k,0))**2) + &
                     M4_VOLEY(i,j,k) / epsinvy(i,j,k) * &
                     (REALPART(diag%buf_E(diag%h_pos_E,i,j,k,1))**2+&
                     IMAGPART(diag%buf_E(diag%h_pos_E,i,j,k,1))**2) + &
                     M4_VOLEZ(i,j,k) / epsinvz(i,j,k) * &
                     (REALPART(diag%buf_E(diag%h_pos_E,i,j,k,2))**2+&
                     IMAGPART(diag%buf_E(diag%h_pos_E,i,j,k,2))**2) + &
                     M4_VOLHX(i,j,k) / M4_MUINVX(i,j,k) * &
                     (REALPART(diag%buf_H(diag%h_pos_H,i,j,k,0))**2+&
                     IMAGPART(diag%buf_H(diag%h_pos_H,i,j,k,0))**2) + &
                     M4_VOLHY(i,j,k) / M4_MUINVY(i,j,k) * &
                     (REALPART(diag%buf_H(diag%h_pos_H,i,j,k,1))**2+&
                     IMAGPART(diag%buf_H(diag%h_pos_H,i,j,k,1))**2) + &
                     M4_VOLHZ(i,j,k) / M4_MUINVZ(i,j,k) * &
                     (REALPART(diag%buf_H(diag%h_pos_H,i,j,k,2))**2+&
                     IMAGPART(diag%buf_H(diag%h_pos_H,i,j,k,2))**2);
             case (2)
                v = diag%sumds
             case (3)                
                v = 0
                m_P = diag%h_pos_P
                if (m_P .eq. 0) then
                   n_P = diag%p - 1
                else
                   n_P = m_P - 1
                endif
                do m4_m = 1, numMATLORENTZobj
                   v = v + (M4_IFELSE_TM({ M4_VOLEX(i,j,k) * abs(diag%buf_E(diag%h_pos_E,i,j,k,0) *&
                        dconjg( diag%buf_P(m_P,i,j,k,0,m4_m-1) - diag%buf_P(n_P,i,j,k,0,m4_m-1) ) ) / DT +},{0. +}) &
                        M4_IFELSE_TM({ M4_VOLEY(i,j,k) * abs(diag%buf_E(diag%h_pos_E,i,j,k,1) *&
                        dconjg( diag%buf_P(m_P,i,j,k,1,m4_m-1) - diag%buf_P(n_P,i,j,k,1,m4_m-1) ) ) / DT +},{0. +}) &
                        M4_IFELSE_TE({ M4_VOLEZ(i,j,k) * abs(diag%buf_E(diag%h_pos_E,i,j,k,2) *&
                        dconjg( diag%buf_P(m_P,i,j,k,2,m4_m-1) - diag%buf_P(n_P,i,j,k,2,m4_m-1) ) ) / DT},{0. }))
                enddo
          end select
          write(out%funit,"(5E15.6E3)",advance = "no") v
       enddo
       write(out%funit,"(5E15.6E3)")
    enddo
  end subroutine WriteMatrixDiagebalMode

  subroutine WriteValuesDiagebalMode(out, mode)
    
    type (T_OUT) :: out
    integer :: mode
    type (T_DIAGEBALMODE) :: diag
    
    diag = diagebalmodeobj(out%objidx)
    
    M4_WRITE_DBG({"WriteValues!"})
    M4_IFELSE_DBG({call EchoRegObj(regobj(out%regidx))})
    
    select case ( mode ) 
       
    case ( 1 )
       write(out%funit,"(5E15.6E3)") diag%dudt, diag%ds, diag%je, diag%res
    case ( 2 )
       write(out%funit,"(5E15.6E3)") diag%sumdudt, diag%sumds,  diag%sumje, diag%sumres
    case ( 3 ) 
       write(out%funit,"(4E15.6E3)") diag%ds, diag%dsx, diag%dsy, diag%dsz
    case ( 4 ) 
       write(out%funit,"(4E15.6E3)") diag%sumds, diag%sumdsx, diag%sumdsy, diag%sumdsz
    end select

  end subroutine WriteValuesDiagebalMode
        
end module diagebalmode_outgpl

!
! Authors:  J.Hamm
! Modified: 1/02/2008
!
!======================================================================
