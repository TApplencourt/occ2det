program occ2det
  implicit none
  print *, 'n_orbital:',  n_orbital
  print *, 'n_int:',  n_int
  

  integer :: i
  print *,  single_index(:)
!  do i=1,n_int
!    write(*,'(B32.32)') occ(i,1)
!    write(*,'(B32.32)') occ(i,2)
!    write(*,*)
!  enddo
  print *, 'n_det:',  n_det
end
