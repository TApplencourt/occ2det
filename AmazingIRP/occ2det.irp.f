program occ2det
  implicit none
  print *, 'n_orbital:',  n_orbital
  print *, 'n_alpha:',  n_alpha
  print *, 'n_int:',  n_int
  print *, 'n_det:',  n_det
  
  integer :: d_init
  if (mode == 0) then
    d_init = n_det+1
  else
    d_init = 1
  endif

  integer :: d,i,j
  do d = d_init, n_det
    do j=1,2
      do i=1,n_int
        write(*, '(B32.32)') gen_dets(i,j,d)
      enddo
    enddo
  enddo

end
