integer, parameter :: bit_kind = 4


BEGIN_PROVIDER [ integer, mode ]
 implicit none
 BEGIN_DOC
! if mode = 0, don't print else print in stdout.
 END_DOC
 character(len=64)              :: buffer
 call get_command_argument(1, buffer)
 read(buffer, '(I4)') mode
END_PROVIDER


BEGIN_PROVIDER [ integer, n_orbital ]
 implicit none
 BEGIN_DOC
! Number of singly occupied mos
 END_DOC
 character(len=64)              :: buffer
 call get_command_argument(2, buffer)
 read(buffer, '(I4)') n_orbital
END_PROVIDER


BEGIN_PROVIDER [ integer, n_alpha ]
 implicit none
 BEGIN_DOC
! Number of alpha electrons 
 END_DOC
 character(len=64)              :: buffer
 call get_command_argument(3, buffer)
 read(buffer, '(I4)') n_alpha
END_PROVIDER


BEGIN_PROVIDER [ integer, size_orbital_bucket ]
 implicit none
 BEGIN_DOC
 ! Number of bits in an integer
 END_DOC
 size_orbital_bucket = 8*bit_kind
END_PROVIDER


BEGIN_PROVIDER [ integer, log_size_orbital_bucket ]
 implicit none
 BEGIN_DOC
 ! Number of bit to shift to divide/multiply by size_orbital_bucket
 END_DOC                          
 log_size_orbital_bucket = trailz(size_orbital_bucket)
END_PROVIDER                      


BEGIN_PROVIDER [ integer, n_int ]
 implicit none
 BEGIN_DOC
 ! Number of integers required for n_orbital bits
 END_DOC
 n_int = n_orbital / size_orbital_bucket + 1
END_PROVIDER


BEGIN_PROVIDER [ integer, occ_int, (n_orbital) ]
 implicit none
 BEGIN_DOC
 ! Occupation pattern, stored as an array of integers
 END_DOC
 integer                        :: i
 character(len=64)              :: buffer
 do i= 1, n_orbital
   call get_command_argument(i+3, buffer)
   read(buffer, '(I1)') occ_int(i)
 enddo
END_PROVIDER


BEGIN_PROVIDER [ integer(bit_kind), occ, (n_int,2) ]
 implicit none
 BEGIN_DOC
 ! Spatial occupation pattern, encoded in bit strings
 END_DOC
 integer :: i, iint, ipos
 occ = 0
 do i=1,n_orbital
   if (occ_int(i) == 0) cycle
   iint = shiftr(i-1,log_size_orbital_bucket) + 1
   ipos = i-shiftl((iint-1),log_size_orbital_bucket)-1
   if (occ_int(i) == 1) then
     occ(iint,1) = ibset( occ(iint,1), ipos )
   else if (occ_int(i) == 2) then
     occ(iint,2) = ibset( occ(iint,2), ipos )
   else
     stop 'bug in occ'
   endif
 enddo
END_PROVIDER


BEGIN_PROVIDER [ integer, n_single_orbital ]
 implicit none
 BEGIN_DOC
 ! Number of singly occupied orbitals
 END_DOC
 integer :: i
 n_single_orbital = popcnt(occ(1,1))
 do i=2,n_int
   n_single_orbital = n_single_orbital + popcnt(occ(i,1))
 enddo
 if (n_single_orbital > 63) then
   stop 'n_single_orbital too large'
 endif
END_PROVIDER


BEGIN_PROVIDER [ integer, single_index, (n_single_orbital) ]
 implicit none
 BEGIN_DOC
 ! Indexes of the singly occupied MOs
 END_DOC
 integer                        :: i, n, ishift
 integer(bit_kind)              :: v

 ishift = 1
 n = 0
 do i=1,n_int
   v = occ(i,1)
   do while (v /= 0_bit_kind)
     n = n+1
     single_index(n) = ishift + trailz(v)
     v = iand(v,v-1)
   enddo
   ishift = ishift + size_orbital_bucket
 enddo
 
END_PROVIDER


BEGIN_PROVIDER [ integer, n_alpha_in_single ]
 implicit none
 BEGIN_DOC
 ! Number of alpha electrons in open shells
 END_DOC
 n_alpha_in_single = n_alpha
 integer :: i
 do i=1,n_orbital
   if (occ_int(i) == 2) then
     n_alpha_in_single -= 1
   endif
 enddo
END_PROVIDER


BEGIN_PROVIDER [ integer, n_det ]
 implicit none
 BEGIN_DOC
 ! Number of generated determinants
 END_DOC
 integer*8, external :: n_combinations
 n_det = n_combinations(n_single_orbital,n_alpha_in_single)
END_PROVIDER


BEGIN_PROVIDER [ integer(bit_kind), gen_dets, (n_int,2,n_det) ]
 implicit none
 BEGIN_DOC
 ! Generated determinants
 END_DOC
 integer                        :: i, k
 integer                        :: ipos(n_single_orbital), iint(n_single_orbital)
 
 integer(bit_kind)              :: v,t,tt,v_prev,diff
 integer                        :: idx
 integer                        :: ispin, ispin2

 ! Compute once the iint and ipos for each singly occupied MO
 do k=1,n_single_orbital
   idx = single_index(k)-1
   ! Find integer
   iint(k) = 1 + shiftr(idx,log_size_orbital_bucket) 
   ! Find position in integer
   ipos(k) = idx-shiftl((iint(k)-1),log_size_orbital_bucket)
 enddo

 v_prev = 0_bit_kind
 v = shiftl(1,n_alpha_in_single) - 1

 if (shiftl(n_alpha_in_single,1) == n_single_orbital) then

    ! Initialize first determinant
    gen_dets(:,1,1) = occ(:,2)
    gen_dets(:,2,1) = occ(:,2)

    do k=1,n_alpha_in_single
        gen_dets(iint(k),1,1) = ibset( gen_dets(iint(k),1,1), ipos(k) )
    enddo

    do k=n_alpha_in_single+1,n_single_orbital
        gen_dets(iint(k),2,1) = ibset( gen_dets(iint(k),2,1), ipos(k) )
    enddo

    ! Time reversal symmetry
    gen_dets(:,1,2) = gen_dets(:,2,1)
    gen_dets(:,2,2) = gen_dets(:,1,1)

    do i=3,n_det,2

      ! Generate next permutation
      v_prev = v
      t = ior(v,v-1)
      tt = t+1
      v = ior(tt, shiftr( and(not(t),tt) - 1, trailz(v)+1) )

      ! Find what has changed between v_prev and v
      diff = ieor(v,v_prev)

      gen_dets(:,1,i) = gen_dets(:,1,i-2)
      gen_dets(:,2,i) = gen_dets(:,2,i-2)

      ! Swap bits only where they have changed from v_prev to v
      do while (diff /= 0_bit_kind)
        k = trailz(diff)+1
        if (btest(v,k-1)) then 
          gen_dets(iint(k),1,i) = ibset( gen_dets(iint(k),1,i), ipos(k) )
          gen_dets(iint(k),2,i) = ibclr( gen_dets(iint(k),2,i), ipos(k) )
        else
          gen_dets(iint(k),1,i) = ibclr( gen_dets(iint(k),1,i), ipos(k) )
          gen_dets(iint(k),2,i) = ibset( gen_dets(iint(k),2,i), ipos(k) )
        endif
        diff = iand(diff,diff-1)
      enddo

      ! Time reversal symmetry
      gen_dets(:,1,i+1) = gen_dets(:,2,i)
      gen_dets(:,2,i+1) = gen_dets(:,1,i)

    enddo

 else

    ! Initialize first determinant
    gen_dets(:,1,1) = occ(:,2)
    gen_dets(:,2,1) = occ(:,2)

    do k=1,n_alpha_in_single
        gen_dets(iint(k),1,1) = ibset( gen_dets(iint(k),1,1), ipos(k) )
    enddo

    do k=n_alpha_in_single+1,n_single_orbital
        gen_dets(iint(k),2,1) = ibset( gen_dets(iint(k),2,1), ipos(k) )
    enddo

    do i=2,n_det

      ! Generate next permutation
      v_prev = v
      t = ior(v,v-1)
      tt = t+1
      v = ior(tt, shiftr( and(not(t),tt) - 1, trailz(v)+1) )

      ! Find what has changed between v_prev and v
      diff = ieor(v,v_prev)

      gen_dets(:,1,i) = gen_dets(:,1,i-1)
      gen_dets(:,2,i) = gen_dets(:,2,i-1)

      ! Swap bits only where they have changed from v_prev to v
      do while (diff /= 0_bit_kind)
        k = trailz(diff)+1
        if (btest(v,k-1)) then 
          gen_dets(iint(k),1,i) = ibset( gen_dets(iint(k),1,i), ipos(k) )
          gen_dets(iint(k),2,i) = ibclr( gen_dets(iint(k),2,i), ipos(k) )
        else
          gen_dets(iint(k),1,i) = ibclr( gen_dets(iint(k),1,i), ipos(k) )
          gen_dets(iint(k),2,i) = ibset( gen_dets(iint(k),2,i), ipos(k) )
        endif
        diff = iand(diff,diff-1)
      enddo

    enddo

 endif

END_PROVIDER


BEGIN_PROVIDER [ integer(bit_kind), bench_v ]
 implicit none
 BEGIN_DOC
 ! Generated determinants
 END_DOC
 integer                        :: i
 integer(bit_kind)              :: v,t,tt

 v = shiftl(1,n_alpha_in_single) - 1
 do i=1,n_det
   t = ior(v,v-1)
   tt = t+1
   v = ior(tt, shiftr( and(not(t),tt) - 1, trailz(v)+1) )
 enddo
 bench_v = v

END_PROVIDER






FUNCTION n_combinations(n,k) RESULT(r)
    IMPLICIT NONE

    INTEGER(4), INTENT(IN) :: n,k
    INTEGER(8)             :: r
    INTEGER(4)             :: d, n0

    IF (k > n) THEN
        r = 0
        RETURN
    ELSE
        r = 1
    ENDIF

    n0 = n
    DO d=1, k
        r = (r*n0) / d
        n0 = n0 - 1
    ENDDO
END FUNCTION n_combinations


