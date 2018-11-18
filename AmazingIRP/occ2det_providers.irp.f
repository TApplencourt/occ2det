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
 integer                        :: ipos, iint
 
 integer(bit_kind)              :: v,t,tt
 integer                        :: idx
 integer                        :: ispin
 integer                        :: single_index_local(n_single_orbital)

 single_index_local(:) = single_index(:)-1

 v = shiftl(1,n_alpha_in_single) - 1

 do i=1,n_det

   ! Initialize gen_dets with occ(:,2)
   gen_dets(:,1,i) = occ(:,2)
   gen_dets(:,2,i) = occ(:,2)

   do k=1,n_single_orbital
      idx = single_index_local(k)

      ispin = 2 - iand(shiftr(v,k-1),1) 

      ! Find integer
      iint = 1 + shiftr(idx,log_size_orbital_bucket) 

      ! Find position in integer
      ipos = idx-shiftl((iint-1),log_size_orbital_bucket)

      ! Set the bit in the temporary determinant
      gen_dets(iint,ispin,i) = ibset( gen_dets(iint,ispin,i), ipos )
   enddo

   ! Generate next permutation
   t = ior(v,v-1)
   tt = t+1
   v = ior(tt, shiftr( and(not(t),tt) - 1, trailz(v)+1) )
 enddo

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


