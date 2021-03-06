PURE FUNCTION n_combination(n,k) RESULT(r)
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
END FUNCTION n_combination

MODULE detMod  
        INTEGER(4), ALLOCATABLE :: l_det_alpha(:,:), l_det_beta(:,:)
        INTEGER(8) :: n_det      
END MODULE detMod


SUBROUTINE gen_permutation(n_orbital, n_int, log_size_orbital_bucket, occ, n_alpha)
      use detMod, ONLY : l_det_alpha, l_det_beta, n_det
      implicit none
      INTEGER(8) :: n_combination

      INTEGER(4), INTENT(IN) :: n_orbital, n_alpha, n_int,log_size_orbital_bucket
      INTEGER(1), DIMENSION(n_orbital), INTENT(IN) :: occ

      INTEGER(4) :: n_alpha_f, n_single_orbital
      INTEGER(4), DIMENSION(n_int) :: det_pattern

      INTEGER(4) :: idx, idx_d, i, det_i, det_p, p,t,tt
      LOGICAL :: o

      INTEGER(4), ALLOCATABLE :: occ_if(:)
      INTEGER(4) :: size_orbital_bucket

      size_orbital_bucket = ibset(0,log_size_orbital_bucket)
      det_pattern(:) = 0
      n_alpha_f = n_alpha
      n_single_orbital = 0
 
      ALLOCATE(occ_if(n_orbital))

      DO IDX = 0, n_orbital-1
        SELECT CASE (occ(idx+1))
            CASE (1)
                n_single_orbital = n_single_orbital + 1
                occ_if(n_single_orbital) = idx
            CASE (2)
                det_i = rshift(idx, log_size_orbital_bucket) + 1
                det_p = and(idx,size_orbital_bucket-1)

                det_pattern(det_i) = ibset(det_pattern(det_i),det_p)
                n_alpha_f = n_alpha_f - 1;
            END SELECT
      END DO

      IF (n_single_orbital > BIT_SIZE(p) * 8) THEN
            print*, 'To many permutation'
            STOP 1
      ENDIF

      n_det = n_combination(n_single_orbital, n_alpha_f);

      IF (.NOT.ALLOCATED(l_det_alpha)) THEN
        ALLOCATE(l_det_alpha(n_int,n_det))
        ALLOCATE(l_det_beta(n_int,n_det))
      ENDIF

      p = lshift(1, n_alpha_f) - 1
        
      DO idx_d = 1, n_det

        l_det_alpha(:,idx_d) = det_pattern(:)
        l_det_beta(:,idx_d) = det_pattern(:)
    
        DO i=1,n_single_orbital
            idx = occ_if(i) 
            det_i = rshift(idx, log_size_orbital_bucket) + 1
            det_p = and(idx,size_orbital_bucket-1)

            IF (btest(p,i-1)) THEN
                l_det_alpha(det_i,idx_d) = ibset(l_det_alpha(det_i,idx_d), det_p)
            ELSE
                l_det_beta(det_i,idx_d) = ibset(l_det_beta(det_i,idx_d), det_p)
            ENDIF
        ENDDO

        ! Compute the permutation
        ! https://graphics.stanford.edu/~seander/bithacks.html#NextBitPermutation
        t = or(p, p-1)
        tt = t+1
        p = or(tt, (rshift(and(not(t),tt) - 1, trailz(p) +1)))

       ENDDO

       DEALLOCATE(occ_if)


END SUBROUTINE


program hello
    use detMod, ONLY : l_det_alpha, l_det_beta, n_det

    implicit none

    INTEGER(1), ALLOCATABLE :: occ(:)
    INTEGER(4) :: n_orbital, n_alpha, n_int
    INTEGER(4), PARAMETER :: size_orbital_bucket = 32
    INTEGER(4), PARAMETER :: log_size_orbital_bucket = trailz(size_orbital_bucket)
   
    CHARACTER(LEN=64) :: buffer
    
    INTEGER :: i, d, d_init, mode

    REAL(8) :: t1, t0

    call get_command_argument(1, buffer)
    read(buffer, '(I4)') mode


    call get_command_argument(2, buffer)
    read(buffer, '(I4)') n_orbital

    call get_command_argument(3, buffer)
    read(buffer, '(I4)') n_alpha

    n_int =  n_orbital  / size_orbital_bucket  + 1 ;
    print*, "n_int:",  n_int

    ALLOCATE(occ(n_orbital))
    DO i= 1, n_orbital
        call get_command_argument(i+3, buffer)
        read(buffer, '(I1)') occ(i)
    ENDDO

    do i=1,100
      call gen_permutation(n_orbital, n_int, log_size_orbital_bucket, occ, n_alpha)
    enddo
    print*, 'n_det:', n_det

    IF (mode.EQ.0) THEN
        d_init = n_det+1
    ELSE
        d_init = 1
    ENDIF
 
    DO D=d_init,n_det
        DO I=1, n_int
             write(*, '(B32.32)') l_det_alpha(I,D)
        ENDDO

        DO I=1, n_int
             write(*, '(B32.32)') l_det_beta(I,D)
        ENDDO

    ENDDO 
end program hello
