! Fortran benchmark driver for SLICOT routines (MA02ED, MA02ES, DLACPY_SLC).
! Times each routine at n = 32, 64, 128, 256, 512, 1024 to match Rust ladder.
! Compile: see Makefile. Run: ./bench_slicot
! Output: one line per (routine, n) with time in microseconds.

program bench_slicot
  implicit none
  integer, parameter :: nsizes = 6
  integer, parameter :: sizes(nsizes) = [32, 64, 128, 256, 512, 1024]
  integer :: n, i, j, iter, niter, lda
  double precision, allocatable :: A(:,:), B(:,:)
  double precision :: t0, t1, t_per_call_us
  integer :: count_rate, count_max, c0, c1
  character(len=1) :: uplo

  call system_clock(count_rate=count_rate, count_max=count_max)

  ! MA02ED: store by symmetry (upper triangle given)
  do i = 1, nsizes
    n = sizes(i)
    if (n > 1024) cycle
    lda = n
    allocate(A(lda,n))
    do j = 1, n
      do iter = 1, n
        A(iter,j) = dble(iter + j) * 0.1d0
      end do
    end do
    uplo = 'U'
    niter = 0
    call system_clock(c0)
    do
      call MA02ED(uplo, n, A, lda)
      niter = niter + 1
      call system_clock(c1)
      if (niter >= 100 .and. (c1 - c0) * 1.0d0 / count_rate >= 1.0d0) exit
      if (niter >= 10000000) exit
    end do
    t_per_call_us = (c1 - c0) * 1.0d6 / (count_rate * niter)
    print '(A,I5,A,F12.3,A)', 'MA02ED  n=', n, '  ', t_per_call_us, ' us/call'
    deallocate(A)
  end do

  ! MA02ES: store by skew-symmetry
  do i = 1, nsizes
    n = sizes(i)
    if (n > 1024) cycle
    lda = n
    allocate(A(lda,n))
    do j = 1, n
      do iter = 1, n
        A(iter,j) = dble(iter + j) * 0.1d0
      end do
    end do
    uplo = 'U'
    niter = 0
    call system_clock(c0)
    do
      call MA02ES(uplo, n, A, lda)
      niter = niter + 1
      call system_clock(c1)
      if (niter >= 100 .and. (c1 - c0) * 1.0d0 / count_rate >= 1.0d0) exit
      if (niter >= 10000000) exit
    end do
    t_per_call_us = (c1 - c0) * 1.0d6 / (count_rate * niter)
    print '(A,I5,A,F12.3,A)', 'MA02ES  n=', n, '  ', t_per_call_us, ' us/call'
    deallocate(A)
  end do

  ! DLACPY_SLC: copy full matrix (UPLO = 'A' or any non-U/L for full)
  do i = 1, nsizes
    n = sizes(i)
    if (n > 1024) cycle
    lda = n
    allocate(A(lda,n), B(lda,n))
    do j = 1, n
      do iter = 1, n
        A(iter,j) = dble(iter + j) * 0.1d0
      end do
    end do
    uplo = 'F'
    niter = 0
    call system_clock(c0)
    do
      call DLACPY_SLC(uplo, n, n, A, lda, B, lda)
      niter = niter + 1
      call system_clock(c1)
      if (niter >= 100 .and. (c1 - c0) * 1.0d0 / count_rate >= 1.0d0) exit
      if (niter >= 10000000) exit
    end do
    t_per_call_us = (c1 - c0) * 1.0d6 / (count_rate * niter)
    print '(A,I5,A,F12.3,A)', 'DLACPY  n=', n, '  ', t_per_call_us, ' us/call'
    deallocate(A, B)
  end do

end program bench_slicot
