//! Minimal BLAS-style helpers for MB04 (column-major layout).
//! Element (i,j) is at index i + j * ld.

/// DAXPY: y := alpha*x + y. incx, incy are strides.
#[inline]
pub(crate) fn daxpy(n: usize, alpha: f64, x: &[f64], incx: usize, y: &mut [f64], incy: usize) {
    if alpha == 0.0 || n == 0 {
        return;
    }
    if incx == 1 && incy == 1 {
        for i in 0..n {
            y[i] += alpha * x[i];
        }
    } else {
        for i in 0..n {
            y[i * incy] += alpha * x[i * incx];
        }
    }
}

/// DCOPY: y := x.
#[inline]
pub(crate) fn dcopy(n: usize, x: &[f64], incx: usize, y: &mut [f64], incy: usize) {
    if n == 0 {
        return;
    }
    if incx == 1 && incy == 1 {
        y[..n].copy_from_slice(&x[..n]);
    } else {
        for i in 0..n {
            y[i * incy] = x[i * incx];
        }
    }
}

/// DGEMV: y := alpha*A*x + beta*y (trans='N') or y := alpha*A'*x + beta*y (trans='T').
/// A is m-by-n, ld is leading dimension of A.
#[inline]
pub(crate) fn dgemv(
    trans: bool,
    m: usize,
    n: usize,
    alpha: f64,
    a: &[f64],
    lda: usize,
    x: &[f64],
    incx: usize,
    beta: f64,
    y: &mut [f64],
    incy: usize,
) {
    if m == 0 || n == 0 {
        if beta != 1.0 {
            for i in 0..(if trans { n } else { m }) {
                y[i * incy] *= beta;
            }
        }
        return;
    }
    if trans {
        // y := alpha*A'*x + beta*y  => y_j += alpha * sum_i A(i,j)*x_i
        if beta != 1.0 {
            for j in 0..n {
                y[j * incy] *= beta;
            }
        }
        for j in 0..n {
            let mut t = 0.0_f64;
            for i in 0..m {
                t += a[i + j * lda] * x[i * incx];
            }
            y[j * incy] += alpha * t;
        }
    } else {
        // y := alpha*A*x + beta*y
        if beta != 1.0 {
            for i in 0..m {
                y[i * incy] *= beta;
            }
        }
        for i in 0..m {
            let mut t = 0.0_f64;
            for j in 0..n {
                t += a[i + j * lda] * x[j * incx];
            }
            y[i * incy] += alpha * t;
        }
    }
}

/// DGER: A := alpha*x*y' + A (rank-1 update). A is m-by-n, lda.
#[inline]
pub(crate) fn dger(
    m: usize,
    n: usize,
    alpha: f64,
    x: &[f64],
    incx: usize,
    y: &[f64],
    incy: usize,
    a: &mut [f64],
    lda: usize,
) {
    if m == 0 || n == 0 || alpha == 0.0 {
        return;
    }
    for j in 0..n {
        let yj = alpha * y[j * incy];
        for i in 0..m {
            a[i + j * lda] += x[i * incx] * yj;
        }
    }
}

/// DDOT: dot product, x'*y.
#[inline]
pub(crate) fn ddot(n: usize, x: &[f64], incx: usize, y: &[f64], incy: usize) -> f64 {
    if n == 0 {
        return 0.0;
    }
    let mut sum = 0.0_f64;
    for i in 0..n {
        sum += x[i * incx] * y[i * incy];
    }
    sum
}

/// DSCAL: x := alpha*x.
#[inline]
pub(crate) fn dscal(n: usize, alpha: f64, x: &mut [f64], incx: usize) {
    if n == 0 || alpha == 1.0 {
        return;
    }
    for i in 0..n {
        x[i * incx] *= alpha;
    }
}

/// DSWAP: swap x and y.
#[inline]
pub(crate) fn dswap(n: usize, x: &mut [f64], incx: usize, y: &mut [f64], incy: usize) {
    if n == 0 {
        return;
    }
    for i in 0..n {
        let xi = x[i * incx];
        let yi = y[i * incy];
        x[i * incx] = yi;
        y[i * incy] = xi;
    }
}

/// DTRMV: x := A*x (trans=false) or A'*x (trans=true). A upper triangular, non-unit, n×n.
#[inline]
pub(crate) fn dtrmv(
    uplo_upper: bool,
    trans: bool,
    n: usize,
    a: &[f64],
    lda: usize,
    x: &mut [f64],
    incx: usize,
) {
    if n == 0 {
        return;
    }
    if trans {
        // x := A'*x, A upper => x_i = sum_{j>=i} A(j,i)*x_j
        for i in (0..n).rev() {
            let mut t = 0.0_f64;
            for j in i..n {
                t += a[j + i * lda] * x[j * incx];
            }
            x[i * incx] = t;
        }
    } else {
        // x := A*x, A upper => x_i = sum_{j<=i} A(i,j)*x_j
        for i in 0..n {
            let mut t = 0.0_f64;
            for j in 0..=i {
                t += a[i + j * lda] * x[j * incx];
            }
            x[i * incx] = t;
        }
    }
}

/// DTRMM side='R': B(m,n) := alpha*B*op(A), A n×n upper.
#[inline]
pub(crate) fn dtrmm_right(
    trans_a: bool,
    m: usize,
    n: usize,
    alpha: f64,
    a: &[f64],
    lda: usize,
    b: &mut [f64],
    ldb: usize,
) {
    if m == 0 || n == 0 || alpha == 0.0 {
        return;
    }
    for j in (0..n).rev() {
        if alpha != 1.0 {
            for i in 0..m {
                b[i + j * ldb] *= alpha;
            }
        }
        for i in (0..j).rev() {
            let aij = if trans_a { a[j + i * lda] } else { a[i + j * lda] };
            if aij != 0.0 {
                for k in 0..m {
                    b[k + j * ldb] += aij * b[k + i * ldb];
                }
            }
        }
    }
}

/// DTRMM side='L': B(m,n) := alpha*op(A)*B, A m×m upper.
#[inline]
pub(crate) fn dtrmm_left(
    trans_a: bool,
    m: usize,
    n: usize,
    alpha: f64,
    a: &[f64],
    lda: usize,
    b: &mut [f64],
    ldb: usize,
) {
    if m == 0 || n == 0 || alpha == 0.0 {
        return;
    }
    for j in 0..n {
        if trans_a {
            // B := alpha*A'*B => B(i,j) = alpha * sum_{k>=i} A(k,i)*B(k,j)
            for i in 0..m {
                let mut t = 0.0_f64;
                for k in i..m {
                    t += a[k + i * lda] * b[k + j * ldb];
                }
                b[i + j * ldb] = alpha * t;
            }
        } else {
            // B := alpha*A*B => B(i,j) = alpha * sum_{k<=i} A(i,k)*B(k,j)
            for i in 0..m {
                let mut t = 0.0_f64;
                for k in 0..=i {
                    t += a[i + k * lda] * b[k + j * ldb];
                }
                b[i + j * ldb] = alpha * t;
            }
        }
    }
}

/// DLARFG: Generate real Householder reflector H = I - tau*u*u', u = (1; v).
/// On entry: alpha is scalar, x[0..n-1] is the rest of the vector.
/// On exit: alpha is overwritten by beta (so H*[alpha;x] = [beta;0]), x is overwritten by v, tau is set.
/// INCX is stride for x.
pub(crate) fn dlarfg(n: usize, alpha: &mut f64, x: &mut [f64], incx: usize, tau: &mut f64) {
    if n <= 1 {
        *tau = 0.0;
        return;
    }
    let mut xnorm = 0.0_f64;
    for i in 0..(n - 1) {
        let xi = x[i * incx];
        xnorm += xi * xi;
    }
    xnorm = xnorm.sqrt();
    if xnorm == 0.0 {
        *tau = 0.0;
        return;
    }
    let sign_alpha = if *alpha >= 0.0 { 1.0_f64 } else { -1.0_f64 };
    let beta = -sign_alpha * (alpha.mul_add(*alpha, xnorm * xnorm).sqrt());
    let alpha_minus_beta = *alpha - beta;
    *tau = (beta - *alpha) / beta;
    *alpha = beta;
    for i in 0..(n - 1) {
        x[i * incx] /= alpha_minus_beta;
    }
}
