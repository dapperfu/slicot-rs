//! NF01AD — Output of Wiener system (SLICOT NF01AD)
//
// x(t+1)=A*x(t)+B*u(t), z(t)=C*x(t)+D*u(t), y(t)=f(z(t),wb). f is NF01AY.

/// X = [wb(1:L); theta]. theta = N*(M+L+1)+L*M (A,B,C,D,x0 in output normal form).
/// For N=0: theta = D (L*M). For N>0 we use direct pack A,B,C,D,x0 for simulation.
pub fn nf01ad(
    nsmp: i32,
    m: i32,
    l: i32,
    ipar: &[i32],
    lipar: i32,
    x: &[f64],
    lx: i32,
    u: &[f64],
    ldu: usize,
    y: &mut [f64],
    ldy: usize,
    dwork: &mut [f64],
    ldwork: i32,
) -> i32 {
    let (nsmp, m, l) = (nsmp as usize, m as usize, l as usize);
    let n = *ipar.get(0).unwrap_or(&0) as usize;
    let nn = *ipar.get(1).unwrap_or(&0) as usize;
    if lipar < 2 {
        return -5;
    }
    let wb_len = (nn * (l + 2) + 1) * l;
    let theta_len = n * (m + l + 1) + l * m;
    if lx < (wb_len + theta_len) as i32 {
        return -7;
    }
    if ldu < nsmp.max(1) || ldy < nsmp.max(1) {
        return -10;
    }
    let wb = &x[0..wb_len];
    let theta = &x[wb_len..wb_len + theta_len];

    if n == 0 {
        for t in 0..nsmp {
            let (z_part, rest) = dwork.split_at_mut(l);
            for j in 0..l {
                let mut zj = 0.0;
                for k in 0..m {
                    zj += theta[j + k * l] * u[t + k * ldu];
                }
                z_part[j] = zj;
            }
            let mut y_t = vec![0.0; l];
            let info = super::nf01ay::nf01ay(
                1,
                l as i32,
                l as i32,
                ipar,
                lipar,
                wb,
                wb_len as i32,
                z_part,
                1,
                &mut y_t,
                1,
                rest,
                ldwork - l as i32,
            );
            if info != 0 {
                return info;
            }
            for j in 0..l {
                y[t + j * ldy] = y_t[j];
            }
        }
        return 0;
    }

    let (n_n, n_m, l_n, l_m) = (n * n, n * m, l * n, l * m);
    if theta.len() < n_n + n_m + l_n + l_m + n {
        return -7;
    }
    let mut state = vec![0.0; n];
    for i in 0..n {
        state[i] = theta[n_n + n_m + l_n + l_m + i];
    }
    let a = &theta[0..n_n];
    let b = &theta[n_n..n_n + n_m];
    let c = &theta[n_n + n_m..n_n + n_m + l_n];
    let d = &theta[n_n + n_m + l_n..n_n + n_m + l_n + l_m];

    let lda = n;
    let ldb = n;
    let ldc = l;
    let ldd = l;

    for t in 0..nsmp {
        let (z_part, rest) = dwork.split_at_mut(l);
        for j in 0..l {
            let mut zj = 0.0;
            for i in 0..n {
                zj += c[j + i * ldc] * state[i];
            }
            for k in 0..m {
                zj += d[j + k * ldd] * u[t + k * ldu];
            }
            z_part[j] = zj;
        }
        let mut y_t = vec![0.0; l];
        let info = super::nf01ay::nf01ay(
            1,
            l as i32,
            l as i32,
            ipar,
            lipar,
            wb,
            wb_len as i32,
            z_part,
            1,
            &mut y_t,
            1,
            rest,
            ldwork - l as i32,
        );
        if info != 0 {
            return info;
        }
        for j in 0..l {
            y[t + j * ldy] = y_t[j];
        }
        let mut new_state = vec![0.0; n];
        for i in 0..n {
            let mut sum = 0.0;
            for ii in 0..n {
                sum += a[i + ii * lda] * state[ii];
            }
            for k in 0..m {
                sum += b[i + k * ldb] * u[t + k * ldu];
            }
            new_state[i] = sum;
        }
        state = new_state;
    }
    0
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_nf01ad_n0() {
        let ipar = [0_i32, 0_i32];
        let x = [0.0; 4];
        let u = [1.0];
        let mut y = [0.0; 2];
        let mut dwork = [0.0; 20];
        assert_eq!(
            nf01ad(1, 1, 1, &ipar, 2, &x, 4, &u, 1, &mut y, 1, &mut dwork, 20),
            0
        );
    }
}
