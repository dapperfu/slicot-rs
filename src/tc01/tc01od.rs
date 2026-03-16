//! TC01OD — Dual of a left/right polynomial matrix representation (SLICOT TC01OD)
//!
//! Given left PMR Q(s)*inv(P(s)) or right inv(P(s))*Q(s), produces the dual right/left PMR
//! by transposing P and Q coefficient arrays (swap roles: P' and Q').

/// Left or right matrix fraction input.
#[derive(Clone, Copy, PartialEq, Eq)]
pub enum Leri {
    /// Left matrix fraction: inv(P)*Q.
    Left,
    /// Right matrix fraction: Q*inv(P).
    Right,
}

/// Computes the dual right (left) PMR of a given left (right) PMR.
/// In-place: PCOEFF becomes P'(s), QCOEFF becomes Q'(s) with dimensions swapped.
///
/// * Left input: P is p×p×indlim, Q is p×m×indlim → dual: P' p×p (transpose each slice), Q' m×p (transpose each slice).
/// * Right input: P is m×m×indlim, Q is m×p×indlim → dual: P' m×m, Q' p×m.
///
/// Storage: pcoeff(i,j,k) = coefficient of s^(indlim-k) for polynomial (i,j). We transpose
/// each k-slice: P' has (i,j) <- P(j,i) for each k; Q' has (i,j) <- Q(j,i), and output Q is m×p (was p×m).
pub fn tc01od(
    leri: Leri,
    m: usize,
    p: usize,
    indlim: usize,
    pcoeff: &mut [f64],
    ldpco1: usize,
    ldpco2: usize,
    qcoeff: &mut [f64],
    ldqco1: usize,
    ldqco2: usize,
) -> i32 {
    if indlim < 1 {
        return -4;
    }
    let (porm, porp) = if leri == Leri::Left {
        (p, m)
    } else {
        (m, p)
    };
    if ldpco1 < porm || ldpco2 < porm {
        return -6;
    }
    if ldqco1 < p.max(m).max(p) || ldqco2 < p.max(m).max(p) {
        return -9;
    }
    // Transpose P: P'(i,j,k) = P(j,i,k). In-place transpose of each slice.
    for k in 0..indlim {
        for i in 0..porm {
            for j in (i + 1)..porm {
                let idx_ij = k * ldpco1 * ldpco2 + i * ldpco2 + j;
                let idx_ji = k * ldpco1 * ldpco2 + j * ldpco2 + i;
                pcoeff.swap(idx_ij, idx_ji);
            }
        }
    }
    // Transpose Q: output Q' is M×P (leading part). Input is porm×porp (left: P×M, right: M×P).
    // Q'(j,i,k) = Q(i,j,k). Copy each k-slice to temp then write transposed.
    let mut tmp = vec![0.0; porm * porp];
    for k in 0..indlim {
        let base = k * ldqco1 * ldqco2;
        for i in 0..porm {
            for j in 0..porp {
                tmp[i * porp + j] = qcoeff[base + i * ldqco2 + j];
            }
        }
        for j in 0..porp {
            for i in 0..porm {
                qcoeff[base + j * ldqco2 + i] = tmp[i * porp + j];
            }
        }
    }
    0
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_tc01od_left_2x2_indlim3() {
        // TC01OD example: M=2, P=2, INDLIM=3, L. P 2×2×3. Fortran column-major (i,j,k)= i + j*LDPCO1 + k*LDPCO1*LDPCO2.
        // P (1,1)=2,3,1; (2,1)=5,7,-6; (1,2)=4,-1,-1; (2,2)=3,2,2 -> flat: k=0: 2,5,4,3; k=1: 3,7,-1,2; k=2: 1,-6,-1,2
        let mut pcoeff = vec![
            2.0, 5.0, 4.0, 3.0, 3.0, 7.0, -1.0, 2.0, 1.0, -6.0, -1.0, 2.0,
        ];
        // Q 2×2: (1,1)=6,-1,5 (1,2)=1,1,1 (2,1)=1,7,5 (2,2)=4,1,-1 -> col-major k=0: 6,1,1,4; k=1: -1,1,7,1; k=2: 5,1,5,-1
        let mut qcoeff = vec![
            6.0, 1.0, 1.0, 4.0, -1.0, 1.0, 7.0, 1.0, 5.0, 1.0, 5.0, -1.0,
        ];
        let ldpco1 = 2;
        let ldpco2 = 2;
        let ldqco1 = 2;
        let ldqco2 = 2;
        let info = tc01od(
            Leri::Left,
            2,
            2,
            3,
            &mut pcoeff,
            ldpco1,
            ldpco2,
            &mut qcoeff,
            ldqco1,
            ldqco2,
        );
        assert_eq!(info, 0);
        // After transpose P': (1,0) gets P(0,1). Col-major (1,0)=index 1. Original (0,1)=index 2 = 4. So pcoeff[1]=4.
        assert!((pcoeff[0] - 2.0).abs() < 1e-10);
        assert!((pcoeff[1] - 4.0).abs() < 1e-10);
        assert!((pcoeff[2] - 5.0).abs() < 1e-10);
        assert!((pcoeff[3] - 3.0).abs() < 1e-10);
        // Q' (transpose): col-major k=0 (0,0)=6,(1,0)=1,(0,1)=1,(1,1)=4
        assert!((qcoeff[0] - 6.0).abs() < 1e-10);
        assert!((qcoeff[1] - 1.0).abs() < 1e-10);
        assert!((qcoeff[2] - 1.0).abs() < 1e-10);
        assert!((qcoeff[3] - 4.0).abs() < 1e-10);
    }
}
