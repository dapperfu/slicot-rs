//! MA02OD — Number of zero rows of a real (skew-)Hamiltonian matrix (SLICOT MA02OD)
//
// H = [A D; E ±A']. DE: lower triangle = E, columns 1..M = upper D. Returns count of zero rows.

#[derive(Clone, Copy, PartialEq, Eq)]
pub enum Ma02OdSkew {
    Hamiltonian,    // E=E', D=D', bottom block = A'
    SkewHamiltonian, // E=-E', D=-D', bottom block = -A'
}

/// Unpack E(i,j) from DE. skew => E(i,j) = -E(j,i) for i < j.
fn get_e(de: &[f64], ldde: usize, _m: usize, skew: bool, i: usize, j: usize) -> f64 {
    if i >= j {
        de[i + j * ldde]
    } else if skew {
        -de[j + i * ldde]
    } else {
        de[j + i * ldde]
    }
}

/// Unpack D(i,j) from DE (upper stored in columns 1..M). skew => D(i,j) = -D(j,i) for i > j.
fn get_d(de: &[f64], ldde: usize, _m: usize, skew: bool, i: usize, j: usize) -> f64 {
    if i <= j {
        de[i + (j + 1) * ldde]
    } else if skew {
        -de[j + (i + 1) * ldde]
    } else {
        de[j + (i + 1) * ldde]
    }
}

/// Returns the number of zero rows of H. A is M×M (LDA×M), DE is M×(M+1) (LDDE×(M+1)).
pub fn ma02od(
    skew: Ma02OdSkew,
    m: usize,
    a: &[f64],
    lda: usize,
    de: &[f64],
    ldde: usize,
) -> usize {
    if m == 0 || lda < m || ldde < m {
        return 0;
    }
    let skew_bool = skew == Ma02OdSkew::SkewHamiltonian;
    let sign = if skew_bool { -1.0 } else { 1.0 };
    let mut count = 0_usize;

    for i in 0..m {
        let mut row_zero = true;
        for j in 0..m {
            if a[i + j * lda].abs() > 0.0 {
                row_zero = false;
                break;
            }
        }
        if row_zero {
            for j in 0..m {
                if get_d(de, ldde, m, skew_bool, i, j).abs() > 0.0 {
                    row_zero = false;
                    break;
                }
            }
        }
        if row_zero {
            count += 1;
        }
    }

    for i in 0..m {
        let mut row_zero = true;
        for j in 0..m {
            if get_e(de, ldde, m, skew_bool, i, j).abs() > 0.0 {
                row_zero = false;
                break;
            }
        }
        if row_zero {
            for j in 0..m {
                if (sign * a[j + i * lda]).abs() > 0.0 {
                    row_zero = false;
                    break;
                }
            }
        }
        if row_zero {
            count += 1;
        }
    }
    count
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_ma02od_zero_dim() {
        assert_eq!(ma02od(Ma02OdSkew::SkewHamiltonian, 0, &[], 0, &[], 0), 0);
    }

    #[test]
    fn test_ma02od_identity_block() {
        // H = [I 0; 0 I] (skew-Hamiltonian with A=I, D=0, E=0). All rows non-zero except we need A=I, D=0, E=0 -> 0 zero rows.
        let m = 2;
        let a = vec![1.0, 0.0, 0.0, 1.0];
        let de = vec![0.0; m * (m + 1)];
        assert_eq!(ma02od(Ma02OdSkew::SkewHamiltonian, m, &a, 2, &de, 2), 0);
    }

    #[test]
    fn test_ma02od_all_zero() {
        let m = 2;
        let a = vec![0.0; 4];
        let de = vec![0.0; m * (m + 1)];
        assert_eq!(ma02od(Ma02OdSkew::SkewHamiltonian, m, &a, 2, &de, 2), 4);
    }
}
