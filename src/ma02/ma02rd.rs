//! MA02RD — Sort D (increasing or decreasing) and rearrange E by same permutation (SLICOT MA02RD)

/// Sorts d in place into increasing ('I') or decreasing ('D') order and applies the same permutation to e.
pub fn ma02rd(id: char, d: &mut [f64], e: &mut [f64]) -> i32 {
    let n = d.len();
    if e.len() != n {
        return -3;
    }
    if n <= 1 {
        return 0;
    }
    let increasing = match id {
        'I' | 'i' => true,
        'D' | 'd' => false,
        _ => return -1,
    };
    let mut indices: Vec<usize> = (0..n).collect();
    if increasing {
        indices.sort_by(|&i, &j| d[i].partial_cmp(&d[j]).unwrap_or(std::cmp::Ordering::Equal));
    } else {
        indices.sort_by(|&i, &j| d[j].partial_cmp(&d[i]).unwrap_or(std::cmp::Ordering::Equal));
    }
    let d_old: Vec<f64> = d.to_vec();
    let e_old: Vec<f64> = e.to_vec();
    for (i, &j) in indices.iter().enumerate() {
        d[i] = d_old[j];
        e[i] = e_old[j];
    }
    0
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_ma02rd_increasing() {
        let mut d = [3.0, 1.0, 2.0];
        let mut e = [30.0, 10.0, 20.0];
        assert_eq!(ma02rd('I', &mut d, &mut e), 0);
        assert_eq!(d, [1.0, 2.0, 3.0]);
        assert_eq!(e, [10.0, 20.0, 30.0]);
    }
}
