//! UE01MD — Default machine-specific parameters (SLICOT extension of ILAENV).
//!
//! Returns default values for blocksize, crossover, etc. Good for a wide range of machines;
//! for optimal performance the user may tune. ISPEC: 1=opt blocksize, 2=min block size,
//! 3=crossover (use unblocked below), 4=number of shifts, 8=crossover multishift QR.

/// Returns the default parameter for the given ISPEC. NAME and OPTS identify the calling
/// routine; N1, N2, N3 are problem dimensions (may be unused). Values aim for good performance
/// on a wide range of computers. Argument order matches Fortran: (ISPEC, NAME, OPTS, N1, N2, N3).
pub fn ue01md(ispec: i32, _name: &str, _opts: &str, _n1: i32, _n2: i32, _n3: i32) -> i32 {
    match ispec {
        1 => 32,   // optimal blocksize
        2 => 2,    // minimum block size for blocked routine
        3 => 128,  // crossover: use unblocked for N < this
        4 => 6,    // number of shifts
        8 => 50,   // crossover for multishift QR
        _ => 1,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_ue01md_ispec_1() {
        assert!(ue01md(1, "MB03WD", "", 0, 0, 0) >= 1);
    }

    #[test]
    fn test_ue01md_ispec_3() {
        let v = ue01md(3, "MB03WD", "S", 100, 0, 0);
        assert!(v >= 0);
    }

    #[test]
    fn test_ue01md_unknown_ispec() {
        assert_eq!(ue01md(99, "X", "", 0, 0, 0), 1);
    }
}
