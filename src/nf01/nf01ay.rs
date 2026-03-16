//! NF01AY — Output of set of neural networks (SLICOT NF01AY)
//
// y = sum_i ws(i)*tanh(wi'*z + bi) + b(n+1) per output. NN neurons, NZ input dim, L outputs.

/// IPAR(1) = NN (neurons). WB = [wb(1)..wb(L)], each wb length NN*(NZ+2)+1:
/// [w1(1:NZ),...,wn(1:NZ), ws(1:n), b(1:n+1)]. Z(NSMP,NZ), Y(NSMP,L) output.
pub fn nf01ay(
    nsmp: i32,
    nz: i32,
    l: i32,
    ipar: &[i32],
    _lipar: i32,
    wb: &[f64],
    _lwb: i32,
    z: &[f64],
    ldz: usize,
    y: &mut [f64],
    ldy: usize,
    _dwork: &mut [f64],
    _ldwork: i32,
) -> i32 {
    let (nsmp, nz, l) = (nsmp as usize, nz as usize, l as usize);
    let nn = *ipar.get(0).unwrap_or(&0) as usize;
    if ldy < nsmp.max(1) || ldz < nsmp.max(1) {
        return -10;
    }
    let wb_len = (nn * (nz + 2) + 1) * l;
    if wb.len() < wb_len {
        return -7;
    }
    if z.len() < nsmp * ldz || y.len() < nsmp * ldy {
        return -9;
    }
    let wb_stride = nn * (nz + 2) + 1;
    for s in 0..nsmp {
        for out in 0..l {
            let wb_off = out * wb_stride;
            let mut sum = 0.0;
            for i in 0..nn {
                let mut dot = 0.0;
                for j in 0..nz {
                    dot += wb[wb_off + i * nz + j] * z[s + j * ldz];
                }
                dot += wb[wb_off + nn * nz + nn + i];
                sum += wb[wb_off + nn * nz + i] * dot.tanh();
            }
            sum += wb[wb_off + nn * nz + nn + nn];
            y[s + out * ldy] = sum;
        }
    }
    0
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_nf01ay_one_sample() {
        let ipar = [1_i32];
        let mut wb = vec![0.0; 10];
        wb[0] = 1.0;
        wb[1] = 0.0;
        wb[2] = 1.0;
        wb[3] = 0.0;
        let z = [0.0];
        let mut y = [0.0];
        let mut dwork = [0.0; 4];
        assert_eq!(
            nf01ay(1, 1, 1, &ipar, 1, &wb, 10, &z, 1, &mut y, 1, &mut dwork, 4),
            0
        );
        assert!((y[0] - 0.0_f64.tanh()).abs() < 1e-10);
    }
}
