# SLICOT → slicot-rs mapping index

One-to-one mapping of SLICOT (Fortran) routines to Pure Rust in the `slicot-rs` crate. Includes all routines from `slicot_module.f90` and `lapack_aux`. Status: `not started` | `in progress` | `done`. Rust function names are lowercase per plan.

| SLICOT | Rust module | Rust function | Status |
|--------|-------------|---------------|--------|
| AB01MD | ab01 | ab01md | done |
| AB01ND | ab01 | ab01nd | done |
| AB01OD | ab01 | ab01od | done |
| AB04MD | ab04 | ab04md | done |
| AB05MD | ab05 | ab05md | done |
| AB05ND | ab05 | ab05nd | done |
| AB05OD | ab05 | ab05od | done |
| AB05PD | ab05 | ab05pd | done |
| AB05QD | ab05 | ab05qd | done |
| AB05RD | ab05 | ab05rd | done |
| AB05SD | ab05 | ab05sd | done |
| AB07MD | ab07 | ab07md | done |
| AB07ND | ab07 | ab07nd | done |
| AB08MD | ab08 | ab08md | done |
| AB08MZ | ab08 | ab08mz | done |
| AB08ND | ab08 | ab08nd | done |
| AB08NW | ab08 | ab08nw | done |
| AB08NX | ab08 | ab08nx | done |
| AB08NY | ab08 | ab08ny | done |
| AB08NZ | ab08 | ab08nz | done |
| AB09AD | ab09 | ab09ad | done |
| AB09AX | ab09 | ab09ax | done |
| AB09BD | ab09 | ab09bd | done |
| AB09BX | ab09 | ab09bx | done |
| AB09CD | ab09 | ab09cd | done |
| AB09CX | ab09 | ab09cx | done |
| AB09DD | ab09 | ab09dd | done |
| AB09ED | ab09 | ab09ed | done |
| AB09FD | ab09 | ab09fd | done |
| AB09GD | ab09 | ab09gd | done |
| AB09HD | ab09 | ab09hd | done |
| AB09HX | ab09 | ab09hx | done |
| AB09HY | ab09 | ab09hy | done |
| AB09ID | ab09 | ab09id | done |
| AB09IX | ab09 | ab09ix | done |
| AB09IY | ab09 | ab09iy | done |
| AB09JD | ab09 | ab09jd | done |
| AB09JV | ab09 | ab09jv | done |
| AB09JW | ab09 | ab09jw | done |
| AB09JX | ab09 | ab09jx | done |
| AB09KD | ab09 | ab09kd | done |
| AB09KX | ab09 | ab09kx | done |
| AB09MD | ab09 | ab09md | done |
| AB09ND | ab09 | ab09nd | done |
| AB13AD | ab13 | ab13ad | done |
| AB13AX | ab13 | ab13ax | done |
| AB13BD | ab13 | ab13bd | done |
| AB13CD | ab13 | ab13cd | done |
| AB13DD | ab13 | ab13dd | done |
| AB13DX | ab13 | ab13dx | done |
| AB13ED | ab13 | ab13ed | done |
| AB13FD | ab13 | ab13fd | done |
| AB13HD | ab13 | ab13hd | done |
| AB13ID | ab13 | ab13id | done |
| AB13MD | ab13 | ab13md | done |
| AB8NXZ | ab8n | ab8nxz | done |
| AG07BD | ag07 | ag07bd | done |
| AG08BD | ag08 | ag08bd | done |
| AG08BY | ag08 | ag08by | done |
| AG08BZ | ag08 | ag08bz | done |
| AG8BYZ | ag8b | ag8byz | done |
| BB01AD | bb01 | bb01ad | done |
| BB02AD | bb02 | bb02ad | done |
| BB03AD | bb03 | bb03ad | done |
| BB04AD | bb04 | bb04ad | done |
| BD01AD | bd01 | bd01ad | done |
| BD02AD | bd02 | bd02ad | done |
| DE01OD | de01 | de01od | done |
| DE01PD | de01 | de01pd | done |
| DF01MD | df01 | df01md | done |
| DG01MD | dg01 | dg01md | done |
| DG01ND | dg01 | dg01nd | done |
| DG01NY | dg01 | dg01ny | done |
| DG01OD | dg01 | dg01od | done |
| DGEGS | dgeg | dgegs | done |
| DGEGV | dgeg | dgegv | done |
| DK01MD | dk01 | dk01md | done |
| DLACPY_SLC | dlac | dlacpy_slc | done |
| DLATZM | dlat | dlatzm | done |
| FB01QD | fb01 | fb01qd | done |
| FB01RD | fb01 | fb01rd | done |
| FB01SD | fb01 | fb01sd | done |
| FB01TD | fb01 | fb01td | done |
| FB01VD | fb01 | fb01vd | done |
| FD01AD | fd01 | fd01ad | done |
| IB01AD | ib01 | ib01ad | done |
| IB01BD | ib01 | ib01bd | done |
| IB01CD | ib01 | ib01cd | done |
| IB01MD | ib01 | ib01md | done |
| IB01MY | ib01 | ib01my | done |
| IB01ND | ib01 | ib01nd | done |
| IB01OD | ib01 | ib01od | done |
| IB01OY | ib01 | ib01oy | done |
| IB01PD | ib01 | ib01pd | done |
| IB01PX | ib01 | ib01px | done |
| IB01PY | ib01 | ib01py | done |
| IB01QD | ib01 | ib01qd | done |
| IB01RD | ib01 | ib01rd | done |
| IB03AD | ib03 | ib03ad | done |
| IB03BD | ib03 | ib03bd | done |
| MA01AD | ma01 | ma01ad | done |
| MA01BD | ma01 | ma01bd | done |
| MA01BZ | ma01 | ma01bz | done |
| MA01CD | ma01 | ma01cd | done |
| MA01DD | ma01 | ma01dd | done |
| MA01DZ | ma01 | ma01dz | done |
| MA02AD | ma02 | ma02ad | done |
| MA02AZ | ma02 | ma02az | done |
| MA02BD | ma02 | ma02bd | done |
| MA02BZ | ma02 | ma02bz | done |
| MA02CD | ma02 | ma02cd | done |
| MA02CZ | ma02 | ma02cz | done |
| MA02DD | ma02 | ma02dd | done |
| MA02ED | ma02 | ma02ed | done |
| MA02ES | ma02 | ma02es | done |
| MA02EZ | ma02 | ma02ez | done |
| MA02FD | ma02 | ma02fd | done |
| MA02GD | ma02 | ma02gd | done |
| MA02GZ | ma02 | ma02gz | done |
| MA02HD | ma02 | ma02hd | done |
| MA02HZ | ma02 | ma02hz | done |
| MA02ID | ma02 | ma02id | done |
| MA02IZ | ma02 | ma02iz | done |
| MA02JD | ma02 | ma02jd | done |
| MA02JZ | ma02 | ma02jz | done |
| MA02MD | ma02 | ma02md | done |
| MA02MZ | ma02 | ma02mz | done |
| MA02NZ | ma02 | ma02nz | done |
| MA02OD | ma02 | ma02od | done |
| MA02OZ | ma02 | ma02oz | done |
| MA02PD | ma02 | ma02pd | done |
| MA02PZ | ma02 | ma02pz | done |
| MA02RD | ma02 | ma02rd | done |
| MA02SD | ma02 | ma02sd | done |
| MB01KD | mb01 | mb01kd | done |
| MB01LD | mb01 | mb01ld | done |
| MB01MD | mb01 | mb01md | done |
| MB01ND | mb01 | mb01nd | done |
| MB01OC | mb01 | mb01oc | done |
| MB01OD | mb01 | mb01od | done |
| MB01OE | mb01 | mb01oe | done |
| MB01OH | mb01 | mb01oh | done |
| MB01OO | mb01 | mb01oo | done |
| MB01OS | mb01 | mb01os | done |
| MB01OT | mb01 | mb01ot | done |
| MB01PD | mb01 | mb01pd | done |
| MB01QD | mb01 | mb01qd | done |
| MB01RB | mb01 | mb01rb | done |
| MB01RD | mb01 | mb01rd | done |
| MB01RH | mb01 | mb01rh | done |
| MB01RT | mb01 | mb01rt | done |
| MB01RU | mb01 | mb01ru | done |
| MB01RW | mb01 | mb01rw | done |
| MB01RX | mb01 | mb01rx | done |
| MB01RY | mb01 | mb01ry | done |
| MB01SD | mb01 | mb01sd | done |
| MB01SS | mb01 | mb01ss | done |
| MB01TD | mb01 | mb01td | done |
| MB01UD | mb01 | mb01ud | done |
| MB01UW | mb01 | mb01uw | done |
| MB01UX | mb01 | mb01ux | done |
| MB01UY | mb01 | mb01uy | done |
| MB01UZ | mb01 | mb01uz | done |
| MB01VD | mb01 | mb01vd | done |
| MB01WD | mb01 | mb01wd | done |
| MB01XD | mb01 | mb01xd | done |
| MB01XY | mb01 | mb01xy | done |
| MB01YD | mb01 | mb01yd | done |
| MB01ZD | mb01 | mb01zd | done |
| MB02CD | mb02 | mb02cd | done |
| MB02CU | mb02 | mb02cu | done |
| MB02CV | mb02 | mb02cv | done |
| MB02CX | mb02 | mb02cx | done |
| MB02CY | mb02 | mb02cy | done |
| MB02DD | mb02 | mb02dd | done |
| MB02ED | mb02 | mb02ed | done |
| MB02FD | mb02 | mb02fd | done |
| MB02GD | mb02 | mb02gd | done |
| MB02HD | mb02 | mb02hd | done |
| MB02ID | mb02 | mb02id | done |
| MB02JD | mb02 | mb02jd | done |
| MB02JX | mb02 | mb02jx | done |
| MB02KD | mb02 | mb02kd | done |
| MB02MD | mb02 | mb02md | done |
| MB02ND | mb02 | mb02nd | done |
| MB02NY | mb02 | mb02ny | done |
| MB02OD | mb02 | mb02od | done |
| MB02PD | mb02 | mb02pd | done |
| MB02QD | mb02 | mb02qd | done |
| MB02QY | mb02 | mb02qy | done |
| MB02RD | mb02 | mb02rd | done |
| MB02RZ | mb02 | mb02rz | done |
| MB02SD | mb02 | mb02sd | done |
| MB02SZ | mb02 | mb02sz | done |
| MB02TD | mb02 | mb02td | done |
| MB02TZ | mb02 | mb02tz | done |
| MB02UD | mb02 | mb02ud | done |
| MB02UU | mb02 | mb02uu | done |
| MB02UV | mb02 | mb02uv | done |
| MB02UW | mb02 | mb02uw | done |
| MB02VD | mb02 | mb02vd | done |
| MB02WD | mb02 | mb02wd | done |
| MB02XD | mb02 | mb02xd | done |
| MB02YD | mb02 | mb02yd | done |
| MB03AB | mb03 | mb03ab | done |
| MB03AD | mb03 | mb03ad | done |
| MB03AE | mb03 | mb03ae | done |
| MB03AF | mb03 | mb03af | done |
| MB03AG | mb03 | mb03ag | done |
| MB03AH | mb03 | mb03ah | done |
| MB03AI | mb03 | mb03ai | done |
| MB03BA | mb03 | mb03ba | done |
| MB03BB | mb03 | mb03bb | done |
| MB03BC | mb03 | mb03bc | done |
| MB03BD | mb03 | mb03bd | done |
| MB03BE | mb03 | mb03be | done |
| MB03BF | mb03 | mb03bf | done |
| MB03BG | mb03 | mb03bg | done |
| MB03BZ | mb03 | mb03bz | done |
| MB03CD | mb03 | mb03cd | done |
| MB03CZ | mb03 | mb03cz | done |
| MB03DD | mb03 | mb03dd | done |
| MB03DZ | mb03 | mb03dz | done |
| MB03ED | mb03 | mb03ed | done |
| MB03FD | mb03 | mb03fd | done |
| MB03FZ | mb03 | mb03fz | done |
| MB03GD | mb03 | mb03gd | done |
| MB03GZ | mb03 | mb03gz | done |
| MB03HD | mb03 | mb03hd | done |
| MB03HZ | mb03 | mb03hz | done |
| MB03ID | mb03 | mb03id | done |
| MB03IZ | mb03 | mb03iz | done |
| MB03JD | mb03 | mb03jd | done |
| MB03JP | mb03 | mb03jp | done |
| MB03JZ | mb03 | mb03jz | done |
| MB03KA | mb03 | mb03ka | done |
| MB03KB | mb03 | mb03kb | done |
| MB03KC | mb03 | mb03kc | done |
| MB03KD | mb03 | mb03kd | done |
| MB03KE | mb03 | mb03ke | done |
| MB03LD | mb03 | mb03ld | done |
| MB03LF | mb03 | mb03lf | done |
| MB03LP | mb03 | mb03lp | done |
| MB03LZ | mb03 | mb03lz | done |
| MB03MD | mb03 | mb03md | done |
| MB03MY | mb03 | mb03my | done |
| MB03ND | mb03 | mb03nd | done |
| MB03NY | mb03 | mb03ny | done |
| MB03OD | mb03 | mb03od | done |
| MB03OY | mb03 | mb03oy | done |
| MB03PD | mb03 | mb03pd | done |
| MB03PY | mb03 | mb03py | done |
| MB03QD | mb03 | mb03qd | done |
| MB03QG | mb03 | mb03qg | done |
| MB03QV | mb03 | mb03qv | done |
| MB03QW | mb03 | mb03qw | done |
| MB03QX | mb03 | mb03qx | done |
| MB03QY | mb03 | mb03qy | done |
| MB03RD | mb03 | mb03rd | done |
| MB03RW | mb03 | mb03rw | done |
| MB03RX | mb03 | mb03rx | done |
| MB03RY | mb03 | mb03ry | done |
| MB03RZ | mb03 | mb03rz | done |
| MB03SD | mb03 | mb03sd | done |
| MB03TD | mb03 | mb03td | done |
| MB03TS | mb03 | mb03ts | done |
| MB03UD | mb03 | mb03ud | done |
| MB03VD | mb03 | mb03vd | done |
| MB03VW | mb03 | mb03vw | done |
| MB03VY | mb03 | mb03vy | done |
| MB03WA | mb03 | mb03wa | done |
| MB03WD | mb03 | mb03wd | done |
| MB03WX | mb03 | mb03wx | done |
| MB03XD | mb03 | mb03xd | done |
| MB03XP | mb03 | mb03xp | done |
| MB03XS | mb03 | mb03xs | done |
| MB03XU | mb03 | mb03xu | done |
| MB03XZ | mb03 | mb03xz | done |
| MB03YA | mb03 | mb03ya | done |
| MB03YD | mb03 | mb03yd | done |
| MB03YT | mb03 | mb03yt | done |
| MB03ZA | mb03 | mb03za | done |
| MB03ZD | mb03 | mb03zd | done |
| MB04AD | mb04 | mb04ad | done |
| MB04AZ | mb04 | mb04az | done |
| MB04BD | mb04 | mb04bd | done |
| MB04BP | mb04 | mb04bp | done |
| MB04BZ | mb04 | mb04bz | done |
| MB04CD | mb04 | mb04cd | done |
| MB04DB | mb04 | mb04db | done |
| MB04DD | mb04 | mb04dd | done |
| MB04DI | mb04 | mb04di | done |
| MB04DL | mb04 | mb04dl | done |
| MB04DP | mb04 | mb04dp | done |
| MB04DS | mb04 | mb04ds | done |
| MB04DY | mb04 | mb04dy | done |
| MB04DZ | mb04 | mb04dz | done |
| MB04ED | mb04 | mb04ed | done |
| MB04FD | mb04 | mb04fd | done |
| MB04FP | mb04 | mb04fp | done |
| MB04GD | mb04 | mb04gd | done |
| MB04HD | mb04 | mb04hd | done |
| MB04ID | mb04 | mb04id | done |
| MB04IY | mb04 | mb04iy | done |
| MB04IZ | mb04 | mb04iz | done |
| MB04JD | mb04 | mb04jd | done |
| MB04KD | mb04 | mb04kd | done |
| MB04LD | mb04 | mb04ld | done |
| MB04MD | mb04 | mb04md | done |
| MB04ND | mb04 | mb04nd | done |
| MB04NY | mb04 | mb04ny | done |
| MB04OD | mb04 | mb04od | done |
| MB04OW | mb04 | mb04ow | done |
| MB04OX | mb04 | mb04ox | done |
| MB04OY | mb04 | mb04oy | done |
| MB04PA | mb04 | mb04pa | done |
| MB04PB | mb04 | mb04pb | done |
| MB04PU | mb04 | mb04pu | done |
| MB04PY | mb04 | mb04py | done |
| MB04QB | mb04 | mb04qb | done |
| MB04QC | mb04 | mb04qc | done |
| MB04QF | mb04 | mb04qf | done |
| MB04QS | mb04 | mb04qs | done |
| MB04QU | mb04 | mb04qu | done |
| MB04RB | mb04 | mb04rb | done |
| MB04RD | mb04 | mb04rd | done |
| MB04RS | mb04 | mb04rs | done |
| MB04RT | mb04 | mb04rt | done |
| MB04RU | mb04 | mb04ru | done |
| MB04RV | mb04 | mb04rv | done |
| MB04RW | mb04 | mb04rw | done |
| MB04RZ | mb04 | mb04rz | done |
| MB04SU | mb04 | mb04su | done |
| MB04TB | mb04 | mb04tb | done |
| MB04TS | mb04 | mb04ts | done |
| MB04TT | mb04 | mb04tt | done |
| MB04TU | mb04 | mb04tu | done |
| MB04TV | mb04 | mb04tv | done |
| MB04TW | mb04 | mb04tw | done |
| MB04TX | mb04 | mb04tx | done |
| MB04TY | mb04 | mb04ty | done |
| MB04UD | mb04 | mb04ud | done |
| MB04VD | mb04 | mb04vd | done |
| MB04VX | mb04 | mb04vx | done |
| MB04WD | mb04 | mb04wd | done |
| MB04WP | mb04 | mb04wp | done |
| MB04WR | mb04 | mb04wr | done |
| MB04WU | mb04 | mb04wu | done |
| MB04XD | mb04 | mb04xd | done |
| MB04XY | mb04 | mb04xy | done |
| MB04YD | mb04 | mb04yd | done |
| MB04YW | mb04 | mb04yw | done |
| MB04ZD | mb04 | mb04zd | done |
| MB05MD | mb05 | mb05md | done |
| MB05MY | mb05 | mb05my | done |
| MB05ND | mb05 | mb05nd | done |
| MB05OD | mb05 | mb05od | done |
| MB05OY | mb05 | mb05oy | done |
| MB3JZP | mb3j | mb3jzp | not started |
| MB3LZP | mb3l | mb3lzp | not started |
| MB3OYZ | mb3o | mb3oyz | not started |
| MB3PYZ | mb3p | mb3pyz | not started |
| MB4DBZ | mb4d | mb4dbz | not started |
| MB4DLZ | mb4d | mb4dlz | not started |
| MB4DPZ | mb4d | mb4dpz | not started |
| MC01MD | mc01 | mc01md | done |
| MC01ND | mc01 | mc01nd | done |
| MC01OD | mc01 | mc01od | done |
| MC01PD | mc01 | mc01pd | done |
| MC01PY | mc01 | mc01py | done |
| MC01QD | mc01 | mc01qd | done |
| MC01RD | mc01 | mc01rd | done |
| MC01SD | mc01 | mc01sd | done |
| MC01SW | mc01 | mc01sw | done |
| MC01SX | mc01 | mc01sx | done |
| MC01SY | mc01 | mc01sy | done |
| MC01TD | mc01 | mc01td | done |
| MC01VD | mc01 | mc01vd | done |
| MC01WD | mc01 | mc01wd | done |
| MC01XD | mc01 | mc01xd | done |
| MC03MD | mc03 | mc03md | done |
| MC03ND | mc03 | mc03nd | done |
| MC03NX | mc03 | mc03nx | done |
| MC03NY | mc03 | mc03ny | done |
| MD03AD | md03 | md03ad | not started |
| MD03BA | md03 | md03ba | not started |
| MD03BB | md03 | md03bb | not started |
| MD03BD | md03 | md03bd | not started |
| MD03BF | md03 | md03bf | not started |
| MD03BX | md03 | md03bx | not started |
| MD03BY | md03 | md03by | not started |
| NF01AD | nf01 | nf01ad | not started |
| NF01AY | nf01 | nf01ay | not started |
| NF01BA | nf01 | nf01ba | not started |
| NF01BB | nf01 | nf01bb | not started |
| NF01BD | nf01 | nf01bd | not started |
| NF01BE | nf01 | nf01be | not started |
| NF01BF | nf01 | nf01bf | not started |
| NF01BP | nf01 | nf01bp | not started |
| NF01BQ | nf01 | nf01bq | not started |
| NF01BR | nf01 | nf01br | not started |
| NF01BS | nf01 | nf01bs | not started |
| NF01BU | nf01 | nf01bu | not started |
| NF01BV | nf01 | nf01bv | not started |
| NF01BW | nf01 | nf01bw | not started |
| NF01BX | nf01 | nf01bx | not started |
| NF01BY | nf01 | nf01by | not started |
| SB01BD | sb01 | sb01bd | done |
| SB01BX | sb01 | sb01bx | done |
| SB01BY | sb01 | sb01by | done |
| SB01DD | sb01 | sb01dd | done |
| SB01FY | sb01 | sb01fy | done |
| SB01MD | sb01 | sb01md | done |
| SB02CX | sb02 | sb02cx | done |
| SB02MD | sb02 | sb02md | done |
| SB02MR | sb02 | sb02mr | done |
| SB02MS | sb02 | sb02ms | done |
| SB02MT | sb02 | sb02mt | done |
| SB02MU | sb02 | sb02mu | done |
| SB02MV | sb02 | sb02mv | done |
| SB02MW | sb02 | sb02mw | done |
| SB02MX | sb02 | sb02mx | done |
| SB02ND | sb02 | sb02nd | done |
| SB02OD | sb02 | sb02od | done |
| SB02OU | sb02 | sb02ou | done |
| SB02OV | sb02 | sb02ov | done |
| SB02OW | sb02 | sb02ow | done |
| SB02OX | sb02 | sb02ox | done |
| SB02OY | sb02 | sb02oy | done |
| SB02PD | sb02 | sb02pd | done |
| SB02QD | sb02 | sb02qd | done |
| SB02RD | sb02 | sb02rd | done |
| SB02RU | sb02 | sb02ru | done |
| SB02SD | sb02 | sb02sd | done |
| SB03MD | sb03 | sb03md | done |
| SB03MU | sb03 | sb03mu | done |
| SB03MV | sb03 | sb03mv | done |
| SB03MW | sb03 | sb03mw | done |
| SB03MX | sb03 | sb03mx | done |
| SB03MY | sb03 | sb03my | done |
| SB03OD | sb03 | sb03od | done |
| SB03OR | sb03 | sb03or | done |
| SB03OS | sb03 | sb03os | done |
| SB03OT | sb03 | sb03ot | done |
| SB03OU | sb03 | sb03ou | done |
| SB03OV | sb03 | sb03ov | done |
| SB03OY | sb03 | sb03oy | done |
| SB03OZ | sb03 | sb03oz | done |
| SB03PD | sb03 | sb03pd | done |
| SB03QD | sb03 | sb03qd | done |
| SB03QX | sb03 | sb03qx | done |
| SB03QY | sb03 | sb03qy | done |
| SB03RD | sb03 | sb03rd | done |
| SB03SD | sb03 | sb03sd | done |
| SB03SX | sb03 | sb03sx | done |
| SB03SY | sb03 | sb03sy | done |
| SB03TD | sb03 | sb03td | done |
| SB03UD | sb03 | sb03ud | done |
| SB04MD | sb04 | sb04md | done |
| SB04MR | sb04 | sb04mr | done |
| SB04MU | sb04 | sb04mu | done |
| SB04MW | sb04 | sb04mw | done |
| SB04MY | sb04 | sb04my | done |
| SB04ND | sb04 | sb04nd | done |
| SB04NV | sb04 | sb04nv | done |
| SB04NW | sb04 | sb04nw | done |
| SB04NX | sb04 | sb04nx | done |
| SB04NY | sb04 | sb04ny | done |
| SB04OD | sb04 | sb04od | done |
| SB04OW | sb04 | sb04ow | done |
| SB04PD | sb04 | sb04pd | done |
| SB04PX | sb04 | sb04px | done |
| SB04PY | sb04 | sb04py | done |
| SB04QD | sb04 | sb04qd | done |
| SB04QR | sb04 | sb04qr | done |
| SB04QU | sb04 | sb04qu | done |
| SB04QY | sb04 | sb04qy | done |
| SB04RD | sb04 | sb04rd | done |
| SB04RV | sb04 | sb04rv | done |
| SB04RW | sb04 | sb04rw | done |
| SB04RX | sb04 | sb04rx | done |
| SB04RY | sb04 | sb04ry | done |
| SB06ND | sb06 | sb06nd | done |
| SB08CD | sb08 | sb08cd | done |
| SB08DD | sb08 | sb08dd | done |
| SB08ED | sb08 | sb08ed | done |
| SB08FD | sb08 | sb08fd | done |
| SB08GD | sb08 | sb08gd | done |
| SB08HD | sb08 | sb08hd | done |
| SB08MD | sb08 | sb08md | done |
| SB08MY | sb08 | sb08my | done |
| SB08ND | sb08 | sb08nd | done |
| SB08NY | sb08 | sb08ny | done |
| SB09MD | sb09 | sb09md | not started |
| SB10AD | sb10 | sb10ad | not started |
| SB10DD | sb10 | sb10dd | not started |
| SB10ED | sb10 | sb10ed | not started |
| SB10FD | sb10 | sb10fd | not started |
| SB10HD | sb10 | sb10hd | not started |
| SB10ID | sb10 | sb10id | not started |
| SB10JD | sb10 | sb10jd | not started |
| SB10KD | sb10 | sb10kd | not started |
| SB10LD | sb10 | sb10ld | not started |
| SB10MD | sb10 | sb10md | not started |
| SB10PD | sb10 | sb10pd | not started |
| SB10QD | sb10 | sb10qd | not started |
| SB10RD | sb10 | sb10rd | not started |
| SB10SD | sb10 | sb10sd | not started |
| SB10TD | sb10 | sb10td | not started |
| SB10UD | sb10 | sb10ud | not started |
| SB10VD | sb10 | sb10vd | not started |
| SB10WD | sb10 | sb10wd | not started |
| SB10YD | sb10 | sb10yd | not started |
| SB10ZD | sb10 | sb10zd | not started |
| SB10ZP | sb10 | sb10zp | not started |
| SB16AD | sb16 | sb16ad | not started |
| SB16AY | sb16 | sb16ay | not started |
| SB16BD | sb16 | sb16bd | not started |
| SB16CD | sb16 | sb16cd | not started |
| SB16CY | sb16 | sb16cy | not started |
| SG02AD | sg02 | sg02ad | not started |
| SG02CV | sg02 | sg02cv | not started |
| SG02CW | sg02 | sg02cw | not started |
| SG02CX | sg02 | sg02cx | not started |
| SG02ND | sg02 | sg02nd | not started |
| SG03AD | sg03 | sg03ad | not started |
| SG03AX | sg03 | sg03ax | not started |
| SG03AY | sg03 | sg03ay | not started |
| SG03BD | sg03 | sg03bd | not started |
| SG03BR | sg03 | sg03br | not started |
| SG03BS | sg03 | sg03bs | not started |
| SG03BT | sg03 | sg03bt | not started |
| SG03BU | sg03 | sg03bu | not started |
| SG03BV | sg03 | sg03bv | not started |
| SG03BW | sg03 | sg03bw | not started |
| SG03BX | sg03 | sg03bx | not started |
| SG03BY | sg03 | sg03by | not started |
| SG03BZ | sg03 | sg03bz | not started |
| TB01ID | tb01 | tb01id | not started |
| TB01IZ | tb01 | tb01iz | not started |
| TB01KD | tb01 | tb01kd | not started |
| TB01KX | tb01 | tb01kx | not started |
| TB01LD | tb01 | tb01ld | not started |
| TB01MD | tb01 | tb01md | done |
| TB01ND | tb01 | tb01nd | not started |
| TB01PD | tb01 | tb01pd | not started |
| TB01PX | tb01 | tb01px | not started |
| TB01TD | tb01 | tb01td | not started |
| TB01TY | tb01 | tb01ty | not started |
| TB01UD | tb01 | tb01ud | not started |
| TB01UX | tb01 | tb01ux | not started |
| TB01UY | tb01 | tb01uy | not started |
| TB01VD | tb01 | tb01vd | not started |
| TB01VY | tb01 | tb01vy | not started |
| TB01WD | tb01 | tb01wd | not started |
| TB01WX | tb01 | tb01wx | not started |
| TB01XD | tb01 | tb01xd | not started |
| TB01XZ | tb01 | tb01xz | not started |
| TB01YD | tb01 | tb01yd | not started |
| TB01ZD | tb01 | tb01zd | not started |
| TB03AD | tb03 | tb03ad | not started |
| TB03AY | tb03 | tb03ay | not started |
| TB04AD | tb04 | tb04ad | not started |
| TB04AY | tb04 | tb04ay | not started |
| TB04BD | tb04 | tb04bd | not started |
| TB04BV | tb04 | tb04bv | not started |
| TB04BW | tb04 | tb04bw | not started |
| TB04BX | tb04 | tb04bx | not started |
| TB04CD | tb04 | tb04cd | not started |
| TB05AD | tb05 | tb05ad | not started |
| TC01OD | tc01 | tc01od | not started |
| TC04AD | tc04 | tc04ad | not started |
| TC05AD | tc05 | tc05ad | not started |
| TD03AD | td03 | td03ad | not started |
| TD03AY | td03 | td03ay | not started |
| TD04AD | td04 | td04ad | not started |
| TD05AD | td05 | td05ad | not started |
| TF01MD | tf01 | tf01md | not started |
| TF01MX | tf01 | tf01mx | not started |
| TF01MY | tf01 | tf01my | not started |
| TF01ND | tf01 | tf01nd | not started |
| TF01OD | tf01 | tf01od | not started |
| TF01PD | tf01 | tf01pd | not started |
| TF01QD | tf01 | tf01qd | not started |
| TF01RD | tf01 | tf01rd | not started |
| TG01AD | tg01 | tg01ad | not started |
| TG01AZ | tg01 | tg01az | not started |
| TG01BD | tg01 | tg01bd | not started |
| TG01CD | tg01 | tg01cd | not started |
| TG01DD | tg01 | tg01dd | not started |
| TG01ED | tg01 | tg01ed | not started |
| TG01FD | tg01 | tg01fd | not started |
| TG01FZ | tg01 | tg01fz | not started |
| TG01GD | tg01 | tg01gd | not started |
| TG01HD | tg01 | tg01hd | not started |
| TG01HU | tg01 | tg01hu | not started |
| TG01HX | tg01 | tg01hx | not started |
| TG01HY | tg01 | tg01hy | not started |
| TG01ID | tg01 | tg01id | not started |
| TG01JD | tg01 | tg01jd | not started |
| TG01JY | tg01 | tg01jy | not started |
| TG01KD | tg01 | tg01kd | not started |
| TG01KZ | tg01 | tg01kz | not started |
| TG01LD | tg01 | tg01ld | not started |
| TG01LY | tg01 | tg01ly | not started |
| TG01MD | tg01 | tg01md | not started |
| TG01ND | tg01 | tg01nd | not started |
| TG01NX | tg01 | tg01nx | not started |
| TG01OA | tg01 | tg01oa | not started |
| TG01OB | tg01 | tg01ob | not started |
| TG01OD | tg01 | tg01od | not started |
| TG01OZ | tg01 | tg01oz | not started |
| TG01PD | tg01 | tg01pd | not started |
| TG01QD | tg01 | tg01qd | not started |
| TG01WD | tg01 | tg01wd | not started |
| UD01BD | ud01 | ud01bd | not started |
| UD01CD | ud01 | ud01cd | not started |
| UD01DD | ud01 | ud01dd | not started |
| UD01MD | ud01 | ud01md | not started |
| UD01MZ | ud01 | ud01mz | not started |
| UD01ND | ud01 | ud01nd | not started |
| UE01MD | ue01 | ue01md | not started |
| ZGEGS | zgeg | zgegs | not started |
| ZGEGV | zgeg | zgegv | not started |
| ZLATZM | zlat | zlatzm | not started |
