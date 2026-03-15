# SLICOT-rs feature status

Implemented vs unimplemented SLICOT routines. The authoritative per-routine list is [docs/SLICOT_MAPPING.md](SLICOT_MAPPING.md).

## Summary by module

| Module | Implemented | Unimplemented | Total | % done |
|--------|-------------|----------------|-------|--------|
| ab01 | 3 | 0 | 3 | 100 |
| ab04 | 1 | 0 | 1 | 100 |
| ab05 | 7 | 0 | 7 | 100 |
| ab07 | 2 | 0 | 2 | 100 |
| ab08 | 7 | 0 | 7 | 100 |
| ab09 | 24 | 0 | 24 | 100 |
| ab13 | 11 | 0 | 11 | 100 |
| ab8n | 1 | 0 | 1 | 100 |
| ag07 | 1 | 0 | 1 | 100 |
| ag08 | 3 | 0 | 3 | 100 |
| ag8b | 1 | 0 | 1 | 100 |
| bb01 | 1 | 0 | 1 | 100 |
| bb02 | 1 | 0 | 1 | 100 |
| bb03 | 1 | 0 | 1 | 100 |
| bb04 | 1 | 0 | 1 | 100 |
| bd01 | 1 | 0 | 1 | 100 |
| bd02 | 1 | 0 | 1 | 100 |
| de01 | 2 | 0 | 2 | 100 |
| df01 | 1 | 0 | 1 | 100 |
| dg01 | 4 | 0 | 4 | 100 |
| dgeg | 2 | 0 | 2 | 100 |
| dk01 | 1 | 0 | 1 | 100 |
| dlac | 1 | 0 | 1 | 100 |
| dlat | 1 | 0 | 1 | 100 |
| fb01 | 5 | 0 | 5 | 100 |
| fd01 | 1 | 0 | 1 | 100 |
| ib01 | 13 | 0 | 13 | 100 |
| ib03 | 2 | 0 | 2 | 100 |
| ma01 | 6 | 0 | 6 | 100 |
| ma02 | 28 | 0 | 28 | 100 |
| mb01 | 35 | 0 | 35 | 100 |
| mb02 | 35 | 0 | 35 | 100 |
| mb03 | 79 | 0 | 79 | 100 |
| mb04 | 70 | 0 | 70 | 100 |
| mb05 | 5 | 0 | 5 | 100 |
| mb3j | 1 | 0 | 1 | 100 |
| mb3l | 1 | 0 | 1 | 100 |
| mb3o | 1 | 0 | 1 | 100 |
| mb3p | 1 | 0 | 1 | 100 |
| mb4d | 3 | 0 | 3 | 100 |
| mc01 | 15 | 0 | 15 | 100 |
| mc03 | 4 | 0 | 4 | 100 |
| md03 | 7 | 0 | 7 | 100 |
| nf01 | 16 | 0 | 16 | 100 |
| sb01 | 6 | 0 | 6 | 100 |
| sb02 | 21 | 0 | 21 | 100 |
| sb03 | 24 | 0 | 24 | 100 |
| sb04 | 24 | 0 | 24 | 100 |
| sb06 | 1 | 0 | 1 | 100 |
| sb08 | 10 | 0 | 10 | 100 |
| sb09 | 1 | 0 | 1 | 100 |
| sb10 | 0 | 21 | 21 | 0 |
| sb16 | 0 | 5 | 5 | 0 |
| sg02 | 0 | 5 | 5 | 0 |
| sg03 | 0 | 13 | 13 | 0 |
| tb01 | 1 | 21 | 22 | 5 |
| tb03 | 0 | 2 | 2 | 0 |
| tb04 | 0 | 7 | 7 | 0 |
| tb05 | 0 | 1 | 1 | 0 |
| tc01 | 0 | 1 | 1 | 0 |
| tc04 | 0 | 1 | 1 | 0 |
| tc05 | 0 | 1 | 1 | 0 |
| td03 | 0 | 2 | 2 | 0 |
| td04 | 0 | 1 | 1 | 0 |
| td05 | 0 | 1 | 1 | 0 |
| tf01 | 0 | 8 | 8 | 0 |
| tg01 | 0 | 30 | 30 | 0 |
| ud01 | 0 | 6 | 6 | 0 |
| ue01 | 0 | 1 | 1 | 0 |
| zgeg | 0 | 2 | 2 | 0 |
| zlat | 0 | 1 | 1 | 0 |

| **Total** | **487** | **138** | **625** | **77** |

## Full table (all routines)

| SLICOT | Rust module | Rust function | Status |
|--------|-------------|---------------|--------|
| AB01MD | ab01 | ab01md | Implemented |
| AB01ND | ab01 | ab01nd | Implemented |
| AB01OD | ab01 | ab01od | Implemented |
| AB04MD | ab04 | ab04md | Implemented |
| AB05MD | ab05 | ab05md | Implemented |
| AB05ND | ab05 | ab05nd | Implemented |
| AB05OD | ab05 | ab05od | Implemented |
| AB05PD | ab05 | ab05pd | Implemented |
| AB05QD | ab05 | ab05qd | Implemented |
| AB05RD | ab05 | ab05rd | Implemented |
| AB05SD | ab05 | ab05sd | Implemented |
| AB07MD | ab07 | ab07md | Implemented |
| AB07ND | ab07 | ab07nd | Implemented |
| AB08MD | ab08 | ab08md | Implemented |
| AB08MZ | ab08 | ab08mz | Implemented |
| AB08ND | ab08 | ab08nd | Implemented |
| AB08NW | ab08 | ab08nw | Implemented |
| AB08NX | ab08 | ab08nx | Implemented |
| AB08NY | ab08 | ab08ny | Implemented |
| AB08NZ | ab08 | ab08nz | Implemented |
| AB09AD | ab09 | ab09ad | Implemented |
| AB09AX | ab09 | ab09ax | Implemented |
| AB09BD | ab09 | ab09bd | Implemented |
| AB09BX | ab09 | ab09bx | Implemented |
| AB09CD | ab09 | ab09cd | Implemented |
| AB09CX | ab09 | ab09cx | Implemented |
| AB09DD | ab09 | ab09dd | Implemented |
| AB09ED | ab09 | ab09ed | Implemented |
| AB09FD | ab09 | ab09fd | Implemented |
| AB09GD | ab09 | ab09gd | Implemented |
| AB09HD | ab09 | ab09hd | Implemented |
| AB09HX | ab09 | ab09hx | Implemented |
| AB09HY | ab09 | ab09hy | Implemented |
| AB09ID | ab09 | ab09id | Implemented |
| AB09IX | ab09 | ab09ix | Implemented |
| AB09IY | ab09 | ab09iy | Implemented |
| AB09JD | ab09 | ab09jd | Implemented |
| AB09JV | ab09 | ab09jv | Implemented |
| AB09JW | ab09 | ab09jw | Implemented |
| AB09JX | ab09 | ab09jx | Implemented |
| AB09KD | ab09 | ab09kd | Implemented |
| AB09KX | ab09 | ab09kx | Implemented |
| AB09MD | ab09 | ab09md | Implemented |
| AB09ND | ab09 | ab09nd | Implemented |
| AB13AD | ab13 | ab13ad | Implemented |
| AB13AX | ab13 | ab13ax | Implemented |
| AB13BD | ab13 | ab13bd | Implemented |
| AB13CD | ab13 | ab13cd | Implemented |
| AB13DD | ab13 | ab13dd | Implemented |
| AB13DX | ab13 | ab13dx | Implemented |
| AB13ED | ab13 | ab13ed | Implemented |
| AB13FD | ab13 | ab13fd | Implemented |
| AB13HD | ab13 | ab13hd | Implemented |
| AB13ID | ab13 | ab13id | Implemented |
| AB13MD | ab13 | ab13md | Implemented |
| AB8NXZ | ab8n | ab8nxz | Implemented |
| AG07BD | ag07 | ag07bd | Implemented |
| AG08BD | ag08 | ag08bd | Implemented |
| AG08BY | ag08 | ag08by | Implemented |
| AG08BZ | ag08 | ag08bz | Implemented |
| AG8BYZ | ag8b | ag8byz | Implemented |
| BB01AD | bb01 | bb01ad | Implemented |
| BB02AD | bb02 | bb02ad | Implemented |
| BB03AD | bb03 | bb03ad | Implemented |
| BB04AD | bb04 | bb04ad | Implemented |
| BD01AD | bd01 | bd01ad | Implemented |
| BD02AD | bd02 | bd02ad | Implemented |
| DE01OD | de01 | de01od | Implemented |
| DE01PD | de01 | de01pd | Implemented |
| DF01MD | df01 | df01md | Implemented |
| DG01MD | dg01 | dg01md | Implemented |
| DG01ND | dg01 | dg01nd | Implemented |
| DG01NY | dg01 | dg01ny | Implemented |
| DG01OD | dg01 | dg01od | Implemented |
| DGEGS | dgeg | dgegs | Implemented |
| DGEGV | dgeg | dgegv | Implemented |
| DK01MD | dk01 | dk01md | Implemented |
| DLACPY_SLC | dlac | dlacpy_slc | Implemented |
| DLATZM | dlat | dlatzm | Implemented |
| FB01QD | fb01 | fb01qd | Implemented |
| FB01RD | fb01 | fb01rd | Implemented |
| FB01SD | fb01 | fb01sd | Implemented |
| FB01TD | fb01 | fb01td | Implemented |
| FB01VD | fb01 | fb01vd | Implemented |
| FD01AD | fd01 | fd01ad | Implemented |
| IB01AD | ib01 | ib01ad | Implemented |
| IB01BD | ib01 | ib01bd | Implemented |
| IB01CD | ib01 | ib01cd | Implemented |
| IB01MD | ib01 | ib01md | Implemented |
| IB01MY | ib01 | ib01my | Implemented |
| IB01ND | ib01 | ib01nd | Implemented |
| IB01OD | ib01 | ib01od | Implemented |
| IB01OY | ib01 | ib01oy | Implemented |
| IB01PD | ib01 | ib01pd | Implemented |
| IB01PX | ib01 | ib01px | Implemented |
| IB01PY | ib01 | ib01py | Implemented |
| IB01QD | ib01 | ib01qd | Implemented |
| IB01RD | ib01 | ib01rd | Implemented |
| IB03AD | ib03 | ib03ad | Implemented |
| IB03BD | ib03 | ib03bd | Implemented |
| MA01AD | ma01 | ma01ad | Implemented |
| MA01BD | ma01 | ma01bd | Implemented |
| MA01BZ | ma01 | ma01bz | Implemented |
| MA01CD | ma01 | ma01cd | Implemented |
| MA01DD | ma01 | ma01dd | Implemented |
| MA01DZ | ma01 | ma01dz | Implemented |
| MA02AD | ma02 | ma02ad | Implemented |
| MA02AZ | ma02 | ma02az | Implemented |
| MA02BD | ma02 | ma02bd | Implemented |
| MA02BZ | ma02 | ma02bz | Implemented |
| MA02CD | ma02 | ma02cd | Implemented |
| MA02CZ | ma02 | ma02cz | Implemented |
| MA02DD | ma02 | ma02dd | Implemented |
| MA02ED | ma02 | ma02ed | Implemented |
| MA02ES | ma02 | ma02es | Implemented |
| MA02EZ | ma02 | ma02ez | Implemented |
| MA02FD | ma02 | ma02fd | Implemented |
| MA02GD | ma02 | ma02gd | Implemented |
| MA02GZ | ma02 | ma02gz | Implemented |
| MA02HD | ma02 | ma02hd | Implemented |
| MA02HZ | ma02 | ma02hz | Implemented |
| MA02ID | ma02 | ma02id | Implemented |
| MA02IZ | ma02 | ma02iz | Implemented |
| MA02JD | ma02 | ma02jd | Implemented |
| MA02JZ | ma02 | ma02jz | Implemented |
| MA02MD | ma02 | ma02md | Implemented |
| MA02MZ | ma02 | ma02mz | Implemented |
| MA02NZ | ma02 | ma02nz | Implemented |
| MA02OD | ma02 | ma02od | Implemented |
| MA02OZ | ma02 | ma02oz | Implemented |
| MA02PD | ma02 | ma02pd | Implemented |
| MA02PZ | ma02 | ma02pz | Implemented |
| MA02RD | ma02 | ma02rd | Implemented |
| MA02SD | ma02 | ma02sd | Implemented |
| MB01KD | mb01 | mb01kd | Implemented |
| MB01LD | mb01 | mb01ld | Implemented |
| MB01MD | mb01 | mb01md | Implemented |
| MB01ND | mb01 | mb01nd | Implemented |
| MB01OC | mb01 | mb01oc | Implemented |
| MB01OD | mb01 | mb01od | Implemented |
| MB01OE | mb01 | mb01oe | Implemented |
| MB01OH | mb01 | mb01oh | Implemented |
| MB01OO | mb01 | mb01oo | Implemented |
| MB01OS | mb01 | mb01os | Implemented |
| MB01OT | mb01 | mb01ot | Implemented |
| MB01PD | mb01 | mb01pd | Implemented |
| MB01QD | mb01 | mb01qd | Implemented |
| MB01RB | mb01 | mb01rb | Implemented |
| MB01RD | mb01 | mb01rd | Implemented |
| MB01RH | mb01 | mb01rh | Implemented |
| MB01RT | mb01 | mb01rt | Implemented |
| MB01RU | mb01 | mb01ru | Implemented |
| MB01RW | mb01 | mb01rw | Implemented |
| MB01RX | mb01 | mb01rx | Implemented |
| MB01RY | mb01 | mb01ry | Implemented |
| MB01SD | mb01 | mb01sd | Implemented |
| MB01SS | mb01 | mb01ss | Implemented |
| MB01TD | mb01 | mb01td | Implemented |
| MB01UD | mb01 | mb01ud | Implemented |
| MB01UW | mb01 | mb01uw | Implemented |
| MB01UX | mb01 | mb01ux | Implemented |
| MB01UY | mb01 | mb01uy | Implemented |
| MB01UZ | mb01 | mb01uz | Implemented |
| MB01VD | mb01 | mb01vd | Implemented |
| MB01WD | mb01 | mb01wd | Implemented |
| MB01XD | mb01 | mb01xd | Implemented |
| MB01XY | mb01 | mb01xy | Implemented |
| MB01YD | mb01 | mb01yd | Implemented |
| MB01ZD | mb01 | mb01zd | Implemented |
| MB02CD | mb02 | mb02cd | Implemented |
| MB02CU | mb02 | mb02cu | Implemented |
| MB02CV | mb02 | mb02cv | Implemented |
| MB02CX | mb02 | mb02cx | Implemented |
| MB02CY | mb02 | mb02cy | Implemented |
| MB02DD | mb02 | mb02dd | Implemented |
| MB02ED | mb02 | mb02ed | Implemented |
| MB02FD | mb02 | mb02fd | Implemented |
| MB02GD | mb02 | mb02gd | Implemented |
| MB02HD | mb02 | mb02hd | Implemented |
| MB02ID | mb02 | mb02id | Implemented |
| MB02JD | mb02 | mb02jd | Implemented |
| MB02JX | mb02 | mb02jx | Implemented |
| MB02KD | mb02 | mb02kd | Implemented |
| MB02MD | mb02 | mb02md | Implemented |
| MB02ND | mb02 | mb02nd | Implemented |
| MB02NY | mb02 | mb02ny | Implemented |
| MB02OD | mb02 | mb02od | Implemented |
| MB02PD | mb02 | mb02pd | Implemented |
| MB02QD | mb02 | mb02qd | Implemented |
| MB02QY | mb02 | mb02qy | Implemented |
| MB02RD | mb02 | mb02rd | Implemented |
| MB02RZ | mb02 | mb02rz | Implemented |
| MB02SD | mb02 | mb02sd | Implemented |
| MB02SZ | mb02 | mb02sz | Implemented |
| MB02TD | mb02 | mb02td | Implemented |
| MB02TZ | mb02 | mb02tz | Implemented |
| MB02UD | mb02 | mb02ud | Implemented |
| MB02UU | mb02 | mb02uu | Implemented |
| MB02UV | mb02 | mb02uv | Implemented |
| MB02UW | mb02 | mb02uw | Implemented |
| MB02VD | mb02 | mb02vd | Implemented |
| MB02WD | mb02 | mb02wd | Implemented |
| MB02XD | mb02 | mb02xd | Implemented |
| MB02YD | mb02 | mb02yd | Implemented |
| MB03AB | mb03 | mb03ab | Implemented |
| MB03AD | mb03 | mb03ad | Implemented |
| MB03AE | mb03 | mb03ae | Implemented |
| MB03AF | mb03 | mb03af | Implemented |
| MB03AG | mb03 | mb03ag | Implemented |
| MB03AH | mb03 | mb03ah | Implemented |
| MB03AI | mb03 | mb03ai | Implemented |
| MB03BA | mb03 | mb03ba | Implemented |
| MB03BB | mb03 | mb03bb | Implemented |
| MB03BC | mb03 | mb03bc | Implemented |
| MB03BD | mb03 | mb03bd | Implemented |
| MB03BE | mb03 | mb03be | Implemented |
| MB03BF | mb03 | mb03bf | Implemented |
| MB03BG | mb03 | mb03bg | Implemented |
| MB03BZ | mb03 | mb03bz | Implemented |
| MB03CD | mb03 | mb03cd | Implemented |
| MB03CZ | mb03 | mb03cz | Implemented |
| MB03DD | mb03 | mb03dd | Implemented |
| MB03DZ | mb03 | mb03dz | Implemented |
| MB03ED | mb03 | mb03ed | Implemented |
| MB03FD | mb03 | mb03fd | Implemented |
| MB03FZ | mb03 | mb03fz | Implemented |
| MB03GD | mb03 | mb03gd | Implemented |
| MB03GZ | mb03 | mb03gz | Implemented |
| MB03HD | mb03 | mb03hd | Implemented |
| MB03HZ | mb03 | mb03hz | Implemented |
| MB03ID | mb03 | mb03id | Implemented |
| MB03IZ | mb03 | mb03iz | Implemented |
| MB03JD | mb03 | mb03jd | Implemented |
| MB03JP | mb03 | mb03jp | Implemented |
| MB03JZ | mb03 | mb03jz | Implemented |
| MB03KA | mb03 | mb03ka | Implemented |
| MB03KB | mb03 | mb03kb | Implemented |
| MB03KC | mb03 | mb03kc | Implemented |
| MB03KD | mb03 | mb03kd | Implemented |
| MB03KE | mb03 | mb03ke | Implemented |
| MB03LD | mb03 | mb03ld | Implemented |
| MB03LF | mb03 | mb03lf | Implemented |
| MB03LP | mb03 | mb03lp | Implemented |
| MB03LZ | mb03 | mb03lz | Implemented |
| MB03MD | mb03 | mb03md | Implemented |
| MB03MY | mb03 | mb03my | Implemented |
| MB03ND | mb03 | mb03nd | Implemented |
| MB03NY | mb03 | mb03ny | Implemented |
| MB03OD | mb03 | mb03od | Implemented |
| MB03OY | mb03 | mb03oy | Implemented |
| MB03PD | mb03 | mb03pd | Implemented |
| MB03PY | mb03 | mb03py | Implemented |
| MB03QD | mb03 | mb03qd | Implemented |
| MB03QG | mb03 | mb03qg | Implemented |
| MB03QV | mb03 | mb03qv | Implemented |
| MB03QW | mb03 | mb03qw | Implemented |
| MB03QX | mb03 | mb03qx | Implemented |
| MB03QY | mb03 | mb03qy | Implemented |
| MB03RD | mb03 | mb03rd | Implemented |
| MB03RW | mb03 | mb03rw | Implemented |
| MB03RX | mb03 | mb03rx | Implemented |
| MB03RY | mb03 | mb03ry | Implemented |
| MB03RZ | mb03 | mb03rz | Implemented |
| MB03SD | mb03 | mb03sd | Implemented |
| MB03TD | mb03 | mb03td | Implemented |
| MB03TS | mb03 | mb03ts | Implemented |
| MB03UD | mb03 | mb03ud | Implemented |
| MB03VD | mb03 | mb03vd | Implemented |
| MB03VW | mb03 | mb03vw | Implemented |
| MB03VY | mb03 | mb03vy | Implemented |
| MB03WA | mb03 | mb03wa | Implemented |
| MB03WD | mb03 | mb03wd | Implemented |
| MB03WX | mb03 | mb03wx | Implemented |
| MB03XD | mb03 | mb03xd | Implemented |
| MB03XP | mb03 | mb03xp | Implemented |
| MB03XS | mb03 | mb03xs | Implemented |
| MB03XU | mb03 | mb03xu | Implemented |
| MB03XZ | mb03 | mb03xz | Implemented |
| MB03YA | mb03 | mb03ya | Implemented |
| MB03YD | mb03 | mb03yd | Implemented |
| MB03YT | mb03 | mb03yt | Implemented |
| MB03ZA | mb03 | mb03za | Implemented |
| MB03ZD | mb03 | mb03zd | Implemented |
| MB04AD | mb04 | mb04ad | Implemented |
| MB04AZ | mb04 | mb04az | Implemented |
| MB04BD | mb04 | mb04bd | Implemented |
| MB04BP | mb04 | mb04bp | Implemented |
| MB04BZ | mb04 | mb04bz | Implemented |
| MB04CD | mb04 | mb04cd | Implemented |
| MB04DB | mb04 | mb04db | Implemented |
| MB04DD | mb04 | mb04dd | Implemented |
| MB04DI | mb04 | mb04di | Implemented |
| MB04DL | mb04 | mb04dl | Implemented |
| MB04DP | mb04 | mb04dp | Implemented |
| MB04DS | mb04 | mb04ds | Implemented |
| MB04DY | mb04 | mb04dy | Implemented |
| MB04DZ | mb04 | mb04dz | Implemented |
| MB04ED | mb04 | mb04ed | Implemented |
| MB04FD | mb04 | mb04fd | Implemented |
| MB04FP | mb04 | mb04fp | Implemented |
| MB04GD | mb04 | mb04gd | Implemented |
| MB04HD | mb04 | mb04hd | Implemented |
| MB04ID | mb04 | mb04id | Implemented |
| MB04IY | mb04 | mb04iy | Implemented |
| MB04IZ | mb04 | mb04iz | Implemented |
| MB04JD | mb04 | mb04jd | Implemented |
| MB04KD | mb04 | mb04kd | Implemented |
| MB04LD | mb04 | mb04ld | Implemented |
| MB04MD | mb04 | mb04md | Implemented |
| MB04ND | mb04 | mb04nd | Implemented |
| MB04NY | mb04 | mb04ny | Implemented |
| MB04OD | mb04 | mb04od | Implemented |
| MB04OW | mb04 | mb04ow | Implemented |
| MB04OX | mb04 | mb04ox | Implemented |
| MB04OY | mb04 | mb04oy | Implemented |
| MB04PA | mb04 | mb04pa | Implemented |
| MB04PB | mb04 | mb04pb | Implemented |
| MB04PU | mb04 | mb04pu | Implemented |
| MB04PY | mb04 | mb04py | Implemented |
| MB04QB | mb04 | mb04qb | Implemented |
| MB04QC | mb04 | mb04qc | Implemented |
| MB04QF | mb04 | mb04qf | Implemented |
| MB04QS | mb04 | mb04qs | Implemented |
| MB04QU | mb04 | mb04qu | Implemented |
| MB04RB | mb04 | mb04rb | Implemented |
| MB04RD | mb04 | mb04rd | Implemented |
| MB04RS | mb04 | mb04rs | Implemented |
| MB04RT | mb04 | mb04rt | Implemented |
| MB04RU | mb04 | mb04ru | Implemented |
| MB04RV | mb04 | mb04rv | Implemented |
| MB04RW | mb04 | mb04rw | Implemented |
| MB04RZ | mb04 | mb04rz | Implemented |
| MB04SU | mb04 | mb04su | Implemented |
| MB04TB | mb04 | mb04tb | Implemented |
| MB04TS | mb04 | mb04ts | Implemented |
| MB04TT | mb04 | mb04tt | Implemented |
| MB04TU | mb04 | mb04tu | Implemented |
| MB04TV | mb04 | mb04tv | Implemented |
| MB04TW | mb04 | mb04tw | Implemented |
| MB04TX | mb04 | mb04tx | Implemented |
| MB04TY | mb04 | mb04ty | Implemented |
| MB04UD | mb04 | mb04ud | Implemented |
| MB04VD | mb04 | mb04vd | Implemented |
| MB04VX | mb04 | mb04vx | Implemented |
| MB04WD | mb04 | mb04wd | Implemented |
| MB04WP | mb04 | mb04wp | Implemented |
| MB04WR | mb04 | mb04wr | Implemented |
| MB04WU | mb04 | mb04wu | Implemented |
| MB04XD | mb04 | mb04xd | Implemented |
| MB04XY | mb04 | mb04xy | Implemented |
| MB04YD | mb04 | mb04yd | Implemented |
| MB04YW | mb04 | mb04yw | Implemented |
| MB04ZD | mb04 | mb04zd | Implemented |
| MB05MD | mb05 | mb05md | Implemented |
| MB05MY | mb05 | mb05my | Implemented |
| MB05ND | mb05 | mb05nd | Implemented |
| MB05OD | mb05 | mb05od | Implemented |
| MB05OY | mb05 | mb05oy | Implemented |
| MB3JZP | mb3j | mb3jzp | Unimplemented |
| MB3LZP | mb3l | mb3lzp | Unimplemented |
| MB3OYZ | mb3o | mb3oyz | Unimplemented |
| MB3PYZ | mb3p | mb3pyz | Unimplemented |
| MB4DBZ | mb4d | mb4dbz | Unimplemented |
| MB4DLZ | mb4d | mb4dlz | Unimplemented |
| MB4DPZ | mb4d | mb4dpz | Unimplemented |
| MC01MD | mc01 | mc01md | Implemented |
| MC01ND | mc01 | mc01nd | Implemented |
| MC01OD | mc01 | mc01od | Implemented |
| MC01PD | mc01 | mc01pd | Implemented |
| MC01PY | mc01 | mc01py | Implemented |
| MC01QD | mc01 | mc01qd | Implemented |
| MC01RD | mc01 | mc01rd | Implemented |
| MC01SD | mc01 | mc01sd | Implemented |
| MC01SW | mc01 | mc01sw | Implemented |
| MC01SX | mc01 | mc01sx | Implemented |
| MC01SY | mc01 | mc01sy | Implemented |
| MC01TD | mc01 | mc01td | Implemented |
| MC01VD | mc01 | mc01vd | Implemented |
| MC01WD | mc01 | mc01wd | Implemented |
| MC01XD | mc01 | mc01xd | Implemented |
| MC03MD | mc03 | mc03md | Implemented |
| MC03ND | mc03 | mc03nd | Implemented |
| MC03NX | mc03 | mc03nx | Implemented |
| MC03NY | mc03 | mc03ny | Implemented |
| MD03AD | md03 | md03ad | Implemented |
| MD03BA | md03 | md03ba | Implemented |
| MD03BB | md03 | md03bb | Implemented |
| MD03BD | md03 | md03bd | Implemented |
| MD03BF | md03 | md03bf | Implemented |
| MD03BX | md03 | md03bx | Implemented |
| MD03BY | md03 | md03by | Implemented |
| NF01AD | nf01 | nf01ad | Implemented |
| NF01AY | nf01 | nf01ay | Implemented |
| NF01BA | nf01 | nf01ba | Implemented |
| NF01BB | nf01 | nf01bb | Implemented |
| NF01BD | nf01 | nf01bd | Implemented |
| NF01BE | nf01 | nf01be | Implemented |
| NF01BF | nf01 | nf01bf | Implemented |
| NF01BP | nf01 | nf01bp | Implemented |
| NF01BQ | nf01 | nf01bq | Implemented |
| NF01BR | nf01 | nf01br | Implemented |
| NF01BS | nf01 | nf01bs | Implemented |
| NF01BU | nf01 | nf01bu | Implemented |
| NF01BV | nf01 | nf01bv | Implemented |
| NF01BW | nf01 | nf01bw | Implemented |
| NF01BX | nf01 | nf01bx | Implemented |
| NF01BY | nf01 | nf01by | Implemented |
| SB01BD | sb01 | sb01bd | Implemented |
| SB01BX | sb01 | sb01bx | Implemented |
| SB01BY | sb01 | sb01by | Implemented |
| SB01DD | sb01 | sb01dd | Implemented |
| SB01FY | sb01 | sb01fy | Implemented |
| SB01MD | sb01 | sb01md | Implemented |
| SB02CX | sb02 | sb02cx | Implemented |
| SB02MD | sb02 | sb02md | Implemented |
| SB02MR | sb02 | sb02mr | Implemented |
| SB02MS | sb02 | sb02ms | Implemented |
| SB02MT | sb02 | sb02mt | Implemented |
| SB02MU | sb02 | sb02mu | Implemented |
| SB02MV | sb02 | sb02mv | Implemented |
| SB02MW | sb02 | sb02mw | Implemented |
| SB02MX | sb02 | sb02mx | Implemented |
| SB02ND | sb02 | sb02nd | Implemented |
| SB02OD | sb02 | sb02od | Implemented |
| SB02OU | sb02 | sb02ou | Implemented |
| SB02OV | sb02 | sb02ov | Implemented |
| SB02OW | sb02 | sb02ow | Implemented |
| SB02OX | sb02 | sb02ox | Implemented |
| SB02OY | sb02 | sb02oy | Implemented |
| SB02PD | sb02 | sb02pd | Implemented |
| SB02QD | sb02 | sb02qd | Implemented |
| SB02RD | sb02 | sb02rd | Implemented |
| SB02RU | sb02 | sb02ru | Implemented |
| SB02SD | sb02 | sb02sd | Implemented |
| SB03MD | sb03 | sb03md | Implemented |
| SB03MU | sb03 | sb03mu | Implemented |
| SB03MV | sb03 | sb03mv | Implemented |
| SB03MW | sb03 | sb03mw | Implemented |
| SB03MX | sb03 | sb03mx | Implemented |
| SB03MY | sb03 | sb03my | Implemented |
| SB03OD | sb03 | sb03od | Implemented |
| SB03OR | sb03 | sb03or | Implemented |
| SB03OS | sb03 | sb03os | Implemented |
| SB03OT | sb03 | sb03ot | Implemented |
| SB03OU | sb03 | sb03ou | Implemented |
| SB03OV | sb03 | sb03ov | Implemented |
| SB03OY | sb03 | sb03oy | Implemented |
| SB03OZ | sb03 | sb03oz | Implemented |
| SB03PD | sb03 | sb03pd | Implemented |
| SB03QD | sb03 | sb03qd | Implemented |
| SB03QX | sb03 | sb03qx | Implemented |
| SB03QY | sb03 | sb03qy | Implemented |
| SB03RD | sb03 | sb03rd | Implemented |
| SB03SD | sb03 | sb03sd | Implemented |
| SB03SX | sb03 | sb03sx | Implemented |
| SB03SY | sb03 | sb03sy | Implemented |
| SB03TD | sb03 | sb03td | Implemented |
| SB03UD | sb03 | sb03ud | Implemented |
| SB04MD | sb04 | sb04md | Implemented |
| SB04MR | sb04 | sb04mr | Implemented |
| SB04MU | sb04 | sb04mu | Implemented |
| SB04MW | sb04 | sb04mw | Implemented |
| SB04MY | sb04 | sb04my | Implemented |
| SB04ND | sb04 | sb04nd | Implemented |
| SB04NV | sb04 | sb04nv | Implemented |
| SB04NW | sb04 | sb04nw | Implemented |
| SB04NX | sb04 | sb04nx | Implemented |
| SB04NY | sb04 | sb04ny | Implemented |
| SB04OD | sb04 | sb04od | Implemented |
| SB04OW | sb04 | sb04ow | Implemented |
| SB04PD | sb04 | sb04pd | Implemented |
| SB04PX | sb04 | sb04px | Implemented |
| SB04PY | sb04 | sb04py | Implemented |
| SB04QD | sb04 | sb04qd | Implemented |
| SB04QR | sb04 | sb04qr | Implemented |
| SB04QU | sb04 | sb04qu | Implemented |
| SB04QY | sb04 | sb04qy | Implemented |
| SB04RD | sb04 | sb04rd | Implemented |
| SB04RV | sb04 | sb04rv | Implemented |
| SB04RW | sb04 | sb04rw | Implemented |
| SB04RX | sb04 | sb04rx | Implemented |
| SB04RY | sb04 | sb04ry | Implemented |
| SB06ND | sb06 | sb06nd | Implemented |
| SB08CD | sb08 | sb08cd | Implemented |
| SB08DD | sb08 | sb08dd | Implemented |
| SB08ED | sb08 | sb08ed | Implemented |
| SB08FD | sb08 | sb08fd | Implemented |
| SB08GD | sb08 | sb08gd | Implemented |
| SB08HD | sb08 | sb08hd | Implemented |
| SB08MD | sb08 | sb08md | Implemented |
| SB08MY | sb08 | sb08my | Implemented |
| SB08ND | sb08 | sb08nd | Implemented |
| SB08NY | sb08 | sb08ny | Implemented |
| SB09MD | sb09 | sb09md | Implemented |
| SB10AD | sb10 | sb10ad | Unimplemented |
| SB10DD | sb10 | sb10dd | Unimplemented |
| SB10ED | sb10 | sb10ed | Unimplemented |
| SB10FD | sb10 | sb10fd | Unimplemented |
| SB10HD | sb10 | sb10hd | Unimplemented |
| SB10ID | sb10 | sb10id | Unimplemented |
| SB10JD | sb10 | sb10jd | Unimplemented |
| SB10KD | sb10 | sb10kd | Unimplemented |
| SB10LD | sb10 | sb10ld | Unimplemented |
| SB10MD | sb10 | sb10md | Unimplemented |
| SB10PD | sb10 | sb10pd | Unimplemented |
| SB10QD | sb10 | sb10qd | Unimplemented |
| SB10RD | sb10 | sb10rd | Unimplemented |
| SB10SD | sb10 | sb10sd | Unimplemented |
| SB10TD | sb10 | sb10td | Unimplemented |
| SB10UD | sb10 | sb10ud | Unimplemented |
| SB10VD | sb10 | sb10vd | Unimplemented |
| SB10WD | sb10 | sb10wd | Unimplemented |
| SB10YD | sb10 | sb10yd | Unimplemented |
| SB10ZD | sb10 | sb10zd | Unimplemented |
| SB10ZP | sb10 | sb10zp | Unimplemented |
| SB16AD | sb16 | sb16ad | Unimplemented |
| SB16AY | sb16 | sb16ay | Unimplemented |
| SB16BD | sb16 | sb16bd | Unimplemented |
| SB16CD | sb16 | sb16cd | Unimplemented |
| SB16CY | sb16 | sb16cy | Unimplemented |
| SG02AD | sg02 | sg02ad | Unimplemented |
| SG02CV | sg02 | sg02cv | Unimplemented |
| SG02CW | sg02 | sg02cw | Unimplemented |
| SG02CX | sg02 | sg02cx | Unimplemented |
| SG02ND | sg02 | sg02nd | Unimplemented |
| SG03AD | sg03 | sg03ad | Unimplemented |
| SG03AX | sg03 | sg03ax | Unimplemented |
| SG03AY | sg03 | sg03ay | Unimplemented |
| SG03BD | sg03 | sg03bd | Unimplemented |
| SG03BR | sg03 | sg03br | Unimplemented |
| SG03BS | sg03 | sg03bs | Unimplemented |
| SG03BT | sg03 | sg03bt | Unimplemented |
| SG03BU | sg03 | sg03bu | Unimplemented |
| SG03BV | sg03 | sg03bv | Unimplemented |
| SG03BW | sg03 | sg03bw | Unimplemented |
| SG03BX | sg03 | sg03bx | Unimplemented |
| SG03BY | sg03 | sg03by | Unimplemented |
| SG03BZ | sg03 | sg03bz | Unimplemented |
| TB01ID | tb01 | tb01id | Unimplemented |
| TB01IZ | tb01 | tb01iz | Unimplemented |
| TB01KD | tb01 | tb01kd | Unimplemented |
| TB01KX | tb01 | tb01kx | Unimplemented |
| TB01LD | tb01 | tb01ld | Unimplemented |
| TB01MD | tb01 | tb01md | Implemented |
| TB01ND | tb01 | tb01nd | Unimplemented |
| TB01PD | tb01 | tb01pd | Unimplemented |
| TB01PX | tb01 | tb01px | Unimplemented |
| TB01TD | tb01 | tb01td | Unimplemented |
| TB01TY | tb01 | tb01ty | Unimplemented |
| TB01UD | tb01 | tb01ud | Unimplemented |
| TB01UX | tb01 | tb01ux | Unimplemented |
| TB01UY | tb01 | tb01uy | Unimplemented |
| TB01VD | tb01 | tb01vd | Unimplemented |
| TB01VY | tb01 | tb01vy | Unimplemented |
| TB01WD | tb01 | tb01wd | Unimplemented |
| TB01WX | tb01 | tb01wx | Unimplemented |
| TB01XD | tb01 | tb01xd | Unimplemented |
| TB01XZ | tb01 | tb01xz | Unimplemented |
| TB01YD | tb01 | tb01yd | Unimplemented |
| TB01ZD | tb01 | tb01zd | Unimplemented |
| TB03AD | tb03 | tb03ad | Unimplemented |
| TB03AY | tb03 | tb03ay | Unimplemented |
| TB04AD | tb04 | tb04ad | Unimplemented |
| TB04AY | tb04 | tb04ay | Unimplemented |
| TB04BD | tb04 | tb04bd | Unimplemented |
| TB04BV | tb04 | tb04bv | Unimplemented |
| TB04BW | tb04 | tb04bw | Unimplemented |
| TB04BX | tb04 | tb04bx | Unimplemented |
| TB04CD | tb04 | tb04cd | Unimplemented |
| TB05AD | tb05 | tb05ad | Unimplemented |
| TC01OD | tc01 | tc01od | Unimplemented |
| TC04AD | tc04 | tc04ad | Unimplemented |
| TC05AD | tc05 | tc05ad | Unimplemented |
| TD03AD | td03 | td03ad | Unimplemented |
| TD03AY | td03 | td03ay | Unimplemented |
| TD04AD | td04 | td04ad | Unimplemented |
| TD05AD | td05 | td05ad | Unimplemented |
| TF01MD | tf01 | tf01md | Unimplemented |
| TF01MX | tf01 | tf01mx | Unimplemented |
| TF01MY | tf01 | tf01my | Unimplemented |
| TF01ND | tf01 | tf01nd | Unimplemented |
| TF01OD | tf01 | tf01od | Unimplemented |
| TF01PD | tf01 | tf01pd | Unimplemented |
| TF01QD | tf01 | tf01qd | Unimplemented |
| TF01RD | tf01 | tf01rd | Unimplemented |
| TG01AD | tg01 | tg01ad | Unimplemented |
| TG01AZ | tg01 | tg01az | Unimplemented |
| TG01BD | tg01 | tg01bd | Unimplemented |
| TG01CD | tg01 | tg01cd | Unimplemented |
| TG01DD | tg01 | tg01dd | Unimplemented |
| TG01ED | tg01 | tg01ed | Unimplemented |
| TG01FD | tg01 | tg01fd | Unimplemented |
| TG01FZ | tg01 | tg01fz | Unimplemented |
| TG01GD | tg01 | tg01gd | Unimplemented |
| TG01HD | tg01 | tg01hd | Unimplemented |
| TG01HU | tg01 | tg01hu | Unimplemented |
| TG01HX | tg01 | tg01hx | Unimplemented |
| TG01HY | tg01 | tg01hy | Unimplemented |
| TG01ID | tg01 | tg01id | Unimplemented |
| TG01JD | tg01 | tg01jd | Unimplemented |
| TG01JY | tg01 | tg01jy | Unimplemented |
| TG01KD | tg01 | tg01kd | Unimplemented |
| TG01KZ | tg01 | tg01kz | Unimplemented |
| TG01LD | tg01 | tg01ld | Unimplemented |
| TG01LY | tg01 | tg01ly | Unimplemented |
| TG01MD | tg01 | tg01md | Unimplemented |
| TG01ND | tg01 | tg01nd | Unimplemented |
| TG01NX | tg01 | tg01nx | Unimplemented |
| TG01OA | tg01 | tg01oa | Unimplemented |
| TG01OB | tg01 | tg01ob | Unimplemented |
| TG01OD | tg01 | tg01od | Unimplemented |
| TG01OZ | tg01 | tg01oz | Unimplemented |
| TG01PD | tg01 | tg01pd | Unimplemented |
| TG01QD | tg01 | tg01qd | Unimplemented |
| TG01WD | tg01 | tg01wd | Unimplemented |
| UD01BD | ud01 | ud01bd | Unimplemented |
| UD01CD | ud01 | ud01cd | Unimplemented |
| UD01DD | ud01 | ud01dd | Unimplemented |
| UD01MD | ud01 | ud01md | Unimplemented |
| UD01MZ | ud01 | ud01mz | Unimplemented |
| UD01ND | ud01 | ud01nd | Unimplemented |
| UE01MD | ue01 | ue01md | Unimplemented |
| ZGEGS | zgeg | zgegs | Unimplemented |
| ZGEGV | zgeg | zgegv | Unimplemented |
| ZLATZM | zlat | zlatzm | Unimplemented |
