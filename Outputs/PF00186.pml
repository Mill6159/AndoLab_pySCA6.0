delete all
load ../Inputs/1RX2.pdb, main
hide all
bg_color white
show cartoon, (chain A)
color white

set_color color1, [0.500,1.000,0.000]
create sector1, (resi 133,55,27,32,25,107,23,153,39,121,18,90,13,38,63,11,44,51,24,46,7,125,111,113,42,49,94,6,22,15,54,122,92,57,96,95,35,43,31,61,50,21,14,100,40,52,47,59,53) & (chain A)
color color1, sector1
show spheres, sector1
show surface, sector1

set_color color_ic1, [1.000,0.000,0.000]
create ic_1, (resi 133,55,27,32,25,107,23,153,39,18,90,13,38,63) & (chain A)
color color_ic1, ic_1
show spheres, ic_1
show surface, ic_1

set_color color_ic2, [0.500,1.000,0.000]
create ic_2, (resi 42,57,31,35,113,43,122,96,54,95,46,44,15,14,94,49,7,59,61) & (chain A)
color color_ic2, ic_2
show spheres, ic_2
show surface, ic_2

set_color color_ic3, [0.000,1.000,1.000]
create ic_3, (resi 22,21,121,52,24) & (chain A)
color color_ic3, ic_3
show spheres, ic_3
show surface, ic_3

set_color color_ic4, [0.500,0.000,1.000]
create ic_4, (resi 47,125,53,40,92,100,111,51,6,11,50) & (chain A)
color color_ic4, ic_4
show spheres, ic_4
show surface, ic_4

zoom
set transparency, 0.4
ray
png PF00186
