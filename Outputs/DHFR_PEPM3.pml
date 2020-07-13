delete all
load ../Inputs/1RX2.pdb, main
hide all
bg_color white
show cartoon, (chain A)
color white

set_color color1, [1.000,0.000,1.000]
create sector1, (resi 27,55,121,32,133,28,23,13,51,71,25,38,39,107,158,63,18,105,90,47,149,11,46,144,40,94,35,57,122,49,7,60,54,42,61,31,92,14,22,50,123,44,96,95,43,21,15,126,147,115,103,113,111,59,52,53,41,5,81,6,45,100,125) & (chain A)
color color1, sector1
show spheres, sector1
show surface, sector1

set_color color_ic1, [1.000,0.000,0.000]
create ic_1, (resi 27,55,121,32,133,28,23,13,51,71,25,38,39,107,158,63,18,105) & (chain A)
color color_ic1, ic_1
show spheres, ic_1
show surface, ic_1

set_color color_ic2, [1.000,1.000,0.000]
create ic_2, (resi 42,57,31,46,35,44,22,43,95,54,96,113,15,14,94,7,49,61) & (chain A)
color color_ic2, ic_2
show spheres, ic_2
show surface, ic_2

set_color color_ic3, [0.000,1.000,0.000]
create ic_3, (resi 47,59,53,40,52,100,50,81,103) & (chain A)
color color_ic3, ic_3
show spheres, ic_3
show surface, ic_3

set_color color_ic4, [0.000,1.000,1.000]
create ic_4, (resi 125,45,111,6,11,60,92,41,90,126) & (chain A)
color color_ic4, ic_4
show spheres, ic_4
show surface, ic_4

set_color color_ic5, [0.000,0.000,1.000]
create ic_5, (resi 122,21,115,123,5,147) & (chain A)
color color_ic5, ic_5
show spheres, ic_5
show surface, ic_5

set_color color_ic6, [1.000,0.000,1.000]
create ic_6, (resi 144,149) & (chain A)
color color_ic6, ic_6
show spheres, ic_6
show surface, ic_6

zoom
set transparency, 0.4
ray
png DHFR_PEPM3
