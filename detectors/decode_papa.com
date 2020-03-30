#!/bin/csh
unalias rm
unalias mv
unalias cp
unalias lpr
alias rm rm -f
alias lpr lpr -Pdec -h 
#
goto skip1 
# Nota: 3800 frames with 100 photons and about 20H cpu, and 380 000 photons.
# so with 200 photons and 30mn cpu, 100 frames and 20 000 photons.  
 runs decode_papa 12 g_sex_4sec_550_35nm.dat 150,4000
 mv long.bdf long_gs.bdf
 mv modsq.bdf modsq_gs.bdf
 mv bisp1.bdf bisp1_gs.bdf
 runs pstcopy0 long_gs /AUTOSCALE/ID=long_gs_25-07-91
 lpr pstcopy.tmp
 goto end
skip0:
# runs pstcopy0 modsq_gs.bdf /BLACK=0.002/WHITE=0./ID=modsq_gs_24-07-91
 runs pstcopy0 modsq_gs.bdf /AUTOSCALE/ID=modsq_gs_25-07-91
 lpr pstcopy.tmp
 goto end
skip1:
 runs inv_bispec 1 12 modsq_gs bisp1_gs
skip2:
 runs fft2 -1 ref n imf testr testi 
 mv testr.bdf testr_gs.bdf
 mv testi.bdf testi_gs.bdf
 runs pstcopy0 testr_gs ID=testr_gs_25-07-91
 lpr pstcopy.tmp
goto skip4
 runs positiv 12,10 testr_gs testi_gs tr ti 
 runs pstcopy0 tr
 lpr pstcopy.tmp 
skip4:
 cp modsq_gs.bdf modsq.bdf
 rm profile1.log
 rm modsq.pro
 runs profile1 N modsq_centred.par
 $EXEC/plot1.exe < plot1.com
end:
