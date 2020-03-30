#!/bin/csh
unalias rm
unalias mv
unalias cp
unalias lpr
alias rm rm -f
alias lpr lpr -Pdec -h 
#
goto skip4 
#
# runs decode_mama 12 m3.data 810,112,1,256,224 10,19000 ffaverag 
# runs decode_mama 12 m3.data 810,128,1,256,256 10,19000 ffaverag 

 mv long.bdf long_03.bdf
 mv modsq.bdf modsq_03.bdf
 mv bisp1.bdf bisp1_03.bdf
 runs pstcopy0 long_03
 lpr pstcopy.tmp
skip0:
 runs pstcopy0 modsq_03.bdf /BLACK=0.002/WHITE=0./ID=modsq_03_22-07-91
 lpr pstcopy.tmp
skip1:
 runs inv_bispec 1 12 modsq_03 bisp1_03
skip2:
 runs fft2 -1 ref n imf testr testi 
 mv testr.bdf testr_03.bdf
 mv testi.bdf testi_03.bdf
 runs pstcopy0 testr_03 ID=testr_03_22-07-91
 lpr pstcopy.tmp
goto skip4
 runs positiv 12,10 testr_03 testi_03 tr ti 
 runs pstcopy0 tr
 lpr pstcopy.tmp 
skip4:
 cp modsq_03.bdf modsq.bdf
 rm profile1.log
 rm modsq.pro
 runs profile1 N modsq_centred.par
 $EXEC/plot1.exe < plot1.com
end:
