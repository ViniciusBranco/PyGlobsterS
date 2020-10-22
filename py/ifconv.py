import os
import subprocess as sp
import numpy as np
import shutil

trossmin, trossmax = (-1, 5)				# T_Ross convergence limits
ifconvlistname 	= 'ifconvlist.txt'			# list of 'ifconv' model files for convergence test
ifconvoutput 	= 'ifconvoutput.txt'		# ifconv output file name for the convegent models
ifconvnotoutput	= 'ifconvnotoutput.txt'		# ifconv output file name for the non-convergent models

os.chdir('atlas12/ifconv/')

### ===============================================================================
sp.call("ls *.log > "+ifconvlistname, shell=True)
ifconvlist = [line.rstrip() for line in open(ifconvlistname)]

notconv = []; conv = []; convnotlog = []
for ifconvfile in ifconvlist:
	convergence = [line.rstrip() for line in open(ifconvfile)]
	
	ifconv = 0; i = 0; j = 1; notconvlayers = []
	if 'NOT' in convergence[4]:
		while i == 0:				
			line = convergence[-j]
			if line is not '':
				k_ross = float(line.split()[2])
				k_ross10 = int(np.round(np.log10(float(k_ross)),0))

				if (k_ross10 > trossmin) & (k_ross10 < trossmax):					
					notconvlayers.append(f'[tROSSmin = {trossmin}] < [log10( K_ROSS = {k_ross} ) = {k_ross10}] < [tROSSmax = {trossmax}] at layer {convergence[-j].split()[0] }')
					ifconv = 1						
				j += 1

			else:
				i = 1

	if ifconv == 1:
		notconv.append(ifconvfile)
		convnotlog.append(ifconvfile)
		convnotlog.append(notconvlayers)
		convnotlog.append('')

	else:
		conv.append(ifconvfile)
		
np.savetxt('convnot.log', convnotlog, fmt='%s')
np.savetxt(ifconvnotoutput, notconv, fmt='%s')
np.savetxt(ifconvoutput, conv, fmt='%s')
### ===============================================================================

os.chdir('../../')

