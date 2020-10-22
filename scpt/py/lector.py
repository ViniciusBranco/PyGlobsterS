import os
import subprocess as sp
import numpy as np


NGCs = ['NGC0104','NGC5927','NGC1904'] # NGC = 'NGC0104'
abunds = ['Alpha','AlphaCNONaHe']
schwag = ['sch', 'wag']
# sw = 'wag'
# for sw in schwag:
for NGC in NGCs:
# 	for abund in abunds:

		# specpath, specfile = ('libraries/synthetic/gc/'+NGC+'/', 'reb'+sw+'synthe_'+NGC+'_'+abund+'.asc')
# specpath, specfile = ('libraries/synthetic/gc/'+NGC+'/', 'synthe_'+NGC+'_'+abund+'.asc')

	specpath = 'libraries/synthetic/gc/'+NGC+'/specs/reb/'
	speclist = np.loadtxt(specpath+'speclist.txt', dtype=str)

	for specfile in speclist:

		os.chdir('shscripts/')
		lines = [line.rstrip() for line in open('lector.sh')]
		lines[3] = 'ln -s ../'+specpath+specfile
		lines[7] = 's'
		lines[8] = specfile

		lines[13] = 'mv output.txt ../'+specpath
		lines[14] = 'mv '+specfile+'_INDICES ../'+specpath+'LECTOR'
		lines[15] = 'mv '+specfile+'_LINE ../'+specpath+'LECTOR'
		lines[16] = 'mv '+specfile+'_ROSE ../'+specpath+'LECTOR'

		np.savetxt('lector.sh', lines, fmt='%s')					# save the edited sh script file
		sp.call('chmod +x lector.sh', shell=True)						# make sh script executable
		sp.call('./lector.sh', shell=True)								# Run sh script
		print(f'\nOutput saved at: {specpath}.')
		os.chdir('../')