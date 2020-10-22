import os
import subprocess as sp
import numpy as np
import sys
import shutil
'''

### --- PySYNTHE 1.1.1.3  --- ###

(last version, adition of functions, changes on functions, small updates). LAST UPDATE: September, 06, 2018.. Oficial date of creation: August 21, 2018.

Basic instructions:

	1. If you want to run one script file on SYNTHE, set the single script file name in 'shscript' then leave single_script() uncommented and multi_scripts() commented.
	2. Otherwise, leave single_script() commented and multi_scripts() uncommented. Check the documentation to learn what all the methods do.
	3. For any run, you must provide the following list of inputs:
		3.1 Model	(usualy a .dat file, from 'models' folder).
		3.2 Atomic lines (usualy an ascii file, from 'GFs' folder).	
		3.3 Wavelength interval and resolution (initial and final values of wavelength of the spectra you want SYNTHE to calculate - MUST BE in nanometers)
		3.4 Number of wavelength interval divisions (how many parts your wavelength range will be cut - must be n >= 1)

	4. If you want to run a full SYNTHE script (Molecular and Atomic lines):
		4.1 leave the following command lines in the multi_run() and single_run() methods UNCOMMENTED:

			lines[27] = 'ln -s '+atomiclines+' fort.11'
			scptname = scpt[:-3]+'_'+model[:-4]+'_'+wlbeg+'_'+wlend+'_'+atomiclines[:-4]+'.sh'

		4.2 and leave COMMENTED the following command line on both methods:

			scptname = scpt[:-3]+'_'+model[:-4]+'_'+str(idx)+'_'+wlbeg+'_'+wlend+'.sh'	


	5. If you want to evaluate ATOMIC LINES ONLY (same instruction as that for the full SYNTHE script):
		5.1 leave the following command lines in the multi_run() and single_run() methods UNCOMMENTED:

			lines[27] = 'ln -s '+atomiclines+' fort.11'
			scptname = scpt[:-3]+'_'+model[:-4]+'_'+wlbeg+'_'+wlend+'_'+atomiclines[:-4]+'.sh'

		5.2 and leave COMMENTED the following command line on both methods:

			scptname = scpt[:-3]+'_'+model[:-4]+'_'+str(idx)+'_'+wlbeg+'_'+wlend+'.sh'

	6. If you want to evaluate MOLECULAR LISTS:
		6.1 Choose among single and multi_run method and leave the following command lines in the multi_run() and single_run() methods COMMENTED:

			lines[27] = 'ln -s '+atomiclines+' fort.11'
			scptname = scpt[:-3]+'_'+model[:-4]+'_'+wlbeg+'_'+wlend+'_'+atomiclines[:-4]+'.sh'

		6.2 and leave UNCOMMENTED the following command line on both methods:

			scptname = scpt[:-3]+'_'+model[:-4]+'_'+str(idx)+'_'+wlbeg+'_'+wlend+'.sh'

Simple order of run on ipython:

	1. Folder: Main. Run: pysynthe.py
	2. Folder: Plotz. Run: wrapspecs.py
	3. Folder: Plotz/Specs. Run: convolve.py
	4. Folder: Plotz. Run: plot.py or plot_mols.py


Lastest updates:
06/09/2018: creation of the Basic Instructions, list of ToDo;
			wraspecs.py should be run before convolve.py;
			wrap_spec() and wrap_mol() methods on wrapspecs.py should be chosen prior the run.
05/09/2018: single_script() and multi_script methods() added (see documentation for further information);
			plot_mols.py now plots each molecular spectrum and all plots in the same plot, i.e., each_mol() and all_mol() methods added there.

ToDo:
- Make wrapspecs.py work before convolve.py so that one may delete fragmented spec files and work with just one. OK
- Describe the code: documentation.
- On multi_run() and single_run(), set a 'if' statement for the use of atomic or molecular lines so that the user does not have to worry about commenting the script.
- Find a way to converge all py files that compose this PySYNTHE solution for the SYNTHE Suite ir order to avoid multiple opened terminals.
- There are too manu "plot something.py" files. Merge all!
- Make the code create some folders in order to acomodate the relative scripts, specs and plots of the same run.

~ V. Branco, 2018. vinicius.branco.silva@usp.br.

'''
### ------------- DO NOT CHANGE THE CODE BELOW UNLESS YOU KNOW WHAT YOU ARE DOING ------------- ###


def comment_specs(outputpath, filename):

	print ("Commenting spectra header for plotting...")

	# os.chdir(outputpath)

	# sp.call("ls synthe*.asc > speclist.txt", shell=True)
	# spec_list = [line.rstrip() for line in open('speclist.txt')]
	# for spec in spec_list:
	lines = [line.rstrip() for line in open(outputpath+filename)]
	lines[0] = '#'+lines[0]
	lines[1] = '#'+lines[1]
	np.savetxt(outputpath+spec, lines, fmt='%s')

	sp.call("ls synthe*.asc > speclist.txt", shell=True)
	# os.chdir('../')

def pysynthe(shscriptpath, shscriptfile, modelpath, modelfile, abundpath, abundfile, gfpath, gffile, temp, wb, we, sr, broadening, header, AV, outputpath, molecule):

	rot = 'Rotate On'; header = header + '.header'
	if (header != 'intensi.header'):			# Se header diff intensi.header (ou seja, se for fluxo), o bloco rotate eh comentado
		rot = 'Rotate Off'

	TiO = 'TiO On'
	if (temp >= 4500):
		TiO = 'TiO Off'

	lines = [line.rstrip() for line in open(modelpath+modelfile)]; title = lines[1]
	print('=============================================================================')
	print(f'Calculating {title[6:]}\n{TiO} {rot} with script {shscriptfile}')	
	print(f'Reference Model: {modelfile}')
	print(f'Reference Abundances: {abundfile}')
	print(f'Reference Molecule Line List: {molecule}') if molecule != 0 else print(f'Reference Atomic Line List: {gffile}')
	print(f'Wavelength Interval: {wb}-{we} (nm) in the {AV}')
	print(f'Sample Resolution: {sr}')
	print(f'Broadening Resolution: {broadening}\n')

	lines = [line.rstrip() for line in open(shscriptpath+shscriptfile)]
	lines[5] = 'ln -s ../../'+modelpath+'* .'	
	lines[7] = 'ln -s ../../'+abundpath+'* .'

	lines[10] = 'cat ../'+header+' > synthe.mod'
	lines[11] = 'head -n 3 '+modelfile+' >> synthe.mod'
	lines[12] = 'cat ../conv.txt '+abundfile+' >> synthe.mod'
	lines[13] = 'tail -n 75 '+modelfile+' >> synthe.mod'
	lines[24] = AV+'       '+wb+'     '+we+'     '+sr+'.   0.     0         -10 .001         0   00'		# SYNBEG parameters

	if shscriptfile == 'synthemolecules.sh':
		lines[30] = 'ln -s '+molecule+' fort.11'
		lines[58] = 'GAUSSIAN  '+broadening+'.    RESOLUTION' 		# Broadening
		scpt = shscriptfile[:-3]+'_'+modelfile[:-4]+'_'+wb+'_'+we+'_'+molecule.split('/')[-1][:-4]+'_'+AV+'.sh'
		lines[68] = 'mv spec.asc ../../'+outputpath+scpt[:-3]+'.asc'		

	elif shscriptfile == 'syntheatoms.sh':
		lines[29] = 'ln -s '+gffile+' fort.11'														
		lines[57] = 'GAUSSIAN  '+broadening+'.    RESOLUTION' 		# Broadening
		scpt = shscriptfile[:-3]+'_'+modelfile[:-4]+'_'+wb+'_'+we+'_'+gffile[:-4]+'_'+AV+'.sh'
		lines[67] = 'mv spec.asc ../../'+outputpath+scpt[:-3]+'.asc'

	else:
		lines[29] = 'ln -s '+gffile+' fort.11'														
		lines[237] = 'GAUSSIAN  '+broadening+'.    RESOLUTION' 		# Broadening
		scpt = shscriptfile[:-3]+'_'+modelfile[:-4]+'_'+wb+'_'+we+'_'+gffile[:-4]+'_'+AV+'.sh'
		lines[247] = 'mv spec.asc ../../'+outputpath+scpt[:-3]+'.asc'

	if (header != 'intensi.header'):						# Se header diff intensi.header (ou seja, se for fluxo), o bloco rotate eh comentado
		for i in range(1, len(lines[1:])):
			if 'rotate.exe' in lines[i]:
				for k in range(-2,5):
					lines[i+k] = '#'+lines[i+k]
				
	if ((temp >= 4500) & (shscriptfile == 'synthe.sh')):	# Se temperatura maior que 4500K, o bloco TiO eh comentado
		for i in range(1, len(lines[1:])):
			if 'tioschwenke' in lines[i]:
				for k in range(0,7):
					lines[i+k] = '#'+lines[i+k]

	os.mkdir(shscriptpath+'temp_'+scpt[:-3]) 		# create the folder where everything will happen
	os.chdir(shscriptpath+'temp_'+scpt[:-3])		# change env to the folder where everything will happen

	np.savetxt(scpt, lines, fmt='%s')				# save the edited sh script file
	sp.call('chmod +x '+scpt, shell=True)			# make sh script executable
	sp.call("./"+scpt, shell=True)					# Run SYNTHE script
	print(f'\nSaved at: {outputpath+scpt[:-3]}.asc\n')	

	os.chdir('../../')								# bring the env back to main app folder
	shutil.rmtree(shscriptpath+'temp_'+scpt[:-3])	# remove temporary folder

	lines = [line.rstrip() for line in open(outputpath+scpt[:-3]+'.asc')]
	lines[0] = '#'+lines[0]
	lines[1] = '#'+lines[1]
	np.savetxt(outputpath+scpt[:-3]+'.asc', lines, fmt='%s')

	sp.call("ls "+outputpath+"*.asc > "+outputpath+"speclist.txt", shell=True)
	# comment_specs(outputpath, scpt[:-3])
