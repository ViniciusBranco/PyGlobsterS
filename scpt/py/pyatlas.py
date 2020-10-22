import os
import subprocess as sp
import numpy as np
import shutil

def ifconv(trossmin, trossmax, ifconvlistname, ifconvoutput, ifconvnotoutput):
	# ifconvlistname 	= 'ifconvlist.txt'			# list of 'ifconv' model files for convergence test
	# ifconvoutput 	= 'convmodels.txt'			# ifconv output file name for the convegent models
	# ifconvnotoutput	= 'convnotmodels.txt'		# ifconv output file name for the non-convergent models
	# trossmin, trossmax = (-1, 5)				# T_Ross convergence limits	

	os.chdir('atlas12/ifconv/')
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
	os.chdir('../../')

def pyatlas(refmodelpath, refmodel, abundpath, abundfile, abundtag, alphaenh, modscriptpath, modscriptfile, idmod, 
	iterations, nlayers, enhHe, FeH, Teff, Logg, vTurb, tRoss, mixlength, overshooting, outputmodelpath, outputrefname):

	if (FeH >= 0):
		mpfeh = 'p'
	else:
		mpfeh = 'm'

	Title = 'TITLE ATLAS12 #'+str(idmod)+' '+abundtag+' Teff '+str(Teff)+' Logg '+str(Logg)+' FeH '+str(FeH)+' Y '+str(enhHe)+' vTurb '+str(vTurb)
	savefilename = outputrefname+'_t'+str(Teff)+'_g'+str(Logg)+'_'+mpfeh+str(np.abs(FeH))+'_p'+str(alphaenh)+'_y'+str(enhHe)+'_'+abundtag+'.mod'
	print('===============================================================================')
	print(f'Modeling {Title[6:]}')
	print(f'Reference Model: {refmodel}')
	print(f'Reference Abundances: {abundfile}')
	print(f'Output file: {savefilename}\n')

	abundlines 			= [line.rstrip() for line in open(abundpath+abundfile)]
	atlaslines 			= [line.rstrip() for line in open(modscriptpath+modscriptfile)]

	atlaslines[4] 		= 'ln -s ../../'+refmodelpath+refmodel+' fort.3'				# changes the reference model symbolic link
	atlaslines[21] 		= 'CONVECTION OVER '+str(mixlength)+' '+str(overshooting)+' 36'
	atlaslines[33] 		= Title
	atlaslines[36] 		= 'CONVECTION OVER '+str(mixlength)+' '+str(overshooting)+' 36'
	atlaslines[37] 		= 'ITERATIONS '+str(iterations)
	atlaslines[40] 		= 'SCALE MODEL '+str(nlayers)+' '+str(tRoss)+' 0.125 '+str(Teff)+'. '+str(Logg)
	atlaslines[41] 		= 'VTURB '+str(vTurb)+'E+05'
	atlaslines[-3] 		= 'awk "NR <= 22 || NR >= 45" fort.7 > '+savefilename
	
	atlaslines 			= atlaslines[:23]+abundlines+atlaslines[24:42]+abundlines+atlaslines[43:]  # add abundance change lines

	tempdirectory = modscriptpath+'temp_'+savefilename[:-3]+'_'+modscriptfile[:-3]+'_'+abundtag	# set temporary dictory path and folder name

	os.mkdir(tempdirectory) # create the folder where everything will happen
	os.chdir(tempdirectory) # change env to the folder where everything will happen

	np.savetxt('temp_'+modscriptfile, atlaslines, fmt='%s')		# save temporary sh script
	sp.call('chmod +x temp_'+modscriptfile, shell=True)			# make temporary sh script executable
	sp.call('./temp_'+modscriptfile, shell=True)			 	# run temporary ATLAS12 sh script

	print(f'Saved at: {outputmodelpath+savefilename}\n')
	
	sp.call('cp run2.log ../../atlas12/ifconv/'+savefilename[:-4]+'.log', shell=True)	# copy run2.log from model calculation to ifconv folder for convergence test
	sp.call('cp '+savefilename+' ../../'+outputmodelpath, shell=True)					# copy output model to the model folder
	os.chdir('../../')																	# bring the env back to main app folder
	sp.call('./atlas12/ifconv.pl atlas12/ifconv/'+savefilename[:-4]+'.log > atlas12/ifconv/ifconv_'+savefilename[:-4]+'.log', shell=True)	# Run convergence test
	sp.call('rm atlas12/ifconv/'+savefilename[:-4]+'.log', shell=True)					# remove run2.log file from the model calculation

	shutil.rmtree(tempdirectory)	# remove temporary folder

	