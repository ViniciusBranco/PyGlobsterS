import sys
import numpy as np
import pyscripts.pyatlas as patl
import pyscripts.pysynthe as psyn
import pyscripts.rebin as rebin
import pyscripts.mad as mad
import pyscripts.gf as gf


# update: 20.10.2020 2

# ======================================================================================================
# To run PyGlobsterS, use the following command on terminal: $ python3 pyglobsters.py x x x x x
# where 'x' is a boolean control flag to execute the internal methods: PyAtlas, PySynthe, rebin, MAD, GF
# 0: do not execute; 1: execute
# ======================================================================================================


# FUTURE UPDATES =======================================================================================
# - use the following command to open a new gnome-terminal and make pseudo-parallel coding
# import os; os.system("gnome-terminal -e 'bash -c \"<command>; exec bash\"'")
# whereas e.g. <command> = 'ipython' or 'sudo apt update' or even 'python3 pyglobsters.py x x x x x' !!!
# ======================================================================================================


##### ========================= PyGlobsterS init ================================================= #####


class bcolors:
    HEADER = '\033[95m'
    OKBLUE = '\033[94m'
    OKGREEN = '\033[92m'
    WARNING = '\033[93m'
    FAIL = '\033[91m'
    ENDC = '\033[0m'
    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'

n = 6
if len(sys.argv) == n:
    patl_bool, psyn_bool, rebin_bool, mad_bool, gf_bool = sys.argv[1:]
elif len(sys.argv) > n:
    print(f'\n {bcolors.FAIL}PyGlobsterS ERROR1:{bcolors.ENDC} too many values to unpack (got {len(sys.argv[1:])}, expected {n} booleans: PyAtlas, PySynthe, rebin, MAD, GF).\n')
    sys.exit()
elif len(sys.argv) < n:
    print(f'\n {bcolors.FAIL}PyGlobsterS ERROR2:{bcolors.ENDC} not enough values to unpack (got {len(sys.argv[1:])}, expected {n} booleans: PyAtlas, PySynthe, rebin, MAD, GF).\n')
    sys.exit()

if patl_bool.isnumeric() & psyn_bool.isnumeric() & rebin_bool.isnumeric() & mad_bool.isnumeric() & gf_bool.isnumeric():
	patl_bool, psyn_bool, rebin_bool, mad_bool, gf_bool = int(patl_bool), int(psyn_bool), int(rebin_bool), int(mad_bool), int(gf_bool)

	if (patl_bool > 1) | (psyn_bool > 1) | (rebin_bool > 1) | (mad_bool > 1) | (gf_bool > 1):
		print(f'\n {bcolors.FAIL}PyGlobsterS ERROR3:{bcolors.ENDC} one or more values are not boolean (0 or 1).\n')
		sys.exit()
	else:
		print('\n##### ========== PyGlobsterS ========== #####')
else:
	print(f'\n {bcolors.FAIL}PyGlobsterS ERROR4:{bcolors.ENDC} one or more values are not digit.\n')
	sys.exit()


##### ========================= Entry parameters for PyATLAS ===================================== #####


if bool(patl_bool):

	typeobj			= 'gc'					# type of the objetc to model: 'gc' or 'stars'
	obj				= 'NGC1904'				# name of the object to model: 'sun', 'arcturus', '47tuc', 'NGCxxxx'
	# NGC0104, NGC1904, NGC5927

	refmodelpath 	= 'models/stars/sun/'			# reference model path
	# refmodels		= [line.rstrip() for line in open(refmodelpath+'refmodellist.txt')]
	convnotmodels 	= [line.rstrip() for line in open('atlas12/ifconv/convnotmodels2.txt')]

	ifconvlistname 	= 'ifconvlist.txt'			# list of 'ifconv' model files for convergence test
	ifconvoutput 	= 'convmodels3.txt'			# ifconv output file name for the convegent models
	ifconvnotoutput	= 'convnotmodels3.txt'		# ifconv output file name for the non-convergent models

	abundpath 		= 'abunds/'+obj+'/'
	abundtags 		= ['AlphaCNONaHe', 'Alpha']	# keep abundtag consistent with abundfiles
	abundlist 		= abundpath+'abundanceslist.txt'
	abundfiles 		= [line.rstrip() for line in open(abundlist)]

	modscriptpath, modscriptfile 	= ('shscripts/','atlas12.sh')	# ATLAS12 script path and script file

	outputmodelpath = 'models/'+typeobj+'/'+obj+'/'; outputrefname 	= obj 	# Output Model path and file name

	Teff 			= 3800
	Logg 			= 1.0
	vTurb 			= 2.0						# Turbulent Velocity in km/s

	idmod 			= 1							# first model ID (number)
	iterations 		= 45						# number of iterations Atlas will do to converge
	nlayers			= 72						# numer of layers to calculate. Max 72. Does not make sense much less (ref atlas cookbook)
	tRoss 			= -6.875					# rosseland optical depth: 10**(tRoss); max tRoss = 2, 10**(2) = 100.
	mixlength		= 1.25						# mixing length parameter
	overshooting 	= 0							# 0 for OFF, 1 for ON
	trossmin, trossmax = (-1, 5)				# T_Ross convergence limits

	metalvar 		= {'C': -0.3, 'N': 1.2, 'O':-0.45, 'Na': 0.6, 'Mg': 0, 'Al': 0}							# chemical variation in dex
	metalidx 		= {'6': metalvar['C'], '7': metalvar['N'], '8': metalvar['O'], '11': metalvar['Na'],	# chemical variation idx
						'12': metalvar['Mg'], '13': metalvar['Al']}	


##### ========================= Entry parameters for PySYNTHE ==================================== #####


if bool(psyn_bool) == 1:

	typeobj			= 'gc'					# type of the objetc to model: 'gc' or 'stars'
	obj				= 'NGC5927'				# name of the object to model: 'sun', 'arcturus', '47tuc', 'NGCxxxx'
	# NGC0104, NGC1904, NGC5927

	shscriptpath 	= 'shscripts/'						# sh script path
	shscriptfile 	= 'synthe.sh'						# sh script name

	modelpath 		= 'models/'+typeobj+'/'+obj+'/'		# atmospheric model path
	modellist 		= modelpath+'modellist.txt'
	modelfiles 		= [line.rstrip() for line in open(modellist)]

	abundpath 		= 'abunds/'+obj+'/'	
	abundfiles 		= [line.rstrip() for line in open(abundpath+'abundanceslist.txt')]
	abundtag 		= ['Alpha', 'AlphaCNONaHe']	# keep abundtag consistent with abundfiles

	gfpath			= 'gfs/'
	gffiles 		= ['gfallvinicius20.dat']	# ['gfallcoelho14.dat', 'gfallcastelli16.dat', 'gfallkurucz17.dat', 'gfallvinicius20.dat']

	wb	 			= '250.0'					# Wavelength Begin - examples of format: 515.0; 1042.0 nm
	we	 			= '1050.0'					# Wavelength End - examples of format: 515.0; 1042.0 nm

	sr 				= '300000'					# Spectrum Resolution - examples of format: 100000; 50000
	broadening		= '150000' 					# broadening function - examples of format: 100000; 50000 - using 'GAUSSIAN  '<value>'    RESOLUTION'	#
	header 			= 'flux'					# Choose among flux or intensi .header for SYNTHE. If flux, PySYNTHE will NOT evaluate rotate.exe.
	AV 				= 'AIR'						# AIR or VAC - if lambda <= 200 nm: VAC; else if lambda >= 200 nm AIR, usually.

	outputpath 	= 'libraries/synthetic/'+typeobj+'/'+obj+'/'	# Output path




##### ========================= Entry parameters for rebin.py ==================================== #####


if bool(rebin_bool) == 1:

	rebinobs_path 	= 'libraries/observed/stars/arcturus/'	# observed spectrum path of reference to rebin
	# rebinobs_file 	= 'sun_300_1000_kurucz.asc'			# observed spectrum file of reference to rebin
	rebinobs_file 	= 'arc_372_930_noao.asc'				# observed spectrum file of reference to rebin
	rebinsyn_path 	= outputpath							# output synthetic spectra to rebin


##### ========================= Entry parameters for mad.py ====================================== #####


if bool(mad_bool) == 1:

	madobs_path 	= rebinobs_path			# observed spectrum path of reference to evaluate the MAD
	madobs_file 	= rebinobs_file			# observed spectrum file of reference to evaluate the MAD
	madsyn_path 	= 'libraries/rebinned/'	# output synthetic spectra to evaluate the MAD
	madtolerance	= 0.99					# only 1% of telluric impact is allowed: madtolerance - 0.99
	maddelta		= 0.02					# delta = 0.02 nm or delta = 0.2 A


##### ========================= Entry parameters for gf.py ======================================= #####


if bool(gf_bool) == 1:

	gfspath			= gfpath
	consolidated 	= 'misc/mads/consolidated_teste.asc'


##### ========================= Modules Execution ================================================ #####


## ===== PyATLAS ===== ##

if bool(patl_bool):
	print('\n\n## ===== PyATLAS ===== ##\n')

	aux = 0
	for convnotmodel in convnotmodels: # convnotmodel = convnotmodels[0]
		refmodel = 'sunasp09t5777logg44377v100.dat' # refmodels[aux]

		FeH = float(convnotmodel.split('_')[4][1:])*(-1)

		if FeH == -0.47:
			obj = 'NGC5927'
		elif FeH == -1.579:
			obj = 'NGC1904'
		elif FeH == -0.768:
			obj = 'NGC0104'	

		Teff, Logg 	= (int(convnotmodel.split('_')[2][1:]), float(convnotmodel.split('_')[3][1:]))
		alphaenh 	= float(convnotmodel.split('_')[5][1:])
		enhHe 		= float(convnotmodel.split('_')[6][1:])

		abundtag 		= convnotmodel.split('_')[-1][:-4]
		abundpath 		= 'abunds/'+obj+'/'
		abundlist 		= abundpath+'abundanceslist.txt'
		abundfiles 		= [line.rstrip() for line in open(abundlist)]

		if abundtag == 'AlphaCNONaHe':
			abundfile = abundfiles[0]
		elif abundtag == 'Alpha':
			abundfile = abundfiles[1]

		outputmodelpath = 'models/'+typeobj+'/'+obj+'/'; outputrefname 	= obj 	# Output Model path and file name

		patl.pyatlas(refmodelpath, refmodel, abundpath, abundfile, abundtag, alphaenh, modscriptpath, modscriptfile, idmod, 
			iterations, nlayers, enhHe, FeH, Teff, Logg, vTurb, tRoss, mixlength, overshooting, outputmodelpath, outputrefname)	
		idmod += 1; aux += 1;

	# for refmodel in refmodels: # refmodel = refmodels[0]

	# 	if obj == 'NGC5927':
	# 		FeH = -0.47
	# 	elif obj == 'NGC1904':
	# 		FeH = -1.579
	# 	elif obj == 'NGC0104':
	# 		FeH = -0.768	

	# 	Teff = int(refmodel.split('_')[0][1:])
	# 	Logg = float(refmodel.split('_')[1][1:])

	# 	for abundfile in abundfiles: # abundfile = abundfiles[0]

	# 		abundtag 	= abundfile.split('_')[3]
	# 		alphaenh 	= float(abundfile.split('_')[5][1:])
	# 		enhHe 		= float(abundfile.split('_')[6][1:-4])

	# 		patl.pyatlas(refmodelpath, refmodel, abundpath, abundfile, abundtag, alphaenh, modscriptpath, modscriptfile, idmod, 
	# 			iterations, nlayers, enhHe, FeH, Teff, Logg, vTurb, tRoss, mixlength, overshooting, outputmodelpath, outputrefname)	
	# 		idmod += 1

	patl.ifconv(trossmin, trossmax, ifconvlistname, ifconvoutput, ifconvnotoutput)

## ===== PySYNTHE ===== ##

if bool(psyn_bool) == 1:
	print('\n\n## ===== PySYNTHE ===== ##\n')
	if shscriptfile == 'synthemolecules.sh':
		for molecule in moleculelist:
			for modelfile in modelfiles:
				for abundfile in abundfiles:
					psyn.pysynthe(shscriptpath, shscriptfile, modelpath, modelfile, abundpath, abundfile, gfpath, 0, temp, wb, we, sr, broadening, 
						header, AV, outputpath, molecule)
	else:					
		for gffile in gffiles: # gffile = gffiles[0]
			for modelfile in modelfiles: # modelfile = modelfiles[0]
				abundtag 	= modelfile.split('_')[-1][:-4]
				if abundtag == 'AlphaCNONaHe':
					abundfile = abundfiles[0]
				elif abundtag == 'Alpha':
					abundfile = abundfiles[1]
				temp = int(modelfile.split('_')[1][1:])
				psyn.pysynthe(shscriptpath, shscriptfile, modelpath, modelfile, abundpath, abundfile, gfpath, gffile, temp, wb, we, sr, broadening, 
					header, AV, outputpath, 0)
		psyn.comment_specs(outputpath)

## ===== Rebin ===== ##

if bool(rebin_bool) == 1:
	print('\n\n## ===== Rebin ===== ##\n')

	outputpath = rebin.rebin(rebinobs_path, rebinobs_file, rebinsyn_path)

## ===== MAD ===== ##

if bool(mad_bool) == 1:
	print('\n\n## ===== MAD ===== ##\n')

	mad.mad(madobs_path, madobs_file, madsyn_path, madtolerance, maddelta)

## ===== GF ===== ##

if bool(gf_bool) == 1:
	print('\n\n## ===== GF ===== ##\n')

	gf.gf(sorted(gffiles), gfspath, consolidated)

print('\n##### ========== Finished ========== #####\n')
sys.exit()
