import numpy as np 
import csv
# import pyphot


def loadfile(file, path=''): # funcao para abrir arquivo tabulado (ex: espectros)
	return np.array([list(i) for i in tuple(map(tuple, np.loadtxt(path+file).T))]) 


def loadlector(arxiv):
	lines = []
	with open(arxiv) as f:
		reader = csv.reader(f); data = list(reader)
		lines.append(data)
	return lines[0]


def getlicks(lines):
	lickname = []; lickvalue = []

	for line in lines: # 
		if len(line[0].split()) == 3:
			lickname.append(line[0].split()[0])
			lickvalue.append(line[0].split()[1])

	return lickname,lickvalue


def createDict(arg):
	name, value = (arg[0], arg[1])
	dictionary = {}
	for i in range(0,len(name)):
		dictionary[name[i]] = float(value[i])
	return dictionary


def calcMgFe(licks):
	Mgb, Fe5270, Fe5335 = (licks['Mgb'], licks['Fe5270'], licks['Fe5335'])
	x = Mgb * ((Fe5270 + Fe5335) /2)
	return np.round(np.float(np.sqrt( x )),3)


def calcdiff(bands, lick1, lick2):
	
	diff = []
	for band in bands: 
		try: 
			diff.append( np.round(lick2[band] - lick1[band],3) )
		except:
			diff.append(0)

	return createDict( (bands, diff) )



deltalickcoenames 	= ['Hdelta_F', 'Hgamma_F', 'CN1', 'CN2', 'Ca4227', 'G4300', 'H_beta', 'Mgb', 'Fe5270', 'Fe5335', 'Fe5406', 'Na_D'] # deltas Coelho et al 2011
deltalickcoevalue	= [  -0.068,    0.274    , 0.084, 0.087,   -0.651,  -0.572,  -0.032, -0.158,    0.004,    0.010,   -0.005,  1.346]

deltalickcoe =  createDict( (deltalickcoenames,deltalickcoevalue) )

lectorbandswag = ['Hdelta/4045','Hdelta/4063','SrII/4045','SrII/4063','Hgamma/Gband','Hgamma/4325',
				'4289/4271','4384/4352','p[Fe/H]','CaII','3888/3859','p4220/4208','FeI_HR','CaI_HR','HgammaHR',
				'Hg_sigma_130','Hg_sigma_125','Hg_sigma_200',
				'Hg_sigma_275','Hdelta_A','Hgamma_A','Hdelta_F','Hgamma_F',
				'CN1','CN2','Ca4227','G4300','Fe4383','H_beta','H_beta_o',
				'Fe5015','Mg1','Mg2','Mgb','Fe5270','Fe5335','Fe5406','Ca4455','Fe4531','Fe4668' ,'Fe5709' ,'Fe5782','Na_D','TiO_1','TiO_2','CaT*','PaT','CaT','[MgFe]']

lectorbandssch = ['Hdelta/4045','Hdelta/4063','SrII/4045','SrII/4063','Hgamma/Gband','Hgamma/4325',
				'4289/4271','4384/4352','p[Fe/H]','CaII','3888/3859','p4220/4208','FeI_HR','CaI_HR','HgammaHR',
				'Hg_sigma_130','Hg_sigma_125','Hg_sigma_200',
				'Hg_sigma_275','Hdelta_A','Hgamma_A','Hdelta_F','Hgamma_F',
				'CN1','CN2','Ca4227','G4300','Fe4383','H_beta','H_beta_o',
				'Fe5015','Mg1','Mg2','Mgb','Fe5270','Fe5335','Fe5406','Ca4455','Fe4531','Fe4668' ,'Fe5709' ,'Fe5782','Na_D','TiO_1','[MgFe]']


NGCs = ['NGC0104', 'NGC1904', 'NGC5927']

for NGC in NGCs: # NGC = 'NGC0104'

	lickspath = 'libraries/synthetic/gc/'+NGC+'/LECTOR/'

	licks1 = createDict(getlicks(loadlector(lickspath+'rebconvsch_synthe_'+NGC+'_Alpha.asc_INDICES'))); licks1['[MgFe]'] = calcMgFe(licks1)
	licks2  = createDict(getlicks(loadlector(lickspath+'rebconvsch_synthe_'+NGC+'_AlphaCNONaHe.asc_INDICES'))); licks2['[MgFe]'] = calcMgFe(licks2)

	lickw1 = createDict(getlicks(loadlector(lickspath+'rebreconvwag_synthe_'+NGC+'_Alpha.asc_INDICES'))); lickw1['[MgFe]'] = calcMgFe(lickw1)
	lickw2 = createDict(getlicks(loadlector(lickspath+'rebreconvwag_synthe_'+NGC+'_AlphaCNONaHe.asc_INDICES'))); lickw2['[MgFe]'] = calcMgFe(lickw2)

	licksch = createDict(getlicks(loadlector('libraries/observed/gc/schiavon/'+NGC+'/'+NGC+'_martins_schiavon.txt_INDICES'))); licksch['[MgFe]'] = calcMgFe(licksch)
	lickwag = createDict(getlicks(loadlector('libraries/observed/gc/waggs/norm_'+NGC+'/corr_waggs_norm_'+NGC+'.asc_INDICES'))); lickwag['[MgFe]'] = calcMgFe(lickwag)

	deltalicksch = calcdiff(lectorbandssch, licks1, licks2)
	deltalickwag = calcdiff(lectorbandswag, lickw1, lickw2)
	
	np.savetxt('libraries/synthetic/gc/'+NGC+'convdeltalicksch.txt', 
		np.column_stack(( lectorbandssch  , list(licksch.values()) , list(licks1.values()) , list(licks2.values()), list(deltalicksch.values()) )), 
		fmt='%13s', header = '       Band      Schiavon          Gen1          Gen2       Delta I (Gen2 - Gen1)'  )

	np.savetxt('libraries/synthetic/gc/'+NGC+'convdeltalickwag.txt', 
		np.column_stack(( lectorbandswag , list(lickwag.values()) , list(lickw1.values()) , list(lickw2.values()), list(deltalickwag.values()) )), 
		fmt='%13s', header = '       Band         WAGGS          Gen1          Gen2       Delta I (Gen2 - Gen1)'  )




	# ### PyPhot ###

	# lib = pyphot.LickLibrary()

	# pathwagg 			= 'libraries/observed/gc/waggs/norm_'+obj+'/'
	# filewagg 			= 'corr_waggs_norm_'+obj+'.asc'
	# wagwave, wagflux 	= loadfile(filewagg, pathwagg)

	# pathschi 			= 'libraries/observed/gc/schiavon/'+obj+'/'
	# fileschi 			= obj+'_martins_schiavon.txt'
	# schwave, schflux 	= loadfile(fileschi, pathschi)

	# pathsyn = 'libraries/synthetic/gc/'+obj+'/'

	# schpop1file = 'reb_convschiavon_synthe_'+obj+'_Alpha_250_1050_gfallvinicius.asc'			# Schiavon Alpha
	# schpop2file	= 'reb_convschiavon_synthe_'+obj+'_AlphaCNONaHe_250_1050_gfallvinicius.asc'		# Schiavon AlphaCNONaHe
	# wagpop1file	= 'reb_convwaggs_synthe_'+obj+'_Alpha_250_1050_gfallvinicius.asc'				# Waggs Alpha
	# wagpop2file	= 'reb_convwaggs_synthe_'+obj+'_AlphaCNONaHe_250_1050_gfallvinicius.asc'		# Waggs AlphaCNONaHe

	# schwavepop1, schfluxpop1 = loadfile(schpop1file, pathsyn)
	# schwavepop2, schfluxpop2 = loadfile(schpop2file, pathsyn)
	# wagwavepop1, wagfluxpop1 = loadfile(wagpop1file, pathsyn)
	# wagwavepop2, wagfluxpop2 = loadfile(wagpop2file, pathsyn)	

	# pyphot_schidxs = []
	# pyphot_wagidxs = []
	# pyphot_schpop1 = []
	# pyphot_schpop2 = []
	# pyphot_wagpop1 = []
	# pyphot_wagpop2 = []

	# for band in lib[pyphotbands]: # band = lib[pyphotbands][0]
	# 	f = lib[pyphotbands]
	# 	try:
	# 		pyphot_schidxs.append(band.get(schwave, schflux, axis=1))
	# 		pyphot_schpop1.append(band.get(schwavepop1, schfluxpop1, axis=1))
	# 		pyphot_schpop2.append(band.get(schwavepop2, schfluxpop2, axis=1))

	# 		pyphot_wagidxs.append(band.get(wagwave, wagflux, axis=1))
	# 		pyphot_wagpop1.append(band.get(wagwavepop1, wagfluxpop1, axis=1))
	# 		pyphot_wagpop2.append(band.get(wagwavepop2, wagfluxpop2, axis=1))
	# 	except:
	# 		pyphot_schidxs.append(0)
	# 		pyphot_schpop1.append(0)
	# 		pyphot_schpop2.append(0)

	# 		pyphot_wagidxs.append(0)
	# 		pyphot_wagpop1.append(0)
	# 		pyphot_wagpop2.append(0)

	# flag = obj+' PyPhot Schiavon'
	# pyphot_schidxs.append(calcMgFe(pyphot_schidxs[7], pyphot_schidxs[8], pyphot_schidxs[9]))
	# flag = obj+' PyPhot Syn Schiavon Alpha'
	# pyphot_schpop1.append(calcMgFe(pyphot_schpop1[7], pyphot_schpop1[8], pyphot_schpop1[9]))
	# flag = obj+' PyPhot Syn Schiavon AlphaCNONaHe'
	# pyphot_schpop2.append(calcMgFe(pyphot_schpop2[7], pyphot_schpop2[8], pyphot_schpop2[9]))

	# flag = obj+' PyPhot Waggs'
	# pyphot_wagidxs.append(calcMgFe(pyphot_wagidxs[7], pyphot_wagidxs[8], pyphot_wagidxs[9]))
	# flag = obj+' PyPhot Syn Waggs Alpha'
	# pyphot_wagpop1.append(calcMgFe(pyphot_wagpop1[7], pyphot_wagpop1[8], pyphot_wagpop1[9]))
	# flag = obj+' PyPhot Syn Waggs AlphaCNONaHe'
	# pyphot_wagpop2.append(calcMgFe(pyphot_wagpop2[7], pyphot_wagpop2[8], pyphot_wagpop2[9])) 

	# pyphot_schIdeltas = (np.array(pyphot_schpop2) - np.array(pyphot_schpop1)).tolist()
	# pyphot_wagIdeltas = (np.array(pyphot_wagpop2) - np.array(pyphot_wagpop1)).tolist()

	# pyphot_schidxs = np.round(np.array(pyphot_schidxs),3)
	# pyphot_wagidxs = np.round(np.array(pyphot_wagidxs),3)
	# pyphot_schpop1 = np.round(np.array(pyphot_schpop1),3)
	# pyphot_schpop2 = np.round(np.array(pyphot_schpop2),3)
	# pyphot_wagpop1 = np.round(np.array(pyphot_wagpop1),3)
	# pyphot_wagpop2 = np.round(np.array(pyphot_wagpop2),3)
	# pyphot_schIdeltas = np.round(np.array(pyphot_schIdeltas),3)
	# pyphot_wagIdeltas = np.round(np.array(pyphot_wagIdeltas),3)

	# if save == 1:

	# 	np.savetxt('libraries/synthetic/gc/'+obj+'/LickIndexes_'+obj+'_Lector.txt', 
	# 		np.column_stack(( bands, lector_schidxs, lector_schpop1, lector_schpop2, lector_schIdeltas, 
	# 			lector_wagidxs, lector_wagpop1, lector_wagpop2, lector_wagIdeltas, CoelhoDeltaLickdxs )), 
	# 		fmt='%12s', header='   Indexes       SchObs     SchAlpha   SchCNONaHe    SchiSynDeltas    WaggsObs  WaggsAlpha  WaggsCNONaHe WaggsSynDeltas  CoelhoIdeltas')

	# 	np.savetxt('libraries/synthetic/gc/'+obj+'/LickIndexes_'+obj+'_PyPhot.txt', 
	# 		np.column_stack(( bands, pyphot_schidxs, pyphot_schpop1, pyphot_schpop2, pyphot_schIdeltas, 
	# 			pyphot_wagidxs, pyphot_wagpop1, pyphot_wagpop2, pyphot_wagIdeltas, CoelhoDeltaLickdxs )), 
	# 		fmt='%12s', header='   Indexes       SchObs     SchAlpha   SchCNONaHe    SchiSynDeltas    WaggsObs  WaggsAlpha  WaggsCNONaHe WaggsSynDeltas  CoelhoIdeltas')


	# 	np.savetxt('libraries/synthetic/gc/'+obj+'/LickIndexes_'+obj+'_Schiavon.txt', 
	# 		np.column_stack(( bands, pyphot_schidxs, lector_schidxs, pyphot_schpop1, lector_schpop1, 
	# 			pyphot_schpop2, lector_schpop2, pyphot_schIdeltas, lector_schIdeltas, CoelhoDeltaLickdxs )), 
	# 		fmt='%12s', header='   Indexes    PyPhotObs    LectorObs  PyPhotAlpha  LectorAlpha   PyPhotCNONaHe LectorCNONaHe PyPhotDeltas LectorDeltas   CoelhoIdeltas')

	# 	np.savetxt('libraries/synthetic/gc/'+obj+'/LickIndexes_'+obj+'_Waggs.txt', 
	# 		np.column_stack(( bands, pyphot_wagidxs, lector_wagidxs, pyphot_wagpop1, lector_wagpop1, 
	# 			pyphot_wagpop2, lector_wagpop2, pyphot_wagIdeltas, lector_wagIdeltas, CoelhoDeltaLickdxs )),
	# 		fmt='%12s', header='   Indexes    PyPhotObs    LectorObs  PyPhotAlpha  LectorAlpha   PyPhotCNONaHe LectorCNONaHe PyPhotDeltas LectorDeltas   CoelhoIdeltas')	

