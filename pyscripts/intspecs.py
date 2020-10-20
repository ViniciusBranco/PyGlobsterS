import numpy as np 
from scipy.integrate import simps
import subprocess as sp
import os


def loadfile(path, file):
	return np.array([list(i) for i in tuple(map(tuple, np.loadtxt(path+file).T))])


def calcstars(specfile):
	
	specwave, specflux, specfcont, specfnorm = loadfile('', specfile); specflux = specflux;
	specteff, speclogg = (float(specfile.split('_')[2][1:]), float(specfile.split('_')[3][1:]))
	
	idx = np.where( (Teff_lib == specteff) & (logg_lib == speclogg) )[0]

	print(f'Calculando (Teff , Logg) = ({specteff} , {speclogg}): {Vfinal[idx].size} estrelas - {abundance}')
	
	passed	 = np.interp(specwave, vfilterwave, vfiltertran2)		# interpola o filtro de transmissao para o vetor de fluxo do especto sintetico
	integral = simps( (passed*specflux), specwave )					# calcula a integral usando metodo de simpson
	soma	 = sum(passed*specflux)									# calcula a soma do resultado da multiplicacao dos bins
	
	if Vfinal[idx].size != 0:
		Ci = ( 10**( -( Vfinal[idx]/2.5 ) ) ) / (integral)					# Peso de cada estrela do par teff e logg 
		
		starintflux = 0
		for c in Ci:
			starintflux += specflux*c 										# weighed flux

		return starintflux
	else:
		return specflux*0


### =================================================================================================== ###
###                                                                                                     ###
### ===== Leitura dos modelos da biblioteca Lucimara e filtro de transmissao V ======================== ###

vfilterwave, vfiltertran2, vfiltertran = loadfile('libraries/synthetic/gc/', 'vfilter_transmission.txt')
bfilterwave, bfiltertran = loadfile('libraries/synthetic/gc/', 'Bessel_V-1.txt')

objs = ['NGC0104', 'NGC1904', 'NGC5927']

for obj in objs: # obj = 'NGC5927'

	print('\nIntegrando espectros de', obj)
	lucipath = 'libraries/synthetic/gc/'+obj+'/'
	d, Teff_cmd, logg_cmd, Teff_lib, logg_lib, libcode, BVfinal, Vfinal = loadfile(lucipath, obj+'_martins_coelhoresult_cmd.txt')
	waveluci, fluxsci, fluxlib = loadfile(lucipath, obj+'_martins_coelhoresult_specs.txt')


	### =================================================================================================== ###
	###                                                                                                     ###
	### ===== Integracao dos espectros ==================================================================== ###
	
	specpath = 'libraries/synthetic/gc/'+obj+'/specs/'; os.chdir(specpath)

	sp.call('ls synthe* > speclist.txt', shell=True); speclist = np.loadtxt('speclist.txt', dtype=str)

	pop1intspecname = 'synthe_'+obj+'_Alpha.asc'
	pop2intspecname = 'synthe_'+obj+'_AlphaCNONaHe.asc'

	intfluxpop1, intfluxpop2, ipop1, ipop2 = (0,0,0,0)

	for specfile in speclist:	# specfile = speclist[0]

		abundance = specfile.split('_')[7]
		if abundance == 'Alpha':
			intfluxpop1 += calcstars(specfile); ipop1 += 1
		elif abundance == 'AlphaCNONaHe':
			intfluxpop2 += calcstars(specfile); ipop2 += 1
		else:
			print('ERROR: Abundance not found!')

	wave = loadfile('', specfile)[0]

	os.chdir('../')

	np.savetxt(pop1intspecname, np.column_stack(( wave, intfluxpop1 )))
	np.savetxt(pop2intspecname, np.column_stack(( wave, intfluxpop2 )))

	os.chdir('../../../../')