
                ***************************************
                               LECTOR
                ***************************************

INTRODUCTION
============

LECTOR is a FORTRAN 77 code for measuring line-strengths in one-dimensional
ascii spectra with either lineal or logarithmic wavelength scale. The code
returns the values of the Lick indices (Worthey et al. 1994 (ApJS, 94, 687);
Worthey & Ottaviani 1997 (ApJS, 111, 377)), the Hgamma age indicators of
Vazdekis & Arimoto 1999 (ApJ, 525, 144) and Vazdekis et al. 2001 (ApJ, 549,
274), the indices of Rose 1994 (AJ, 107, 206) and Jones & Worthey 1995 (ApJ,
446, L31), and the CaII triplet and Paschen indices of Cenarro et al. 2001
(MNRAS, 326, 959) (on the basis of a subroutine written by Cardiel & Gorgas).
The code allows you to measure as many indices as you wish if the limits of two
pseudocontinua (at each side of the feature) and the feature itself (i.e.
Lick-style index definition) are provided in a separate file called BANDS. If
requested the program provides index error estimates on the basis of photon
statistics following the formulae given in Cardiel et al. 1998 (A&AS, 127,
597), Vazdekis & Arimoto 1999 (ApJ, 525, 144) and Cenarro et al. 2001 (MNRAS,
326, 959).


INSTALLATION
============

All the required files are contained in the file "LECTOR.tar.gz" and the
installation is straightforward:

> gunzip LECTOR.tar.gz
> tar xvf LECTOR.tar
> f77 -o LECTOR LECTOR.f  (or g77 -o LECTOR LECTOR.f)

to execute the program type

> LECTOR 


INDICES
=======  

The LECTOR program measures the Rose 1994 (AJ,107,206), Jones & Worthey 1995
(ApJ,446,L31) and Cenarro et al. 2001 (MNRAS, 326, 959) set of indices
automatically if they fall into the spectral range covered by the provided
spectrum, otherwise their values will be "99.999". Lick-style indices can be
measured if their two pseudocontinua (at each side of the feature) and the
feature itself are provided in a file called BANDS. The Lick-style indices can
be either expressed in pseudo-equivalent widths or in magnitudes. In this way
we can measure the popular indices of the Lick system (Worthey et al. 1994
(ApJS,94,687); Worthey & Ottaviani 1997 (ApJS,111,377). 

The file BANDS should contain the definitions of the indices. For example, the
BANDS file provided in this distribution include the Lick indices that fall
in the two optical spectral ranges covered by our model SEDs (Vazdekis 1999,
ApJ,513,224), as well as the age indicators of Vazdekis & Arimoto 1999
(Hg_sigma_125, Hg_sigma_200, Hg_sigma_275) and Vazdekis et al. 2001
(Hg_sigma_130). More indices can be measured if appropriate band definitions
are written in the BANDS file. Index definitions should be given after the
second line!. The first lines of the BANDS file look like:


Please write/delete indices after the second horizontal line
_______________________________________________________________________
Blue pseudocont.        Feature       Red pseudocont.  A/M Index name
_______________________________________________________________________
4329.000 4340.468  4333.250 4363.000  4352.500 4368.250 A  Hg_sigma_130


Column 7 tells the program whether the index is expressed as pseudo-equivalent
width in Angstroms (A) or in magnitudes (M).


USE
===

Only one-dimensional ascii spectra are accepted. The spectra might be in lineal
or logarithmic wavelength scale and you do not need to specify this extent. For
example see the /noao/onedspec/wspectext IRAF task to convert fits to ascii
files. The program will first ask you for the following question:

"Single spectrum or list of spectra(s/l)?". 

If a single spectrum were to be measured you should answer "s", otherwise you
should type "l". If the option is "s" the code asks for the name of the
spectrum. If the answer is "l", a list of spectra should be given. The list
should contain at least one column with the names of the spectra. A very
useful option would be to write in a second column the redshift values to be
applied. Otherwise the program will ask for a redshift value to be applied to
the whole set of spectra in the list. An example of a list of spectra is given
below:

spectrum1.ascii  0.
spectrum2.ascii  1027.
spectrum3.ascii  975.

LECTOR does not shift the spectrum but the wavelengths limits defined in
BANDS. Finally the program asks whether the redshift is given in velocity
(Km/s) or z.

SUGESTION: you can correct your spectrum for the REDSHIFT, ROTATION CURVE and
from ANY WAVELENGTH SHIFT by cross-correlation (e.g., using /noao/rv/fxcor
IRAF task) with at least a number of single-age (e.g. 3, 6, 12Gyr) single
metallicity (e.g. [M/H]=0.0) models, and choose the model that provides the
largest cross-correlation peak height.

Finally, the program asks whether you wish to estimate errors. This is
performed on the basis of photon statistics following the formulae given in
Cardiel et al. 1998 (A&AS, 127, 597), Vazdekis & Arimoto 1999 (ApJ, 525, 144)
and Cenarro et al. 2001 (MNRAS, 326, 959). Should calculate errors the code
asks for the conversion electrons/ADU and Readout Noise (RN, in electrons).
The RN should correspond to the total number of pixels in your selected
aperture, including the slit, width. The following formula should work out:

RNtot=RNpix*SQRT(Npix)

where RNpix is the readout noise corresponding to a single pixel and Npix is
the total number of pixels. 

Note: in the current version no error estimates are provided for the Jones &
Worthey 1995 and Cenarro et al. 2001 indices. 

OUTPUT FILES
============

There are three sets of output files. The first two sets are meant to be used
for working and plotting purposes:

1.- Index measurements can be found in the following tables:

"givenlist_LINE" 

where "givenlist" is either the name of the spectrum or the list of spectra.
This table include the obtained values for the indices written in BANDS and
those of Cenarro et al 2001 (in the last three columns). The Cenarro et al
indices are calculated by default. If the BANDS file of the current
distribution is used (with no modification at all) the resulting files will
have the following order:

Column   Index            Source
------ -----------   --------------------
1      Hg_sigma_130  Vazdekis et al 2001
2      Hg_sigma_125  Vazdekis & Arimoto 2001
3      Hg_sigma_200          "           
4      Hg_sigma_275          "
5      Hdelta_A      Worthey & Ottaviani 1997
6      Hgamma_A              "
7      Hdelta_F              "
8      Hgamma_F              "
9      CN1           Worthey et al 1994
10     CN2                   "
11     Ca4227                "
12     G4300                 "
13     Fe4383                "
14     H_beta                "
15     Fe5015                "
16     Mg1                   "
17     Mg2                   "
18     Mgb                   "
19     Fe5270                "
20     Fe5335                "
21     Fe5406                "
22     CaT*          Cenarro et al 2001
23     PaT                   "
24     CaT                   "

"givenlist_ROSE" 

which include the values for the indices of Rose 1994 and Jones & Worthey 1995
according to this order:

Column   Index            Source
------ -----------   --------------------
1      Hdelta/4045   Rose 1994
2      Hdelta/4063	   "
3      SrII/4045	   "
4      SrII/4063	   "
5      Hgamma/Gband	   "
6      Hgamma/4325	   "
7      4289/4271	   "
8      4384/4352	   "
9      p[Fe/H]  	   "
10     CaII		   "
11     3888/3859	   "
12     4220/4208	   "
13     FeI_HR	     Jones & Worthey 1995
14     CaI_HR		   "
15     Hgamma_HR	   "


2.- Index error estimates can be found in:

"givenlist_LINE_ERR"

"givenlist_ROSE_ERR"


3.- Finally, an easy reading version of these tables is provided in:

"givenlist_INDICES"



