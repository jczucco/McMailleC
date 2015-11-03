===========================================================================
            McMaille is a program for indexing by Monte Carlo
                (Maille in french = cell in english)
===========================================================================

The 2-theta peak positions extracted from a peak hunting program are
used together with the intensities in order to build a pseudo powder
pattern to which are compared patterns calculated from the cell parameters 
proposed by a Monte Carlo process. Peak shapes are columnar in version 3.0.
The best cells are refined, more or less. This is similar to the (still 
unavailable ?) software by B.M. Karuki et al., J. Synchrotron Rad. 6. (1999) 
87-92, though the latter uses a genetic algorithm and the raw data.
                                                 Armel Le Bail
                                                 alb@cristal.org
                                                 September 2002
                                                 November  2006
===========================================================================

The package MvMaille-V4.zip contains mainly :

  McMaille.exe      : executable for MS Windows, single processor
  McMaille.for      : the FORTRAN source code
  license.html      : the GNU Public License (GPL)  

  PMcMaille.zip     : contains the executable for the multi core processors 
                      (named McMaille.exe) and a DLL named libguide40.dll
                      which has to be installed in the same directory.
                      McMaille.for with OpenMP directives is also there.

  McMaille-v4.html  : the complete manual
  short-manual.html : a short manual version for the automated "black-
                          box mode

  McMaille.pdf      : copy of the McMaille publication
  benchmarks.pdf    : copy of the benchmarks publication

  tests.zip         : test files
  benchmarks.zip    : more test files (benchmarks)
  sdpdrr2.zip       : more test files (from the 'Structure Determination
                      by Powder Diffractometry Round Robin 2')

  cub.hkl, hex.hkl, rho.hkl, tet.hkl, ort.hkl, mon.hkl, tri.hkl
                    : the prepared lists of hkl Miller indices which have
                      to be installed in the same directory as the
                      McMaille.exe program 


===========================================================================

           The manual : McMaille-v4.html included in the package is
             also available at http://www.cristal.org/mcmaille/
                or at http://sdpd.univ-lemans.fr/mcmaille/

===========================================================================

            In case of successful use, please cite the paper :
             A. Le Bail, Powder Diffraction 19 (2004) 249-254.
                       (included into the package)

===========================================================================

           See comparisons of indexing software including McMaille :
           J. Bergmann, A. Le Bail, R. Shirley and V. Zlokazov, 
                     Z. Kristallogr. 219 (2004) 783-790.
                       (included into the package)

                 and visit the Indexing Benchmarks Web page at :
                  http://sdpd.univ-lemans.fr/uppw/benchmarks/

===========================================================================
