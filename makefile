range_spectrum:	range_spectrum.cpp
	`root-config --cxx --cflags` -o range_spectrum range_spectrum.cpp `root-config --glibs`

transmission_spectrum:	transmission_spectrum.cpp
	`root-config --cxx --cflags` -o transmission_spectrum transmission_spectrum.cpp `root-config --glibs`

energy_degradation:	energy_degradation.cpp
	`root-config --cxx --cflags` -o energy_degradation energy_degradation.cpp `root-config --glibs`

dEdX_plotter:	dEdX_plotter.cpp
	`root-config --cxx --cflags` -o dEdX_plotter dEdX_plotter.cpp `root-config --glibs`

espectrum:	espectrum.cpp
	`root-config --cxx --cflags` -o espectrum espectrum.cpp `root-config --glibs`
