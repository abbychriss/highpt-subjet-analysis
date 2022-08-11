# highpt-jet-analysis
This analysis takes the approach of a two point energy correlation function, using CMB techniques to study high-pt jet substructures. We plot angular power spectra (Cls) to look for angular correlations between the distribution of particles within a jet.
We take advantage of the healpy package in HEALPix (Hierarchical Equal Area isoLatitude Pixelation) scheme to produce sky maps of the particles within a jet.
We used CMSSW version 12.4.0 to simulate, digitize, and reconstruct our Monte Carlo events and PYTHIA version 8.3 to simulate collisions.
Use PYTHIA generated anti-kt jets with an R parameter of 0.8 to produce our Monte Carlo events. Started with 10 generated events, then moved to 100.
First we boost the generated particles in each jet into the CM frame of the jet, then plot the angular distributions of the particles in a single jet in a HEALPix Mollweide view sky map. Sky maps are plotted on a log scale weighted by energy. Size of particles are determined by radius of pixels around the particle position vector in polar coordinates.
Next we extract a Cl from the map of a single jet using the healpy anafast function. Plot in multipole moments.
We chose nside = 128 to tesselate our spheres and lmax = 150. Later we increase nside to 512 and lmax to 1440 to study smaller angular scales.
Loop through all the jets in a given event, and every generated event in the ntuple, concatenating the angular power spectra from each jet into a large list of Cls.
Then we plot the average Cl of all gen jets in all gen events and use healpy synfast function to reproduce a statistical sampling of a jet.
These statistical samplings represent a Poisson mean per pixel. We would throw a Poisson trial per pixel extrapolate the probability that a particle occurs in a given location.
The particle flux per pixel in the map made from Cls of many gen jets is determined by Poisson statistics per pixel. We are looking for a correlation between pixels separated by a given angular distance weighted by energy.
We are looking for the equivalent of CMB acoustic peaks in our angular power spectra which may appear as peaks in the average Cl over many events.
If we find this, it may represent a quantum mechanical fluctuation in the color field that occurs during hadronization, in a highly non-perturbative and strongly coupled region of QCD.
