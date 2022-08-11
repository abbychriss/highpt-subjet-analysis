#!/usr/bin/env python
# coding: utf-8

# In[41]:


import ROOT
import healpy as hp
import sys
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import gridspec
from matplotlib.legend import Legend
import array
from math import sqrt, fabs, sin, exp, log10
from ROOT import TFile, TTree, TChain, TBranch, TH1D, TH1I, TH1F, TH2F, Math
from ROOT import TLorentzVector, TVector2
from ROOT.Math import LorentzVector, VectorUtil


# In[42]:


tfile = TFile.Open('ntuple_jet_tracks.root')
print(tfile)
tree = tfile.Get("qcdtree")
print(tree)


# In[43]:


branches = tree.GetListOfBranches()
leaves = tree.GetListOfLeaves()


# In[44]:


for lv in leaves:
    lvname = lv.GetName()
    print(lv)


# In[45]:


nevt = 100
nentries = tree.GetEntries()
print("Number of events: ", nentries, "  Printout every: ", 1)


# In[46]:


for iev in range(nentries):
    if iev%nevt==0:
        print("Processing event: ", iev)
    tree.GetEntry(iev)
    if iev%nevt==0:
        print("Genjets: ", tree.nGenJets)
    for ijet in range(tree.nGenJets):
        if iev%nevt==0:
            print(tree.nGenJetParticles[ijet], tree.genJetM[ijet], tree.genJetPx[ijet], tree.genJetPy[ijet], tree.genJetPz[ijet])
print("      JetM                 Px                     Py                Pz ")


# In[47]:


tree.GetEntry(0)
ntracks = tree.nGenTracks
print("Number of tracks: ", ntracks)
print("          P                Pt")
for itrack in range(ntracks):
    print(itrack, tree.trackP[itrack], tree.trackPt[itrack])


# In[48]:


print("          Eta                Phi")
for itrack in range(ntracks):
    print(itrack, tree.trackEta[itrack], tree.trackPhi[itrack])


# In[49]:


print("         Pt                  Eta                 Phi                 E")
tracketaphi_ls = []
for itrack in range(ntracks):
    pt = tree.trackPt[itrack]
    eta = tree.trackEta[itrack]
    phi = tree.trackPhi[itrack]
    E = (tree.trackP[itrack])**2
    tracketaphi = LorentzVector('ROOT::Math::PtEtaPhiE4D<double>')(pt,eta,phi,E)
    tracketaphi_ls.append(tracketaphi)
    print(itrack, tracketaphi.pt(), tracketaphi.eta(), tracketaphi.phi(), tracketaphi.E())
    # tracketaphi are the eta, phi vectors of each track in genTracks


# In[50]:


print("         Pt                  Eta                 Phi                 E")
partetaphi_test = []
matches = []
nmatches = []
for ijet in range(1): #tree.nGenJets):
    for ipart in range(tree.nGenJetParticles[ijet]):
        pdgId = tree.genJetParticlePdgId[ipart]
        m = tree.genJetParticleMvec[ipart]
        px = tree.genJetParticlePxvec[ipart]
        py = tree.genJetParticlePyvec[ipart]
        pz = tree.genJetParticlePzvec[ipart]
        E = np.sqrt(m*m + px*px + py*py + pz*pz)
        p4part = LorentzVector('ROOT::Math::PxPyPzE4D<double>')(px,py,pz,E)

        partetaphi = LorentzVector('ROOT::Math::PtEtaPhiE4D<double>')(p4part.pt(), p4part.eta(), p4part.phi() , p4part.E())
        if (abs(pdgId)==211 or abs(pdgId)==321 or abs(pdgId)==2212 or abs(pdgId)==11 or abs(pdgId)==13 or abs(pdgId)==15):
            partetaphi_test += (partetaphi,)
            print(ipart, partetaphi.pt(), partetaphi.eta(), partetaphi.phi(), partetaphi.E())
            for track in tracketaphi_ls:
                deltaR = VectorUtil.DeltaR(partetaphi, track)
                if deltaR<=0.001:
                    print("matched")
                    matches += [[tracketaphi_ls.index(track), deltaR]]
                    nmatches.append(tracketaphi_ls.index(track))

print("\n Out of", len(tracketaphi_ls), "tracks,", len(set(nmatches)), "have matches to charged gen jet particles. \n")
print("[Track index, deltaR]: \n")
print(matches)
        # partetaphi are the lorentz 4 vectors that contain eta and phi components for the particles in the gen jet


# In[55]:


allgenpart = False #True
event_test = []
for iev in range(nentries):
    jet_test = []
    if iev%nevt==0:
        print("Processing event: ",iev)
    tree.GetEntry(iev)
    print("Genjets: ",tree.nGenJets)
    iptot = 0
    for ijet in range(tree.nGenJets): 
        partetaphi_ls = []
        p4jet = LorentzVector('ROOT::Math::PxPyPzE4D<double>')(0.0,0.0,0.0,0.0)
        for ipart in range(tree.nGenJetParticles[ijet]):
            status = tree.genJetParticleStatus[iptot+ipart]
            pdgId = tree.genJetParticlePdgId[iptot+ipart]
            m = tree.genJetParticleMvec[iptot+ipart]
            px = tree.genJetParticlePxvec[iptot+ipart]
            py = tree.genJetParticlePyvec[iptot+ipart]
            pz = tree.genJetParticlePzvec[iptot+ipart]
            E = np.sqrt(m*m + px*px + py*py + pz*pz)
            p4part = LorentzVector('ROOT::Math::PxPyPzE4D<double>')(px,py,pz,E)
            partetaphi = LorentzVector('ROOT::Math::PtEtaPhiE4D<double>')(p4part.pt(), p4part.eta(), p4part.phi(), p4part.E())
            # select particles
            if(not allgenpart):
                if (abs(pdgId)==211 or abs(pdgId)==321 or abs(pdgId)==2212 or abs(pdgId)==11 or abs(pdgId)==13 or abs(pdgId)==15):
                    for itrack in range(tree.nGenTracks):
                        pt = tree.trackPt[itrack]
                        eta = tree.trackEta[itrack]
                        phi = tree.trackPhi[itrack]
                        E = (tree.trackP[itrack])**2
                        tracketaphi = LorentzVector('ROOT::Math::PtEtaPhiE4D<double>')(pt,eta,phi,E)
                        deltaR = VectorUtil.DeltaR(partetaphi, tracketaphi)
                        if deltaR<=0.001:
                            print("Event",iev,"jet",ijet,"particle",ipart,"matched")
                            print(p4part.pt(), p4part.eta(), p4part.phi() , p4part.E())
                            partetaphi_ls.append(partetaphi)
            else:
                partetaphi_ls.append(partetaphi)
            p4jet += p4part 
        cmJet = p4jet.BoostToCM()
        print("cmJet Boost:",cmJet.x(),cmJet.y(),cmJet.z())
        p4jetcmjet = VectorUtil.boost(p4jet, cmJet)
        print("cmJet:",p4jetcmjet.mass(),p4jetcmjet.px(),p4jetcmjet.py(),p4jetcmjet.pz())
        jet_particle_test = []
        sump2cmjet = 0.0
        for partetaphi in partetaphi_ls:
            p4partetaphi = LorentzVector('ROOT::Math::PxPyPzE4D<double>')(partetaphi.px(), partetaphi.py(), partetaphi.pz(), partetaphi.energy())
            p4partcmjet = VectorUtil.boost(p4partetaphi, cmJet)
            print("Boosted particle:", p4partcmjet.px(), p4partcmjet.py(), p4partcmjet.pz(), p4partcmjet.energy())
            jet_particle_test.append(p4partcmjet) # p4part p4partcmjet for lab or rest frame
            sump2cmjet+=p4partcmjet.px()*p4partcmjet.px()+p4partcmjet.py()*p4partcmjet.py()+p4partcmjet.pz()*p4partcmjet.pz()
        print("number of particles = ",len(jet_particle_test)," sump2cmjet = ",sump2cmjet)
        if (sump2cmjet>5.0*5.0):
            print("new jet sump2 = ",sump2cmjet)
            jet_test.append(jet_particle_test)
        iptot+=tree.nGenJetParticles[ijet]
    event_test.append(jet_test)


# In[60]:


#Sky map of particles in jet that have matches to tracks
singlepixel = False

nside = 128*4

iev = 0
NPIX = hp.nside2npix(nside)
cls_list = []
print("number of events = ",len(event_test))
for jet_test in event_test:
    print("number of jets = ",len(jet_test))
    iev+=1
    ijet = 0
    for jet_particle_test in jet_test:
        print("number of particles = ",len(jet_particle_test))
        ijet+=1
        m = np.zeros(NPIX) # blank map
        for part in jet_particle_test:
            if(singlepixel):
                ipix_disc = hp.pixelfunc.ang2pix(nside=nside, theta=part.theta(), phi=part.phi())
            else: # choose radius
                vec = hp.ang2vec(part.theta(), part.phi())
                ipix_disc = hp.query_disc(nside=nside, vec=vec, radius=np.radians(1)) # all pixels within 1 degree of vector
            m[ipix_disc] += part.energy() # energy of particle
        cls_hp = hp.sphtfunc.anafast(m)
        cls_list.append(cls_hp)
        if iev in range(1,10):
            hp.visufunc.mollview(m, title=r'Single jet - DeltaR ≤ 0.001', unit=r'GeV', min=0.01, norm='log', badcolor='black')
            hp.graticule(color='grey') # draw grid lines


# In[64]:


singlepixel = False
lmax = 360*4 #150
ls = np.array(range(1,lmax+1))

#Sky map of particles in jet that have matches to tracks
nside = 128*4
iev = 0
NPIX = hp.nside2npix(nside)
cls_list = []
print("number of events = ",len(event_test))
for jet_test in event_test:
    print("number of jets = ",len(jet_test))
    iev+=1
    ijet = 0
    for jet_particle_test in jet_test:
        print("number of particles = ",len(jet_particle_test))
        ijet+=1
        m = np.zeros(NPIX) # blank map
        for part in jet_particle_test:
            if(singlepixel):
                ipix_disc = hp.pixelfunc.ang2pix(nside=nside, theta=part.theta(), phi=part.phi())
                energy = part.energy()
            else: # choose radius
                vec = hp.ang2vec(part.theta(), part.phi())
                ipix_disc = hp.query_disc(nside=nside, vec=vec, radius=np.radians(0.1)) # all pixels within 0.1 degree of vector
                energy = part.energy()/len(ipix_disc)
            m[ipix_disc] += energy # energy of particle
        cls_hp = hp.sphtfunc.anafast(m)
        cls_list.append(cls_hp)
print("number of cls = ",len(cls_list))
plt.figure(figsize=(8,8))
icls = 0
for cls_test in cls_list:
    icls+=1
    if (icls == 1):
        cls_hp = cls_test
    else:
        cls_hp += cls_test 
plt.semilogy(ls,cls_hp[ls],label='nside = {}'.format(nside), color='#1ABC9C')
#plt.ylim((5e-8,None))
plt.legend()


# In[67]:


nside=128*4
lmax=360*4
mapredo = hp.sphtfunc.synfast(cls_hp, nside, lmax=80)
hp.visufunc.mollview(mapredo, title=r'Jet from Cls', unit=r'GeV', min=0.01, norm='log', badcolor='black')
hp.graticule()


# In[ ]:




