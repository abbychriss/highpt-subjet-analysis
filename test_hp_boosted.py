#!/usr/bin/env python
# coding: utf-8

# In[1]:


import ROOT
import healpy as hp


# In[4]:


import sys
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import gridspec
from matplotlib.legend import Legend
import array
from math import sqrt, fabs, sin, exp, log10
from ROOT import TFile, TTree, TChain, TBranch, TH1D, TH1I, TH1F, TH2F, Math
from ROOT import TLorentzVector
from ROOT.Math import LorentzVector, VectorUtil


# In[5]:


lmax = 150
ls = np.array(range(1,lmax+1))
print(ls)


# In[6]:


cl_test = np.ones(lmax)
cl_maps = np.concatenate((cl_test, [0]))
cl_maps = np.roll(cl_maps, 1)

nside = 128
nside_coarse = 4
maps = hp.sphtfunc.synfast(cl_maps, nside=nside, lmax=None, pol=False)
hp.mollview(maps, title=r'Test map', unit=r'arb')
cls_hp = hp.sphtfunc.anafast(maps)


# In[7]:


plt.figure(figsize=(8,8))
np.random.seed(5)

plt.semilogy(ls,cls_hp[ls],label='nside = {}'.format(nside))
plt.legend()


# In[20]:


NPIX = hp.nside2npix(nside)
print(NPIX)
ipix = hp.ang2pix(nside, np.pi/2, np.pi*3/4) # nside, theta, phi in radians and RING pixel ordering
print(ipix)
vec = hp.ang2vec(np.pi/2, np.pi*3/4)
print(vec)
ipix_disc = hp.query_disc(nside=nside, vec=vec, radius=np.radians(10)) # all pixels within 10 degrees of vector
m = np.arange(NPIX) # top/bottom gradiant map
m[ipix_disc] = m.max() # max value of map
hp.mollview(m, title=r'Test particle', unit=r'GeV')
hp.graticule() # draw grid lines


# In[21]:


tfile = TFile.Open('ntuple_cern.root')
#tfile = TFile.Open('ntuple.root')
#tfile = TFile.Open('ntuple_1k.root')
print(tfile)
tree = tfile.Get("qcdtree")
print(tree)


# In[22]:


branches = tree.GetListOfBranches()
leaves = tree.GetListOfLeaves()


# In[23]:


for lv in leaves:
    lvname = lv.GetName()
    print(lv)


# In[24]:


nevt = 10
nentries = tree.GetEntries()
print("Number of events: ",nentries, "  printout every ",nevt)


# In[25]:


for iev in range(nentries):
    if iev%nevt==0:
        print("Processing event: ",iev)
    tree.GetEntry(iev)
    
    if iev%nevt==0:
        print("Genjets: ",tree.nGenJets)
    iptot = 0
    for ijet in range(tree.nGenJets):
        if iev%nevt==0:
            print(tree.nGenJetParticles[ijet],tree.genJetM[ijet],tree.genJetPx[ijet],tree.genJetPy[ijet],tree.genJetPz[ijet])
        #for ipart in range(tree.nGenJetParticles[ijet]):
            #print(tree.genJetParticleMvec[iptot+ipart],tree.genJetParticlePxvec[iptot+ipart],tree.genJetParticlePyvec[iptot+ipart],tree.genJetParticlePzvec[iptot+ipart])
        iptot+=tree.nGenJetParticles[ijet]


# In[26]:


for iev in range(nentries):
    if iev%nevt==0:
        print("Processing event: ",iev)
    tree.GetEntry(iev)
    
    if iev%nevt==0:
        print("Genjets: ",tree.nGenJets)
    iptot = 0
    for ijet in range(tree.nGenJets):
        if iev%nevt==0:
            print(tree.nGenJetParticles[ijet],tree.genJetM[ijet],tree.genJetPx[ijet],tree.genJetPy[ijet],tree.genJetPz[ijet])
        p4jet = LorentzVector('ROOT::Math::PxPyPzE4D<double>')(0.0,0.0,0.0,0.0)
        for ipart in range(tree.nGenJetParticles[ijet]):
            m = tree.genJetParticleMvec[iptot+ipart]
            px = tree.genJetParticlePxvec[iptot+ipart]
            py = tree.genJetParticlePyvec[iptot+ipart]
            pz = tree.genJetParticlePzvec[iptot+ipart]
            E = numpy.sqrt(m*m + px*px + py*py + pz*pz)
            p4part = LorentzVector('ROOT::Math::PxPyPzE4D<double>')(px,py,pz,E)
            p4jet += p4part
            #print(tree.genJetParticleMvec[iptot+ipart],tree.genJetParticlePxvec[iptot+ipart],tree.genJetParticlePyvec[iptot+ipart],tree.genJetParticlePzvec[iptot+ipart])
        iptot+=tree.nGenJetParticles[ijet]
        if iev%nevt==0:
            print(" redo:",p4jet.mass(),p4jet.px(),p4jet.py(),p4jet.pz())


# In[31]:


jet_test = []
for iev in range(nentries):
    if iev%nevt==0:
        print("Processing event: ",iev)
    tree.GetEntry(iev)
    
    if iev%nevt==0:
        print("Genjets: ",tree.nGenJets)
    iptot = 0
    for ijet in range(tree.nGenJets):
        #print(tree.nGenJetParticles[ijet],tree.genJetM[ijet],tree.genJetPx[ijet],tree.genJetPy[ijet],tree.genJetPz[ijet])
        p4jet = LorentzVector('ROOT::Math::PxPyPzE4D<double>')(0.0,0.0,0.0,0.0)
        for ipart in range(tree.nGenJetParticles[ijet]):
            m = tree.genJetParticleMvec[iptot+ipart]
            px = tree.genJetParticlePxvec[iptot+ipart]
            py = tree.genJetParticlePyvec[iptot+ipart]
            pz = tree.genJetParticlePzvec[iptot+ipart]
            E = numpy.sqrt(m*m + px*px + py*py + pz*pz)
            p4part = LorentzVector('ROOT::Math::PxPyPzE4D<double>')(px,py,pz,E)
            # p4part are the jet generator particle 4-momenta in the lab frame
            p4jet += p4part
            #print(tree.genJetParticleMvec[iptot+ipart],tree.genJetParticlePxvec[iptot+ipart],tree.genJetParticlePyvec[iptot+ipart],tree.genJetParticlePzvec[iptot+ipart])
        #print(" redo:",p4jet.mass(),p4jet.px(),p4jet.py(),p4jet.pz())
        #cmJet = p4jet.BoostToCM()
        #if iev%nevt==0:
            #print(" cmJet Boost:",cmJet.x(),cmJet.y(),cmJet.z())
        #p4jetcmjet = VectorUtil.boost(p4jet, cmJet)
        #if iev%nevt==0:
            #print(" cmJet:",p4jetcmjet.mass(),p4jetcmjet.px(),p4jetcmjet.py(),p4jetcmjet.pz())
        #p4jetcm = LorentzVector('ROOT::Math::PxPyPzE4D<double>')(0.0,0.0,0.0,0.0)
        jet_particle_test = []
        sump2 = 0.0
        for ipart in range(tree.nGenJetParticles[ijet]):
            m = tree.genJetParticleMvec[iptot+ipart]
            px = tree.genJetParticlePxvec[iptot+ipart]
            py = tree.genJetParticlePyvec[iptot+ipart]
            pz = tree.genJetParticlePzvec[iptot+ipart]
            E = numpy.sqrt(m*m + px*px + py*py + pz*pz)
            p4part = LorentzVector('ROOT::Math::PxPyPzE4D<double>')(px,py,pz,E)
            #p4partcmjet = VectorUtil.boost(p4part, cmJet)
            # p4partcmjet are the jet generator particle 4-momenta in the jet center-of-mass
            #p4jetcm += p4partcmjet
#            if (iev==0 and ijet==1): # test on 2nd jet in 1st event
            jet_particle_test.append(p4part) # p4part p4partcmjet for lab or rest frame
            sump2+=p4part.px()*p4part.px()+p4part.py()*p4part.py()+p4part.pz()*p4part.pz()
        if (sump2>100.0*100.0):
            print("new jet sump2 = ",sump2)
            jet_test.append(jet_particle_test)
        if iev%nevt==0:
            print(" redo:",p4jet.mass(),p4jet.px(),p4jet.py(),p4jet.pz())
        iptot+=tree.nGenJetParticles[ijet]


# In[32]:


ijet = 0
cls_list = []
for jet_particle_test in jet_test:
    m = np.zeros(NPIX) # blank map
    ijet+=1
    for part in jet_particle_test:
        vec = hp.ang2vec(part.theta(), part.phi())
        ipix_disc = hp.query_disc(nside=nside, vec=vec, radius=np.radians(1)) # all pixels within 1 degrees of vector
        m[ipix_disc] = 1 #part.energy() # energy of particle
    cls_hp = hp.sphtfunc.anafast(m)
    cls_list.append(cls_hp)
    if (ijet == 3):
        hp.mollview(m, title=r'Test particle', unit=r'GeV')
        hp.graticule() # draw grid lines


# In[29]:


plt.figure(figsize=(8,8))

icls = 0
for cls_test in cls_list:
    icls+=1
    if (icls == 1):
        cls_hp = cls_test
    else:
        cls_hp+= cls_test
plt.semilogy(ls,cls_hp[ls],label='nside = {}'.format(nside))
plt.legend()


# In[ ]:




