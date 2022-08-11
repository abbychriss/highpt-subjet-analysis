#!/usr/bin/env python
# coding: utf-8

# In[1]:


import ROOT
import healpy as hp


# In[2]:


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


# In[3]:


tfile = TFile.Open('ntuple_cern2.root')
#tfile = TFile.Open('ntuple.root')
#tfile = TFile.Open('ntuple_1k.root')
print(tfile)
tree = tfile.Get("qcdtree")
print(tree)


# In[4]:


branches = tree.GetListOfBranches()
leaves = tree.GetListOfLeaves()


# In[5]:


for lv in leaves:
    lvname = lv.GetName()
    print(lv)


# In[6]:


nevt = 10
nentries = tree.GetEntries()
print("Number of events: ",nentries, "  printout every ",nevt)


# In[7]:


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


# In[8]:


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
            status = tree.genJetParticleStatus[iptot+ipart]
            pdgId = tree.genJetParticlePdgId[iptot+ipart]
            m = tree.genJetParticleMvec[iptot+ipart]
            px = tree.genJetParticlePxvec[iptot+ipart]
            py = tree.genJetParticlePyvec[iptot+ipart]
            pz = tree.genJetParticlePzvec[iptot+ipart]
            E = np.sqrt(m*m + px*px + py*py + pz*pz)
            p4part = LorentzVector('ROOT::Math::PxPyPzE4D<double>')(px,py,pz,E)
            p4jet += p4part
            #print(tree.genJetParticleStatus[iptot+ipart],tree.genJetParticlePdgId[iptot+ipart],tree.genJetParticleMvec[iptot+ipart],tree.genJetParticlePxvec[iptot+ipart],tree.genJetParticlePyvec[iptot+ipart],tree.genJetParticlePzvec[iptot+ipart])
        iptot+=tree.nGenJetParticles[ijet]
        if iev%nevt==0:
            print(" redo:",p4jet.mass(),p4jet.px(),p4jet.py(),p4jet.pz())


# In[9]:


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
            status = tree.genJetParticleStatus[iptot+ipart]
            pdgId = tree.genJetParticlePdgId[iptot+ipart]
            m = tree.genJetParticleMvec[iptot+ipart]
            px = tree.genJetParticlePxvec[iptot+ipart]
            py = tree.genJetParticlePyvec[iptot+ipart]
            pz = tree.genJetParticlePzvec[iptot+ipart]
            E = np.sqrt(m*m + px*px + py*py + pz*pz)
            p4part = LorentzVector('ROOT::Math::PxPyPzE4D<double>')(px,py,pz,E)
            # p4part are the jet generator particle 4-momenta in the lab frame
            p4jet += p4part
            #print(tree.genJetParticleMvec[iptot+ipart],tree.genJetParticlePxvec[iptot+ipart],tree.genJetParticlePyvec[iptot+ipart],tree.genJetParticlePzvec[iptot+ipart])
        #print(" redo:",p4jet.mass(),p4jet.px(),p4jet.py(),p4jet.pz())
        cmJet = p4jet.BoostToCM()
        if iev%nevt==0:
            print(" cmJet Boost:",cmJet.x(),cmJet.y(),cmJet.z())
        p4jetcmjet = VectorUtil.boost(p4jet, cmJet)
        if iev%nevt==0:
            print(" cmJet:",p4jetcmjet.mass(),p4jetcmjet.px(),p4jetcmjet.py(),p4jetcmjet.pz())
        p4jetcm = LorentzVector('ROOT::Math::PxPyPzE4D<double>')(0.0,0.0,0.0,0.0)
        jet_particle_test = []
        sump2cmjet = 0.0
        for ipart in range(tree.nGenJetParticles[ijet]):
            status = tree.genJetParticleStatus[iptot+ipart]
            pdgId = tree.genJetParticlePdgId[iptot+ipart]
            m = tree.genJetParticleMvec[iptot+ipart]
            px = tree.genJetParticlePxvec[iptot+ipart]
            py = tree.genJetParticlePyvec[iptot+ipart]
            pz = tree.genJetParticlePzvec[iptot+ipart]
            E = np.sqrt(m*m + px*px + py*py + pz*pz)
            p4part = LorentzVector('ROOT::Math::PxPyPzE4D<double>')(px,py,pz,E)
            p4partcmjet = VectorUtil.boost(p4part, cmJet)
            # p4partcmjet are the jet generator particle 4-momenta in the jet center-of-mass
            p4jetcm += p4partcmjet
#            if (iev==0 and ijet==1): # test on 2nd jet in 1st event
#restrict to charged particles charged pions abs(pdgId)=211, charged kaons abs(pdgId)=321 and protons abd(pdgId)=2212, electrons abs(pdgId)=11, muons abs(pdgId)=13, tau abs(pdgId)=15
            if (abs(pdgId)==211 or abs(pdgId)==321 or abs(pdgId)==2212 or abs(pdgId)==11 or abs(pdgId)==13 or abs(pdgId)==15):
                jet_particle_test.append(p4partcmjet) # p4part p4partcmjet for lab or rest frame
                sump2cmjet+=p4partcmjet.px()*p4partcmjet.px()+p4partcmjet.py()*p4partcmjet.py()+p4partcmjet.pz()*p4partcmjet.pz()
        if (sump2cmjet>50.0*50.0):
            print("new jet sump2 = ",sump2cmjet)
            jet_test.append(jet_particle_test)
        if iev%nevt==0:
            print(" cmredo:",p4jetcm.mass(),p4jetcm.px(),p4jetcm.py(),p4jetcm.pz())
        iptot+=tree.nGenJetParticles[ijet]


# In[10]:


lmax = 150
ls = np.array(range(1,lmax+1))
nside = 128
NPIX = hp.nside2npix(nside)
for n in range(1,10,1):
    ijet = 0
    cls_list = []
    for jet_particle_test in jet_test:
        m = np.zeros(NPIX) # blank map
        ijet+=1
        for part in jet_particle_test:
            vec = hp.ang2vec(part.theta(), part.phi())
            ipix_disc = hp.query_disc(nside=nside, vec=vec, radius=np.radians(n)) # all pixels within 1 degree of vector
            m[ipix_disc] += part.energy()/(np.pi*n*n) # energy of particle
        cls_hp = hp.sphtfunc.anafast(m)
        cls_list.append(cls_hp)
        if ijet==3:
            hp.visufunc.mollview(m, title=r'Single jet - radius '+ str(n), unit=r'GeV', min=0.01, norm='log', badcolor='#FFD8CD')
            hp.graticule() # draw grid lines


# In[11]:


lmax = 150
ls = np.array(range(1,lmax+1))
nside = 128
NPIX = hp.nside2npix(nside)
radii=[]
plt.figure(figsize=(8,8))

for n in range(1,10,1):
    ijet = 0
    cls_list = []
    for jet_particle_test in jet_test:
        m = np.zeros(NPIX) # blank map
        ijet+=1
        for part in jet_particle_test:
            vec = hp.ang2vec(part.theta(), part.phi())
            ipix_disc = hp.query_disc(nside=nside, vec=vec, radius=np.radians(n)) # all pixels within 1 degree of vector
            m[ipix_disc] += part.energy()/(n*n) # energy of particle
        cls_hp = hp.sphtfunc.anafast(m)
        cls_list.append(cls_hp)
        #if ijet==3:
            #hp.visufunc.mollview(m, title=r'Single jet - radius '+ str(n), unit=r'GeV', min=0.01, norm='log', badcolor='#FFD8CD')
            #hp.graticule() # draw grid lines
            
#Plot angular power spectrum for each jet
    icls = 0
    for cls_test in cls_list:
        icls+=1
        if (icls == 1):
            cls_hp = cls_test
        else:
            cls_hp+= cls_test
    #plt.semilogy(ls,cls_hp[ls],label='radius = ' + str(n))
    plt.plot(ls,cls_hp[ls],label='radius = ' + str(n) + '°')

    radii.append(n)
print('Radii: ' + str(radii))
plt.legend()


# In[ ]:





# In[12]:


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
            status = tree.genJetParticleStatus[iptot+ipart]
            pdgId = tree.genJetParticlePdgId[iptot+ipart]
            m = tree.genJetParticleMvec[iptot+ipart]
            px = tree.genJetParticlePxvec[iptot+ipart]
            py = tree.genJetParticlePyvec[iptot+ipart]
            pz = tree.genJetParticlePzvec[iptot+ipart]
            E = np.sqrt(m*m + px*px + py*py + pz*pz)
            p4part = LorentzVector('ROOT::Math::PxPyPzE4D<double>')(px,py,pz,E)
            p4jet += p4part
            #print(tree.genJetParticleStatus[iptot+ipart],tree.genJetParticlePdgId[iptot+ipart],tree.genJetParticleMvec[iptot+ipart],tree.genJetParticlePxvec[iptot+ipart],tree.genJetParticlePyvec[iptot+ipart],tree.genJetParticlePzvec[iptot+ipart])
        iptot+=tree.nGenJetParticles[ijet]
        if iev%nevt==0:
            print(" redo:",p4jet.mass(),p4jet.px(),p4jet.py(),p4jet.pz())
            
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
            status = tree.genJetParticleStatus[iptot+ipart]
            pdgId = tree.genJetParticlePdgId[iptot+ipart]
            m = tree.genJetParticleMvec[iptot+ipart]
            px = tree.genJetParticlePxvec[iptot+ipart]
            py = tree.genJetParticlePyvec[iptot+ipart]
            pz = tree.genJetParticlePzvec[iptot+ipart]
            E = np.sqrt(m*m + px*px + py*py + pz*pz)
            p4part = LorentzVector('ROOT::Math::PxPyPzE4D<double>')(px,py,pz,E)
            # p4part are the jet generator particle 4-momenta in the lab frame
            p4jet += p4part
            #print(tree.genJetParticleMvec[iptot+ipart],tree.genJetParticlePxvec[iptot+ipart],tree.genJetParticlePyvec[iptot+ipart],tree.genJetParticlePzvec[iptot+ipart])
        #print(" redo:",p4jet.mass(),p4jet.px(),p4jet.py(),p4jet.pz())
        cmJet = p4jet.BoostToCM()
        if iev%nevt==0:
            print(" cmJet Boost:",cmJet.x(),cmJet.y(),cmJet.z())
        p4jetcmjet = VectorUtil.boost(p4jet, cmJet)
        if iev%nevt==0:
            print(" cmJet:",p4jetcmjet.mass(),p4jetcmjet.px(),p4jetcmjet.py(),p4jetcmjet.pz())
        p4jetcm = LorentzVector('ROOT::Math::PxPyPzE4D<double>')(0.0,0.0,0.0,0.0)
        jet_particle_test = []
        sump2cmjet = 0.0
        for ipart in range(tree.nGenJetParticles[ijet]):
            status = tree.genJetParticleStatus[iptot+ipart]
            pdgId = tree.genJetParticlePdgId[iptot+ipart]
            m = tree.genJetParticleMvec[iptot+ipart]
            px = tree.genJetParticlePxvec[iptot+ipart]
            py = tree.genJetParticlePyvec[iptot+ipart]
            pz = tree.genJetParticlePzvec[iptot+ipart]
            E = np.sqrt(m*m + px*px + py*py + pz*pz)
            p4part = LorentzVector('ROOT::Math::PxPyPzE4D<double>')(px,py,pz,E)
            p4partcmjet = VectorUtil.boost(p4part, cmJet)
            # p4partcmjet are the jet generator particle 4-momenta in the jet center-of-mass
            p4jetcm += p4partcmjet
#            if (iev==0 and ijet==1): # test on 2nd jet in 1st event
#restrict to charged particles charged pions abs(pdgId)=211, charged kaons abs(pdgId)=321 and protons abd(pdgId)=2212, electrons abs(pdgId)=11, muons abs(pdgId)=13, tau abs(pdgId)=15
            if (abs(pdgId)==211 or abs(pdgId)==321 or abs(pdgId)==2212 or abs(pdgId)==11 or abs(pdgId)==13 or abs(pdgId)==15):
                jet_particle_test.append(p4part) #cmjet) # p4part p4partcmjet for lab or rest frame
                sump2cmjet+=p4partcmjet.px()*p4partcmjet.px()+p4partcmjet.py()*p4partcmjet.py()+p4partcmjet.pz()*p4partcmjet.pz()
        if (sump2cmjet>50.0*50.0):
            print("new jet sump2 = ",sump2cmjet)
            jet_test.append(jet_particle_test)
        if iev%nevt==0:
            print(" cmredo:",p4jetcm.mass(),p4jetcm.px(),p4jetcm.py(),p4jetcm.pz())
        iptot+=tree.nGenJetParticles[ijet]


# In[13]:


lmax = 150
ls = np.array(range(1,lmax+1))
nside = 128
NPIX = hp.nside2npix(nside)
for n in range(1,10,1):
    ijet = 0
    cls_list = []
    for jet_particle_test in jet_test:
        m = np.zeros(NPIX) # blank map
        ijet+=1
        for part in jet_particle_test:
            vec = hp.ang2vec(part.theta(), part.phi())
            ipix_disc = hp.query_disc(nside=nside, vec=vec, radius=np.radians(n)) # all pixels within 1 degree of vector
            m[ipix_disc] += part.energy()/(n*n) # energy of particle
        cls_hp = hp.sphtfunc.anafast(m)
        cls_list.append(cls_hp)
        if ijet==3:
            hp.visufunc.mollview(m, title=r'Single jet - radius '+ str(n), unit=r'GeV', min=0.01, norm='log', badcolor='#FFD8CD')
            hp.graticule() # draw grid lines


# In[15]:


lmax = 150
ls = np.array(range(1,lmax+1))
nside = 128
NPIX = hp.nside2npix(nside)
radii=[]
plt.figure(figsize=(8,8))

for n in range(1,10,1):
    ijet = 0
    cls_list = []
    for jet_particle_test in jet_test:
        m = np.zeros(NPIX) # blank map
        ijet+=1
        for part in jet_particle_test:
            vec = hp.ang2vec(part.theta(), part.phi())
            ipix_disc = hp.query_disc(nside=nside, vec=vec, radius=np.radians(n)) # all pixels within 1 degree of vector
            m[ipix_disc] += part.energy()/(n*n) # energy of particle
        cls_hp = hp.sphtfunc.anafast(m)
        cls_list.append(cls_hp)
        #if ijet==3:
            #hp.visufunc.mollview(m, title=r'Single jet - radius '+ str(n), unit=r'GeV', min=0.01, norm='log', badcolor='#FFD8CD')
            #hp.graticule() # draw grid lines
            
#Plot angular power spectrum for each jet
    icls = 0
    for cls_test in cls_list:
        icls+=1
        if (icls == 1):
            cls_hp = cls_test
        else:
            cls_hp+= cls_test
    #plt.semilogy(ls,cls_hp[ls],label='radius = ' + str(n))
    plt.plot(ls,cls_hp[ls],label='radius = ' + str(n) + '°')

    radii.append(n)
print('Radii: ' + str(radii))
plt.legend()


# In[ ]:




