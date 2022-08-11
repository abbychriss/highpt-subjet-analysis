#!/usr/bin/env python
# coding: utf-8

# In[16]:


import ROOT
import healpy


# In[18]:


import sys
import numpy
import array
from math import sqrt, fabs, sin, exp, log10
from ROOT import TFile, TTree, TChain, TBranch, TH1D, TH1I, TH1F, TH2F, Math
from ROOT import TLorentzVector
from ROOT.Math import LorentzVector, VectorUtil


# In[19]:


tfile = TFile.Open('ntuple.root')
print(tfile)
tree = tfile.Get("qcdtree")
print(tree)


# In[20]:


branches = tree.GetListOfBranches()
leaves = tree.GetListOfLeaves()


# In[21]:


for lv in leaves:
    lvname = lv.GetName()
    print(lv)


# In[22]:


nentries = tree.GetEntries()
print("Number of events: ",nentries)


# In[23]:


for iev in range(nentries):
    if iev%1==0:
        print("Processing event: ",iev)
    tree.GetEntry(iev)
    
    print(tree.nGenJets)
    iptot = 0
    for ijet in range(tree.nGenJets):
        print(tree.nGenJetParticles[ijet],tree.genJetM[ijet],tree.genJetPx[ijet],tree.genJetPy[ijet],tree.genJetPz[ijet])
        for ipart in range(tree.nGenJetParticles[ijet]):
            print(tree.genJetParticleMvec[iptot+ipart],tree.genJetParticlePxvec[iptot+ipart],tree.genJetParticlePyvec[iptot+ipart],tree.genJetParticlePzvec[iptot+ipart])
        iptot+=tree.nGenJetParticles[ijet]


# In[24]:


for iev in range(nentries):
    if iev%1==0:
        print("Processing event: ",iev)
    tree.GetEntry(iev)
    
    print(tree.nGenJets)
    iptot = 0
    for ijet in range(tree.nGenJets):
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
        print(" redo:",p4jet.mass(),p4jet.px(),p4jet.py(),p4jet.pz())


# In[31]:


for iev in range(nentries):
    if iev%1==0:
        print("Processing event: ",iev)
    tree.GetEntry(iev)
    
    print(tree.nGenJets)
    iptot = 0
    for ijet in range(tree.nGenJets):
        print(tree.nGenJetParticles[ijet],tree.genJetM[ijet],tree.genJetPx[ijet],tree.genJetPy[ijet],tree.genJetPz[ijet])
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
        print(" redo:",p4jet.mass(),p4jet.px(),p4jet.py(),p4jet.pz())
        cmJet = p4jet.BoostToCM()
        print(" cmJet Boost:",cmJet.x(),cmJet.y(),cmJet.z())
        p4jetcmjet = VectorUtil.boost(p4jet, cmJet)
        print(" cmJet:",p4jetcmjet.mass(),p4jetcmjet.px(),p4jetcmjet.py(),p4jetcmjet.pz())
        p4jetcm = LorentzVector('ROOT::Math::PxPyPzE4D<double>')(0.0,0.0,0.0,0.0)
        for ipart in range(tree.nGenJetParticles[ijet]):
            m = tree.genJetParticleMvec[iptot+ipart]
            px = tree.genJetParticlePxvec[iptot+ipart]
            py = tree.genJetParticlePyvec[iptot+ipart]
            pz = tree.genJetParticlePzvec[iptot+ipart]
            E = numpy.sqrt(m*m + px*px + py*py + pz*pz)
            p4part = LorentzVector('ROOT::Math::PxPyPzE4D<double>')(px,py,pz,E)
            p4partcmjet = VectorUtil.boost(p4part, cmJet)
            # p4partcmjet are the jet generator particle 4-momenta in the jet center-of-mass
            p4jetcm += p4partcmjet
        print(" cmredo:",p4jetcm.mass(),p4jetcm.px(),p4jetcm.py(),p4jetcm.pz())
        iptot+=tree.nGenJetParticles[ijet]


# In[ ]:




