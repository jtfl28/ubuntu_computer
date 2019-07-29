# -*- coding: utf-8 -*-
"""
Created on Fri Jul  5 12:39:06 2019

@author: jtfl2
"""

from __future__ import print_function
#%matplotlib inline
import pdb_parser as pp
import pandas as pd
from torch import tensor
import torch
import os
import numpy as np
import random
import torch.nn as nn
import torch.nn.functional as F
from torch.autograd import Variable
import torch.optim as optim


# Set random seem for reproducibility
manualSeed = 999
#manualSeed = random.randint(1, 10000) # use if you want new results
print("Random Seed: ", manualSeed)
random.seed(manualSeed)
torch.manual_seed(manualSeed)

d = os.getcwd()
dataroot = d + '\\made_adj'
# Number of workers for dataloader
workers = 2

# Batch size during training
batch_size = 128

# Spatial size of training images. All images will be resized to this
#   size using a transformer.
done = []
pep_length = 12
maxV = pp.find_max(pep_length)
image_size = maxV[0]
r = torch.randn(12, 12)

#def Encoder(N, P, pep_length):
#    l1 = Variable(torch.randn(len(N), int(len(N)/2))
#    l1 = Variable(torch.randn(int(len(N)/2), int(len(N)/4))
#    l1 = Variable(torch.randn(int(len(N)/4), pep_length)
#    l1 = Variable(torch.randn(int(len(N)/2), int(len(N)/2))
#    l1 = Variable(torch.randn(int(len(N)/2), int(len(N)/2))
#    l1 = Variable(torch.randn(lint(len(N)/2), ilen(N))

    

class Encoder(nn.Module):
    S_weights = nn.ParameterList([])
    N_weights = nn.ParameterList([])
    Zs = None
    Zn = None
    def __init__(self, image_size, pep_length):
        super(Encoder, self).__init__()
        self.S_weights.append(Variable(torch.randn(image_size, int(image_size/2))))
        self.S_weights.append(Variable(torch.randn(int(image_size/2), int(image_size/4))))
        self.S_weights.append(Variable(torch.randn(int(image_size/4), pep_length)))
        self.S_weights.append(Variable(torch.randn(pep_length, int(image_size/4))))
        self.S_weights.append(Variable(torch.randn(int(image_size/4), int(image_size/2))))
        self.S_weights.append(Variable(torch.randn(int(image_size/2), image_size)))
        self.N_weights.append(Variable(torch.randn(int(image_size/2), int(image_size))))
        self.N_weights.append(Variable(torch.randn(int(image_size/4), int(image_size/2))))
        self.N_weights.append(Variable(torch.randn(pep_length, int(image_size/4))))
        self.N_weights.append(Variable(torch.randn(int(image_size/4), pep_length)))
        self.N_weights.append(Variable(torch.randn(int(image_size/2), int(image_size/4))))
        self.N_weights.append(Variable(torch.randn(int(image_size/2), image_size)))
        nn.Parameter(self.S_weights)
        nn.Parameter(self.N_weights)
    def forward(self, N, S):
        H = N
        P = S
        for i in range(len(self.S_weights)):
            D = torch.diag(torch.sum(H, dim = 1))
            temp = torch.mm(torch.pow(D, (-0.5)), S)
            temp = torch.mm(temp, torch.pow(D, (-0.5)))
            temp = torch.mm(N_weights[i], temp)
            H = F.relu(torch.mm(temp, H)
            temp = torch.mm(torch.pow(torch.t(S_weights[i]), 0.5), P)
            P = F.relu(torch.mm(temp, torch.pow(S_weights[i], 0.5))

   


#for i in range(100):
#    
#    pdb = random.choice(os.listdir(d + '\\pdb'))
#    if len(done) != len(os.listdir(d + '\\pdb')):
#        while pdb in done:
#            pdb = random.choice(os.listdir(d + '\\pdb'))
#    else:
#        done = []
#    done.append(pdb)
#    found = pp.make_structure(d + '\\pdb\\' + i, 12, maxValues = maxV)
#    for dic in found:  
#        
#        N = torch.autograd.variable(tensor(dic['index'] , requires_grad = False))
#        A1 = torch.autograd.variable(tensor(np.multiply(dic['primary'],dic['secondary']), requires_grad = False))
#        A2 = torch.autograd.variable(tensor(dic['secondary'], requires_grad = False))
#
#
#class Generator(nn.Module):
#    def __init__(self):
#        super(Generator, self).__init__()
#        w1 = torch.random(())
#        
#    def forward(self, input):
#        return self.main(input)

        
