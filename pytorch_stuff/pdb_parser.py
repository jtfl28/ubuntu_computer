# -*- coding: utf-8 -*-
"""
Created on Tue Jun 25 14:54:46 2019

@author: jtfl2
"""

import prody as pd
import os
import numpy as np
import cleaning
import motif as m
from functools import partial
import PeptideBuilder as pb
import Bio.PDB 

d = os.getcwd()
ElementSymbols = 'H, He, Li, Be, B, C, N, O, F, Ne, Na, Mg, Al, Si, P, S, Cl, Ar, K, Ca, Sc, Ti, V, Cr, Mn, Fe, Co, Ni, Cu, Zn, Ga, Ge, As, Se, Br, Kr, Rb, Sr, Y, Zr, Nb, Mo, Tc, Ru, Rh, Pd, Ag, Cd, In, Sn, Sb, Te, I, Xe, Cs, Ba, La, Ce, Pr, Nd, Pm, Sm, Eu, Gd, Tb, Dy, Ho, Er, Tm, Yb, Lu, Hf, Ta, W, Re, Os, Ir, Pt, Au, Hg, Tl, Pb, Bi, Po, At, Rn, Fr, Ra, Ac, Th, Pa, U, Np, Pu, Am, Cm, Bk, Cf, Es, Fm, Md, No, Lr, Rf, Db, Sg, Bh, Hs, Mt, Uun, Uuu, Uub, Uut, Uuq, Uup, Uuh, Uus, Uuo'.split(', ')
aa = 'A,N,D,C,Q,E,G,H,I,L,K,M,F,P,R,S,T,W,Y,V'.split(',')
seq = ''
file = None
seq_to_pdb = {}
pdb = ''
mV = [1,1]

def find_max(length):
    global seq, file, seq_to_pdb, pdb_to_seq
    seq = 'W'*length
    out = Bio.PDB.PDBIO()
    i = pb.make_structure(seq, [180]*len(seq),[180]*len(seq))
    out.set_structure(i)
    out.save("ignore.pdb")
    cleaning.cleanATOM("ignore.pdb", out_file= d + "\\ignore", ext = '.pdb')
    file = pd.parsePDB("ignore.pdb")
    seq, seq_to_pdb, pdb_to_seq = find_seq(file)
    final = main_calc(seq)
    atom_length = len(final['primary'])
    eucl_length = final['secondary'].max()
    return atom_length, eucl_length

def find_index(full_pep, pep):
    for i in range(len(full_pep)):
        if pep == full_pep[i:i+len(pep)]:
            return [i, i+len(pep)]
    return -1

def find_seq(pdb):
    seq = ''
    seq_to_pdb = {}
    pdb_to_seq = {}
    count = 0
    for i in pdb.iterAtoms():
        if str(i.getName()) == 'N':
            seq += i.getSequence()
            seq_to_pdb[len(seq)-1] = [count]
            pdb_to_seq[count] = len(seq)-1
        else:
            seq_to_pdb[len(seq)-1].append(count)
            pdb_to_seq[count] = len(seq)-1
        count += 1
    return seq, seq_to_pdb, pdb_to_seq


def possible_bonds(i, pep, ind):
    possible = []
    possible.append(pep[i-ind[0]])
    if i != ind[0]:
        possible.append(pep[i -ind[0] -1])
    if i != ind[1] - 1:
        possible.append(pep[i - ind[0] +1])
    return possible
    
class memoize(object):
    def __init__(self, func):
        self.func = func
        self.cache = {}

    def __call__(self, *args):
        if args in self.cache:
            return self.cache[args]
        else:
            result = self.func(*args)
            self.cache[args] = result
            return result

    def __get__(self, obj, objtype):
        return partial(self.__call__, obj)
    

@memoize
def make_structure(pdb1, length, maxValues):
    global seq, file, seq_to_pdb, pdb, mV
    mV = maxValues
    pdb = pdb1
    cleaning.cleanATOM(d + '\\pdb\\' + pdb, out_file= d + '\\cleaned\\' + pdb, ext = '.pdb')
    file = pd.parsePDB(d + '\\cleaned\\' + pdb)
    seq, seq_to_pdb, pdb_to_seq = find_seq(file)
    info = m.pep_parser(seq, length)
    peps = []
    if not os.path.exists(d + '\\made_adj\\' + pdb.replace('.pdb','')):
        os.mkdir(d + '\\made_adj\\' + pdb.replace('.pdb',''))

    completed = []
    for motif in info:
        if motif == 'full':
            continue
        for pep in info[motif]:
            peps.append(pep)
            completed.append(main_calc(pep))
    return completed


def main_calc(pep):
    print(pep)
    Nterm = False
    ind = find_index(seq, pep)
    new_pep = []
    for i in range(ind[0],ind[1]):
        for j in seq_to_pdb[i]:
            new_pep.append([file[j], i])

    lengths = np.zeros((mV[0], mV[0]))
    index = np.zeros((mV[0],2))
    bonded = np.zeros((mV[0], mV[0]))

#    lengths = np.zeros((len(new_pep), len(new_pep)))
#    index = np.zeros((len(new_pep),2))
#    bonded = np.zeros((len(new_pep), len(new_pep)))
    for i in range(len(new_pep)):
        ele = float(ElementSymbols.index(str(new_pep[i][0].getName()[0]))) + 1
        if ele == 6.0:
            index[i][1] = 0.0
        elif ele == 7.0 and Nterm:
            index[i][1] = 1.0
        elif ele == 8.0:
            index[i][1] = 2.0
        if Nterm == False and new_pep[i][0].getName() == 'N':
            Nterm = True
        index[i][0] = ele#/len(ElementSymbols)
        for j in range(len(new_pep)):
            if i != j and lengths[i][j] == 0:
                lengths[i][j] = (float(pd.calcDistance(new_pep[i][0],new_pep[j][0])))
                lengths[j][i] = lengths[i][j]
            if (lengths[i][j] < 1.9) and (new_pep[i][0] != new_pep[j][0]):
                if new_pep[j][0].getSequence() in possible_bonds(new_pep[i][1], pep, ind):
                    bonded[i][j] = lengths[i][j]
                    bonded[j][i] = bonded[i][j]
    
    final = {}
    final['sequence'] = pep
    final['index'] = index
    final['secondary'] = lengths
    final['primary'] = bonded
    return final



#done = make_structure('1akj.pdb',12)
                
                
                
                
                
                
                
                
                
                
                
                
                