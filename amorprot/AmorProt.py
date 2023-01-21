import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import MACCSkeys


class AmorProt:
    
    def __init__(self, maccs=True, ecfp4=True, ecfp6=True, rdkit=True):
        
        self.AA_dict = {'G':'NCC(=O)O', 
                        'A':'N[C@@]([H])(C)C(=O)O', 
                        'R':'N[C@@]([H])(CCCNC(=N)N)C(=O)O', 
                        'N':'N[C@@]([H])(CC(=O)N)C(=O)O', 
                        'D':'N[C@@]([H])(CC(=O)O)C(=O)O', 
                        'C':'N[C@@]([H])(CS)C(=O)O', 
                        'E':'N[C@@]([H])(CCC(=O)O)C(=O)O', 
                        'Q':'N[C@@]([H])(CCC(=O)N)C(=O)O', 
                        'H':'N[C@@]([H])(CC1=CN=C-N1)C(=O)O', 
                        'I':'N[C@@]([H])(C(CC)C)C(=O)O', 
                        'L':'N[C@@]([H])(CC(C)C)C(=O)O', 
                        'K':'N[C@@]([H])(CCCCN)C(=O)O', 
                        'M':'N[C@@]([H])(CCSC)C(=O)O', 
                        'F':'N[C@@]([H])(Cc1ccccc1)C(=O)O', 
                        'P':'N1[C@@]([H])(CCC1)C(=O)O', 
                        'S':'N[C@@]([H])(CO)C(=O)O', 
                        'T':'N[C@@]([H])(C(O)C)C(=O)O', 
                        'W':'N[C@@]([H])(CC(=CN2)C1=C2C=CC=C1)C(=O)O', 
                        'Y':'N[C@@]([H])(Cc1ccc(O)cc1)C(=O)O', 
                        'V':'N[C@@]([H])(C(C)C)C(=O)O'}
        
        self.maccs = maccs
        self.ecfp4 = ecfp4
        self.ecfp6 = ecfp6
        self.rdkit = rdkit
        
        if self.maccs:
            self.maccs_dict = {}
            for aa in self.AA_dict.keys():
                mol = Chem.MolFromSmiles(self.AA_dict[aa])
                self.maccs_dict[aa] = np.array(MACCSkeys.GenMACCSKeys(mol)).tolist()
        
        if self.ecfp4:
            self.ecfp4_dict = {}
            for aa in self.AA_dict.keys():
                mol = Chem.MolFromSmiles(self.AA_dict[aa])
                self.ecfp4_dict[aa] = np.array(AllChem.GetMorganFingerprintAsBitVect(mol, radius=2, nBits=1024)).tolist()
        
        if self.ecfp6:
            self.ecfp6_dict = {}
            for aa in self.AA_dict.keys():
                mol = Chem.MolFromSmiles(self.AA_dict[aa])
                self.ecfp6_dict[aa] = np.array(AllChem.GetMorganFingerprintAsBitVect(mol, radius=3, nBits=1024)).tolist()
        
        if self.rdkit:
            self.rdkit_dict = {}
            for aa in self.AA_dict.keys():
                mol = Chem.MolFromSmiles(self.AA_dict[aa])
                self.rdkit_dict[aa] = np.array(AllChem.RDKFingerprint(mol)).tolist()
    
    def fingerprint(self, seq):
        pos = np.arange(1, len(seq)+1)/(len(seq))
        
        arrays = []
        
        if self.maccs:
            maccs_list = []
            for i in range(len(seq)):
                aa = seq[i]
                maccs_list.append((((np.sin(pos[i]/10))/10)+0.85)*np.array(self.maccs_dict[aa]))
            maccs_array = np.array(maccs_list, dtype=np.float32)
            arrays.append(maccs_array)
        
        if self.ecfp4:
            ecfp4_list = []
            for i in range(len(seq)):
                aa = seq[i]
                ecfp4_list.append((((np.sin(pos[i]/10))/10)+0.85)*np.array(self.ecfp4_dict[aa]))
            ecfp4_array = np.array(ecfp4_list, dtype=np.float32)
            arrays.append(ecfp4_array)
        
        if self.ecfp6:
            ecfp6_list = []
            for i in range(len(seq)):
                aa = seq[i]
                ecfp6_list.append((((np.sin(pos[i]/10))/10)+0.85)*np.array(self.ecfp6_dict[aa]))
            ecfp6_array = np.array(ecfp6_list, dtype=np.float32)
            arrays.append(ecfp6_array)
        
        if self.rdkit:
            rdkit_list = []
            for i in range(len(seq)):
                aa = seq[i]
                rdkit_list.append((((np.sin(pos[i]/10))/10)+0.85)*np.array(self.rdkit_dict[aa]))
            rdkit_array = np.array(rdkit_list, dtype=np.float32)
            arrays.append(rdkit_array)
        
        profp = np.concatenate(arrays, axis=1)
        
        profp = profp.sum(axis=0)
        profp = profp/profp.max()
        
        return profp
