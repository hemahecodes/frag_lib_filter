import argparse
import glob
import os
from rdkit import Chem
from rdkit.Chem import QED
from yaml import load, Loader
from collections import ChainMap

class Library:

    def __init__(self, path, filters=None, save=False): 
        self.path = path # not using absolute path since we defined path on parse_args()
        self.filters = filters 
        self.sd_files = self._retrieve_files() # not passing self.path since we access it within method
        # deleted self.fragments since we can omit it
        self.fragments_dum = [Fragment(f) for f in self._read_mols()] 
        self.mols = self._read_mols()        

    def filters_f(self): # transformed condition into method
        self.parsed_filters = self._load_filters()
        self.filtered_fragments = self._apply_filters()

    def save(self, output='filtered_molecules.sdf'): # transformed condition into method
        writer = Chem.SDWriter(output)
        for mol in self.filtered_fragments:
            writer.write(mol)
        return output
   

    def _retrieve_files(self): # deleted argument path since we can access it within method
        sdf_path = os.path.join(self.path, "*.sdf")
        sd_files = glob.glob(sdf_path)
        return sd_files

    
    def _read_mols(self): # deleted argument files since we can access it within method
        for file in self.sd_files:
            mols = Chem.SDMolSupplier(file, removeHs=False)
            output = [mol for mol in mols if mol] # loop into a single-line loop
            return output

    def _load_filters(self):
        
        with open(self.filters, "r") as f:
            yaml = load(f, Loader=Loader)                                                                                                                                               
            filters = yaml["filters"]

        for f in filters:
            """make a comment explaining this"""
            for key, value in f.items():
                try:
                    f[key] = [v.strip() for v in value.split("-")]
                except:
                    f[key] = value
        filters = dict(ChainMap(*filters))
     
        return filters

    def _apply_filters(self):
       
        filtered = []        
        i=0  
        for frag in self.fragments_dum:
            if float(self.parsed_filters['mw'][0])< frag.mw < float(self.parsed_filters['mw'][1]):
                if float(self.parsed_filters['logP'][0]) < frag.logP < float(self.parsed_filters['logP'][1]):
                    if int(self.parsed_filters['hbd'][0]) < frag.hbd < int(self.parsed_filters['hbd'][1]): 
                        if int(self.parsed_filters['hba'][0]) < frag.hba < int(self.parsed_filters['hba'][1]):
                            if int(self.parsed_filters['psa'][0]) < frag.psa < int(self.parsed_filters['psa'][1]):
                                if int(self.parsed_filters['rotb'][0]) <=  frag.rotb <=  int(self.parsed_filters['rotb'][1]):
                                    if bool(self.parsed_filters['arom'][0]) ==  frag.arom:
                                        filtered.append(self.mols[i])
            i+=1
        
             
        return filtered 

     # deleted _save_sd method 

class Fragment():

    def __init__(self, mol=None, mw=None, logP=None, hba=None, hbd=None, psa=None, rotb=None, arom=None):
        self.molecule = mol if mol else None
        self.molecule_name = mol.GetProp('_Name') if mol else None
        self._qed = QED.properties(self.molecule) if mol else None
        self.mw = mw if not mol else self._qed.MW
        self.logP = logP if not mol else self._qed.ALOGP
        self.hba = hba if not mol else self._qed.HBA
        self.hbd = hbd if not mol else self._qed.HBD
        self.psa = psa if not mol else self._qed.PSA
        self.rotb = rotb if not mol else self._qed.ROTB
        self.arom = arom if not mol else self._qed.AROM
        
    def __str__(self): #changed to str to ensure readability
        return "Molecule {name}\nMW = {mw}\nlogP = {logp}\n".format(name=self.molecule_name, mw=self.mw, logp=self.logP)


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--dir", type=str, required=True, help="Directory with fragment SD files.")
    parser.add_argument("--save", type=str, required=False, help="SD file name, if you want to save the filtered output.")
    parser.add_argument("--filters", type=str, required=False, help="YML file with filters")
    args = parser.parse_args()
    
    return os.path.abspath(args.dir), args.save, args.filters

def main():
    path, save, filters = parse_args()
    lib = Library(path, filters, save)
    if filters:
        lib.filters_f() 
    if save:
        lib.save() 
    

if __name__ == "__main__":
    main()

