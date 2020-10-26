import argparse
import glob
import os
from rdkit import Chem
from rdkit.Chem import QED
from yaml import load, Loader
from collections import ChainMap

class Library:

    def __init__(self, path, filters=None, save=False):
        #Do we really need absolute path, plus it was done before on the parse_args()
        self.path = os.path.abspath(path) 
        #self.path = path
        self.sd_files = self._retrieve_files(self.path) #Do not pass self.pass u can access it within the method 
        #self.sd_files = self._retrieve_files()
        self.fragments = self._read_mols(self.sd_files) #We really need this attribute self.fragments, and no need to pass self.sd_files u can accessit within method
        self.fragments_dum = [Fragment(f) for f in self.fragments]  
        #self.fragments_dum = [Fragment(f) for f in self._read_mols()] 
        self.filters = filters #That should not be here
        self.save = True #That should not be here
        #self._parsed_filters = None
        #self._filtered_fragments = None
        
        ###Better a method to be used externally -> adds more flexibility####
        if self.filters:
            self.parsed_filters = self._load_filters()
            self.filtered_fragments = self._apply_filters()
        """
        def filter(self):
            self.parsed_filters = self._load_filters()
            self.filtered_fragments = self._apply_filters()
        """
    

            
        ###Better a method to be used externally -> adds more flexibility#####
        if self.save:
            self._save_sd(self.filtered_fragments)
            
        """
        def save(self, output='filtered_molecules.sdf'):
            writer = Chem.SDWriter(output)
            for mol in fragments:
                writer.write(mol)
            return output
        """

    #def _retrieve_files(self):
    def _retrieve_files(self, path):
        #sdf_path = os.path.join(self.path, "*.sdf")
        sdf_path = os.path.join(path, "*.sdf")
        sd_files = glob.glob(sdf_path)
        return sd_files

    #def _read_mols(self):
    def _read_mols(self, files):
        """All this function could be a nice one liner!! (line below)"""
        #output = [mol for mol in Chem.SDMolSupplier(file, removeHs=False) for file in self.sd_files if mol]
        #return output
        output = []
        for file in files:
            mols = Chem.SDMolSupplier(file, removeHs=False)
            #output = [mol for mol in mols for file in files if mol]
            for mol in mols:
                if mol:
                    output.append(mol)
        return output

    def _load_filters(self):
        #Why do we need to do the absolute path???
        filters_file = os.path.abspath(self.filters)

        with open(filters_file, "r") as f:
            yaml = load(f, Loader=Loader)                                                                                                                                               
            filters = yaml["filters"]

        for f in filters:
            """make a comment explaining this"""
            for key, value in f.items():
                try:
                    f[key] = [v.strip() for v in value.split("-")]
                except:
                    f[key] = value
            #print(f)
        filters = dict(ChainMap(*filters))
        #print(filters)
        ####Do we need this???#####
        for value in filters.values():
            for v in value:
                #print(v)
                pass
                
        #print('filters: ',filters)
        return filters

    def _apply_filters(self): # TO DO
       
        filtered = []
        dummy = Fragment(**self.parsed_filters)
        #print('dummy: ',dummy.mw)
        i=0
        #print(self.fragments_dum)  
        for frag in self.fragments_dum:
            
            #print(frag.mw, frag.logP, frag.hbd,frag.hba,frag.psa,frag.rotb,frag.arom)
            if float(dummy.mw[0])< frag.mw < float(dummy.mw[1]):
                if float(dummy.logP[0]) < frag.logP < float(dummy.logP[1]):
                    if int(dummy.hbd[0]) < frag.hbd < int(dummy.hbd[1]):
                        if int(dummy.hba[0]) < frag.hba < int(dummy.hba[1]):
                            if int(dummy.psa[0]) < frag.psa < int(dummy.psa[1]):
                                if int(dummy.rotb[0]) <=  frag.rotb <=  int(dummy.rotb[1]):
                                    if bool(dummy.arom[0]) ==  frag.arom:
                                        filtered.append(self.fragments[i])
            i+=1
        
        #print(filtered)
        return filtered # temporary
    
    ###No need of this method if u got a public one above
    def _save_sd(self, fragments):
        writer = Chem.SDWriter('filtered_molecules.sdf')
        #print(fragments)
        for mol in fragments:
            #print(type(mol)) 
            writer.write(mol)
 
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
        
    #Should be method __str__ have a look at: https://www.youtube.com/watch?v=5cvM-crlDvg&ab_channel=CoreySchafer
    def __repr__(self):
        return "Molecule {name}\nMW = {mw}\nlogP = {logp}\n".format(name=self.molecule_name, mw=self.mw, logp=self.logP)


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--dir", type=str, required=True, help="Directory with fragment SD files.")
    parser.add_argument("--save", type=str, required=False, help="SD file name, if you want to save the filtered output.")
    parser.add_argument("--filters", type=str, required=False, help="YML file with filters")
    args = parser.parse_args()
    #Why absolute path????
    return os.path.abspath(args.dir), args.save, args.filters

####START CHECKING FROM HERE#####
def main():
    path, save, filters = parse_args()
    #I do not see a point of having a class that runs everything
    # when is instanciated, let's do some methods and decide the user whether
    #or not to run them later
    lib = Library(path)
    #lib = Library(Path, filters, save)
    #if filters:
        # lib.filters
    #if save:
        #lib.save   
    for f in lib.filtered_fragments:
        #print(f)
        pass
        

if __name__ == "__main__":
    main()

