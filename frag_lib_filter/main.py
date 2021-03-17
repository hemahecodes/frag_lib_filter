import argparse
import glob
import os
import yaml
from rdkit import Chem
from rdkit.Chem import QED
from tqdm import tqdm
from yaml import load, Loader
from collections import ChainMap
import time
from multiprocessing import Pool
from functools import partial
class Library:

    def __init__(self, path, filters=None,n_workers=100, save=False): 
        self.path = path # not using absolute path since we defined path on parse_args()
        self.filters = filters
        self.n_workers=n_workers
        self.errors=0 
        self.sd_files = self._retrieve_files() # not passing self.path since we access it within method
        #self.splitted=[]
        #self.sd_files= self.split_in_chunks(self.sd_files)
        #tqdm(self.parallelize(self.main,self.sd_files,self.n_workers))
        self.main(self.sd_files)
    def parallelize(self,func,iterable,n_workers, **kwargs):
        f=partial(func,**kwargs)
        if  n_workers> 1:
            with Pool(n_workers) as pool:
                r = list(tqdm(pool.map(func,iterable)))
                pool.close()
                pool.join()
            return 0
        else:
            return list(map(f,iterable))

    def split_in_chunks(self,sdf_files,n_chunks=25000):
        output_folder=f"split_sdfs_{n_chunks}"
        output_files=[]
        for sdf_file in sdf_files:
            #self.start_time2=time.time()
            print("Processing file", sdf_file)
            #print(" --- %s seconds ---" % (time.time()-start_time))
            name, ext = os.path.splitext(sdf_file)
            name = name.replace("/","_")
            mols = 0
            splits = 1
            try:
                if not os.path.exists(output_folder):
                    os.mkdir(output_folder)
            except:
                pass
            fw = open(os.path.join(output_folder,name+f"_{splits}"+ext),'w')
            output_files.append(os.path.join(output_folder,name+f"_{splits}"+ext))
            with open(sdf_file) as f:
                for line in f:
                    fw.write(line)
                    if line.startswith("$$$$"):
                        mols+=1
                    if mols == n_chunks:
                        mols=0
                        fw.close()
                        splits+=1
                        fw=open(os.path.join(output_folder, name+f"_{splits}"+ext),'w')
                        output_files.append(os.path.join(output_folder, name+f"_{splits}"+ext))
        return output_files
               
    def main(self, files_sdf):
        for file_sdf in files_sdf:
            # deleted self.fragments since we can omit it
            mols= Chem.SDMolSupplier(file_sdf,removeHs=False)
            name = os.path.basename(file_sdf).rsplit(".")[0]
            final_name="%s_filtered.sdf" % name
            for mol in mols:
                self.molecule = mol
                self.fragments_dum = [Fragment(mol)]
                self.filters_f()
                #print('Saving molecule:', file_sdf)
                #print(" --- %s seconds ---" % (time.time()-self.start_time2))
            
                self.save(output='frag_lib_filtered.sdf') 

    def filters_f(self): # transformed condition into method
        self.parsed_filters = self._load_filters()
        print('parsed filters:', self.parsed_filters)
        self.filtered_fragments = self._apply_filters()

    def save(self, output='filtered_molecules.sdf'): # transformed condition into method
       output= open(output,'a')
       writer = Chem.SDWriter(output)
       for mol in self.filtered_fragments:
           writer.write(mol)
        
   

    def _retrieve_files(self): # deleted argument path since we can access it within method
        sdf_path = os.path.join(self.path, "*.sdf")
        sd_files = glob.glob(sdf_path)
        return sd_files

    def _load_filters(self):
        with open(self.filters, 'r') as filters_file: 
            yaml = load(filters_file, Loader=Loader)
            filters=yaml["filters"]
        filters_dict={}
        for k in filters:
            print(k)
            for key, value in filters[k].items():
                try:
                    filters_dict[k].append(value)
                except:
                    filters_dict[k] = []
                    filters_dict[k].append(value)       
        #filters = dict(ChainMap(*filters))
        
        return filters_dict

    def _apply_filters(self):
        filtered = []          
        for frag in self.fragments_dum:
            try:
                if float(self.parsed_filters['mw'][0])< frag.mw < float(self.parsed_filters['mw'][1]):
                    if float(self.parsed_filters['logP'][0]) < frag.logP < float(self.parsed_filters['logP'][1]):
                        if int(self.parsed_filters['hbd'][0]) < frag.hbd < int(self.parsed_filters['hbd'][1]):
                            if int(self.parsed_filters['hba'][0]) < frag.hba < int(self.parsed_filters['hba'][1]):
                                if int(self.parsed_filters['psa'][0]) < frag.psa < int(self.parsed_filters['psa'][1]):
                                    if int(self.parsed_filters['rotb'][0]) <=  frag.rotb <=  int(self.parsed_filters['rotb'][1]):
                                        if bool(self.parsed_filters['arom'][0]) ==  frag.arom:
                                            filtered.append(self.molecule)
                                            self.errors+=1
            except:               
               self.errors+=1
        
             
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
    parser.add_argument("--filters", type=str, required=False, help="YML file with filters")
    parser.add_argument("--n_workers", type = int, required = False, help = "Number of CPUs")
    args = parser.parse_args()
    
    return os.path.abspath(args.dir), args.filters, args.n_workers

def main():
    path, filters, n_workers= parse_args()
    lib = Library(path, filters,n_workers)
    errors = lib.errors
    print("Number of wrong molecules: %s" % (lib.errors))


if __name__ == "__main__":
    start_time=time.time()
    main()
    print(" --- %s seconds ---" % (time.time()-start_time))
