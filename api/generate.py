import os
import sys

from rdkit import Chem

sys.path.append('../experiments/')
sys.path.append('../')
sys.path.append('../src/')
#from do_data_generation import generate

from experiments import do_data_generation
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('-n','--number', type=int, help='Number Gen', required=True)

if __name__ == '__main__':
    
    args = vars(parser.parse_args())
        
    num = args['number']
    
    
    saved_path ='/Users/dml/My Projects/first_web_app/backend/app/storage/molgen/data/temp'
    
    model_path = '/Users/dml/My Projects/virtual_libraries/experiments/results/Dimers_TEMP/models/40.h5'
    
    temp = 0.7
    max_len = 140
    
    
    pad_char = do_data_generation.FP.PROCESSING_FIXED['pad_char']
    temp=0.7
    start_char = do_data_generation.FP.PROCESSING_FIXED['start_char']
    end_char = do_data_generation.FP.PROCESSING_FIXED['end_char']
    indices_token = do_data_generation.FP.INDICES_TOKEN
    token_indices = do_data_generation.FP.TOKEN_INDICES
    
    model = do_data_generation.load_model(model_path)
    
    sampled =[]
    for i in range(num):
        sampled.append(do_data_generation.sample(model, temp, start_char, end_char, max_len, indices_token, token_indices))
    
    
    valids = []
    n_valid = 0
    
    for gen_smile in sampled:
        if len(gen_smile)!=0 and isinstance(gen_smile, str):
            gen_smile = gen_smile.replace(pad_char,'')
            gen_smile = gen_smile.replace(end_char,'')
            gen_smile = gen_smile.replace(start_char,'')
            
            mol = Chem.MolFromSmiles(gen_smile)
            if mol is not None: 
                cans = Chem.MolToSmiles(mol, isomericSmiles=True, canonical=True)
                if len(cans)>=1:
                    n_valid+=1
                    valids.append(cans)
                    
                    
    if n_valid!=0:
        # Now let's pruned our valid guys
        unique_set = set(valids)
        n_unique = len(unique_set)
        novo_tr = list(unique_set - set(data_training))
        n_novo_tr = len(novo_tr)
        novo_val = list(unique_set - set(data_validation))
        n_novo_val = len(novo_val)
        novo_analysis = {'n_valid': n_valid,
                            'n_unique': n_unique,
                            'n_novo_tr': n_novo_tr,
                            'n_novo_val': n_novo_val,
                            'novo_tr': novo_tr}
        
    
    
    
    temp_txt_file = os.path.join(saved_path,'novo_mol.txt')
    
    with open(temp_txt_file, 'w') as file:
        for item in novo_tr:
            file.write("%s\n" % item)
    
    print(sampled)