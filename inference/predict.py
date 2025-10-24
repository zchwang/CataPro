import torch as th
import torch.nn as nn
import pandas as pd
import numpy as np
from utils import *
from model import *
from act_model import KcatModel as _KcatModel
from act_model import KmModel as _KmModel
from act_model import ActivityModel
from torch.utils.data import DataLoader, Dataset
from argparse import RawDescriptionHelpFormatter
import argparse
import os

class EnzymeDatasets(Dataset):
    def __init__(self, values):
        self.values = values

    def __getitem__(self, idx):
        return self.values[idx]

    def __len__(self):
        return len(self.values)

def get_datasets(inp_fpath, ProtT5_model, MolT5_model):
    inp_df = pd.read_csv(inp_fpath, index_col=0)
    ezy_ids = inp_df["Enzyme_id"].values
    ezy_type = inp_df["type"].values
    ezy_keys = [f"{_id}_{t}" for _id, t in zip(ezy_ids, ezy_type)]
    sequences = inp_df["sequence"].values 
    smiles = inp_df["smiles"].values
    
    seq_ProtT5 = Seq_to_vec(sequences, ProtT5_model)
    smi_molT5 = get_molT5_embed(smiles, MolT5_model)
    smi_macc = GetMACCSKeys(smiles)
    
    feats = th.from_numpy(np.concatenate([seq_ProtT5, smi_molT5, smi_macc], axis=1)).to(th.float32)
    datasets = EnzymeDatasets(feats)
    dataloader = DataLoader(datasets)
    
    return ezy_keys, smiles, dataloader

def inference(kcat_model, Km_model, act_model, dataloader, device="cuda:0"):
    kcat_model.eval()
    Km_model.eval()
    act_model.eval()
    with th.no_grad():
        pred_list = []
        for step, data in enumerate(dataloader):
            data = data.to(device)
            ezy_feats = data[:, :1024]
            sbt_feats = data[:, 1024:]
            pred_kcat = kcat_model(ezy_feats, sbt_feats).cpu().numpy()
            pred_Km = Km_model(ezy_feats, sbt_feats).cpu().numpy()
            pred_act = act_model(ezy_feats, sbt_feats)[-1].cpu().numpy()
            pred_list.append(np.concatenate([pred_kcat, pred_Km, pred_act], axis=1))
        
        return np.concatenate(pred_list, axis=0)
            
if __name__ == "__main__":
    d = "RUN CATAPRO ..."
    parser = argparse.ArgumentParser(description=d, formatter_class=RawDescriptionHelpFormatter)
    parser.add_argument("-inp_fpath", type=str, default="sample_inp.csv",
                         help="Input (.fasta). The path of enzyme file.")
    parser.add_argument("-model_dpath", type=str, default="../models",
                         help="Input. The path of saved models.")
    parser.add_argument("-batch_size", type=int, default=64,
                        help="Input. Batch size")
    parser.add_argument("-device", type=str, default="cuda",
                        help="Input. The device: cuda or cpu.")
    parser.add_argument("-out_fpath", type=str, default="catapro_predict_score.csv",
                        help="Input. Store the predicted kinetic parameters in this file..")
    args = parser.parse_args()

    inp_fpath = args.inp_fpath
    model_dpath = args.model_dpath
    batch_size = args.batch_size 
    device = args.device
    out_fpath = args.out_fpath
    kcat_model_dpath = f"{model_dpath}/kcat_models"
    Km_model_dpath = f"{model_dpath}/Km_models"
    act_model_dpath = f"{model_dpath}/act_models"
    ProtT5_model = f"{model_dpath}/prot_t5_xl_uniref50/"
    MolT5_model = f"{model_dpath}/molt5-base-smiles2caption"

    ezy_ids, smiles_list, dataloader = get_datasets(inp_fpath, ProtT5_model, MolT5_model)
    
    pred_kcat_list = []
    pred_Km_list = []
    pred_act_list = []
    for fold in range(10):    
        kcat_model = KcatModel(device=device)
        kcat_model.load_state_dict(th.load(f"{kcat_model_dpath}/{fold}_bestmodel.pth", map_location=device))
        Km_model = KmModel(device=device)
        Km_model.load_state_dict(th.load(f"{Km_model_dpath}/{fold}_bestmodel.pth", map_location=device))
        act_model = ActivityModel(device=device)
        act_model.load_state_dict(th.load(f"{act_model_dpath}/{fold}_bestmodel.pth", map_location=device))

        pred_score = inference(kcat_model, Km_model, act_model, dataloader, device)
        pred_kcat_list.append(pred_score[:, :1])
        pred_Km_list.append(pred_score[:, 1:2])
        pred_act_list.append(pred_score[:, -1:])
    
    pred_kcat = np.mean(np.concatenate(pred_kcat_list, axis=1), axis=1, keepdims=True)
    pred_Km = np.mean(np.concatenate(pred_Km_list, axis=1), axis=1, keepdims=True)
    pred_act = np.mean(np.concatenate(pred_act_list, axis=1), axis=1, keepdims=True)
    
    final_score = np.concatenate([np.array(ezy_ids).reshape(-1, 1), np.array(smiles_list).reshape(-1, 1), pred_kcat, pred_Km, pred_act], axis=1)
    final_df = pd.DataFrame(final_score, columns=["fasta_id", "smiles", "pred_log10[kcat(s^-1)]", "pred_log10[Km(mM)]", "pred_log10[kcat/Km(s^-1mM^-1)]"])
    final_df.to_csv(out_fpath)
