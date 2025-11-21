import torch as th
import torch.nn as nn

class KcatModel(nn.Module):
    def __init__(self, rate=0.0, device="cuda:0"):
        super(KcatModel, self).__init__()

        self.prot_norm = nn.BatchNorm1d(1024).to(device)
        self.molt5_norm = nn.BatchNorm1d(768).to(device)
        self.decoder = nn.Sequential(nn.Linear(1959, 256), nn.BatchNorm1d(256), nn.Dropout(p=rate), nn.ReLU(),).to(device)

        self.out = nn.Sequential(nn.Linear(256, 1)).to(device)

    def forward(self, ezy_feats, sbt_feats):

        prot_feats = self.prot_norm(ezy_feats)
        molt5_feats = self.molt5_norm(sbt_feats[:, :768])
        macc_feats = sbt_feats[:, 768:]
        cplx_feats = th.cat([prot_feats, molt5_feats, macc_feats], axis=1)
        feats = self.decoder(cplx_feats) # 512

        out = self.out(feats)

        return out, feats
    
class KmModel(nn.Module):
    def __init__(self, rate=0.0, device="cuda:0"):
        super(KmModel, self).__init__()

        self.prot_norm = nn.BatchNorm1d(1024).to(device)
        self.molt5_norm = nn.BatchNorm1d(768).to(device)
        self.decoder = nn.Sequential(nn.Linear(1959, 256), nn.BatchNorm1d(256), nn.Dropout(p=rate), nn.ReLU(),).to(device)

        self.out = nn.Sequential(nn.Linear(256, 1)).to(device)

    def forward(self, ezy_feats, sbt_feats):

        prot_feats = self.prot_norm(ezy_feats)
        molt5_feats = self.molt5_norm(sbt_feats[:, :768])
        macc_feats = sbt_feats[:, 768:]
        cplx_feats = th.cat([prot_feats, molt5_feats, macc_feats], axis=1)
        feats = self.decoder(cplx_feats) # 512

        out = self.out(feats)

        return out, feats

class ActivityModel(nn.Module):
    def __init__(self, rate=0.0, alpha=0.4, device="cuda:0"):
        super(ActivityModel, self).__init__()
        self.alpha = alpha

        self.kcat_model = KcatModel(device=device).to(device)
        self.Km_model = KmModel(device=device).to(device)

        self.prot_norm = nn.BatchNorm1d(1024).to(device)
        self.molt5_norm = nn.BatchNorm1d(768).to(device)

        self.decoder = nn.Sequential(nn.Linear(1959, 256), nn.BatchNorm1d(256), nn.Dropout(p=rate), nn.ReLU()).to(device)
        self.attn = nn.Sequential(nn.Linear(256, 256), nn.Softmax(dim=1)).to(device)
        self.out = nn.Linear(256, 1).to(device)

    def forward(self, ezy_feats, sbt_feats):

        pred_kcat, _ = self.kcat_model(ezy_feats, sbt_feats)
        pred_Km, _ = self.Km_model(ezy_feats, sbt_feats)

        #out = th.log10(th.pow(10, pred_kcat) / th.pow(10, pred_Km))
        pred_activity_1 = pred_kcat - pred_Km  # log(kcat/Km) = log(kcat) - log(Km)

        ezy_feats = self.prot_norm(ezy_feats)
        molt5_feats = self.molt5_norm(sbt_feats[:, :768])
        macc_feats = sbt_feats[:, 768:]
        cplx_feats = th.cat([ezy_feats, molt5_feats, macc_feats], axis=1)
        feats = self.decoder(cplx_feats)
        attn_score = self.attn(feats)
        attn_feats = attn_score * feats

        pred_activity_2 = self.out(attn_feats)
        pred_activity = pred_activity_1.detach() * (1 - self.alpha) + pred_activity_2 * self.alpha

        return pred_kcat, pred_Km, pred_activity
