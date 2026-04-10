import os
import pandas as pd
from sklearn.decomposition import PCA
import torch
import torch.nn.functional as F
import numpy as np
from torch_geometric.data import HeteroData
from torch_geometric.nn import SAGEConv, HeteroConv, Linear
from torch_geometric.transforms import ToUndirected
from sklearn.metrics import roc_auc_score, average_precision_score
import torch.optim as optim
from torch.optim.lr_scheduler import CosineAnnealingWarmRestarts
from torch.amp import autocast, GradScaler
import argparse
import warnings
import random
warnings.filterwarnings('ignore')

parser = argparse.ArgumentParser(description="GRN Application: Predict All TF-TG Pairs")
#
parser.add_argument('--weight_decay', type=float, default=2e-3)
parser.add_argument('--dropout_rate', type=float, default=0.45, help="Dropout")
parser.add_argument('--edge_dropout_rate', type=float, default=0.55, help="edge drop")
parser.add_argument('--learning_rate', type=float, default=3e-3)
parser.add_argument('--gnn_hidden', type=int, default=512)
parser.add_argument('--gnn_layers', type=int, default=3)
parser.add_argument('--gnn_out', type=int, default=32)
parser.add_argument('--use_gradient_clip', action='store_true')
parser.add_argument('--max_grad_norm', type=float, default=1)
parser.add_argument('--max_epochs', type=int, default=400)
parser.add_argument('--patience', type=int, default=40)
# pca and top k
parser.add_argument('--pca_ratio', type=float, default=0.15, help="PCA ratio")
parser.add_argument('--pca_max_components', type=int, default=32, help="PCA max")
parser.add_argument('--pca_min_components', type=int, default=8, help="PCA min")
parser.add_argument('--topk_ratio', type=float, default=0.02, help="TopK ratio")
parser.add_argument('--topk_max', type=int, default=6, help="TopK max")
parser.add_argument('--topk_min', type=int, default=2, help="TopK min")
# apply
parser.add_argument('--application_data_path', type=str, default="./", help="data path")
parser.add_argument('--output_path', type=str, default="./application_results", help="path to save result")
args = parser.parse_args()

os.makedirs(args.output_path, exist_ok=True)
# ========== 
torch.manual_seed(42)
np.random.seed(42)
random.seed(42)
torch.cuda.manual_seed_all(42)
device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')

tfsim1 = pd.read_csv(f"./tf1000gocos.txt", sep="\t", index_col=0)
tfsim2 = pd.read_csv(f"./tf1000pearson.txt", sep="\t", index_col=0)
tfsim3 = pd.read_csv(f"./tf1000_max.txt", sep="\t", index_col=0)
tfsim = pd.concat([tfsim1,tfsim2,tfsim3], axis=1)
tgsim1 = pd.read_csv(f"./tg1000gocos.txt", sep="\t", index_col=0)
tgsim2 = pd.read_csv(f"./tg1000pearson.txt", sep="\t", index_col=0)
tgsim = pd.concat([tgsim1,tgsim2], axis=1)
#
trainedge = pd.read_csv(
    os.path.join(args.application_data_path, "train_application_1000.txt"),
    sep="\t", index_col=None, header=None, names=['source', 'target', 'value']
)
valedge = pd.read_csv(
    os.path.join(args.application_data_path, "val_application_1000.txt"),
    sep="\t", index_col=None, header=None, names=['source', 'target', 'value']
)
predict_pairs = pd.read_csv(
    os.path.join(args.application_data_path, "all_pairs_to_predict_1000.txt"),
    sep="\t", index_col=None, header=None, names=['source', 'target']
)

tfs = tfsim1.index.tolist()  
tgs = tgsim1.index.tolist()
tf_to_idx = {f"tf_{tf}": idx for idx, tf in enumerate(tfs)}
tg_to_idx = {f"tg_{tg}": idx for idx, tg in enumerate(tgs)}

# ===================
def similarity_to_features(sim_matrix):
    sim_matrix = sim_matrix.values
    n_samples = sim_matrix.shape[0]
    n_components = min(int(n_samples * args.pca_ratio), args.pca_max_components)
    n_components = max(n_components, args.pca_min_components)
    return torch.tensor(PCA(n_components).fit_transform(sim_matrix), dtype=torch.float)

def create_topk_edges(sim_matrix):
    n_nodes = sim_matrix.shape[0]
    k = min(max(args.topk_min, int(args.topk_ratio * n_nodes)), args.topk_max)
    sim = sim_matrix.values.copy()
    np.fill_diagonal(sim, -np.inf)
    edges = []
    for i in range(n_nodes):
        topk_indices = np.argpartition(-sim[i], k)[:k]
        edges.extend([(i, int(j)) for j in topk_indices])
    return torch.tensor(edges, dtype=torch.long).t().contiguous()

def create_edges(edge_df):
    pos_edges, neg_edges = [], []
    for _, row in edge_df.iterrows():
        src_idx = tf_to_idx.get(f"tf_{row['source']}")
        dst_idx = tg_to_idx.get(f"tg_{row['target']}")
        if src_idx is not None and dst_idx is not None:
            if row['value'] == 1:
                pos_edges.append((src_idx, dst_idx))
            else:
                neg_edges.append((src_idx, dst_idx))
    
    pos_tensor = torch.tensor(pos_edges, dtype=torch.long).t().contiguous() if pos_edges else torch.empty((2, 0), dtype=torch.long)
    neg_tensor = torch.tensor(neg_edges, dtype=torch.long).t().contiguous() if neg_edges else torch.empty((2, 0), dtype=torch.long)
    return pos_tensor, neg_tensor

def create_edge_label_index_and_labels(pos_edges, neg_edges):
    if pos_edges.numel() == 0:
        return neg_edges, torch.zeros(neg_edges.size(1))
    if neg_edges.numel() == 0:
        return pos_edges, torch.ones(pos_edges.size(1))
    edge_label_index = torch.cat([pos_edges, neg_edges], dim=1)
    edge_label = torch.cat([torch.ones(pos_edges.size(1)), torch.zeros(neg_edges.size(1))])
    num_edges = edge_label_index.size(1)  
    shuffle_idx = torch.randperm(num_edges)  
    shuffled_edge_label_index = edge_label_index[:, shuffle_idx]
    shuffled_edge_label = edge_label[shuffle_idx]
    return shuffled_edge_label_index, shuffled_edge_label

#
tf_features = similarity_to_features(tfsim).to(device)
tg_features = similarity_to_features(tgsim).to(device)
#
tf_edges_gocos = create_topk_edges(tfsim1).to(device)
tf_edges_pearson = create_topk_edges(tfsim2).to(device)
tf_edges_tmscore = create_topk_edges(tfsim3).to(device)
tg_edges_gocos = create_topk_edges(tgsim1).to(device)
tg_edges_pearson = create_topk_edges(tgsim2).to(device)

tf_sim_edges = torch.cat([tf_edges_gocos, tf_edges_pearson,tf_edges_tmscore], dim=1)
tf_sim_edges = torch.unique(tf_sim_edges, dim=1)
tg_sim_edges = torch.cat([tg_edges_gocos, tg_edges_pearson], dim=1)
tg_sim_edges = torch.unique(tg_sim_edges, dim=1)

# 
train_pos_edges, train_neg_edges = create_edges(trainedge)
val_pos_edges, val_neg_edges = create_edges(valedge)
train_edge_label_index, train_edge_label = create_edge_label_index_and_labels(train_pos_edges, train_neg_edges)
val_edge_label_index, val_edge_label = create_edge_label_index_and_labels(val_pos_edges, val_neg_edges)

train_edge_label_index, train_edge_label = train_edge_label_index.to(device), train_edge_label.to(device)
val_edge_label_index, val_edge_label = val_edge_label_index.to(device), val_edge_label.to(device)

# 
data = HeteroData()
data['tf'].x = tf_features
data['tg'].x = tg_features
data['tf', 'tf_sim', 'tf'].edge_index = tf_sim_edges
data['tg', 'tg_sim', 'tg'].edge_index = tg_sim_edges
data['tf', 'regulates', 'tg'].edge_index = train_pos_edges.to(device)
data = ToUndirected()(data).to(device)

# ==========
class HeteroGNN(torch.nn.Module):
    def __init__(self, hidden_channels, out_channels, num_layers=3, 
                 dropout_rate=0.45, edge_dropout_rate=0.55):
        super().__init__()
        self.dropout_rate = dropout_rate
        self.edge_dropout_rate = edge_dropout_rate
        
        self.convs = torch.nn.ModuleList()
        self.convs.append(self._create_conv_layer((-1, -1), hidden_channels))
        for _ in range(num_layers - 1):
            self.convs.append(self._create_conv_layer(hidden_channels, hidden_channels))
        self.lin = torch.nn.ModuleDict({
            'tf': Linear(hidden_channels, out_channels),
            'tg': Linear(hidden_channels, out_channels)
        })

    def _create_conv_layer(self, in_channels, out_channels):
        conv_dict = {
            ('tf', 'tf_sim', 'tf'): SAGEConv(in_channels, out_channels),
            ('tg', 'tg_sim', 'tg'): SAGEConv(in_channels, out_channels),
            ('tf', 'regulates', 'tg'): SAGEConv(in_channels, out_channels),
            ('tg', 'rev_regulates', 'tf'): SAGEConv(in_channels, out_channels)
        }
        return HeteroConv(conv_dict, aggr='sum')

    def forward(self, x_dict, edge_index_dict):
        # edge drop when train
        if self.training and self.edge_dropout_rate > 0:
            edge_index_dict = {
                key: self.dropout_edges(edge_index, self.edge_dropout_rate)
                for key, edge_index in edge_index_dict.items()
            }

        for i, conv in enumerate(self.convs):
            out_dict = conv(x_dict, edge_index_dict)
            x_dict = {key: F.elu(out_dict[key]) for key in x_dict.keys()}
            x_dict = {key: F.dropout(x, p=self.dropout_rate, training=self.training) for key, x in x_dict.items()}

        # project
        final_dict = {
            node_type: self.lin[node_type](x_dict[node_type])
            for node_type in ['tf', 'tg']
        }
        return final_dict

    def dropout_edges(self, edge_index, p):
        if p <= 0 or edge_index.numel() == 0:
            return edge_index
        mask = torch.rand(edge_index.size(1), device=edge_index.device) > p
        return edge_index[:, mask]
class MLPLinkPredictor(torch.nn.Module):
    def __init__(self, in_channels, mlp_hidden=128, dropout_rate=0.4):
        super().__init__()
        self.dropout_rate = dropout_rate  
        self.lin1 = Linear(3 * in_channels, mlp_hidden)
        self.lin2 = Linear(mlp_hidden, mlp_hidden // 2)
        self.lin3 = Linear(mlp_hidden // 2, 1)
    def forward(self, z_dict, edge_label_index):
        row, col = edge_label_index
        z_tf = z_dict['tf'][row]
        z_tg = z_dict['tg'][col]
        hadamard = z_tf * z_tg
        z = torch.cat([
            z_tf, z_tg,         
            hadamard
        ], dim=-1)
        #MLP
        z = F.dropout(F.elu(self.lin1(z)), p=self.dropout_rate, training=self.training)
        z = F.dropout(F.elu(self.lin2(z)), p=self.dropout_rate, training=self.training)
        return self.lin3(z).view(-1)
#
model = HeteroGNN(
    hidden_channels=args.gnn_hidden,
    out_channels=args.gnn_out,
    num_layers=args.gnn_layers,
    dropout_rate=args.dropout_rate,
    edge_dropout_rate=args.edge_dropout_rate
).to(device)

predictor = MLPLinkPredictor(
    in_channels=args.gnn_out,
    mlp_hidden=args.gnn_hidden//2,
    dropout_rate=args.dropout_rate
).to(device)

optimizer = optim.AdamW(
    list(model.parameters()) + list(predictor.parameters()),
    lr=args.learning_rate,
    weight_decay=args.weight_decay
)
scheduler = CosineAnnealingWarmRestarts(optimizer, T_0=20, T_mult=1, eta_min=args.learning_rate * 0.05)
scaler = GradScaler() if device.type == 'cuda' else None

# ========== 
best_val_auc = 0.0
best_val_aupr = 0.0
counter = 0
best_model_state = None
best_predictor_state = None

print("begin train (application phase)...")
for epoch in range(1, args.max_epochs + 1):
    # train
    model.train()
    predictor.train()
    optimizer.zero_grad()
    
    if scaler:
        with autocast("cuda"):
            z_dict = model(data.x_dict, data.edge_index_dict)
            pred = predictor(z_dict, train_edge_label_index)
            loss = F.binary_cross_entropy_with_logits(pred, train_edge_label)
        scaler.scale(loss).backward()
        if args.use_gradient_clip:
            scaler.unscale_(optimizer)
            torch.nn.utils.clip_grad_norm_(
                list(model.parameters()) + list(predictor.parameters()),
                max_norm=args.max_grad_norm
            )
        scaler.step(optimizer)
        scaler.update()
    else:
        z_dict = model(data.x_dict, data.edge_index_dict)
        pred = predictor(z_dict, train_edge_label_index)
        loss = F.binary_cross_entropy_with_logits(pred, train_edge_label)
        loss.backward()
        if args.use_gradient_clip:
            torch.nn.utils.clip_grad_norm_(
                list(model.parameters()) + list(predictor.parameters()),
                max_norm=args.max_grad_norm
            )
        optimizer.step()
    
    scheduler.step()
    
    # val
    model.eval()
    predictor.eval()
    with torch.no_grad():
        z_dict = model(data.x_dict, data.edge_index_dict)
        val_pred = predictor(z_dict, val_edge_label_index)
        val_pred_probs = torch.sigmoid(val_pred).cpu().numpy()
        val_labels = val_edge_label.cpu().numpy()
        val_auc = roc_auc_score(val_labels, val_pred_probs)
        val_aupr = average_precision_score(val_labels, val_pred_probs)
    
    # 
    if val_auc > best_val_auc:
        best_val_auc = val_auc
        best_val_aupr = val_aupr
        counter = 0
        best_model_state = model.state_dict()
        best_predictor_state = predictor.state_dict()
    else:
        counter += 1
        if counter >= args.patience:
            print(f"Early stop at epoch {epoch} | best val AUC: {best_val_auc:.3f}")
            break
    
    #
    if epoch % 10 == 0:
        print(f"Epoch {epoch:03d} | Loss: {loss.item():.3f} | Val AUC: {val_auc:.3f} | Val AUPR: {val_aupr:.3f}")

#
model_filename = f"application_model.pth"
model_save_path = os.path.join(args.output_path, model_filename)
torch.save({
    'model_state_dict': best_model_state,
    'predictor_state_dict': best_predictor_state,
    'tf_to_idx': tf_to_idx,
    'tg_to_idx': tg_to_idx
}, model_save_path)
print(f"Best model saved to: {model_save_path}")


checkpoint = torch.load(model_save_path, map_location=device, weights_only=True)
model.load_state_dict(checkpoint['model_state_dict'])
predictor.load_state_dict(checkpoint['predictor_state_dict'])


def create_predict_edge_index(predict_pairs, tf_to_idx, tg_to_idx,device):
    predict_edges = []
    valid_pairs = []  
    for _, row in predict_pairs.iterrows():
        src_idx = tf_to_idx.get(f"tf_{row['source']}")
        dst_idx = tg_to_idx.get(f"tg_{row['target']}")
        if src_idx is not None and dst_idx is not None:
            predict_edges.append((src_idx, dst_idx))
            valid_pairs.append((row['source'], row['target']))
    edge_index = torch.tensor(predict_edges, dtype=torch.long).t().contiguous().to(device)
    return edge_index, valid_pairs

#
predict_edge_index, valid_predict_pairs = create_predict_edge_index(predict_pairs, tf_to_idx, tg_to_idx,device)
print(f"Total predict pairs: {len(predict_pairs)} | Valid pairs: {len(valid_predict_pairs)}")

#
model.eval()
predictor.eval()
with torch.no_grad():
    z_dict = model(data.x_dict, data.edge_index_dict)
    pred_logits = predictor(z_dict, predict_edge_index)
    pred_probs = torch.sigmoid(pred_logits).cpu().numpy()  

pred_result_df = pd.DataFrame({
    'TF': [pair[0] for pair in valid_predict_pairs],
    'Target': [pair[1] for pair in valid_predict_pairs],
    'pred_prob': pred_probs
})
pred_filename = f"all_predictions.csv"
pred_save_path = os.path.join(args.output_path, pred_filename)
pred_result_df.to_csv(pred_save_path, index=False, sep=",")
print(f"Full prediction results saved to: {pred_save_path}")

