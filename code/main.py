import os
import pandas as pd
from sklearn.decomposition import PCA
import torch
import torch.nn.functional as F
import numpy as np
from torch_geometric.data import HeteroData
from torch_geometric.nn import SAGEConv, HeteroConv, Linear
from torch_geometric.transforms import ToUndirected
from sklearn.metrics import roc_auc_score, average_precision_score, accuracy_score, recall_score, precision_score
import torch.optim as optim
from torch.optim.lr_scheduler import CosineAnnealingWarmRestarts
from torch.amp import autocast, GradScaler
import argparse
import warnings
import random
import time
warnings.filterwarnings('ignore')
parser = argparse.ArgumentParser(description="script")
# 
parser.add_argument('--cell_type', type=str, required=True)
parser.add_argument('--net_type', type=str, required=True)
parser.add_argument('--gene_num', type=int, required=True)
parser.add_argument('--fold', type=int, required=True)
parser.add_argument('--base_path', type=str, default="/home/bsb/mywork/data/")
# 
parser.add_argument('--weight_decay', type=float, default=2e-3)
parser.add_argument('--dropout_rate', type=float, default=0.5, help="Dropout")
parser.add_argument('--edge_dropout_rate', type=float, default=0.45, help="edge drop")
parser.add_argument('--learning_rate', type=float, default=3e-3)
parser.add_argument('--gnn_hidden', type=int, default=256)
parser.add_argument('--gnn_layers', type=int, default=3)
parser.add_argument('--gnn_out', type=int, default=16)
parser.add_argument('--use_gradient_clip', action='store_true')
parser.add_argument('--max_grad_norm', type=float, default=1)
parser.add_argument('--max_epochs', type=int, default=400)
parser.add_argument('--patience', type=int, default=30)
# 
parser.add_argument('--pca_ratio', type=float, default=0.15, help="PCA ratio")
parser.add_argument('--pca_max_components', type=int, default=32, help="PCA max")
parser.add_argument('--pca_min_components', type=int, default=8, help="PCA min")
parser.add_argument('--topk_ratio', type=float, default=0.02, help="TopK ratio")
parser.add_argument('--topk_max', type=int, default=6, help="TopK max")
parser.add_argument('--topk_min', type=int, default=2, help="TopK min")
args = parser.parse_args()

# ========== seed==========
torch.manual_seed(42)
np.random.seed(42)
random.seed(42)
torch.cuda.manual_seed_all(42)
device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
#result path
script_dir = os.path.dirname(os.path.abspath(__file__))
model_dir = os.path.join(script_dir, "../example/model")
val_dir = os.path.join(script_dir, "../example/val")
test_dir = os.path.join(script_dir, "../example/test")
os.makedirs(model_dir, exist_ok=True)
os.makedirs(val_dir, exist_ok=True)
os.makedirs(test_dir, exist_ok=True)
#dataload
task_dir = os.path.join(args.base_path, args.cell_type, args.net_type)
os.chdir(task_dir)  
start_time = time.time()
tfsim1 = pd.read_csv(f"similarity/tf{args.gene_num}gocos.txt", sep="\t", index_col=0)
tfsim2 = pd.read_csv(f"similarity/tf{args.gene_num}pearson.txt", sep="\t", index_col=0)
tfsim3 = pd.read_csv(f"similarity/tf{args.gene_num}_max.txt", sep="\t", index_col=0)
tfsim = pd.concat([tfsim1,tfsim2,tfsim3], axis=1)
tgsim1 = pd.read_csv(f"similarity/tg{args.gene_num}gocos.txt", sep="\t", index_col=0)
tgsim2 = pd.read_csv(f"similarity/tg{args.gene_num}pearson.txt", sep="\t", index_col=0)
tgsim = pd.concat([tgsim1,tgsim2], axis=1)
train_adj = pd.read_csv(f"groundtruth/adj_matrix{args.gene_num}_fold_{args.fold}.txt", sep="\t", index_col=0)
trainedge = pd.read_csv(f"groundtruth/train{args.gene_num}_fold_{args.fold}.txt", sep="\t", 
                        index_col=None, header=None, names=['source', 'target', 'value'])
valedge = pd.read_csv(f"groundtruth/val{args.gene_num}_fold_{args.fold}.txt", sep="\t", 
                      index_col=None, header=None, names=['source', 'target', 'value'])
#
tfs = train_adj.index.tolist()
tgs = train_adj.columns.tolist()
tf_to_idx = {f"tf_{tf}": idx for idx, tf in enumerate(tfs)}
tg_to_idx = {f"tg_{tg}": idx for idx, tg in enumerate(tgs)}
#data preprocess function
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

#feature
tf_features = similarity_to_features(tfsim).to(device)
tg_features = similarity_to_features(tgsim).to(device)
#sim edge
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
# ========== create heterodata ==========
data = HeteroData()
data['tf'].x = tf_features
data['tg'].x = tg_features
data['tf', 'tf_sim', 'tf'].edge_index = tf_sim_edges
data['tg', 'tg_sim', 'tg'].edge_index = tg_sim_edges
data['tf', 'regulates', 'tg'].edge_index = train_pos_edges.to(device)
data = ToUndirected()(data).to(device)

# 
class HeteroGNN(torch.nn.Module):
    def __init__(self, hidden_channels, out_channels, num_layers=2, 
                 dropout_rate=0.5, edge_dropout_rate=0.45):
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

# ========== train and val ==========
best_val_auc = 0.0
best_val_aupr = 0.0
counter = 0
best_model_state = None
best_predictor_state = None

print("begin train...")
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
        best_val_aupr=val_aupr
        counter = 0
        best_model_state = model.state_dict()
        best_predictor_state = predictor.state_dict()
    else:
        counter += 1
        if counter >= args.patience:
            print(f"stop {epoch} epochs，best val AUC: {best_val_auc:.3f}")
            break
    # print epoch
    if epoch % 10 == 0:
        print(f"epochs {epoch:03d} | loss: {loss.item():.3f} | val AUC: {val_auc:.3f} | val AUPR: {val_aupr:.3f}")
        
end_time = time.time()  # time end
total_time = end_time - start_time
print(total_time)
# ========== save model ==========
model_filename = f"v1model_{args.cell_type}_{args.net_type}_{args.gene_num}_fold{args.fold}.pth"
model_save_path = os.path.join(model_dir, model_filename)
torch.save({
    'model_state_dict': best_model_state,
    'predictor_state_dict': best_predictor_state
}, model_save_path)
print(f"model save to: {model_save_path}")

# ========== save val result ==========
val_result_filename = f"v1val_{args.cell_type}_{args.net_type}_{args.gene_num}_fold{args.fold}.csv"
val_result_path = os.path.join(val_dir, val_result_filename)
val_result_df = pd.DataFrame({
    **{k: [v] for k, v in vars(args).items()},
    'best_val_auc': [round(best_val_auc, 3)],
    'best_val_aupr': [round(best_val_aupr,3)]
})
val_result_df.to_csv(val_result_path, index=False)
print(f"val result to: {val_result_path}")

#load test
testedge = pd.read_csv(f"groundtruth/test{args.gene_num}_fold_{args.fold}.txt", sep="\t", 
                       index_col=None, header=None, names=['source', 'target', 'value'])
test_pos_edges, test_neg_edges = create_edges(testedge)
test_edge_label_index, test_edge_label = create_edge_label_index_and_labels(test_pos_edges, test_neg_edges)
test_edge_label_index, test_edge_label = test_edge_label_index.to(device), test_edge_label.to(device)

#get param
checkpoint = torch.load(model_save_path, map_location=device, weights_only=True)
model.load_state_dict(checkpoint['model_state_dict'])
predictor.load_state_dict(checkpoint['predictor_state_dict'])
# predict test
model.eval()
predictor.eval()
with torch.no_grad():
    z_dict_test = model(data.x_dict, data.edge_index_dict)
    test_pred = predictor(z_dict_test, test_edge_label_index)
    test_pred_probs = torch.sigmoid(test_pred).cpu().numpy()
    test_labels = test_edge_label.cpu().numpy()
    test_auc = roc_auc_score(test_labels, test_pred_probs)
    test_aupr = average_precision_score(test_labels, test_pred_probs)
    test_pred_binary = (test_pred_probs > 0.5).astype(int)
    test_accuracy = accuracy_score(test_labels, test_pred_binary)
    test_recall = recall_score(test_labels, test_pred_binary)
    test_precision = precision_score(test_labels, test_pred_binary) 
print(f"test AUC: {test_auc:.3f} | test AUPR: {test_aupr:.3f}")
# ========== save test ==========
test_result_filename = f"v1test_{args.cell_type}_{args.net_type}_{args.gene_num}_fold{args.fold}.csv"
test_result_path = os.path.join(test_dir, test_result_filename)
test_result_df = pd.DataFrame({
    **{k: [v] for k, v in vars(args).items()},
    'best_val_auc': [round(best_val_auc, 3)],
    'best_val_aupr': [round(best_val_aupr,3)],
    'test_auc': [round(test_auc, 3)],
    'test_aupr': [round(test_aupr, 3)],
    'test_accuracy': [round(test_accuracy, 3)],
    'test_recall': [round(test_recall, 3)],
    'test_precision': [round(test_precision, 3)],
    'time':[round(total_time,2)]
})
test_result_df.to_csv(test_result_path, index=False)
print(f"test result to: {test_result_path}")
