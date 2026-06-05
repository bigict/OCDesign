'''
method1: using the identity matrix to reduce 
'''
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.manifold import MDS
from sklearn.preprocessing import StandardScaler
from scipy.spatial.distance import squareform
from sklearn.manifold import TSNE
from sklearn.decomposition import PCA  # 用于初始化t-SNE
import seaborn as sns
from scipy.spatial import ConvexHull
from scipy.stats import gaussian_kde

def load_sequences_list(fasta_file):
    with open(fasta_file) as f:
        sequences_list = [line.strip() for line in f.readlines()]
    return sequences_list
def get_identity_matrix(sequences_list):
    def cal_identity(seq1, seq2):
        same = 0
        for i,j in zip(seq1, seq2):
            if i == j:same += 1
            else: pass
        return same/len(seq1)
        
    length = len(sequences_list)
    identity_matrix = np.zeros((length, length))
    for i in range(length):
        seq_i = sequences_list[i]
        for j in range(i, length):
            seq_j = sequences_list[j]
            identity = cal_identity(seq_i, seq_j)
            identity_matrix[i,j] = identity
            identity_matrix[j,i] = identity
    return identity_matrix
    
def reduce_symmetric_matrix(matrix, seed = 42, labels = None, perplexity = None):
    import numpy as np
    import matplotlib.pyplot as plt
    from sklearn.manifold import MDS
    from sklearn.manifold import TSNE
    
    # MDS 降维到 2 维
    mds = MDS(n_components=2, random_state=seed, dissimilarity='precomputed')
    tsne = TSNE(n_components=2, random_state=seed, perplexity=perplexity)
    
    coords = mds.fit_transform(matrix)   # 形状 (29, 2)
    coords = tsne.fit_transform(matrix)
    
    theta = np.radians(135)
    rot_matrix = np.array([[np.cos(theta), -np.sin(theta)],
                           [np.sin(theta),  np.cos(theta)]])
    
    new_coords =  coords @ rot_matrix
    print(new_coords.shape)
    return coords
    # 绘制二维散点图
    
    # plt.scatter(coords[:, 0], coords[:, 1], label = labels, cmap='viridis', s=80, edgecolors='k')

    # for i, (xi, yi, label) in enumerate(zip(coords[:, 0], coords[:, 1], labels)):
    #     plt.text(xi, yi, label, fontsize=12, ha='center', va='bottom')
    
    def plot_part(start, end, save_file, color):
        print(start, end, save_file, color)
        fig = plt.figure(figsize=fig_size)
        # df = pd.DataFrame({'x':coords[start:end, 0], 'y':coords[start:end, 1], 'label':labels[start:end]})
        df = pd.DataFrame({'x':new_coords[start:end,0], 'y':new_coords[start:end,1], 'label':labels[start:end]})
        sns.scatterplot(data = df, x = 'x', y = 'y', color = color)
        plt.legend()
        plt.xlabel('MDS dimension 1')
        plt.ylabel('MDS dimension 2')
        plt.xlim([-350,400])
        plt.ylim([-400,850])
        # plt.title('MDS projection of distance matrix')
        plt.axis('off')
        plt.tight_layout()
        plt.savefig(save_file + '.svg', dpi = 300)
        plt.show()
    zips = [
        (0, -1, 'all', '#B6B3D6'),
        (0,0 + 1,'native', '#D5D3DE'),
        (1,10 + 1, 'round1', '#B6B3D6'),
        (11, 11 + 14,'round2', '#F6DFD6'),
        (25, 25 + 14, 'round3', '#E9687A')]
    for a_zip in zips:plot_part(a_zip[0],a_zip[1],a_zip[2], a_zip[3])
    
    
    return coords

def reduce_esm2_embedding(X):
    pca = PCA(n_components=2)
    X_pca = pca.fit_transform(X)   # 形状 (29, 2)
    
    # 绘制二维散点图
    plt.figure(figsize=(8, 6))
    plt.scatter(X_pca[:, 0], X_pca[:, 1], c=np.arange(29), cmap='viridis', s=80, edgecolors='k')
    plt.xlabel('PC1 (%.2f%%)' % (pca.explained_variance_ratio_[0]*100))
    plt.ylabel('PC2 (%.2f%%)' % (pca.explained_variance_ratio_[1]*100))
    plt.title('PCA projection of 29 samples')
    plt.colorbar(label='Sample index')
    plt.grid(True, linestyle='--', alpha=0.6)
    plt.tight_layout()
    plt.show()
    return None 

def plot_entropy_scatter():

    return None

def reduce_matrix(matrix, seed = 42, labels = None, perplexity = None):
    import numpy as np
    import matplotlib.pyplot as plt
    from sklearn.manifold import MDS
    from sklearn.manifold import TSNE
    
    # MDS 降维到 2 维
    mds = MDS(n_components=2, random_state=seed, dissimilarity='precomputed')
    tsne = TSNE(n_components=2, random_state=seed, perplexity=perplexity)
    
    coords = mds.fit_transform(matrix)   # 形状 (29, 2)
    # coords = tsne.fit_transform(matrix)
    
    theta = np.radians(135)
    rot_matrix = np.array([[np.cos(theta), -np.sin(theta)],
                           [np.sin(theta),  np.cos(theta)]])
    
    new_coords =  coords @ rot_matrix
    print(new_coords.shape)
    return coords
    
def preprocess(csv_file, perplexity = 30):
    source_df = pd.read_csv(csv_file)
    
    with35 = True
    if not with35:df = source_df[~source_df['Name'].str.contains('PDN35', na=False)]
    else:df = source_df
        
    labels = df["Name"]
    fasta_list = list(df['Sequence'])
   
    matrix = 1 - get_identity_matrix(fasta_list)
    # coords = reduce_symmetric_matrix(matrix, seed = 742, labels = labels,perplexity = 10)
    coords = reduce_matrix(matrix, seed = 742, labels = labels,perplexity = perplexity)
    return coords, source_df