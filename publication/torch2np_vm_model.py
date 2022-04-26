# %% [markdown] tags=[]
# ### Conversion of Retro* value model to Numpy
#
# This notebook shows how the value model trained for Retro* is converted to a format supported by AiZynthFinder
#
# To use this notebook run the following commands first and execute this notebook in the `retro_star_repo` folder
#
# ```
# git clone https://github.com/binghong-ml/retro_star.git retro_star_repo
# cd retro_star_repo
# wget https://www.dropbox.com/s/ar9cupb18hv96gj/retro_data.zip?dl=0 -O retro_data.zip
# unzip retro_data.zip
# python -m pip install graphviz
# ```
#
# we assume that you have a Python environment installed wiht AiZynthFinder. The extraction of weights is also possible to do in an environment 
# setup for Retro*, and then you don't need to install `graphviz`

# %%
import sys
sys.path.append("retro_star/packages/mlp_retrosyn/")

import torch
import pickle
import numpy as np
from retro_star.model import ValueMLP

# %% [markdown]
# ## Extraction of weights from PyTorch model

# %%
device = "cpu"
fp_dim = 2048
n_layers = 1
latent_dim = 128
dropout_rate = 0.0
chkpt = torch.load("saved_models/best_epoch_final_4.pt")
model = ValueMLP(
        n_layers=n_layers,
        fp_dim=fp_dim,
        latent_dim=latent_dim,
        dropout_rate=dropout_rate,
        device=device,
).to(device)
model.load_state_dict(chkpt)

# %%
weights0 = chkpt['layers.0.weight'].numpy()
weights1 = chkpt['layers.3.weight'].numpy()
biases0 = chkpt['layers.0.bias'].numpy()
biases1 = chkpt['layers.3.bias'].numpy()

out = (
    [weights0.T.tolist(), weights1.T.tolist()], 
    [biases0.tolist(), biases1.tolist()]
)
with open("retrostar_value_model.pickle", "wb") as fileobj:
    pickle.dump(out, fileobj)


# %% [markdown]
# ## Validation

# %%
# This is taken from retro_star.common.smiles_to_fp

import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem

def smiles_to_fp(s, fp_dim=2048, pack=False):
    mol = Chem.MolFromSmiles(s)
    fp = AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits=fp_dim)
    onbits = list(fp.GetOnBits())
    arr = np.zeros(fp.GetNumBits(), dtype=np.bool)
    arr[onbits] = 1

    if pack:
        arr = np.packbits(arr)

    return arr


# %%
smi = "Cc1ccc2nc3ccccc3c(Nc3ccc(NC(=S)Nc4ccccc4)cc3)c2c1"
fp = smiles_to_fp(smi, fp_dim=fp_dim).reshape(1, -1)
torch_fp = torch.FloatTensor(fp).to(device)

# %%
# Result with original PyTorch model
model(torch_fp)

# %%
# This is the implementation in the aizynthfinder.search.retrostar.cost.RetroStarCost class, excluding the drop-out
vec = np.matmul(fp, weights0.T) +  biases0
vec = vec * (vec > 0)
vec = np.matmul(vec, weights1.T) + biases1
np.log(1 + np.exp(vec))

# %%
# Only possible to execute in an environment with AiZynthFinder installed
from aizynthfinder.search.retrostar.cost import RetroStarCost
from aizynthfinder.chem import Molecule
RetroStarCost("retrostar_value_model.pickle", dropout_rate=0)(Molecule(smiles=smi))
