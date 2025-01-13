import numpy as np

#Loading npz files
data = np.load('scores.model_idx_0.npz')

#List all components in the file
print(data.files)

