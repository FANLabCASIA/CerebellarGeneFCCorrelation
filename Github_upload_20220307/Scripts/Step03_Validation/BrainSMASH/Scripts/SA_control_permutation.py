import pandas as pd
import numpy as np
import os

sample_coordinate = pd.read_csv("your_sample_coordinate_filename.csv")
sample_brainmap = pd.read_csv("your_brainmap_filename.csv")

# calculate distance matrix
from brainsmash.workbench.geo import volume
coord_file = sample_coordinate.values
output_dir = "your_outputpath"
# if not os.path.exists(output_dir):
#     os.makedirs(output_dir)
filenames = volume(coord_file, output_dir)

# variogram fitting
import matplotlib as plt 
from brainsmash.mapgen.eval import sampled_fit
%matplotlib inline
# These are three of the key parameters affecting the variogram fit, keep "knn" roughly consistent with the number of samples
kwargs = {'ns': 50,
        'knn': 100,¡¡
        'pv': 70}
for i in range(sample_brainmap.shape[1]):
    brain_map = sample_brainmap.iloc[:, i].values
    # Running this command will generate a matplotlib figure
    sampled_fit(brain_map, filenames['D'], filenames['index'], nsurr=10, **kwargs)
    # plt.savefig(output_dir + "/SAfitting_plot/" + str(i) + ".jpg")

# randomly generate surrogate brain maps with SA matched to target brain map
from brainsmash.mapgen.sampled import Sampled
import time
from tqdm import tqdm
start = time.time()
permutation_times = 10000
ans = np.zeros([permutation_times, sample_brainmap.shape[0], sample_brainmap.shape[1]])

for i in tqdm(range(sample_brainmap.shape[1])):
    brain_map = sample_brainmap.iloc[:, i].values
    gen = Sampled(x=brain_map, D=filenames['D'], index=filenames['index'], **kwargs)
    surrogate_maps = gen(n = permutation_times)
    ans[:, :, i] = surrogate_maps.reshape(permutation_times, sample_brainmap.shape[0])
np.save("your_output_file.csv", ans)
    
# assign sample_points to different ROI 
label_id = np.unique(sample_coordinate[["label"]].values)
network_index_1 = []
network_index_2 = []

for i in label_id:
    curr_index = np.argwhere((sample_coordinate[["label"]].values.squeeze() == i) & (sample_coordinate[["donor"]].values.squeeze() == 9861)).squeeze()
    network_index_1.append(curr_index)
    curr_index = np.argwhere((sample_coordinate[["label"]].values.squeeze() == i) & (sample_coordinate[["donor"]].values.squeeze() == 10021)).squeeze()
    network_index_2.append(curr_index)

# save permutation results
permutation_times = 10000
permutationmap_dir = output_dir + "Permutationmap_10000/"

for times in tqdm(range(permutation_times)):
    curr_brainmap = ans[times, :, :]
    mean_brainmap = np.zeros([2, len(network_index_1), curr_brainmap.shape[1]])
    for i in range(label_id.shape[0]):
        if network_index_1[i].shape == ():
            mean_brainmap[0, i, :] = np.mean(curr_brainmap[network_index_1[i], :][None, :], axis = 0)
        else:
            mean_brainmap[0, i, :] = np.mean(curr_brainmap[network_index_1[i], :], axis = 0)

        if network_index_2[i].shape == ():
            mean_brainmap[1, i, :] = np.mean(curr_brainmap[network_index_2[i], :][None, :], axis = 0)
        else:
            mean_brainmap[1, i, :] = np.mean(curr_brainmap[network_index_2[i], :], axis = 0)
           
    np.savetxt("your_output_file_1.csv", mean_brainmap[0], delimiter = ',')
    np.savetxt("your_output_file_2.csv", mean_brainmap[1], delimiter = ',')
