import pickle
import matplotlib.pyplot as plt
import numpy as np
import scipy.io as sio
from mpl_toolkits.mplot3d import Axes3D
results="/home/jindong/machineLearning/Human-Pose-Estimation/logs/eval_ibhgc_vol_softmax_VolumetricTriangulationNet@19.04.2021-17:33:12/checkpoints/0000/results.pkl"
with open(results, 'rb') as f:
    data = pickle.load(f)
#data keypoints_3d [frames,joints,positionxyz]
keypoints_position=data['keypoints_3d']
