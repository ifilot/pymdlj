#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import numpy as np
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt

from pymdlj import PyMDLJ

sim = PyMDLJ()
results = sim.simulate(os.path.join(os.path.dirname(__file__), '..', 'tests', 'settings', 'default.param.in'))
p0 = results['positions_initial']
v0 = results['velocities_initial']
p = results['positions']
v = results['velocities']
vm0 = [np.linalg.norm(vv) for vv in v0]
vm = [np.linalg.norm(vv) for vv in v]

fig = plt.figure(figsize=(12,4))
ax = fig.add_subplot(121, projection='3d')
ax.scatter(p0[:,0], p0[:,1], p0[:,2], c=vm0, marker='o')

ax = fig.add_subplot(122, projection='3d')
ax.scatter(p[:,0], p[:,1], p[:,2], c=vm, marker='o')
plt.show()