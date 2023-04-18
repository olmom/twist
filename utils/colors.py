import numpy as np
from matplotlib.colors import ListedColormap
import matplotlib.cm as cm

# the colormap to create
BlueBlackRed = None

# create listedColormap
bottom = cm.get_cmap('Blues', 256)
top = cm.get_cmap('Reds_r', 256)
mycolormap = np.vstack((bottom(np.linspace(0.25, 1, 64)),
                        np.array([
                        [0.03137255, 0.08823529, 0.41960784, 1.],
                        [0.02137255, 0.04823529, 0.21960784, 1.],
                        [0.01137255, 0.02823529, 0.11960784, 1.],
                        [0.00037255, 0.00823529, 0.00960784, 1.],
                        ])
                       ))
mycolormap = np.vstack((mycolormap,
                        np.array([
                        [0.00960784, 0.00823529, 0.00037255, 1.],
                        [0.11960784, 0.02823529, 0.01137255, 1.],
                        [0.21960784, 0.04823529, 0.02137255, 1.],
                        [0.41960784, 0.08823529, 0.03137255, 1.],
                        ])
                       ))
mycolormap = np.vstack((mycolormap,
                        top(np.linspace(0, 0.75, 64)),
                       ))

BlueBlackRed = ListedColormap(mycolormap, name='BlueBlackRed')