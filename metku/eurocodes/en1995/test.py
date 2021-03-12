from matplotlib import cbook
from matplotlib import cm
from matplotlib.colors import LightSource
import matplotlib.pyplot as plt
import numpy as np
import math

# # Load and format data
# dem = cbook.get_sample_data('jacksboro_fault_dem.npz', np_load=True)
# z = dem['elevation']
# nrows, ncols = z.shape
# x = np.linspace(dem['xmin'], dem['xmax'], ncols)
# y = np.linspace(dem['ymin'], dem['ymax'], nrows)
# x, y = np.meshgrid(x, y)
#
# region = np.s_[5:50, 5:50]
# x, y, z = x[region], y[region], z[region]
#
# # Set up plot
# fig, ax = plt.subplots(subplot_kw=dict(projection='3d'))
#
# ls = LightSource(270, 45)
# # To use a custom hillshading mode, override the built-in shading and pass
# # in the rgb colors of the shaded surface calculated from "shade".
# rgb = ls.shade(z, cmap=cm.gist_earth, vert_exag=0.1, blend_mode='soft')
# surf = ax.plot_surface(x, y, z, rstride=1, cstride=1, facecolors=rgb,
#                        linewidth=0, antialiased=False, shade=False)

#plt.show()



# plt.plot([1, 1.5, 2, 2.5, 3, 4, 5, 6, 8], [0.141, 0.196, 0.229, 0.249, 0.263, 0.281, 0.299, 0.307, 0.313], 'go-', label='kerroin alpha', linewidth=2)
# #plt.show()
#
#
# print(405 * 90 ** 3 * (0.08 * math.log(405 / 90, 2.3) + 0.14))
# print(495 * 105 ** 3 * (0.08 * math.log(495 / 105, 2.3) + 0.14))
# print(405 * 140 ** 3 * (0.08 * math.log(405 / 140, 2.1) + 0.14))
# print(450 * 165 ** 3 * (0.08 * math.log(450 / 165, 2) + 0.14))
#
# print(360 * 42 ** 3 * (0.08 * math.log(360 / 42, 2.73) + 0.14))
# print(405 * 66 ** 3 * (0.08 * math.log(405 / 66, 2.6) + 0.14))
# print(40 * 10 ** 3 * (0.08 * math.log(40 / 10, 2.5) + 0.14))



a = np.arange(15).reshape(3, 5)
print(a)

ixgrid = np.ix_([0, 1], [2, 3, 2, 0])


print(a[ixgrid])
