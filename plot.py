import numpy as np
import matplotlib.pyplot as plt 

x1 = np.array([1, 2, 3, 4, 5, 6, 7, 8, 9, 10])
y1 = np.array([44.82, 22.729, 15.836, 12.434, 19.706, 16.107, 11.154, 9.940, 16.768, 12.009])
y11 = np.array([20.257, 11.719, 9.369, 8.140, 9.245, 8.398, 7.203, 6.704, 7.756, 7.76])
y2 = np.array([2.436, 1.436, 1.167, 1.045, 1.421, 1.082, 0.998, 1.001, 1.094, 1.1 ])

y_omp1 = np.array([17.423, 9.162, 6.331, 5.295, 5.456, 4.595, 4.167, 4.146, 4.944, 4.813])
y_omp2 = np.array([19.951, 11.687, 9.374, 8.329, 8.754, 7.667, 7.000, 6.492, 7.635, 7.63])
y_omp3 = np.array([43, 22.48, 15.817, 8.68, 9.51, 8.28, 7.59, 7.74, 7.9, 7.81])
y_omp4 = np.array([2.386, 1.422, 1.148, 1.007, 1.065, 0.943, 0.858, 0.823, 1.082, 1.076])

plt.plot(x1, y11, label="pthreads")
plt.plot(x1, y_omp2, label="OpenMP")
plt.legend()
plt.xlabel("Number of threads")
plt.ylabel("Time [s]")
plt.title("Timer versus number of threads")
plt.show()