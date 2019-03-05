import numpy as np
import matplotlib.pyplot as plt 

x1 = np.array([1, 2, 3, 4, 5, 6, 7, 8, 9, 10])
y1 = np.array([44.82, 22.729, 15.836, 12.434, 19.706, 16.107, 11.154, 9.940, 16.768, 12.009])
y2 = np.array([5.226, 2.79, 2.03, 1.595, 1.65, 1.55, 1.4, 1,36, 1.65, 1.592])

y_omp1 = np.array([17.423, 9.162, 6.331, 5.295, 5.456, 4.595, 4.167, 4.146, 4.944, 4.813])
y_omp2 = np.array([106.610, 81.042, 58.761, 45.089, 48.680, 42.267, 36.625, 36.090, 39.051, 36.653])
y_omp3 = np.array([43, 22.48, 15.817, 8.68, 9.51, 8.28, 7.59, 7.74, 7.9, 7.81])
y_omp4 = np.array([5.21, 2.77, 1.94, 1.45, 1.62, 1.4, 1.3, 1.25, 1.565, 1.50])

plt.plot(x1, y1, label="pthreads")
plt.plot(x1, y_omp3, label="OpenMP")
plt.legend()
plt.xlabel("Number of threads")
plt.ylabel("Time [s]")
plt.title("Timer versus number of threads")
plt.show()