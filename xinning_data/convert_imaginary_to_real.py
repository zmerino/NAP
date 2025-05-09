import numpy as np
import os

# Load the data
# load_filename = 'f_sa0.05_sb0.005.txt'
# save_filename = 'f_sa0.05_sb0.005_real.txt'
load_filename = 'f_sa0.05_sb0.002.txt'
save_filename = 'f_sa0.05_sb0.002_real.txt'
base_dir = os.path.join(os.getcwd(), 'xinning_data')
load_path = os.path.join(base_dir, load_filename)
save_path = os.path.join(base_dir, save_filename)
print(load_path)
data = np.loadtxt(load_path, dtype=complex)

# Extract the real part
real_data = data.real

# Save to a new file
np.savetxt(save_path, real_data)
print("Converted file saved as real_values.txt")
