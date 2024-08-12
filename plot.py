import matplotlib.pyplot as plt

# Example data (replace these with your actual data)
tolerance_limits = [1e-2, 1e-3, 1e-4, 1e-5, 1e-6, 1e-7, 1e-8, 1e-9]


## torus
error_norms= [2.066,2.0726, 2.1388, 2.4592, 2.1829 ,2.3259, 2.3259, 2.3142]
time_taken=[57.7082, 55.7021, 55.501, 54.3027, 54.9259 , 54.2747, 54.2535, 54.3963]

# Plotting norm of the error
plt.figure(figsize=(10, 6))
plt.plot(tolerance_limits, error_norms, marker='o', linestyle='-', color='blue')
plt.xlabel('Tolerance Limit')
plt.ylabel('Norm of the Error')
plt.title('Norm of the Error as a Function of Tolerance Limit')
plt.grid(True)
plt.show()

# Plotting time taken
plt.figure(figsize=(10, 6))
plt.plot(tolerance_limits, time_taken, marker='o', linestyle='-', color='green')
plt.xlabel('Tolerance Limit')
plt.ylabel('Time Taken (seconds)')
plt.title('Time Taken as a Function of Tolerance Limit')
plt.grid(True)
plt.show()