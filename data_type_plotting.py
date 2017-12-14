import matplotlib.pyplot as plt
import numpy as np

input1= 'arange1_relative_errors.txt'
input2='arange2_relative_errors.txt'


all1= np.genfromtxt(input1, names= True)
all2= np.genfromtxt(input2, names= True)

plt.scatter(all1['time'], all1['poisson'], label= "Ely expected poisson noise", color = 'green', alpha = 0.5)
plt.scatter(all2['time'], all2['poisson'], label= "Ben expected poisson noise", color = 'cyan', alpha = 0.5)
plt.scatter(all1['time'], all1['std_dev'], label = "Ely's rel_std_dev", color= 'red')
plt.scatter(all2['time'], all2['std_dev'], label = "Ben's  rel_std_dev", color = 'blue')
plt.legend(numpoints=1, fontsize=14, loc='best' )
plt.show()

plt.scatter(all1['time'], all1['std_dev']/all2['std_dev']-1)
plt.xlabel('time bin (s)')
plt.axhline(y=0 ,linestyle = '-', color = 'r' , xmin = 0, xmax = 100000, linewidth = 1)
plt.title("Residuals of Ely's standard deviation divided by Ben's")
plt.show()
