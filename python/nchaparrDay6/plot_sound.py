from netCDF3 import Dataset
from scipy.interpolate import UnivariateSpline
import matplotlib.pyplot as plt
import numpy as np

filename='soundings.nc';
nc_file=Dataset(filename)
var_names=nc_file.variables.keys()
print "variable names: ",var_names
print "global attributes: ",nc_file.ncattrs()
#print "col_names: ",nc_file.col_names  This line gives me an error

fig=plt.figure(1)
fig.clf()
ax1=fig.add_subplot(111)
z_interp=np.arange(2000.,25000.,100.)
Temp_array = np.zeros(len(z_interp)) #initial array to start stacking Temperatures
for var_name,the_var in nc_file.variables.items():
     print var_name, the_var
     interp_temp=UnivariateSpline(the_var[:,1],the_var[:,2]) #note that Feb-28-2
                                                             # 0012-00z only goe
                                                             #s up to 11820
     Temp_array = np.vstack((Temp_array, interp_temp(z_interp))) #stacking Temperatures
     ax1.plot(interp_temp(z_interp),z_interp)
     fig.canvas.draw()
     plt.title('Interpolated Temperatures vs Height')
     plt.xlabel('Temperature(C)')
     plt.ylabel('Height(m)')
plt.show()

Temp_array = np.delete(Temp_array, 0, 0) #deleting initial row of zeros

#get averages accross soundings
Mean_Temp = np.mean(Temp_array, 0)
Std_Dev_Temp = np.std(Temp_array, 0)
Meanplus = np.add(Mean_Temp, Std_Dev_Temp) #mean plus one stdev
Meanminus = np.add(Mean_Temp, -Std_Dev_Temp)#mean minus one stdev

fig = plt.figure(2)
fig.clf()
ax2 = fig.add_subplot(111)
ax2.plot(Mean_Temp, z_interp, 'g', Meanplus, z_interp, 'r', Meanminus, z_interp, 'r')
plt.title('Mean Temperature plus and  minus a Standard Deviation')
plt.xlabel('Temperature(C)')
plt.ylabel('Height(m)')
fig.canvas.draw()
plt.show()
