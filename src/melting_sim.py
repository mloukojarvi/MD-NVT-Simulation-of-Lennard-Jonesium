from matplotlib import pyplot as plt
from matplotlib.colors import ListedColormap, Normalize
import numpy as np
import os

## Set parameters for the simulations ##
N = 500
total_steps = 20000
equilib_steps = 15000
dt = 0.005
gamma_ = 1.0
		
min_dens = 0.86
max_dens = 1.05
dens_step = 0.01

min_temp = 0.5
max_temp = 2.5
temp_step = 0.02

def run_melting_simulations():
	# Runs the MD simulations in a region close to the melting line given by the literature
	
	# Simulated densities
	densities = list(np.arange(min_dens, max_dens+dens_step/2, dens_step))
	
	def temp_limit(density, limit):
		# Finds the temperature limits for a given density
		if limit.lower() == 'lower':
			m, c = 6.6, -5.5
		elif limit.lower() == 'upper':
			m, c = 6.6, -4.4 
		return m*density+c
	
	# Compile the simulation code
	os.system("cd ../src/ && make")
	os.system("cd ../run/")
	
	for dens in densities:
		min_temp_ = max(temp_limit(dens, 'lower'), min_temp)
		max_temp_ = min(temp_limit(dens, 'upper'), max_temp)
		
		# Simulated temperatures for the given density
		temps = list(np.arange(min_temp_, max_temp_, temp_step))
		
		for temp in temps:
			# Run the simulation for the given density and temperature
			os.system(f"../src/main {N} {temp} {dens} {total_steps} {dt} {gamma_}")

def sort_by_dens():
	# Sorts the data into a dictionary where the keys are the simulated densities and the items are
	# dictionaries whose keys are the simulated temperatures for the given density and items
	# the average potential energy, radial distribution function histogram and the corresponding
	# distances for the histogram's x-axis
	
	data = {}
	for filename in os.listdir('./'):
		data_type = filename[:3]
		if data_type in ['mea', 'rdf']:
			with open(filename, 'r') as f:
				for idx, line in enumerate(f):
					line = list(map(float, line.split()))
					if idx == 0:
						N, temp, dens, tot_steps, dt, spd_bins, rdf_bins, sc, rc = line
						dens = float(f'{dens:4.3f}')
						temp = float(f'{temp:4.3f}')
						avg_pot = 0
						rdf = np.array([0]*int(rdf_bins))
						dr = (rc-0.8)/rdf_bins
						r = [0.8+i*dr for i in range(int(rdf_bins))]
						if dens not in data:
							data[dens] = {}
						if temp not in data[dens]:
							data[dens][temp] = {'avg_pot': 0, 'rdf': [], 'r':r}
					elif idx > equilib_steps:
						if data_type == 'mea':
							time, inst_temp, PE, KE, msd = line
							avg_pot = (PE+avg_pot)/2
						elif data_type == 'rdf':
							line = np.array(line)
							rdf = rdf + line
				if data_type == 'mea':
					data[dens][temp]['avg_pot'] = avg_pot
				elif data_type == 'rdf':
					data[dens][temp]['rdf'] = rdf/(total_steps-equilib_steps)
	return data
	
def differentiate(x, y):
	# Calculates the numerical derivative of y with respect to x
	
	dy = np.diff(y)/np.diff(x)
	dx = []
	for i in range(len(dy)):
		temp_x = (x[i+1]+x[i])/2
		dx = np.append(dx, temp_x)
	return dx, dy
	
def find_melting_temps(data):
	# Finds the melting temperature for each density
	
	melting_temps = {}
	for dens in data:
		temps = []
		avg_pots = []
		for temp in data[dens]:
			temps.append(temp)
			avg_pots.append(data[dens][temp]['avg_pot'])
		
		
		# The melting temperature is given by the temperature corresponding to the maximum of 
		# the derivative of the average potential energy with respect to the temperature
		dx, dy = differentiate(temps, avg_pots)
		idx = np.argmax(dy)
		melting_temp = dx[idx]
		melting_temps[dens] = melting_temp
	return melting_temps
	
	
def get_color_gradient(c1, c2, n):
	# Creates a color gradient between the colors
	
    mixer = [x/(n-1) for x in range(n)]
    colors = [((1-mix)*np.array(c1) + (mix*np.array(c2))) for mix in mixer]
    return colors
			
def plot_pot(data):
	# Plots the potential energy-temperature curves
	
	plt.rcParams.update({'font.size': 15})
	fig = plt.figure(figsize=(12,6))
	
	colors = get_color_gradient([0,0,1], [1,0,0], len(data))
	norm = Normalize(vmin=min_dens, vmax=max_dens)
	cmap = ListedColormap(colors)
	scalar_map = plt.cm.ScalarMappable(norm=norm, cmap=cmap)
	colorbar = plt.colorbar(scalar_map)
	colorbar.set_label(r'$\rho^*$',  rotation=90)
	
	min_temp_ = float('inf')
	max_temp_ = -float('inf')
	min_pot = float('inf')
	max_pot = -float('inf')
	
	for dens in data:
		data_points = []
		for temp in data[dens]:
			data_points.append((temp, data[dens][temp]['avg_pot']))
		data_points.sort(key=lambda x: x[0])
		x, y = np.array(data_points).T
		min_temp_ = min(min_temp_, min(x))
		max_temp_ = max(max_temp_, max(x))
		min_pot = min(min_pot, min(y))
		max_pot = max(max_pot, max(y))
		idx = int((dens-min_dens)/dens_step)
		c = colors[idx]
		plt.plot(x, y, '-', linewidth=1.5, color=c)
		
	plt.title(f'Average potential energy as a function of temperature')
	plt.xlabel(r'$T^*$')
	plt.ylabel(r"$\langle U^*\rangle$")
	plt.xlim(min_temp_, max_temp_)
	plt.ylim(min_pot, max_pot)
	plt.grid()
	plt.show()
	
def plot_melting_line(melting_temps):
	# Plots the melting line
	
	x = []
	y = []
	for dens in melting_temps:
		x.append(dens)
		y.append(melting_temps[dens])
	
	fig = plt.figure(figsize=(12,15))
	
	## Plot the phase diagram from the literature ##
	tx1 = [0.871, 0.890, 0.924, 0.971, 1.000, 1.050, 1.062]
	ty1 = [0.681, 0.778, 0.997, 1.352, 1.590, 2.060, 2.197]
	fit1 = np.polyfit(tx1, ty1, 4)
	p1 = np.poly1d(fit1)
	fit_x1 = np.arange(0.865, 1.07, 0.001)
	fit_y1 = [p1(x) for x in fit_x1]
	
	tx2 = [0.728, 0.751, 0.778, 0.802, 0.833, 0.864]
	ty2 = [0.955, 0.910, 0.854, 0.801, 0.730, 0.653]
	fit2 = np.polyfit(tx2, ty2, 4)
	p2 = np.poly1d(fit2)
	fit_x2 = np.arange(0.7, 0.865, 0.001)
	fit_y2 = [p2(x) for x in fit_x2]
	
	tx3 = [0.951, 0.967, 0.982, 0.998, 1.012, 1.034, 1.058, 1.079, 1.099]
	ty3 = [0.667, 0.715, 0.789, 0.908, 1.030, 1.232, 1.451, 1.667, 1.860]
	fit3 = np.polyfit(tx3, ty3, 4)
	p3 = np.poly1d(fit3)
	fit_x3 = np.arange(0.951, 1.15, 0.001)
	fit_y3 = [p3(x) for x in fit_x3]
	
	c = '#acacac'
	plt.plot(fit_x1, fit_y1, linewidth=2, color='#404040', linestyle = '-', zorder=1, label='Literature')
	plt.plot(fit_x2, fit_y2, linewidth=2, color=c, linestyle = '-', zorder=2)
	plt.plot(fit_x3, fit_y3, linewidth=2, color=c, linestyle = '-', zorder=3)
	plt.axhline(y = 0.663, color = c, linestyle = '--', linewidth=2, zorder=4)
	plt.scatter(x=[0.863], y=[0.663], s=300, color=c, zorder=5)
	
	plt.text(0.91, 1.84, 'SC', fontdict=None, size='xx-large')
	plt.text(0.85, 0.87, 'L', fontdict=None, size='xx-large')
	plt.text(1.05, 0.875, 'S', fontdict=None, size='xx-large')
	plt.text(0.99, 1.3, 'S+L', fontdict=None, size='xx-large')
	
	## Plot the simulated melting line ##
	ce = 'red'
	plt.errorbar(x, y, yerr=0.040, marker='o', ls='none', mfc='blue', \
				 mec=ce, ms=12, mew=2, ecolor=ce, elinewidth=2.5, capsize=5, zorder=6, label='Simulation')
	
	## Name the axes and set plot parameters ##
	plt.title(r'Melting curve in the $(\rho^*, T^*)$-plane')
	plt.xlabel(r'$\rho^*$')
	plt.ylabel(r"$T^*$")
	plt.legend()
	plt.grid()
	plt.xlim(0.8, 1.1)
	plt.ylim(0.6, 2.2)
	
	plt.show()
	
def plot_rdf(data, melting_temps, dens):
	# Plots the radial distribution function 
	
	min_temp_ = max(melting_temps[dens] - 10*temp_step, min(data[dens].keys()))
	max_temp_ = min(melting_temps[dens] + 11*temp_step, max(data[dens].keys()))
	nc = int((max_temp_-min_temp_)/temp_step)+1
	
	fig = plt.figure(figsize=(11,6))
	colors = get_color_gradient([0,0,1], [1,0,0], nc)
	norm = Normalize(vmin=min_temp, vmax=max_temp)
	cmap = ListedColormap(colors)
	scalar_map = plt.cm.ScalarMappable(norm=norm, cmap=cmap)
	colorbar = plt.colorbar(scalar_map)
	colorbar.set_label(r'$T^*$',  rotation=90)
	
	for temp in data[dens]:
		if min_temp_ <= temp <= max_temp_:
			rdf = data[dens][temp]['rdf']
			r = data[dens][temp]['r']
			idx = int((temp-min_temp_)/temp_step)
			c = colors[idx]
			plt.plot(r, rdf, '-', linewidth=1, color=c)
			
	plt.plot([min(r), max(r)], [1, 1], color='k', linestyle=':', linewidth=2, label=r'$g^*=1$')
	plt.title(r'Radial distribution function, $\rho^*=$' + f'{dens:4.3f}')
	plt.xlabel(r"$r^*$")
	plt.ylabel(r"$g^*$")
	plt.xlim(0.8, 3.0)
	plt.grid()
	plt.legend()
	plt.show()

def melting_simulation():
	## Run the simulations and sort the data ##
	run = input('Run the simulations? (Y/N): ')
	if run.lower() == 'y':
		run_melting_simulations()
	data = sort_by_dens()
	melting_temps = find_melting_temps(data)
	
	## Plot the simulation data ##
	rdf_dens = [0.86, 0.95]
	for dens in rdf_dens:
		plot_rdf(data, melting_temps, dens)
	plot_pot(data)
	plot_melting_line(melting_temps)

if __name__ == '__main__':
	melting_simulation()	    
