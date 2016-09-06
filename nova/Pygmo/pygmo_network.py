from PyGMO import *
import pylab as pl
from time import time
from itertools import cycle

import seaborn as sns
rc = {'figure.figsize':[4,4],'savefig.dpi':250, # 
      'savefig.jpeg_quality':100,'savefig.pad_inches':0.1,
      'lines.linewidth':2}
color = cycle(sns.color_palette('Set2',5))

'''
prob = problem.zdt(2)
alg = algorithm.nsga_II(gen=100)

archi = archipelago(alg,prob,4,40,topology=topology.ring())
 
pos = archi.draw(scale_by_degree=True,n_size=4,e_alpha=0.03)
  
tic = time()
for i in range(100):
    archi.evolve(1)
archi.join();
print(time()-tic)
pos = archi.draw(layout=pos,scale_by_degree=True,n_size=4,e_alpha=0.03)

    
pop = population(prob)
for isl in archi:
    for ind in isl.population:
        pop.push_back(ind.cur_x)

pop.plot_pareto_fronts()

'''
#pop.plot_pareto_fronts(fronts=[-1],rgb=[0,1,0])

#pop.plot_pareto_fronts(rgb=[0,0,1])

'''
tic = time()
prob = problem.zdt(4)
alg = algorithm.nsga_II(gen=1)
pop = population(prob,200)
for i in range(500):
    pop = alg.evolve(pop)
print(time()-tic)
   
pop.plot_pareto_fronts(rgb=[1,0,0])
'''

prob = problem.schwefel(dim = 50)
algo = []
for i in range(1,9):
    algo.append(algorithm.de(gen=500,variant=i))
archi = archipelago(topology=topology.ring())
for i in range(0,8):
    archi.push_back(island(algo[i],prob,20))
print(min([isl.population.champion.f for isl in archi]))
archi.evolve(20)
print(min([isl.population.champion.f for isl in archi]))