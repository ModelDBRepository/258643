"""Differential Evolution Variator"""

import numpy
import random
from deap import tools

def varDE(population, toolbox, cxpb=1.0, mutpb=1.0, jitter=0, low=None, up=None):
    """Part of an evolutionary algorithm applying the variation part
    of the differential evolution algorithm. The modified individuals have their
    fitness invalidated. The individuals are cloned so returned population is
    independent of the input population.

    :param population: A list of individuals to vary.
    :param toolbox: A :class:`~deap.base.Toolbox` that contains the evolution
                    operators.
    :param cxpb: The crossover probability CR: probability of a 'base' parameter
                 value being replaced with a differentially evolved parameter
                 value
    :param mutpb: The mutation rate F: a scaling of the differential vector that
                  gets added to the base parameter
    :param low: The low bound of each parameter
    :param up: The upper bound of each parameter
    :returns: A list of varied individuals that are independent of their
              parents.

    The variation goes as follows. For each member of the population, we will 
    generate a new offspring. For each offspring, we select 3 random members of 
    the parent population: a base, and 2 parents for the differential. Then we 
    select a random parameter 'index' that we will always change. Then, for each 
    parameter, if a random number [0:1] is < Cx, or if the parameter is 'index', 
    we change the parameter from base by adding the differential of that 
    parameter in each of the second two parents (multiplied by rate factor F),
    (jittered by +/- jitter/2).
    That's it! 
    """

    offspring = [toolbox.clone(ind) for ind in population]

    # Differentially evolve each base offspring:
    for individual in offspring:
        # Keep mutating until all parameters are within all boundaries:
        inside = False
        while inside == False:
            base,a,b = tools.selRandom(population, 3)
            while (base == a or base == b or a == b):
                base,a,b = tools.selRandom(population, 3)
            index = random.randrange(len(individual))
            for j, parameter in enumerate(individual):
                if j == index or random.random() < cxpb:
                    diff = (a[j]-b[j])
                    individual[j] = base[j] + (mutpb*diff) + (((random.random()*jitter)-(jitter/2))*diff)
                    # perform 'bounce' away from boundaries to maintain parameter diversity
                    if individual[j] < low[j]:
                        individual[j] = low[j]+(low[j]-individual[j])
                    elif individual[j] > up[j]:
                        individual[j] = up[j]-(individual[j]-up[j])
            # Check boundary
            inside = True
            for j, parameter in enumerate(individual):
                if individual[j] < low[j]:
                    inside = False
                elif individual[j] > up[j]:
                    inside = False
        del individual.fitness.values

    return offspring

__all__ = ['varDE']
