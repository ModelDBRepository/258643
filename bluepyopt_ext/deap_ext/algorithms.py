"""
Copyright (c) 2016, EPFL/Blue Brain Project

 This file is part of BluePyOpt <https://github.com/BlueBrain/BluePyOpt>

 This library is free software; you can redistribute it and/or modify it under
 the terms of the GNU Lesser General Public License version 3.0 as published
 by the Free Software Foundation.

 This library is distributed in the hope that it will be useful, but WITHOUT
 ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more
 details.

 You should have received a copy of the GNU Lesser General Public License
 along with this library; if not, write to the Free Software Foundation, Inc.,
 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.


NOTE: This is the one that does the loops across generations
"""

# pylint: disable=R0914, R0912


import random
import logging
import pickle

import deap.algorithms
import deap.tools
import os
import shutil

logger = logging.getLogger('__main__')


def _evaluate_invalid_fitness(toolbox, population, gen_index=None):
    '''Evaluate the individuals with an invalid fitness

    Returns the count of individuals with invalid fitness
    '''
    invalid_ind = [ind for ind in population if not ind.fitness.valid]
    logger.info(" -- total individuals that need to be assessed (for this generation) is %s (check if it's the same as the population's size)", len(invalid_ind))
    logger.info(" ---- this is the list of all parameters for the 1st individual (check if the number is the same as IND_SIZE): %s " % (invalid_ind[0]))
    from functools import partial
    execute = partial(toolbox.evaluate, gen_index=gen_index)
    fitnesses = toolbox.map(execute, invalid_ind)
    # fitnesses = toolbox.map(toolbox.evaluate, invalid_ind, gen_index)
    for ind, fit in zip(invalid_ind, fitnesses):
        ind.fitness.values = fit

    return len(invalid_ind)


def _evaluate_invalid_fitness_with_features(toolbox, population, gen_index=None):
    '''Evaluate the individuals with an invalid fitness

    Returns the count of individuals with invalid fitness

    Not only saves fitness values to ind.fitnes.values, also
    saves raw feature scores to ind.fitness.features
    '''
    invalid_ind = [ind for ind in population if not ind.fitness.valid]
    # fitnesses, feature_values = toolbox.map(toolbox.evaluate_with_features, invalid_ind)
    from functools import partial
    execute = partial(toolbox.evaluate_with_features, gen_index=gen_index)
    fitness_dicts = toolbox.map(execute, invalid_ind)
    # fitness_dicts = toolbox.map(toolbox.evaluate_with_features, invalid_ind)
    for ind, fit in zip(invalid_ind, fitness_dicts):
        ind.fitness.values = fit['scores'].values()
        ind.fitness.feature_values = fit['feature_values']

    return len(invalid_ind)


def _evaluate_invalid_fitness_minus_SDs_with_features(toolbox, population,
                                                      gen_index=None, numSDs=0,
                                                      savefile_base=None):
    '''Evaluate the individuals with an invalid fitness

    Returns the count of individuals with invalid fitness

    Before saving the fitness values, this method subtracts
    numSDs from each fitness value returned by efel, with a
    minimum fitness score of 0 - so it flattens the error
    function within an 'acceptable' number of SDs

    Not only saves fitness values to ind.fitnes.values, also
    saves raw feature scores to ind.fitness.features
    '''
    invalid_ind = [ind for ind in population if not ind.fitness.valid]
    from functools import partial
    execute = partial(toolbox.evaluate_with_features, gen_index=gen_index, savefile_base=savefile_base)
    # third argument (ind_index) is not accepted by pool.map
    fitness_dicts = toolbox.map(execute, invalid_ind)
    for ind, fit in zip(invalid_ind, fitness_dicts):
        errlst = list(fit['scores'].values())
        for i, item in enumerate(errlst):
            errlst[i] -= numSDs
            if errlst[i] < 0.0:
                errlst[i] = 0.0
        ind.fitness.values = tuple(errlst)
        ind.fitness.feature_values = fit['feature_values']

    return len(invalid_ind)


def _update_history_and_hof(halloffame, history, population):
    '''Update the hall of fame with the generated individuals

    Note: History and Hall-of-Fame behave like dictionaries
    '''
    if halloffame is not None:
        halloffame.update(population)

    history.update(population)


def _update_history(history, population):
    ''' Update history with generated individuals'''
    history.update(population)


def _record_stats(stats, logbook, gen, population, invalid_count):
    '''Update the statistics with the new population'''
    record = stats.compile(population) if stats is not None else {}
    logbook.record(gen=gen, nevals=invalid_count, **record)


def _get_offspring(parents, toolbox, cxpb, mutpb):
    '''return the offsprint, use toolbox.variate if possible'''
    if hasattr(toolbox, 'variate'):
        return toolbox.variate(parents, toolbox, cxpb, mutpb)
    return deap.algorithms.varAnd(parents, toolbox, cxpb, mutpb)


def _write_h5_file(population, filename, featsToUse):
    '''Write parameters, errors and features to a matrix in an h5 file'''
    import numpy
    parlist = []
    errlist = []
    fealist = []
    for ind in population:
        for param in ind:
            parlist.append(param)
        for value in ind.fitness.values:
            errlist.append(value)
        for featname in featsToUse:
            for feat in ind.fitness.feature_values.values():
                if featname in feat:
                    fealist.append(feat[featname])
    pararr = numpy.asarray(parlist)
    errarr = numpy.asarray(errlist)
    feaarr = numpy.asarray(fealist)
    for num, element in enumerate(pararr):
        if element == None:
            pararr[num] = numpy.float64(1e9)
    for num, element in enumerate(errarr):
        if element == None:
            errarr[num] = numpy.float64(1e9)
    for num, element in enumerate(feaarr):
        if element == None:
            feaarr[num] = numpy.float64(1e9)
    parmat = pararr.reshape(len(population),len(pararr)/len(population))
    errmat = errarr.reshape(len(population),len(errarr)/len(population))
    feamat = feaarr.reshape(len(population),len(feaarr)/len(population))
    feamat = feamat.astype(numpy.float64)
    import h5py
    with h5py.File(filename,'a') as f:
        if 'parameters' in f.keys():
            f['parameters'].resize((f['parameters'].shape[0] + parmat.shape[0]), axis = 0)        
            f['parameters'][-parmat.shape[0]:] = parmat
        else:
            f.create_dataset('parameters',data=parmat,maxshape=(None,None))
        if 'errors' in f.keys():
            f['errors'].resize((f['errors'].shape[0] + errmat.shape[0]), axis = 0)        
            f['errors'][-errmat.shape[0]:] = errmat
        else:
            f.create_dataset('errors',data=errmat,maxshape=(None,None))
        if 'features' in f.keys():
            f['features'].resize((f['features'].shape[0] + feamat.shape[0]), axis = 0)        
            f['features'][-feamat.shape[0]:] = feamat
        else:
            f.create_dataset('features',data=feamat,maxshape=(None,None))


def NonDominatedSortingDifferentialEvolutionFeatureCrowdingCheckpoint(
        population,
        toolbox,
        mu,
        cxpb,
        mutpb,
        ngen,
        stats=None,
        halloffame=None,
        cp_frequency=1,
        cp_filename=None,
        h5_filename=None,
        featsToUse=None,
        continue_cp=False,
        numSDs=0,
        savefile_base=None,
        global_info=None):
    """This is the Non-Dominated Sorting Differential Evolution (NSDE) evolutionary algorithm
        Using features for the crowding distance, and subtracting SDs from the error values,
        so all error values are 0 with a number of SDs of the mean

    Args:
        population(list of deap Individuals)
        toolbox(deap Toolbox)
        mu(int): Total parent population size of EA
        cxpb(float): Crossover probability
        mutpb(float): Mutation probability
        ngen(int): Total number of generation to run
        stats(deap.tools.Statistics): generation of statistics
        halloffame(deap.tools.HallOfFame): hall of fame
        cp_frequency(int): generations between checkpoints
        cp_filename(string): path to checkpoint filename
        continue_cp(bool): whether to continue
    """

    if continue_cp:
        # A file name has been given, then load the data from the file
        cp = pickle.load(open(cp_filename, "r"))
        population = cp["population"]
        parents = cp["parents"]
        start_gen = cp["generation"]
        logbook = cp["logbook"]
        random.setstate(cp["rndstate"])
        halloffame = cp["halloffame"]
        history = cp["history"]
    else:
        # Start a new evolution
        start_gen = 1
        parents = population[:]
        logbook = deap.tools.Logbook()
        logbook.header = ['gen', 'nevals'] + (stats.fields if stats else [])
        history = deap.tools.History()

        invalid_count = _evaluate_invalid_fitness_minus_SDs_with_features(
            toolbox, population, start_gen, numSDs, savefile_base)
        _update_history_and_hof(halloffame, history, population)
        _record_stats(stats, logbook, start_gen, population, invalid_count)
        if h5_filename is not None:
            _write_h5_file(population, h5_filename, featsToUse)

    # Begin the generational process
    for gen in range(start_gen + 1, ngen + 1):
        offspring = _get_offspring(parents, toolbox, cxpb, mutpb)
        population = parents + offspring

        invalid_count = _evaluate_invalid_fitness_minus_SDs_with_features(toolbox, offspring, gen, numSDs, savefile_base)
        _update_history_and_hof(halloffame, history, population)
        _record_stats(stats, logbook, gen, population, invalid_count)

        if h5_filename is not None:
            _write_h5_file(offspring, h5_filename, featsToUse)
            shutil.copy(h5_filename, h5_filename+'.bak')
            _write_h5_file(parents, h5_filename+'.parents', featsToUse)
            shutil.copy(h5_filename, h5_filename+'.parents.bak')

        # Select the next generation parents
        parents = toolbox.select(population, mu)

        logger.info(logbook.stream)

        _record_stats(stats, logbook, gen, offspring, invalid_count)
        logger.info(logbook.stream)
        _record_stats(stats, logbook, gen, parents, invalid_count)
        logger.info(logbook.stream)

        # Output some feature information for the 'best' individuals:
        for ind in parents[0:4]:
            outlist = []
            for featname in toolbox.select.keywords['featsToUse']:
                for feat in ind.fitness.feature_values.values():
                    if featname in feat:
                        outlist.append(feat[featname])
            logger.info(outlist)

        if(cp_filename and cp_frequency and
           gen % cp_frequency == 0):
            cp = dict(population=population,
                      generation=gen,
                      parents=parents,
                      halloffame=halloffame,
                      history=history,
                      logbook=logbook,
                      rndstate=random.getstate())
            pickle.dump(cp, open(cp_filename, "wb"))
            logger.debug('Wrote checkpoint to %s', cp_filename)
            shutil.copy(cp_filename, cp_filename+'.bak')

    return population, logbook, history
