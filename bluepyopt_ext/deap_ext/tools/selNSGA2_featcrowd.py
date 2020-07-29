from __future__ import division
import bisect
import math
import random
import logging
import numpy

from scipy.spatial import distance
from itertools import chain
from operator import attrgetter, itemgetter
from collections import defaultdict

logger = logging.getLogger('__main__')

######################################
# Non-Dominated Sorting   (NSGA-II)  #
# With addition of using raw feature #
# values to determine crowdedness    #
######################################

def selNSGA2_featcrowd(individuals, k, nd='standard', numSDs=0, featsToUse=[]):
    """Apply NSGA-II selection operator on the *individuals*. Usually, the
    size of *individuals* will be larger than *k* because any individual
    present in *individuals* will appear in the returned list at most once.
    Having the size of *individuals* equals to *k* will have no effect other
    than sorting the population according to their front rank. The
    list returned contains references to the input *individuals*. For more
    details on the NSGA-II operator see [Deb2002]_.
    
    :param individuals: A list of individuals to select from.
    :param k: The number of individuals to select.
    :param nd: Specify the non-dominated algorithm to use: 'standard' or 'log'.
    :param numSDs: The number of SDs around the mean for which to zero the error
    :param featsToUse: a list containing indices for which features to use for
                       crowdedness values - an empty list means all features will
                       be used for crowdedness
    :returns: A list of selected individuals.

    .. [Deb2002] Deb, Pratab, Agarwal, and Meyarivan, "A fast elitist
       non-dominated sorting genetic algorithm for multi-objective
       optimization: NSGA-II", 2002.
    """

    """shuffle individuals first so that individuals with the same rank
       and crowdedness rating are randomly selected between, instead of
       being selected by age within the population i.e. are they at a
       lower index within the previous generation.
       This mostly helps with the situation where all fitness values
       come back as NAN for the initial generations, then the population
       gets shuffled every generation so more variety of offspring can be
       tried out, instead of using the same population as the base each time"""
    
    random.shuffle(individuals)
    
    if nd == 'standard':
        pareto_fronts = _sortNondominated(individuals, k)
    else:
        raise Exception('selNSGA2: The choice of non-dominated sorting '
                        'method "{0}" is invalid.'.format(nd))

    # logger.info('Non-dominated sorting')
    for front in pareto_fronts:
        _assignCrowdingDistFeat(front, featsToUse)
    # logger.info('Number of fronts: %d', len(pareto_fronts))    

    chosen = list(chain(*pareto_fronts[:-1]))
    k = k - len(chosen)
    logger.info('Fronts:%d, Chosen size: %d; %d slots left to fill', len(pareto_fronts), len(chosen), k)

    if k > 0:
        # NOTE: Other options to consider:
        # - indicator-based hypercube volume maximization
        if False:
            # sort first by total fitness, then by crowding distance
            # Use this when there are many local minima, then by applying
            # 'crowind distance' it helps to spread them out
            # No known criteria yet, just trial and error to choose which one to use
            sorted_front = sorted(pareto_fronts[-1], key=lambda x: sum(x.fitness.values))
            sorted_front = sorted(sorted_front, key=attrgetter("fitness.crowding_dist"), reverse=True)
        else:
            # sort first by crowding distance, then by total fitness
            # Should use this as 'error' is the primary criteria
            #  but it may not work if errors are not close to zero
            # sorted_front = sorted(pareto_fronts[-1], key=attrgetter("fitness.crowding_dist"), reverse=True)
            sorted_front = sorted(pareto_fronts[-1], key=attrgetter("fitness.crowding_dist"), reverse=True)
            sorted_front = sorted(sorted_front, key=lambda x: sum(x.fitness.values))
        chosen.extend(sorted_front[:k])

    return chosen

def _sortNondominated(individuals, k, first_front_only=False):
    """Sort the first *k* *individuals* into different nondomination levels
    using the "Fast Nondominated Sorting Approach" proposed by Deb et al.,
    see [Deb2002]_. This algorithm has a time complexity of :math:`O(MN^2)`,
    where :math:`M` is the number of objectives and :math:`N` the number of
    individuals.

    :param individuals: A list of individuals to select from.
    :param k: The number of individuals to select.
    :param first_front_only: If :obj:`True` sort only the first front and
                             exit.
    :returns: A list of Pareto fronts (lists), the first list includes
              nondominated individuals.

    .. [Deb2002] Deb, Pratab, Agarwal, and Meyarivan, "A fast elitist
       non-dominated sorting genetic algorithm for multi-objective
       optimization: NSGA-II", 2002.
    """
    if k == 0:
        return []

    map_fit_ind = defaultdict(list)
    for ind in individuals:
        map_fit_ind[ind.fitness].append(ind)
    fits = list(map_fit_ind.keys())

    current_front = []
    next_front = []
    dominating_fits = defaultdict(int)
    dominated_fits = defaultdict(list)

    # Rank first Pareto front
    for i, fit_i in enumerate(fits):
        for fit_j in fits[i+1:]:
            if fit_i.dominates(fit_j):
                dominating_fits[fit_j] += 1
                dominated_fits[fit_i].append(fit_j)
            elif fit_j.dominates(fit_i):
                dominating_fits[fit_i] += 1
                dominated_fits[fit_j].append(fit_i)
        if dominating_fits[fit_i] == 0:
            current_front.append(fit_i)

    fronts = [[]]
    for fit in current_front:
        fronts[-1].extend(map_fit_ind[fit])
    pareto_sorted = len(fronts[-1])

    # Rank the next front until all individuals are sorted or
    # the given number of individual are sorted.
    if not first_front_only:
        N = min(len(individuals), k)
        while pareto_sorted < N:
            fronts.append([])
            for fit_p in current_front:
                for fit_d in dominated_fits[fit_p]:
                    dominating_fits[fit_d] -= 1
                    if dominating_fits[fit_d] == 0:
                        next_front.append(fit_d)
                        pareto_sorted += len(map_fit_ind[fit_d])
                        fronts[-1].extend(map_fit_ind[fit_d])
            current_front = next_front
            next_front = []

    return fronts

def _assignCrowdingDistFeat(individuals, featsToUse):
    """Assign a crowding distance to each individual's fitness. The
    crowding distance can be retrieved via the :attr:`crowding_dist`
    attribute of each individual's fitness.
    """
    if len(individuals) == 0:
        return

    distances = [0.0] * len(individuals)
    feats = []
    for i, ind in enumerate(individuals):
        feats.append([])
        for featname in featsToUse:
            for feat in ind.fitness.feature_values.values():
                if featsToUse == [] or featname in feat:
                    if feat[featname] == None:
                        feats[i].append(0)
                    else:
                        feats[i].append(feat[featname])
    
    distances = [0.0] * len(individuals)
    crowd = [(ind, i) for i, ind in enumerate(feats)]
    
    nobj = len(feats[0])
    
    # assign a crowding rank determined by removing most
    # crowded point and giving that rank 1, then re-finding nearest neighbour
    # for all points, taking out the next most crowded and giving that rank 2
    # and so on until all points are ranked. Then give larger ranks to points
    # with feature values at the extremities, so keeping them in the population

    inds_left = numpy.arange(0,len(individuals))
    distance_rank = 0

    while len(feats) > 1:
        feats_arr = numpy.asarray(feats)
        pairwise_dists = distance.pdist(feats_arr,'seuclidean')
        pairwise_dists_sq = distance.squareform(pairwise_dists)
        numpy.fill_diagonal(pairwise_dists_sq,1e9)
        min_dists = numpy.min(pairwise_dists_sq,axis=1)
        min_index = numpy.argmin(min_dists)
        distances[inds_left[min_index]] = distance_rank # min_dists[min_index]
        feats.pop(min_index)
        inds_left = numpy.delete(inds_left,min_index)
        distance_rank += 1

    distance_rank += 1
    distances[inds_left[0]] = distance_rank

    # # Set distance of any points with max or min for a particular feature to infinity
    for i in range(nobj):
        crowd.sort(key=lambda element: element[0][i])
        distances[crowd[0][1]] = float("inf")
        distances[crowd[-1][1]] = float("inf")

    for i, dist in enumerate(distances):
        individuals[i].fitness.crowding_dist = dist


