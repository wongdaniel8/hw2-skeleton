from .utils import Atom, Residue, ActiveSite
import random

def compute_similarity(site_a, site_b):
    """
    Compute the similarity between two given ActiveSite instances.

    Input: two ActiveSite instances
    Output: the similarity between them (a floating point number)
    """
    a_list = [] #list of only residue type and no number
    for el in site_a.residues:
        a_list.append(el.type)
    b_list = []
    for el in site_b.residues:
        b_list.append(el.type)

    ##simple intersection
    a_list = set(a_list)
    b_list = set(b_list)
    overlaps = [x for x in a_list if x in b_list]
    similarity = len(overlaps) / float(min(len(a_list), len(b_list)))
    return similarity


    ##find common subset of both lists, with replacement
    # overlaps = []
    # for el in a_list:
    #     if el in b_list: #list removal issue, not iterating enough 
    #         # print("GGG", a_list, b_list)
    #         overlaps.append(el)
    #         b_list.remove(el)
    # similarity = len(overlaps) / float(max(len(site_a.residues), len(site_b.residues)))
    # return similarity


def cluster_by_partitioning(active_sites):
    """
    Cluster a given set of ActiveSite instances using a partitioning method (k means).

    Input: a list of ActiveSite instances
    Output: a clustering of ActiveSite instances
            (this is really a list of clusters, each of which is list of
            ActiveSite instances)
    """

    def getClusters(assignments, centerIndices):
        clusters = []
        for center in centerIndices:
            centerList = []
            for i in range(0, len(assignments)):
                if assignments[i] == center:
                    centerList.append(active_sites[i])
            clusters.append(centerList)
        return clusters

    def assignToCenters(centerIndices):
        """
        generate a list of length |number of active_sites| with each element of the list being that index's assigned center index
        """
        assignments = []
        for i in range(0, len(active_sites)):
            closestCenter = -100
            shortestDistance = 200
            for j in centerIndices:
                if i == j:
                    closestCenter = j
                    break
                else:
                    distance = 1 - compute_similarity(active_sites[i], active_sites[j])
                    if distance < shortestDistance:
                        closestCenter = j
                        shortestDistance = distance
            assignments.append(closestCenter)
        return assignments

    #initialize:
    k = 3
    iterations = 100
    centerIndices = []
    assignments = []
    centerIndices = random.sample(range(0, len(active_sites)), k)
    assignments = assignToCenters(centerIndices)
    clusters = getClusters(assignments, centerIndices)
    # print("initial score", qualityMetric(clusters))


    #repeat for # of iterations
    for i in range(0, iterations):
        centerIndices = []
        for c in clusters:
            clusterCenter = -1
            shortestCumulativeDistance = 1000
            for candidate in c: #find new center
                cumulativeDistance = 0
                for other in c:
                    cumulativeDistance += compute_similarity(candidate, other)
                if cumulativeDistance <= shortestCumulativeDistance:
                    clusterCenter = candidate
                    shortestCumulativeDistance = cumulativeDistance
            centerIndices.append(active_sites.index(clusterCenter))
        
        #assign each active site to closest center
        assignments = assignToCenters(centerIndices)
        clusters = getClusters(assignments, centerIndices)

    # print("final clusters", clusters)
    # print("final score", qualityMetric(clusters))
    return clusters


def cluster_hierarchically(active_sites):
    """
    Cluster the given set of ActiveSite instances using a hierarchical algorithm.                                                                  #

    Input: a list of ActiveSite instances
    Output: a list of clusterings
            (each clustering is a list of lists of Sequence objects)
    """
    def findNearestClusters(clusters):
        x, y = -1, -1
        shortestDistance = 1000000
        for i in range(0, len(clusters)):

            for j in range(0, len(clusters)):
                minDistance = 100000000    
                if i != j:
                    #find shortest distance between clusters at index i and index j
                    for elementi in clusters[i]:
                        for elementj in clusters[j]:
                            distance = 1 - compute_similarity(elementi, elementj)
                            if distance < minDistance:
                                minDistance = distance 
                    if minDistance < shortestDistance:
                        shortestDistance = minDistance
                        x , y = i, j
        return x, y
    def join(clusters, i, j):
        """
        joins clusters at index i and j
        """
        newClusters = []
        for r in range(0, len(clusters)):
            if r != i and r != j:
                newClusters.append(clusters[r])
        new = clusters[i] + clusters[j]
        newClusters.append(new)
        return newClusters

    k = 3
    clusters = []
    for act in active_sites:
        clusters.append([act])
    # print(clusters)

    while len(clusters) != k:
        i, j = findNearestClusters(clusters)
        # print("i j ", i, j)
        clusters = join(clusters, i, j)
        # print("ccc", clusters)

    # print("final", clusters)
    return clusters

    return []


def qualityMetric(clustering):
    """
    clustering is a list of clusters, each of which is a list of ActiveSite instances
    returns a measure of how well the clustering worked, with 0 being the lowest and larger return values indicating 
    a better clustering. Within each cluster I will compute the pairwise similarity between each active site. I will then 
    add up each sum from each cluster to get a total score.
    """
    score = 0
    for cluster in clustering:
        for element in cluster:
            for other in cluster:
                score += compute_similarity(element, other)
    return score












