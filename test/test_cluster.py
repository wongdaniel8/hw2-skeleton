from hw2skeleton import cluster
from hw2skeleton import io
import os
import random

def test_similarity():
    filename_a = os.path.join("data", "276.pdb")
    filename_b = os.path.join("data", "4629.pdb")
    activesite_a = io.read_active_site(filename_a)
    activesite_b = io.read_active_site(filename_b)
    assert cluster.compute_similarity(activesite_a, activesite_b) == cluster.compute_similarity(activesite_b, activesite_a)
    assert cluster.compute_similarity(activesite_a, activesite_a) == 1
    assert cluster.compute_similarity(activesite_b, activesite_b) == 1

    
def test_partition_clustering():
    # tractable subset
    # pdb_ids = [276, 4629, 10701]
    pdb_ids = [276, 1806, 3458, 3733, 10814, 4629, 10701]
    active_sites = []
    for id in pdb_ids:
        filepath = os.path.join("data", "%i.pdb"%id)
        active_sites.append(io.read_active_site(filepath))
    assert [] not in cluster.cluster_by_partitioning(active_sites)
    assert len(cluster.cluster_by_partitioning(active_sites)) == 3

    pdb_ids = [276]
    active_sites = []
    for id in pdb_ids:
        filepath = os.path.join("data", "%i.pdb"%id)
        active_sites.append(io.read_active_site(filepath))
    assert len(cluster.cluster_by_partitioning(active_sites)) == 1


def test_hierarchical_clustering():
    # tractable subset
    # pdb_ids = [276, 4629, 10701]
    pdb_ids = [276, 1806, 3458, 3733, 10814, 4629, 10701]
    active_sites = []
    for id in pdb_ids:
        filepath = os.path.join("data", "%i.pdb"%id)
        active_sites.append(io.read_active_site(filepath))
    assert [] not in cluster.cluster_by_partitioning(active_sites)
    assert len(cluster.cluster_by_partitioning(active_sites)) == 3

    pdb_ids = [276]
    active_sites = []
    for id in pdb_ids:
        filepath = os.path.join("data", "%i.pdb"%id)
        active_sites.append(io.read_active_site(filepath))
    assert len(cluster.cluster_by_partitioning(active_sites)) == 1

def test_compareTwoMethods():
    allPossible = []
    for filename in os.listdir("data"):
        allPossible.append(int(filename.split(".")[0]))
    hierTotal = 0
    partTotal = 0
    iterations = 100
    numPDBs = 15
    for i in range(0, iterations):
        indices = random.sample(range(0, len(allPossible)), numPDBs)
        pdb_ids = []
        for j in indices:
            pdb_ids.append(allPossible[j])
        active_sites = []
        for id in pdb_ids:
            filepath = os.path.join("data", "%i.pdb"%id)
            active_sites.append(io.read_active_site(filepath))

        hierScore = cluster.qualityMetric(cluster.cluster_hierarchically(active_sites))
        partScore = cluster.qualityMetric(cluster.cluster_by_partitioning(active_sites))

        hierTotal += hierScore
        partTotal += partScore
    print("hierScoreAverage: ", hierTotal / float(iterations))
    print("partScoreAverage: ", partTotal / float(iterations))        


