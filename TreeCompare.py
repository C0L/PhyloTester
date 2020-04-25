import dendropy
import pyvolve
import os
import shutil
import numpy
import time
from subprocess import call


#Using 100 repetitons to reduce variability

#Results
FastTreeURFD = []
RaxMLURFD = []
HPAURFD = []
#Time results
TFastTreeURFD = []
TRaxMLURFD = []
THPAURFD = []


def NormalizeMaxTreeDepth(t, targettreelen):
    l=GetLongestDendropyRootToTipLength(t.seed_node)
    print "depth from root to tip is", l
    print "will normalize total tree length to", targettreelen
    ScaleDendropyBranchLengths(t.seed_node, targettreelen/l)
    l=GetLongestDendropyRootToTipLength(t.seed_node)
    print "depth from root to tip is now", l
    
def AllPathsToLeaf(t, debug=False):
    if t.taxon is not None:
        return [t.edge_length]
    lens=[]
    for c in t.child_nodes():
        #print "on child node", c, c.taxon, c.edge_length
        #print "all paths to this node:", aplc
        for pathlen in AllPathsToLeaf(c):
            if t.edge_length is None:
                print "WARNING: AllPathsToLeaf got None edge length; setting to 0.0"
                print "probably first split (seed node) and may be no problem"
                el=0.0
            else: el=t.edge_length
            lens.append(el+pathlen)
    print "returning lens:", lens, len(lens)
    return lens
    
def NormalizeAverageTreeDepth(t, targettreelen):

    l=AverageTreeDepth(t.seed_node)
    print "average depth from root to tip is", l
    print "will normalize average tree length to", targettreelen
    ScaleDendropyBranchLengths(t.seed_node, targettreelen/l)
    l=AverageTreeDepth(t.seed_node)
    print "average depth from root to tip is now", l
    
def AverageTreeDepth(tseed, debug=False):
    #paths=[]
    # make a list of every distinct path of nodes from root to tips
    paths=AllPathsToLeaf(tseed, debug=debug)
    #print "paths:"
    #for p in paths:
    #    print p
    return sum(paths)/len(paths)


    
def GetLongestDendropyRootToTipLength(n):
    """
        on a rooted tree, with all branch lengths present, find longest depth of one path from root to tip of tree.
    """
    #print "LEN: nodes: ", len(n.child_nodes())
    ls=[GetLongestDendropyRootToTipLength(cn) for cn in n.child_nodes()]
    ls.append(0.0)      # prevent error max() on an empty sequence
    #print "LEN:", ls, max(ls)
    if n.edge.length is None:
        print "GetLongestDendropyRootToTipLength: WARNING edge length is None. Probably was the branch to first split on rooted tree."
        return max(ls)
    else:
        return n.edge.length + max(ls)


def ScaleDendropyBranchLengths(n, bl):
    # on a rooted tree, with all branch lengths present,
    # adjust branch lengths by factor bl
    #print "adjusting n edge from", n.edge.length
    if n.edge.length==None:
        print "ScaleDendropyBranchLengths: WARNING 0 length edge. Probably was the branch to first split on rooted tree."
    else:
        n.edge.length*=bl
    #print "to", n.edge.length
    for c in n.child_nodes():
        #c.edge.length*=bl      will be modified when we call on child
        ScaleDendropyBranchLengths(c, bl)


def generateTree(tns, ntaxa, seqlen):
    #Construct the tree and save as newick file
    t = dendropy.simulate.treesim.birth_death_tree(birth_rate=1.0, death_rate=0, taxon_namespace=tns, num_extant_tips=ntaxa)
    t.write(path='/tmp/pyvt', schema='newick', suppress_rooting=True, suppress_internal_node_labels=True)
    
    #Set pyvolve data type
    m1 = pyvolve.Model("nucleotide")
    p1 = pyvolve.Partition(models=m1, size=seqlen)
    
    #Read tree from dendropy
    pot = pyvolve.read_tree(file='/tmp/pyvt')
    
    #Simulate evolution with no save file
    e1 = pyvolve.Evolver(tree=pot, partitions=p1)
    e1(seqfile=None)
    
    seqs = e1.get_sequences()
    
    ds=dendropy.DnaCharacterMatrix.from_dict(seqs, taxon_namespace=tns)
    ds.write(path="evolvedsequences.fasta", schema="fasta")
    #print ds
    return t

def compareAlg(ntaxa, seqlen):
    reps = 0
    while (reps < 3):
        #Change the names of the taxa in this tree
        tns = dendropy.TaxonNamespace(["T%d" % i for i in range(1, ntaxa+1)])
        tt = generateTree(tns, ntaxa, seqlen)
        tt.print_plot()
        #FastTree
        #WARNING! 100.0% NUCLEOTIDE CHARACTERS -- IS THIS REALLY A PROTEIN ALIGNMENT?
        start_time = time.time()
        call(["./FastTreeMP"], stdin=file("evolvedsequences.fasta"), stdout=file("FastTreeResult.newick", "w"))
        TFastTreeURFD.append(time.time() - start_time)
        nt = dendropy.Tree.get_from_path("FastTreeResult.newick", "newick", taxon_namespace=tns)
        sim = dendropy.calculate.treecompare.weighted_robinson_foulds_distance(tt, nt)
        print "FT", sim
        FastTreeURFD.append(sim)
        
        #RaxML
        start_time = time.time()
        call(["./raxmlHPC-PTHREADS", "-m", "GTRCATX", "-V", "-s", "evolvedsequences.fasta", "-p", "12345", "-n", "T1"])
        TRaxMLURFD.append(time.time() - start_time)
        nt = dendropy.Tree.get_from_path("RAxML_bestTree.T1", "newick", taxon_namespace=tns)
        sim = dendropy.calculate.treecompare.weighted_robinson_foulds_distance(tt, nt)
        print "RX", sim
        RaxMLURFD.append(sim)
       
        #HPA
        start_time = time.time()
        call(["python", "hpa.py", "-f", "evolvedsequences.fasta", "-t", "HPAResults.newick"])
        THPAURFD.append(time.time() - start_time)
        nt = dendropy.Tree.get_from_path("HPAResults.newick", "newick", taxon_namespace=tns)
        print nt.as_string(schema="newick",)
        NormalizeAverageTreeDepth(nt, 1)
        sim = dendropy.calculate.treecompare.weighted_robinson_foulds_distance(tt, nt)
        print "HPA", sim
        HPAURFD.append(sim)

        #Cleanup
        os.remove("evolvedsequences.fasta")
        os.remove("evolvedsequences.fasta.reduced")
        os.remove("HPAResults.newick")
        os.remove("FastTreeResult.newick")
        os.remove("RAxML_bestTree.T1")
        os.remove("RAxML_info.T1")
        os.remove("RAxML_log.T1")
        os.remove("RAxML_parsimonyTree.T1")
        os.remove("RAxML_result.T1")
        os.remove("site_rates_info.txt")
        os.remove("site_rates.txt")
        
        reps+=1

compareAlg(32, 10000)

print "FASTREE:", FastTreeURFD, "Average data", numpy.mean(FastTreeURFD), "Average time", numpy.mean(TFastTreeURFD)
print "RAXML:", RaxMLURFD, "Average data", numpy.mean(RaxMLURFD), "Average time", numpy.mean(TRaxMLURFD)
print "HPA:", HPAURFD, "Average data", numpy.mean(HPAURFD), "Average time", numpy.mean(THPAURFD)

