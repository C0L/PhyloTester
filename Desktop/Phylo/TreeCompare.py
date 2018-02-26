import dendropy
import pyvolve
import os
import shutil
from subprocess import call


#Using 100 repetitons to reduce variability

#Results
FastTreeURFD = []
RaxMLURFD = []
HPAURFD = []

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
    while (reps < 10):
        #Change the names of the taxa in this tree
        tns = dendropy.TaxonNamespace(["T%d" % i for i in range(1, ntaxa+1)])
        tt = generateTree(tns, ntaxa, seqlen)
        
        #FastTree
        #WARNING! 100.0% NUCLEOTIDE CHARACTERS -- IS THIS REALLY A PROTEIN ALIGNMENT?
        call(["./FastTreeMP"], stdin=file("evolvedsequences.fasta"), stdout=file("FastTreeResult.newick", "w"))
        nt = dendropy.Tree.get_from_path("FastTreeResult.newick", "newick", taxon_namespace=tns)
        sim = dendropy.calculate.treecompare.weighted_robinson_foulds_distance(tt, nt)
        print "FT", sim
        FastTreeURFD.append(sim)
        
        #RaxML
        call(["./raxmlHPC", "-m", "GTRGAMMA", "-s", "evolvedsequences.fasta", "-p", "12345", "-n", "T1"])
        nt = dendropy.Tree.get_from_path("RAxML_bestTree.T1", "newick", taxon_namespace=tns)
        sim = dendropy.calculate.treecompare.weighted_robinson_foulds_distance(tt, nt)
        print "RX", sim
        RaxMLURFD.append(sim)
        
        #Python
        call(["python", "hpa.py", "-f", "evolvedsequences.fasta", "-t", "HPAResults.newick"])
        nt = dendropy.Tree.get_from_path("HPAResults.newick", "newick", taxon_namespace=tns)
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

compareAlg(8, 10)

print "FastTree", FastTreeURFD
print "RaxML", RaxMLURFD
print "HPA", HPAURFD

