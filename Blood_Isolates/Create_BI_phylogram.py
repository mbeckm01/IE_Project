#!/usr/bin/python3

import Bio as Bio
from Bio import SeqIO
from Bio import AlignIO
from Bio import Phylo
from Bio.Phylo.TreeConstruction import DistanceCalculator
from Bio.Phylo.TreeConstruction import DistanceTreeConstructor
import matplotlib
import matplotlib.pyplot as plt 


t2= SeqIO.read('P001.TXT','fasta')
t3= SeqIO.read('P002-1.TXT','fasta')
t4= SeqIO.read('P003.TXT','fasta')
t5= SeqIO.read('P004.TXT','fasta')
t6=SeqIO.read('P005.txt','fasta')
t7=SeqIO.read('P006.txt','fasta')
t8=SeqIO.read('P007.txt','fasta')
t9=SeqIO.read('P008.txt','fasta')
t10=SeqIO.read('P009.txt','fasta')


t2.id='CM002_Serratia_marcescens'
t3.id='CM006_Staphylococcus_saprophyticus'
t4.id='CM007_Cardiobacterium_hominis'
t5.id='CM008_Streptococcus_oralis'
t6.id='CM010_Enterococcus_faecalis'
t7.id='CM013_Staphylococcus_aureus'
t8.id='CM017_Staphylococcus_aureus'
t9.id='CM019_Staphylococcus_agalactiae'
t10.id='CM020_Streptococcus_mutans'


#combine all individual sequences
blood_isolates = SeqIO.write([t2,t3,t4,t5,t6,t7,t8,t9,t10], 'blood_isolates2.fasta', 'fasta')

#alignment completed using MUSCLE Clustal: https://www.ebi.ac.uk/Tools/msa/muscle/

with open('/users/micaelabeckman/blood_alignment2.aln') as aln:
	alignment=AlignIO.read(aln,'clustal')
print(type(alignment))

calculator=DistanceCalculator('identity')
distance_matrix=calculator.get_distance(alignment)
print(distance_matrix)

constructor=DistanceTreeConstructor(calculator)

blood_tree=constructor.build_tree(alignment)
blood_tree.rooted=True 
print(blood_tree)

Phylo.write(blood_tree,"blood_isolate_tree2.xml",'phyloxml')

fig=Phylo.draw(blood_tree)

fig=plt.figure(figsize=(13,5),dpi=100)
matplotlib.rc('font',size=10)
matplotlib.rc('xtick',labelsize=10)
matplotlib.rc('ytick',labelsize=10)
axes=fig.add_subplot(1,1,1)
Phylo.draw(blood_tree,axes=axes)
fig.savefig('Blood_isolate_cladogram')

Phylo.convert('blood_isolate_tree2.xml','phyloxml','blood_nexus2.nex', 'nexus')
blood_nexus=Phylo.read('blood_nexus2.nex', 'nexus')
fig=plt.figure(figsize=(13,5),dpi=100)
matplotlib.rc('font',size=9)
matplotlib.rc('xtick',labelsize=10)
matplotlib.rc('ytick',labelsize=10)
axes=fig.add_subplot(1,1,1)
Phylo.draw(blood_nexus,axes=axes)
fig.savefig('blood_isolates_cladogram2')