from __future__ import unicode_literals
import os
import sys
import re
#import module for dumping the trained probability matrices to file
import pickle
from itertools import chain
from math import log
from nltk.probability import (ConditionalProbDist, ConditionalFreqDist, SimpleGoodTuringProbDist)
from nltk.util import ngrams
from nltk.model.api import ModelI
from nltk.corpus import brown
from nltk.probability import LidstoneProbDist
from nltk.model import NgramModel
import nltk

#class to hold taxonomic information and n-gram models
class Taxon:
	#instantiates a taxonomic group
	def __init__(self, name, file_names = [], parent=None, children=[], training_data=None):
		self.name = name
		self.parent = parent
		self.children = children
		self.file_names = file_names
		if training_data == None and len(self.file_names) != 0:
			self.training_data = self.trainingData()
		else:
			self.training_data = training_data
		#only do the above if the training data does not exist
		#i.e. not when merging child classes into a parent class
		'''self.trigram_model = Taxon.makeTrigramModel(training_data)'''
		self.bigram_model = self.makeBigramModel()
	
	#uncomment if you wish to use trigram models
	'''def makeTrigramModel(self, training_corpus):
		self.trigram_model = NgramModel(3, training_corpus)'''
	
	#sets this group's bigram model
	def makeBigramModel(self):
		#print training_data
		if self.training_data is not None and len(self.training_data)!= 0 and self.training_data[0] is not None:
			if len(self.training_data[0])>1:
				#merge the datasets if there is a multidimensional list
				self.training_data = sum(self.training_data,[])
				estimator = lambda fdist, bins: LidstoneProbDist(fdist, 0.2) 
				return NgramModel(2, self.training_data, estimator=estimator)
			else:
				estimator = lambda fdist, bins: LidstoneProbDist(fdist, 0.2)
				return NgramModel(2, self.training_data, estimator=estimator)
	
	#reads all the fasta files containing training data and updates training_data
	def trainingData(self):
		data_set = []
		#print self.file_names
		for i in range(0, len(self.file_names)):
			#print data_set
			data_set.append(DataHandler.openFile(self.file_names[i]))
		return data_set
	
	#sets the parent group of this taxonomic group
	def setParent(self, parent):
		self.parent = parent
	
	#sets child groups of this taxonomic group
	def setChildren(self, children):
		self.children = children
	
	#add files to the training data
	#recommended to add all files before generating n-grams for time's sake
	def setFiles(self, file_names):
		self.file_names = file_names
	
	#returns true if this group has child groups; false otherwise
	def hasChildren(self):
		return self.children != None and len(self.children)!=0
	
	'''def setGroupProbability(probability):
		self.groupProbability = probability'''
	
	#creates a user procedure for adding a taxonomic group
	#this adds flexibility
	@staticmethod
	def makeTaxonomicGroup():
		print "TIP: The easiest way to input your taxonomic data is to begin with "\
		+" the lowest-level taxonomic groups and then merge them!"
		merge_class = raw_input("Is this a merger of lower-level classes?")
		if(merge_class):
			createParentClass()
		name = raw_input("Input the name of this taxonomic group: ")
		parent = raw_input("Input the parent group of this taxonomic group: ")
		print "Input child groups. When you are done, type \"Done\""
		children = []
		child = raw.input("Enter Child: ")
		while(child != "Done" and child != "done"):
			children.append(child)
			child = raw.input("Enter Child: ")
		files = []
		print("Enter filenames containing this group's genetic training set. "/
		+"When all files are entered, type \"Done\". Ensure that all files"/
		+" are located in the working directory.")
		file = raw_input("File Name: ")
		while file!= "Done" and file!= "done":
			files.append(file)
			file = raw_input("File Name: ")
		return Taxon(name, parent, children, files)
	
	#allows the user to input a parent taxon that merges child taxa without
	#re-parsing files
	@staticmethod
	def createParentTaxon():
		name = raw_input("Input the name of this taxonomic group:")
		parent = raw_input("Input the parent group of this taxonomic group: ")
		print "Input child groups. When you are done, type \"Done\""
		children = []
		child = raw.input("Enter Child: ")
		while(child != "Done" and child != "done"):
			children.append(child)
			child = raw.input("Enter Child: ")
		return parentTaxon(name, children, parent=parent)
		
	#merges child classes into parent class, returning a new class
	@staticmethod
	def parentTaxon(name, child_classes, parent=None):
		child_data = []
		files = []
		#print child_classes
		for child in child_classes:

			child_data+=child.training_data
		#return a new parent class created from the childrens' training data
		#don't need the files because they are already parsed
		return Taxon(name, file_names = [], parent=parent, children=child_classes, training_data=child_data)

#handles fasta files
class DataHandler:
	#checks whether the file is ok to open; if so, calls parseFile
	@staticmethod
	def openFile(fileName):
		try:
			with open(fileName):
				return DataHandler.parseFile(fileName)
		except IOError:
  			print 'Problem with file. Please check the file and try again'
	
	#parses a fasta file
	@staticmethod
	def parseFile(fileName):
		genome = []
		is_recording = False
		with open(fileName, "r") as the_file:
			for line in the_file:
				if line.find('>')!=-1:
					is_recording = False
				else:
					is_recording = True
				if is_recording:
					genome.append(line)
		sequence = "".join(genome)
		nucleotides = [c.encode('ascii') for c in sequence if c!='\n']	
		return nucleotides

#class to represent a test sequence that must be categorized
class UnknownSequence:
	#reads data from the sequence
	def __init__(self, file_name):
		self.data = DataHandler.openFile(file_name)
		#this is going to be a 2-dimensional array
		self.taxonomy = []
		
	#calculates probability of corpus under model given perplexity and number of words
	def probFromPerplexity(self, perplexity, numWords):
		return numWords*math.log((1/perplexity),2)

	#calculates perplexity of corpus given an n-gram model
	def perplexity(self, n_gram_model):
		#print self.data
		return n_gram_model.perplexity(self.data)
	
	#determines the percent chance of being in a given taxonomic group from n-grams
	def probability(self, n_gram_model):
		x = self.perplexity(n_gram_model)
		return probFromPerplexity(x, len(self.data))
	
	#generates the most likely taxonomic hierarchy for the given sequence
	#default starts with life but you can start lower down
	#if you know that, for example, this organism is a vertebrate
	def categorize(self, taxonomic_group):
		#start with kingdom (it's basically a tree) as a default
		#if there are no more child classes available, break, you are done
		if not taxonomic_group.hasChildren():
			#exit the method
			return None
		else: 
			#get child with max percent, initialize to the first child
			most_likely_child = taxonomic_group.children[0]
			max_probability = self.perplexity(taxonomic_group.children[0].bigram_model)
			#check if other children give lower perplexity
			for i in range(0,len(taxonomic_group.children)):
				this_probability = self.perplexity(taxonomic_group.children[i].bigram_model)
				if this_probability<max_probability:
					most_likely_child = taxonomic_group.children[i]
					max_probability = max_probability
			self.taxonomy.append([most_likely_child, max_probability])
			#need a way to save these data to the object; choose data structure
			return self.categorize(most_likely_child)
			
	#prints out the taxonomy sequentially
	def print_taxonomy(self):
		for group in self.taxonomy:
			print "The perplexity of this organism's genome when compared to "+\
			str(group[0].name)+" is "+str(group[1])


#a class to drive the entire program		
class Controller:
	#initiates a run of the categorizing program
	def __init__(self):
		print "Welcome to my sequence classifier! For each organism you want "\
			+ "to classify, enter the name of the file containing the "\
			+ "organism's genetic information. Files should be in FASTA format. "\
			+ "If you wish to quit the program, enter /'q/'"
		prompt()
	
	#asks the user for filenames of organisms until the user wishes to quit
	@staticmethod
	def prompt():
		input = raw_input("Enter name of file to categorize: ")
		while(input != "q"):
			#create a new unknown sequence object
			this_organism = UnknownSequence(input)
			#classify the organism
			this_organism.classify()
			#print the organism's taxonomy
			this_organism.print_taxonomy()
		print "Thanks for using the classifier!"
	
def main():
	'''Note: main is meant to be edited if the user wishes to have different functionality 
	I have implemented this main so that every time I run my dataset I do not
	have to re-enter the information about each class. The driver class is meant to
	run the program for a user who does not wish to write his/her own main'''
	
	os.chdir("/Users/ameliaantrim/Documents/CS72/ProjectFiles")
	
	#all the animal taxa
	#def __init__(self, name, file_names = [], parent=None, children=None, training_data=None):
	mammals = Taxon("Mammalia", file_names = ["Mammalia.fasta"], parent = "Metazoa", children = [])
	ray_finned_fish = Taxon("Actinopterygii", file_names = ["Oreochromis niloticus.fasta", "Danio rerio.fasta"], \
	parent = "Metazoa", children = [])
	arthropods = Taxon("Arthropoda", ["Anopheles gambiae.fasta", "Drosophila melanogaster.fasta", "Apis mellifera.fasta"], \
	"Metazoa", None)
	amoebae = Taxon("Amoeba", ["Amoeba proteus.txt", "Dictyostelium discoideum.fasta"], "Metazoa", None)
	
	#animal kingdom
	animals = Taxon.parentTaxon("Metazoa", [mammals, ray_finned_fish, arthropods, amoebae], "Eukarya")
	#note this is just a representative few groups of the animal kingdom
	#the user may add any taxa he or she feels appropriate or suspects the test organism may belong to
	
	#fungi
	fungi = Taxon("Fungi", ["Saccharomyces pastorianus.fasta", "Smittium.fasta"], "Eukarya", None)
	
	#plants
	plants = Taxon("Embryophyta", file_names=["Arabidopsis thaliana.fasta", "Populus trichocarpa.fasta"],\
	parent="Eukarya", children=None)
	
	#create eukaryote domain
	eukaryotes = Taxon.parentTaxon("Eukarya", [animals, plants, fungi], parent="Life")

		
	#change the above to a merger class
	bacteria = Taxon("Bacteria", ["Bacillus thuringiensis.fasta"], "Life", None)
	
	#tie them all together at the base of the tree of life
	life = Taxon.parentTaxon("Life", [eukaryotes, bacteria], parent=None)
	
	human_test = UnknownSequence("sequence-11.fasta")
	staph_test = UnknownSequence("Staphylococcus aureus.fasta")
	fish_test = UnknownSequence("Cyprinus carpio.fasta")
	amoeba_test = UnknownSequence("Chaos carolinisus.fasta")
	
	print "The taxonomy of a human (Mammalia) is estimated as follows: "
	human_test.categorize(life)
	human_test.print_taxonomy()
	
	print "\nThe taxonomy of Staphylococcus aureus (Bacteria) is estimated as follows: "
	staph_test.categorize(life)
	staph_test.print_taxonomy()
	
	print "\nThe taxonomy of the common carp (ray-finned fish) is estimated as follows: "
	fish_test.categorize(life)
	fish_test.print_taxonomy()
	
	print "\nThe taxonomy of the giant amoeba is estimated as follows: "
	amoeba_test.categorize(life)
	amoeba_test.print_taxonomy()
	
if __name__ == "__main__":
    main()
	
