#!/usr/bin/python3

import sys, getopt, random
import time
import subprocess
import errno
import numpy as np
import os.path
import matplotlib.pyplot as plt
from multiprocessing import Process, freeze_support, set_start_method

def AccuracyToQScore(accuracyPC):
	return -10 * np.log10(1 - (accuracyPC/100.))

def ParseQScore(QScoreString):
	qstring = QScoreString.strip()
	totalError = 0
	for char in qstring:
		score = ord(char) - 33
		error = pow(10, -score/10.0)
		totalError += error

	return totalError / len(qstring)

class AlignmentStats:
	def __init__(self, name, refName, readLength):
		self.name = name
		self.refName = refName
		self.readLength = readLength
		self.alignmentLength = 0
		self.matches = 0
		self.mismatches = 0
		self.insertions = 0
		self.deletions = 0
		self.QScore = 0
		self.QScoreAccuracy = 0

	def ParseCFString(self, cfString):
		i = 0
		while i < len(cfString):
			if cfString[i] == ":":
				j = i + 1
				while j < len(cfString) and cfString[j].isdigit():
					j += 1
				self.matches += int(cfString[i+1:j])
				self.alignmentLength += int(cfString[i+1:j])
				i = j
			elif cfString[i] == "*":
				self.mismatches += 1
				self.alignmentLength += 1
				i += 3
			elif cfString[i] == "-":
				j = i + 1
				while j < len(cfString) and cfString[j] in {"a", "c", "g", "t", "n"}:
					j += 1
				self.deletions += (j-i)
				self.alignmentLength += (j-i)
				i = j
			elif cfString[i] == "+":
				j = i + 1
				while j < len(cfString) and cfString[j] in {"a", "c", "g", "t"}:
					j += 1
				self.insertions += (j-i)
				self.alignmentLength += (j-i)
				i = j
			else:
				print("Error: Could not parse character " + cfString[i] + " in cf string " + cfString)
				sys.exit(2)

	def PrintStats(self):
		print("ReadName: " + self.name)
		print("Accuracy: " + str(100 * self.matches / self.alignmentLength) + "%")

	def GetAccuracy(self):
		return 100 * self.matches / self.alignmentLength

	def GetAlignmentLength(self):
		return self.alignmentLength

def runMinimap2(refFilename, readFilename, PAFFilename, outfileName):
	try:
		with open(PAFFilename, "w") as paffile:
			with open(outfileName, "w") as outfile:
				subprocess.run(["minimap2", "-c", "--cs", refFilename, readFilename], stdout=paffile, stderr=outfile ,check=True)
	except subprocess.CalledProcessError as e:
	    print(e.returncode)
	    print(e.output)	

if __name__ == '__main__':
	readFilenames = list()
	runNames = list()
	refFilename = ''
	showPlots = False
	minReadLength = 0
	minReadQuality = 0
	plt.style.use('ggplot')

	try:
		opts, args = getopt.getopt(sys.argv[1:],"i:r:l:sn:q:")
	except getopt.GetoptError:
		print("Option not recognised.")
		print("QualityByAlignment.py -i <readfile1.fq,readfile2.fq,...,readfileN.fq> -n <run1name,run2name,...,runNname> -r <reference fasta>")
	for opt, arg in opts:
		if opt == "-h":
			print("QualityByAlignment.py -i <readfile1.fq,readfile2.fq,...,readfileN.fq> -n <run1name,run2name,...,runNname> -r <reference fasta>") 
			sys.exit()
		elif opt in ("-i"):
			readFilenames = arg.split(",")
		elif opt in ("-s"):
			showPlots = True
		elif opt in ("-r"):
			refFilename = arg
		elif opt in ("-l"):
			minReadLength = int(arg)
		elif opt in ("-n"):
			runNames = arg.split(",")
		elif opt in ("-q"):
			minReadQuality = int(arg)

	if readFilenames == '' or refFilename =='':
		print("Error: Must provide (PAF filename) or (read filename and reference filename)")
		sys.exit(2)

	if len(runNames) != len(readFilenames):
		print("Error: You must specify a run name for each read file, in the same order.")
		sys.exit(2)

	if(len(set(runNames)) != len(runNames)):
		print("Error: Run names must be unique.")
		sys.exit(2)

	processes = list()
	PAFFilenames = dict()
	i = 0
	for readFilename in readFilenames:
		prefix = runNames[i]
		PAFFilename = prefix + "_alignment.paf"
		outfileName = prefix + "_minimap2.out"
		PAFFilenames[prefix] = PAFFilename
		if not os.path.exists(PAFFilename):
			p = Process(target=runMinimap2, args=(refFilename, readFilename,PAFFilename, outfileName,))
			p.start()
			processes.append(p)
		i += 1

	for p in processes:
	    p.join()

	allStatsDict = dict()
	allRefSequenceList = dict()
	minAccuracy = 100

	i = 0
	for readFilename in readFilenames:
		prefix = runNames[i]
		PAFFilename = PAFFilenames[prefix]
		allStatsDict[prefix] = dict()
		allRefSequenceList[prefix] = list()

		alignmentStatsDict = allStatsDict[prefix]
		refSequenceList = allRefSequenceList[prefix]
		try: 
			with open(PAFFilename, 'r') as PAFFile:
				for line in PAFFile:
					fields = line.split()
					readName = fields[0]
					readLength = int(fields[1])
					alignmentLength = int(fields[3]) - int(fields[2])
					refName = fields[5]
					if float(alignmentLength) / readLength >= 0.9:
						rs = AlignmentStats(readName, refName, readLength)
						rs.ParseCFString(fields[-1][5:])
						if rs.GetAccuracy() < minAccuracy:
							minAccuracy = rs.GetAccuracy()
						if readName not in alignmentStatsDict.keys():
							alignmentStatsDict[readName] = rs
						else:
							if rs.GetAlignmentLength() > alignmentStatsDict[readName].GetAlignmentLength():
								alignmentStatsDict[readName] = rs

						if refName not in refSequenceList:
							refSequenceList.append(refName)

		except (OSError, IOError) as e: 
			if getattr(e, 'errno', 0) == errno.ENOENT:
				print("Could not open file " + PAFFilename)
				sys.exit(2)

		refSequenceList.sort()

		totalReads = 0
		try:
			with open(readFilename, 'r') as readFile:
				count = 0
				readName = ""
				for line in readFile:
					if count % 4 == 0:
						fields = line.split()
						readName = fields[0][1:]
					if count % 4 == 3:
						if readName in alignmentStatsDict.keys():
							AvgReadError = ParseQScore(line.strip())
							rs = alignmentStatsDict[readName]
							rs.QScoreAccuracy = 100*(1 - AvgReadError)
							rs.QScore = -10 * np.log10(AvgReadError)
					count += 1
				totalReads = count / 4
		except (OSError, IOError) as e: 
			if getattr(e, 'errno', 0) == errno.ENOENT:
				print("Could not open file " + readFilename)
				sys.exit(2)
		i += 1


	accuracies = list()
	lengths = list()
	QScoreAccuracies = list()

	accuraciesByRef = dict()
	lengthsByRef = dict()
	QScoreAccuraciesByRef = dict()

	i = 0
	for readFilename in readFilenames:
		prefix = runNames[i]
		alignmentStatsDict = allStatsDict[prefix]
		accuraciesByRef[prefix] = dict()
		lengthsByRef[prefix] = dict()
		QScoreAccuraciesByRef[prefix] = dict()

		accuracyList = list()
		lengthList = list()
		QScoreList = list()

		refSequenceList = allRefSequenceList[prefix]
		alignmentStatsDict = allStatsDict[prefix]
		meanAccuracy = 0
		meanRefAccuracy = dict()

		for ref in refSequenceList:
			meanRefAccuracy[ref] = [0,0]
			lengthsByRef[prefix][ref] = list()
			accuraciesByRef[prefix][ref] = list()
			QScoreAccuraciesByRef[prefix][ref] = list()


		for key in alignmentStatsDict.keys():
			rs = alignmentStatsDict[key]
			accuracy = rs.GetAccuracy()
			length = rs.GetAlignmentLength()
			qscoreAcc = rs.QScoreAccuracy
			ref = rs.refName

			if rs.readLength >= minReadLength and rs.QScore >= minReadQuality:
				meanAccuracy += accuracy
				meanRefAccuracy[rs.refName][0] += accuracy
				meanRefAccuracy[rs.refName][1] += 1

				accuracyList.append(accuracy)
				accuraciesByRef[prefix][ref].append(accuracy)
				lengthList.append(length)
				lengthsByRef[prefix][ref].append(length)
				QScoreList.append(qscoreAcc)
				QScoreAccuraciesByRef[prefix][ref].append(qscoreAcc)

		accuracies.append(accuracyList)
		lengths.append(lengthList)
		QScoreAccuracies.append(QScoreList)

		print("Stats for run " + prefix + ": ")
		print("Mean alignment accuracy: " + str(meanAccuracy / len(alignmentStatsDict.keys())) + "\t= Q" + str(AccuracyToQScore(meanAccuracy / len(alignmentStatsDict.keys()))))
		for ref in refSequenceList:
			print("Mean alignment accuracy for " + ref + ": " + str(meanRefAccuracy[ref][0]/meanRefAccuracy[ref][1]))

		i += 1

	#----- Plots --------------------------------------
	#----- Accuracy Box Plot per Run ------------------

	minY = max(minAccuracy - 5, 0)
	plt.boxplot(accuracies, labels=runNames, showfliers=False)
	plt.title("Alignment Accuracies")
	plt.ylabel("Accuracy (%)")
	plt.savefig("./accuracy_all_runs.pdf", bbox_inches='tight')
	if showPlots:
		plt.show()
	plt.clf()

	plt.boxplot(QScoreAccuracies, labels=runNames, showfliers=False)
	plt.title("QScore Accuracies")
	plt.ylabel("Accuracy (%)")
	plt.savefig("./qscore_all_runs.pdf", bbox_inches='tight')
	if showPlots:
		plt.show()
	plt.clf()

	#----- Accuracy Box Plot per Ref ----
	fig, axes = plt.subplots(1, len(runNames), figsize=(20,10))
	i = 0
	for ax in axes:
		prefix = runNames[i]
		accuraciesForPlot = list()
		refSequenceList = allRefSequenceList[prefix]
		for ref in refSequenceList:
			accuraciesForPlot.append(accuraciesByRef[prefix][ref])
		ax.boxplot(accuraciesForPlot, labels=refSequenceList, showfliers=False)
		ax.set_title(runNames[i])
		ax.set_ylabel("Accuracy (%)")
		ax.set_ylim(minY, 102)
		i += 1

	fig.savefig("./accuracy_by_ref.pdf")

	if showPlots:
		plt.show()
	plt.clf()

	#----- Length vs Accuracy Plot ------
	fig, axes = plt.subplots(1, len(runNames), figsize=(15,5))
	i = 0
	for ax in axes:
		lengthList = lengths[i]
		accuracyList = accuracies[i]
		h2d = ax.hist2d(lengthList, accuracyList, bins=100, cmap="Reds", density=True)
		ax.set_xlabel("Alignment Length")
		ax.set_ylabel("Alignment Accuracy")
		ax.set_title(runNames[i])
		ax.set_ylim(minY, 100)
		i += 1
		fig.colorbar(h2d[3], ax=ax)

	if showPlots:
		plt.show()
	fig.savefig("./length_accuracy.pdf")
	plt.clf()

	#----- QScore Accuracy vs Alignment Accuracy Plot -----
	fig, axes = plt.subplots(1, len(runNames), figsize=(15,5))
	i = 0
	for ax in axes:
		h2d = ax.hist2d(QScoreAccuracies[i], accuracies[i], bins=100, cmap="Reds", density=True, range=[[minY,100],[minY,100]])
		ax.grid(True, color = "grey", linewidth = "1.4", linestyle = "-.")
		ax.set_xlabel("Qscore accuracy")
		ax.set_ylabel("Alignment Accuracy")
		ax.set_title(runNames[i])
		ax.set_ylim(minY, 100)
		i += 1 
		fig.colorbar(h2d[3], ax=ax)

	if showPlots:
		plt.show()
	fig.savefig("./qscore_accuracy.pdf")
	plt.clf()

	#------ Accuracy against filtered Read Quality --------
	fig, axes = plt.subplots(1, len(runNames), figsize=(15,5))
	i = 0
	for readFilename in readFilenames:
		ax = axes[i]
		prefix = runNames[i]
		alignmentStatsDict = allStatsDict[prefix]
		accuracyByFilteredQ = list()
		for j in range(20):
			accuracyByFilteredQ.append(list())

		for key in alignmentStatsDict.keys():
			rs = alignmentStatsDict[key]
			if rs.readLength > minReadLength:
				accuracy = rs.GetAccuracy()	
				for j in range(20):
					if rs.QScore >= j:
						accuracyByFilteredQ[j].append(accuracy)

		ax.boxplot(accuracyByFilteredQ, labels=range(20), showfliers=False)
		ax.set_title(runNames[i])
		ax.set_ylabel("Alignment Accuracy (%)")
		ax.set_xlabel("Min Read Quality")
		ax.set_ylim(minY, 102)
		i += 1

	if showPlots:
		plt.show()
	fig.savefig("./accuracy_filteredQ.pdf")
	plt.clf()














