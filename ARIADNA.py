import csv
import numpy as np
import math
from random import randint
from sklearn.externals import joblib
import os
import subprocess
import sys

if len(sys.argv) != 4 or sys.argv[1]=='h' or sys.argv[1]=='help' or sys.argv[1]=='-h' or sys.argv[1]=='-help':
	print "\nARIADNA (ARtificial Intelligence for Ancient DNA) is a machine learning based approach for detecting single-nucleotide variants in read alignments of ancient DNA samples.\n"
	print "Usage: ARIADNA_MP.py input.bam reference.fasta outfile\n"
	print "bam index and reference index required in same directory as bam and reference input\n"
	print "GROM-ANC_v1_0_2 and Tree_ARIADNA.pkl must be in the same directory as ARIADNA.py\n"
	print "See README.txt\n"
	quit()


bamloc=sys.argv[1]
refloc=sys.argv[2]
outloc=sys.argv[3]

try:
	print "Identifying PSNVs"
	subprocess.call(["./GROM-ANC_v1_0_2 -i "+bamloc+" -r "+refloc+" -o GRO1TD.tmp.txt -f"],shell=True)
	print "Collecting features"
	subprocess.call(["./GROM-ANC_v1_0_2 -t GRO1TD.tmp.txt -i "+bamloc+" -r "+refloc+" -o GRI1TD.tmp.txt -f"],shell=True)
	min_dist_between_snv = 0
	snv_prob_index = 5
	ref_base_index = 81
	snv_base_index = 82
	snv_ratio_index = 83
	snv_count_init = 84
	snv_t_index = [84,85,86,87]
	snv_normal_index = [88,89,90,91]
	snv_normal_indel_index1 = 37
	snv_normal_indel_index2 = 96
	dna_upper = ["A","C","G","T"]
	dna_lower = ["a","c","g","t"]
	sheet=open(outloc,"w")
	ow=csv.writer(sheet,delimiter='\t')
	ow.writerow(["#CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO"])
	trainingsetfinal=[]
	trainingcall=[]
	predict17=[]
	TPs={}
	PScall=[]
	clf = joblib.load('Tree_ARIADNA.pkl')
	count=0	
	last_chr = ""
	last_pos = 0
	full=[]
	print "Calling variants"
	genomes=["GRI1TD.tmp.txt"]
	for gen in genomes:
		NN={}
		read=csv.reader(open(gen,"r"),delimiter='\t')	
		for i in read:
			if i[0]=="SNV":
				NN[i[1].upper(),int(i[2])]=1
		count=0
		read=csv.reader(open(gen,"r"),delimiter='\t')
		for i in read:
			nf=0
			if i[0]=="SNV":
				x=i[5:9]+i[10:26]+i[30:59]+i[84:110]+i[111:138]
				for j in range(0,len(x)):
					try:
						x[j]=float(x[j])
					except ValueError:
						break
					if math.isnan(x[j])==True:
						nf=1
				x=x+[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]
				
				for rf in [1,7,8,9,10,11,12,13,14,15,16,17,18,19,32,33,34,35,36,45,46]:
					x[rf]=0
				
				if nf==0:
					for nn in range(int(i[2])-100,int(i[2])+100):
						try:
							x[122]+=NN[i[1].upper(),nn]
						except KeyError:
							x[122]+=0				
					if i[138][0]=="A" or i[138][0]=="a":
						x[112]=1
					if i[138][0]=="T" or i[138][0]=="t":
						x[113]=1	
					if i[138][0]=="C" or i[138][0]=="c":
						x[114]=1
					if i[138][0]=="G" or i[138][0]=="g":
						x[115]=1
					if i[138][2]=="A" or i[138][2]=="a":
						x[116]=1
					if i[138][2]=="T" or i[138][2]=="t":
						x[117]=1	
					if i[138][2]=="C" or i[138][2]=="c":
						x[118]=1
					if i[138][2]=="G" or i[138][2]=="g":
						x[119]=1				
					if i[ref_base_index] in ["A","T","C","G"]:
						if i[81]=="A" or i[81]=="a":
							x[-9]=1
						if i[81]=="T" or i[81]=="t":
							x[-8]=1	
						if i[81]=="C" or i[81]=="c":
							x[-7]=1
						if i[81]=="G" or i[81]=="g":
							x[-6]=1
						if i[81] in ['A','T','C','G']:
							x[-1]=1
						predict17.append(x)
						count+=1										
					if i[ref_base_index] not in ["A","T","C","G"]:
						if i[82]=="A" or i[82]=="a":
							x[-5]=1
						if i[82]=="T" or i[82]=="t":
							x[-4]=1
						if i[82]=="C" or i[82]=="c":
							x[-3]=1
						if i[82]=="G" or i[82]=="g":
							x[-2]=1
						if i[81] in ['A','T','C','G']:
							x[-1]=1
						predict17.append(x)
						count+=1
					full.append(i)				
	for k in range(0,len(predict17)):
		for l in [2,3,5,6,7,8,9,12,13,16,17,18,19,20,21,22,23,27,28,29,32,33,37,38,42,43,45,46,49,50,51,52,53,54,55,56,57,58,59,60,61,62,63,67,68,69,70,71,72,73,74,77,78,81,82,89,90,96,97]:
			predict17[k][l]=predict17[k][l]/predict17[k][4]
	for i in range(0,len(predict17)):
		z=clf.predict(predict17[i])[0]
		if z > .05 and (clf.predict_proba(predict17[i])[0][1] > 0.05):
			if full[i][82].upper() == 'A':
				ow.writerow([full[i][1],full[i][2],'.',full[i][81].upper(),full[i][82].upper(),full[i][120],'.','DP='+str(full[i][10])+";"+"AF="+str(round(float(full[i][84])/float(full[i][10]),4))+";"+"AC="+str(full[i][84])])
			if full[i][82].upper() == 'C':
				ow.writerow([full[i][1],full[i][2],'.',full[i][81].upper(),full[i][82].upper(),full[i][120],'.','DP='+str(full[i][10])+";"+"AF="+str(round(float(full[i][85])/float(full[i][10]),4))+";"+"AC="+str(full[i][85])])				
			if full[i][82].upper() == 'G':
				ow.writerow([full[i][1],full[i][2],'.',full[i][81].upper(),full[i][82].upper(),full[i][120],'.','DP='+str(full[i][10])+";"+"AF="+str(round(float(full[i][86])/float(full[i][10]),4))+";"+"AC="+str(full[i][86])])				
			if full[i][82].upper() == 'T':
				ow.writerow([full[i][1],full[i][2],'.',full[i][81].upper(),full[i][82].upper(),full[i][120],'.','DP='+str(full[i][10])+";"+"AF="+str(round(float(full[i][87])/float(full[i][10]),4))+";"+"AC="+str(full[i][87])])				
				
	sheet.close()
except:
	print "Error: please check bam reference and output"	

subprocess.call(["rm GRO1TD.tmp.txt"],shell=True)
subprocess.call(["rm GRO1TD.tmp.txt.ctx"],shell=True)
subprocess.call(["rm GRI1TD.tmp.txt"],shell=True)
subprocess.call(["rm GRI1TD.tmp.txt.ctx"],shell=True)
