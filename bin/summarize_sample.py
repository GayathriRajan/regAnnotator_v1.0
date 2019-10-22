import argparse, glob, os, re, csv

parser = argparse.ArgumentParser()
parser.add_argument('-d', action="store", dest="sampleDir", help="Sample directory")

#parser.add_argument('-b', action="store", dest="sumBrief", help="Brief Summary File")
#parser.add_argument('-r1', action="store", dest="regionsFile", help=".gt4c File")
#parser.add_argument('-r2', action="store", dest="regionsList", help="Region List File")

args = parser.parse_args()

sumBrief = glob.glob(os.path.join(args.sampleDir,"*.b.sum"))[0]
covFile = glob.glob(os.path.join(args.sampleDir,"*.cov"))[0]
regionsList = glob.glob(os.path.join(args.sampleDir,"*regions.list"))[0]
regionsFile = glob.glob(os.path.join(args.sampleDir,"*.gt4c"))[0]

resFile = re.sub(".b.sum",".sum",sumBrief)

if __name__ == "__main__":
	summaryDict,regionsArr = {},[]
	try:
		with open(regionsFile,'rU') as reg_file:
			for line in reg_file:
				myline=line.replace('"','').replace('\r','').replace('\n','')
				if myline not in '':
					mysplitline=myline.split("\t")
					if line.startswith(">"):
						summaryDict.setdefault(mysplitline[0],{}).setdefault("genLenCov","%s (%s)" % (mysplitline[1],mysplitline[2]))
						summaryDict[mysplitline[0]]["regions"]=regionsArr
						summaryDict[mysplitline[0]]["regionTags"]=[]
						regionsArr=[]
					else:
						regionsArr.append("%s-%s" % (mysplitline[1],mysplitline[2]))

		with open(regionsList,'rU') as reg_list:
			for line in reg_list:
				myline=line.replace('"','').replace('\r','').replace('\n','')
				mysplitline=myline.split("\t")
				if 'NA' not in mysplitline[2]:
					printText = '-'.join([mysplitline[2].replace("Protein:",""),mysplitline[3].replace("Tag:","")])
				elif 'NA' not in mysplitline[1]:
					printText = '-'.join([mysplitline[1].replace("Gene:",""),mysplitline[3].replace("Tag:","")])
				else:
					printText = mysplitline[3].strip("Tag:")

				splitText = printText.split("-")
				myindex = [i for i, e in enumerate(summaryDict[mysplitline[0]]["regionTags"]) if splitText[0] in e]

				if not myindex:
					summaryDict[mysplitline[0]]["regionTags"].append(printText)
				else:
					summaryDict[mysplitline[0]]["regionTags"][myindex[0]]=summaryDict[mysplitline[0]]["regionTags"][myindex[0]]+",%s" % (splitText[1])

		with open(sumBrief,'rU') as sum_brief:
			reader=csv.DictReader(sum_brief,delimiter="\t")
			for line in reader:
				summaryDict[line["Virus"]]["depth"]=int(line["Hits"])
				summaryDict[line["Virus"]]["normDepth"]=float(line["Normalized hits"])

		with open(covFile,'rU') as cov_file:
			reader=csv.reader(cov_file,delimiter="\t")
			for line in reader:
				summaryDict[line[0]]["normGenLenCov"]=line[3]

		with open(resFile, 'w') as res_file:
			res_file.write("Virus\tGenomeLengthCovered\tNormalizedGenomeLengthCovered\tReads\tNormalizedReads\tRegions\tRegionTags\n")
			for k,v in sorted(summaryDict.items(),key=lambda x:x[1]["normDepth"],reverse=True):
				res_file.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (k,v["genLenCov"],v["normGenLenCov"],v["depth"],v["normDepth"],';'.join(v["regions"]),';'.join(v["regionTags"])))

	except Exception, e:
		print "Error: %s" % e