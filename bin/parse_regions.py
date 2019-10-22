import argparse,csv,re,os
from collections import OrderedDict
from itertools import izip_longest
import urllib,urllib2

parser = argparse.ArgumentParser()
parser.add_argument('-i', action="store", dest="resBed", help="List of Genes in Viral Regions")
parser.add_argument('-c', action="store", dest="covFile", help="Coverage File for Viral Regions")
args = parser.parse_args()

def getUniprotInfo(npIDs):
    print "npIDs=%s" % npIDs
    url = 'https://www.uniprot.org/uploadlists/'
    resDict = OrderedDict()

    params = {
        'from': 'P_REFSEQ_AC',
        'to': 'ACC',
        'format': 'tab',
        'query': npIDs,
        'columns': 'reviewed,genes(PREFERRED),protein names,organism'
    }

    data = urllib.urlencode(params)
    request = urllib2.Request(url, data)
    contact = ""  # Please set a contact email address here to help us debug in case of problems (see https://www.uniprot.org/help/privacy).
    request.add_header('User-Agent', 'Python %s' % contact)

    try:
        response = urllib2.urlopen(request, timeout=60)
        page = response.read()

        for item in page.split("\n"):
            if (item not in '') and (not item.startswith(('your'))):
                item=item.split("\t")
                item.insert(0,item[-2])
                del item[-2:]
                item=[re.sub(r"\((strain.*?)\)",r"\1",x) for x in item]
                item=[re.sub(' \(.*\)', '', x) if x else 'NA' for x in item]
                for protID in item[0].split(','):
                    if protID not in resDict or ('reviewed' in item):
                        resDict.setdefault(protID,{}).update({"Prot":{"Protein":item[2]},"Gene":{"Gene":item[3]}})

    except urllib2.URLError as e:
        print "getUniprotID Error: %s" % e

    return resDict


if __name__ == "__main__":

    covDict,RegDict,GeneDict,gffIDs={},OrderedDict(),OrderedDict(),{}
    npList=[]

    annoDir = os.path.join(os.path.dirname(os.path.dirname(args.resBed)), "NCBI-GFF")

    with open(args.covFile,'rU') as cov_file:
        reader1=csv.reader(cov_file,delimiter="\t")
        for line in reader1:
            if line:
                orgID = line[0].split("|")[3]
                if line[0].startswith(">"):
                    covDict.setdefault(orgID,{}).update({"Name":line[0]})
                else:
                    covDict.setdefault(orgID,{}).update({"%s-%s" % (line[1],line[2]):int(line[3])})

    with open(args.resBed,'rU') as resBed_file:
        reader2=csv.reader(resBed_file,delimiter="\t")
        for line in reader2:
            virus=line[0]
            region= "%s-%s" % (line[1],line[2])
            parts=line[3:]

            for part in parts:
                #print part
                RegDict.setdefault(virus, OrderedDict()).setdefault(region, {})
                GeneDict.setdefault(virus,{})

                if part not in '':
                    partArr=re.split(":|;",part)
                    partArr.insert(2,'GenePart')
                    partArr=[x.replace('@',' ') for x in partArr]
                    geneID=partArr[1]
                    partID=partArr[3]
                    #print virus,partArr, len(partArr)

                    partDict = dict(izip_longest(*[iter(partArr[4:])] * 2, fillvalue=""))

                    #print virus,geneID,partID,partDict
                    protein = partDict.get("Protein")
                    if protein:
                        #print protein
                        if len(protein.split(" ")) == 1 or 'hypothetical' in protein:
                            #print protein,partDict["ProteinID"]
                            if partDict.get("ProteinID") and (not partDict.get("ProteinID") in npList):
                                npList.append(partDict["ProteinID"])

                    #print npList

                    if geneID not in RegDict[virus][region]:
                        RegDict[virus][region].setdefault(geneID, [])

                    if partID not in RegDict[virus][region][geneID]:
                        RegDict[virus][region][geneID].append(partID)

                    RegDict[virus][region][geneID] = sorted(RegDict[virus][region][geneID])
                    GeneDict[virus].setdefault(geneID, {}).setdefault(partID, partDict)

    #print RegDict
    #print GeneDict

    #print npList,"\n",len(npList)
    uniprotRes=getUniprotInfo(' '.join(npList))
    #print uniprotRes

    with open(args.resBed.replace("_res_bed.txt", "_region.bed"), "w") as outfile1:
        for k1,v1 in GeneDict.iteritems():
            #print k1,v1
            gffIDs.setdefault(k1,{})
            for k2,v2 in v1.iteritems():
                if v2.keys()==['Gene']:
                    gffIDs[k1].setdefault(k2,GeneDict[k1][k2]["Gene"])
                    gffIDs[k1][k2].update({"GeneID":k2,"Protein":"NA","ProteinID":"NA","Tag":"Gene","Misc":"NA"})
                    gffIDs[k1][k2].pop("Biotype",None)
                if not v2.get("Gene"):
                    v2.setdefault("Gene",{"Gene":"NA"})
                for k3,v3 in v2.iteritems():
                    if 'CDS' in k3:
                        temp1=uniprotRes.get(v3.get("ProteinID"),{})
                        v3.update(temp1.get("Prot",{}))
                        if v2.get("Gene").get("Gene")=="NA":
                            v2.get("Gene").update(temp1.get("Gene",{}))

                    #temp2=k3.split("-")[1] if "-" in k3 else k3
                    temp2=k2
                    temp2=temp2.upper()

                    temp3 = {"Protein": "NA", "ProteinID": "NA"}
                    temp3.update(GeneDict[k1][k2]["Gene"])
                    temp3.update({"GeneID": k2 if k2.isdigit() else "NA"})
                    temp3.pop('Biotype', None)
                    temp3.update({"Tag": k3})
                    if 'CDS' not in k3:
                        temp3.update({"TagDesc": ','.join('{}:{}'.format(key, val) for key, val in sorted(GeneDict[k1][k2][k3].iteritems())).replace("Note:","")})
                    else:
                        temp3.update(GeneDict[k1][k2][k3])
                        temp3.update({"TagDesc":"NA"})

                    if (temp2 not in gffIDs[k1]) and (temp2.lower() not in gffIDs[k1]) and ("Gene" not in k3):
                        gffIDs[k1].setdefault(temp2,temp3)

            #print k1, gffIDs[k1]


            with open(os.path.join(annoDir,k1+"_parsed.gff"),'rU') as annot_file:
                reader3=csv.reader(annot_file,delimiter="\t")
                for line in reader3:
                    if not line[0].startswith("#"):
                        searchID=line[8].split(";")[0].strip("GeneID:").upper()
                        searchTag=line[8].split(";")[1]
                        #print searchID,searchTag
                        if (searchID in gffIDs[k1]) and (gffIDs[k1][searchID]["Tag"]==searchTag):
                            # and ('repeat_region' not in searchTag)
                            #outfile1.write('\t'.join([covDict[k1]["Name"],'\t'.join(line[3:5]),';'.join('{}:{}'.format(key,val) for key,val in sorted(gffIDs[k1][searchID].iteritems())).replace(" ","@")+";Strand:"+line[6]])+"\n")
                            outfile1.write('\t'.join([covDict[k1]["Name"], '\t'.join('{}:{}'.format(key, gffIDs[k1][searchID].get(key,"NA")) for key in ['Gene','Protein','Tag','TagDesc']) + "\tStrand:" + line[6]]) + "\n")



    with open(args.resBed.replace("_res_bed.txt","_region.sum"),"w") as outfile2:
        outfile2.write("Virus\tRegion\tLength(Region)\tMax Coverage(Region)\tGene-Protein\tGeneID-ProteinID\tCDS#\tOther Parts of Gene")

        orderedSam=OrderedDict(sorted(covDict.items(), key=lambda x:x[1]['Name'].split("|")[4])).keys()
        orderedSamMap={k:i for i,k in enumerate(orderedSam)}
        RegDict = OrderedDict(sorted(RegDict.items(), key=lambda x:orderedSamMap[x[0]]))

        #print orderedSam

        for k1,v1 in RegDict.iteritems():
            for k2,v2 in v1.iteritems():
                 regLen=int(k2.split("-")[1])-int(k2.split("-")[0])
                 #prefix="\n%s\t%s\t%s\t%s\t%.2f\t" % (covDict[k1].get("Name"),k2,regLen,covDict[k1].get(k2),covDict[k1].get(k2,0)*150.00/regLen)
                 prefix = "\n%s\t%s\t%s\t%s\t" % (covDict[k1].get("Name"), k2, regLen, covDict[k1].get(k2)) #1.00 * covDict[k1].get(k2, 0) / regLen)
                 outfile2.write("%s%s" % (prefix,'\n\t\t\t\t'.join(["Gene:%s;Protein:%s\tGeneID:%s;ProteinID:%s\t%s\t%s" %
                         (GeneDict[k1][k3].get("Gene",{}).get("Gene","NA"),
                          GeneDict[k1][k3].get(''.join([x for x in v3 if 'CDS' in x]),{}).get("Protein","NA"),
                          ''.join([k3 if k3.isdigit() else "NA"]),GeneDict[k1][k3].get(''.join([x for x in v3 if 'CDS' in x]),{}).get("ProteinID","NA"),
                          (','.join([x for x in v3 if ('CDS' in x) or ('Gene' in x)]) or "NA"),
                          (';'.join([x.split("-")[0] + " ("+','.join('{}:{}'.format(key,val) for key,val in sorted(GeneDict[k1][k3].get(x,"NA").iteritems())) + ")" for x in v3 if ('Gene' not in x) and ('CDS' not in x)]) or "NA"))
                           for (k3,v3) in sorted(v2.iteritems())])))

        #print sorted(RegDict["NC_001699.1"]['398-556'])
        #print GeneDict["NC_001699.1"]["1489517"]



