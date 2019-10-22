import argparse, csv, os, re, shutil
from subprocess import call
from Bio import Entrez, SeqIO
from BCBio import GFF
import ftputil
from itertools import islice,izip_longest

parser = argparse.ArgumentParser()
parser.add_argument('-i', action="store", dest="idFile", help="RefSeq IDs list")
args = parser.parse_args()


def parse_GFF_features(gffFile):
    gffOutFile = gffFile.replace('.gff', '_parsed.gff')
    outFile = gffFile.replace('.gff','_parser.log')
    IDdict, parentList = {}, []

    with open(outFile, 'w') as out_file:
        try:
            with open(gffFile) as gff_file, open(gffOutFile, 'w') as gff_out_file:
                prevfeature = ""
                for line in gff_file:
                    if not line.startswith("#"):
                        splitline = line.strip("\r\n").split("\t")
                        myfeature = re.split(';|=', splitline[8])
                        if (myfeature[1] != 'id0') and ("source" not in splitline[2]):
                            IDdict[myfeature[1]] = dict(izip_longest(*[iter(myfeature[2:])] * 2, fillvalue=""))
                            # out_file.write("%s" % myfeature)
                            # out_file.write("%s\n" % IDdict[myfeature[1]])

                            # if 'gene' in prevfeature and 'cds' not in myfeature[1]:
                            # gff_out_file.write("\n")

                            outfeature = ""

                            if 'id' in myfeature[1]:
                                if "Dbxref" in IDdict[myfeature[1]]:
                                    outfeature = "%(Dbxref)s;" % (IDdict[myfeature[1]])
                                else:
                                    outfeature = outfeature + "GeneID:%s;" % myfeature[1]

                                outfeature = outfeature + "%s-%s;" % (IDdict[myfeature[1]]["gbkey"], myfeature[1])

                                if "rpt_type" in IDdict[myfeature[1]]:
                                    outfeature = outfeature + "Type:%(rpt_type)s;" % (IDdict[myfeature[1]])

                                if "rpt_family" in IDdict[myfeature[1]]:
                                    outfeature = outfeature + "Family:%(rpt_family)s;" % (IDdict[myfeature[1]])

                                if "rpt_unit_seq" in IDdict[myfeature[1]]:
                                    outfeature = outfeature + "Unit_seq:%(rpt_unit_seq)s;" % (IDdict[myfeature[1]])

                                if "regulatory_class" in IDdict[myfeature[1]]:
                                    outfeature = outfeature + "Type:%(regulatory_class)s;" % (IDdict[myfeature[1]])

                                if "number" in IDdict[myfeature[1]]:
                                    outfeature = outfeature + "Exon_#:%(number)s;" % (IDdict[myfeature[1]])

                                if "product" in IDdict[myfeature[1]]:
                                    outfeature = outfeature + "Protein:%(product)s;" % (IDdict[myfeature[1]])

                                outfeature = outfeature + "Note:%s\n" % IDdict[myfeature[1]].get("Note", "NA")

                            elif 'gene' in myfeature[1]:
                                parentList.append(myfeature[1])
                                if "Dbxref" not in IDdict[myfeature[1]]:
                                    IDdict[myfeature[1]]["Dbxref"] = "GeneID:%s" % myfeature[1]

                                outfeature = "%(Dbxref)s;%(gbkey)s;Gene:%(Name)s;Biotype:%(gene_biotype)s\n" % IDdict[
                                    myfeature[1]]


                            elif 'rna' in myfeature[1]:
                                parentList.append(myfeature[1])
                                if not "Dbxref" in IDdict[myfeature[1]]:
                                    IDdict[myfeature[1]]["Dbxref"] = "GeneID:%s" % IDdict[myfeature[1]].get("Parent",
                                                                                                            myfeature[
                                                                                                                1])

                                outfeature = "%(Dbxref)s;" % (IDdict[myfeature[1]])

                                outfeature = outfeature + "%(gbkey)s;" % IDdict[myfeature[1]]

                                # if "gene" in IDdict[myfeature[1]]:
                                # outfeature = outfeature + "Gene:%(gene)s;" % (IDdict[myfeature[1]])

                                if "product" in IDdict[myfeature[1]]:
                                    outfeature = outfeature + "Protein:%(product)s;" % (IDdict[myfeature[1]])

                                outfeature = outfeature + "Note:%s\n" % IDdict[myfeature[1]].get("Note", "NA")


                            elif 'cds' in myfeature[1]:
                                splitDbxref = IDdict[myfeature[1]]["Dbxref"].replace('Genbank', 'ProteinID').split(",")
                                pID = ''.join([x for x in splitDbxref if "ProteinID" in x])
                                gID = ''.join([x for x in splitDbxref if "GeneID" in x])

                                if not gID:
                                    gID = "GeneID:%s" % IDdict[myfeature[1]].get("Parent", myfeature[1])

                                # outfeature = "%s;%s;" % (gID, IDdict[myfeature[1]]["gbkey"])
                                outfeature = "%s;%s;" % (gID, myfeature[1].upper())

                                # if "gene" in IDdict[myfeature[1]]:
                                # outfeature = outfeature + "Gene:%(gene)s;" % (IDdict[myfeature[1]])

                                outfeature = outfeature + "Protein:%s;" % IDdict[myfeature[1]].get("product", "NA")
                                outfeature = outfeature + "%s\n" % pID

                            else:
                                parentList.append(myfeature[1])
                                if "Dbxref" in IDdict[myfeature[1]]:
                                    outfeature = "%(Dbxref)s;" % (IDdict[myfeature[1]])

                                outfeature = outfeature + "%(gbkey)s;" % IDdict[myfeature[1]]

                                # if "gene" in IDdict[myfeature[1]]:
                                # outfeature = outfeature + "Gene:%(gene)s;" % (IDdict[myfeature[1]])

                                if "product" in IDdict[myfeature[1]]:
                                    outfeature = outfeature + "Protein:%(product)s;" % (IDdict[myfeature[1]])

                                outfeature = outfeature + "Note:%s\n" % IDdict[myfeature[1]].get("Note", "NA")

                            prevfeature = myfeature[1]
                            splitline[8] = outfeature
                            line = '\t'.join(splitline)
                            gff_out_file.write(line)
                    else:
                        gff_out_file.write(line.strip("\r"))

        except Exception, e:
            out_file.write("%s\n" % e)


def convertNCBIgb2gff(sID):

    #sID = "LK931492.1"
    gb_file = os.path.join(annoDir, sID + ".gb")
    gff_file = os.path.join(annoDir, sID + ".gff")
    corr_gff_file = os.path.join(annoDir, sID + "_corr.gff")

    Entrez.email = "gayathri_rajan@rush.edu"

    handle = Entrez.efetch(db="nucleotide", id=sID, rettype="gb", retmode="text")
    seq_record = SeqIO.read(handle, "genbank")
    handle.close()

    gb_out_handle = open(gb_file, "w")
    SeqIO.write(seq_record, gb_out_handle, "genbank")
    gb_out_handle.close()

    gff_out_handle = open(gff_file, "w")
    gb_in_handle = open(gb_file)
    GFF.write(SeqIO.parse(gb_in_handle, "genbank"), gff_out_handle)
    gb_in_handle.close()
    gff_out_handle.close()

    gff_in_handle = open(gff_file)
    corr_out_handle = open(corr_gff_file, "w")

    idct, genect, cdsct, rnact = 0, 0, 0, 0
    geneDict, geneIDs = [], []

    for line in gff_in_handle:

        if not line.startswith("#"):
            line = line.replace("db_xref", "Dbxref")
            line = line.replace("note", "Note")
            splitline = line.strip().split("\t")

            if len(splitline)>8:
                myfeature = re.split(';|=', splitline[8])
                lineDict = dict(izip_longest(*[iter(myfeature)] * 2, fillvalue=""))

                if "remark" in line:
                    splitline = []

                elif 'gene' in splitline[2]:
                    splitline[8] = "ID=gene%s;gbkey=%s;Name=%s;gene_biotype=protein_coding;%s" % (
                        genect, splitline[2].title(), lineDict.get("gene", "NA"), splitline[8])
                    genect += 1
                    geneDict.append(lineDict.get("gene", "NA"))
                    geneIDs.append(lineDict.get("Dbxref", "GeneID:gene%s" % (genect - 1)))

                elif 'CDS' in splitline[2]:
                    if "Dbxref" in splitline[8]:
                        splitline[8] = splitline[8].replace("Dbxref=", "Dbxref=Genbank:%s," % lineDict["protein_id"])
                        splitline[8] = "ID=cds%s;gbkey=%s;Parent=%s;%s" % (
                            cdsct, splitline[2], "gene%s" % (genect - 1), splitline[8].split(";translation")[0])
                    else:
                        splitline[8] = "ID=cds%s;gbkey=%s;Parent=%s;Dbxref=Genbank:%s;%s" % (
                            cdsct, splitline[2], "gene%s" % (genect - 1), lineDict["protein_id"],
                            splitline[8].split(";translation")[0])
                    cdsct += 1

                elif 'RNA' in splitline[2]:
                    splitline[8] = "ID=rna%s;gbkey=%s;Parent=%s;%s" % (
                        rnact, splitline[2], "gene%s" % (genect - 1), splitline[8])
                    rnact += 1

                elif (genect > 0) and (geneDict[genect - 1] in lineDict.get("gene", "NA")):
                    splitline[8] = "ID=id%s;gbkey=%s;Dbxref=%s;%s" % (
                    idct, splitline[2], geneIDs[genect - 1], splitline[8])
                    idct += 1

                else:
                    splitline[8] = "ID=id%s;gbkey=%s;%s" % (idct, splitline[2], splitline[8])
                    idct += 1

                line = "\t".join(splitline)
                if not re.match(r'^\s*$', line):
                    corr_out_handle.write(line + "\n")

        else:
            corr_out_handle.write(line)

    gff_in_handle.close()
    corr_out_handle.close()

    os.remove(gb_file)
    os.rename(corr_gff_file, gff_file)
    return gff_file


def getViralProjID(sID):
    handle1 = Entrez.esearch(db="assembly", term=sID)
    seq_record1 = Entrez.read(handle1)
    handle1.close()

    handle2 = Entrez.esummary(db="assembly", id=','.join(seq_record1[u'IdList']))
    seq_record2 = Entrez.read(handle2)
    handle2.close()
    myObj = seq_record2["DocumentSummarySet"]["DocumentSummary"][0]

    if myObj["AssemblyName"]:
        return (myObj["Organism"],myObj["AssemblyName"].split('Proj',1)[-1])

if __name__=="__main__":

    annoDir = os.path.join(os.path.dirname(os.path.dirname(args.idFile)),"NCBI-GFF")
    if not os.path.exists(annoDir):
        os.makedirs(annoDir)

    host = ftputil.FTPHost('ftp.ncbi.nlm.nih.gov', 'anonymous', 'password')
    host.chdir('/genomes/Viruses/')
    dir_list = host.listdir(host.curdir)

    with open(args.idFile) as id_file:
        id_file = csv.reader(id_file, delimiter="\t")
        for line in id_file:
            ncid = line[0]
            print ncid
            if not ncid.startswith("#"):
                gff_file = os.path.join(annoDir, ncid + ".gff")
                parsedGff = os.path.join(annoDir, ncid + "_parsed.gff")
                parsedGffBed = os.path.join(annoDir, ncid + "_parsed.sort.bed")
                if not os.path.isfile(gff_file):
                    try:
                        (name,uid)=getViralProjID(ncid.split(".")[0])
                        sfile = os.path.join(host.getcwd(), ''.join([x for x in dir_list if x.endswith("uid" + uid)]),
                                             ncid.split(".")[0] + ".gff")
                        host.download(sfile, gff_file, None)
                        with open(gff_file) as gff_file_in:
                            fileid = list(islice(gff_file_in, 6))[-1].strip().split("\t")[0]
                            if ncid in fileid:
                                print "Match %s" % fileid
                                parse_GFF_features(gff_file)
                                call("grep -v \"^#\" %s | awk '$9 != \"\" {print}' | awk -F $'\\t' '{OFS=\"\\t\"; print $1,$4,$5,$9}' | awk -F $'\\t' -v OFS=\"\\t\" '{if ($4 ~ /^*[[:space:]]*$/){gsub(\"[[:space:]]\",\"@\",$4); print; next}{print;}}' | sortBed -i - > %s" % (parsedGff,parsedGffBed),shell=True)
                            else:
                                print "Mismatch %s" % fileid

                    except:
                            srcfile = convertNCBIgb2gff(ncid)
                            shutil.move(srcfile,gff_file)
                            parse_GFF_features(gff_file)
                            call("grep -v \"^#\" %s | awk '$9 != \"\" {print}' | awk -F $'\\t' '{OFS=\"\\t\"; print $1,$4,$5,$9}' | awk -F $'\\t' -v OFS=\"\\t\" '{if ($4 ~ /^*[[:space:]]*$/){gsub(\"[[:space:]]\",\"@\",$4); print; next}{print;}}' | sortBed -i - > %s" % (
                                parsedGff, parsedGffBed), shell=True)
                            print "Match %s" % ncid