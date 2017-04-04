from Bio import Entrez

# genes without info --> 19 y 18
listaPrueba = ["18", "100128607", "100129363", "19"]

Entrez.email = "sample@example.org"


# option 1 -- two functions --> need to fix genes without sequence info

def getGeneLocusAndIntervalId(geneId):
    handle = Entrez.efetch(db="gene", id=geneId, retmode="xml")
    record = Entrez.read(handle)
    geneLoci = record[0]["Entrezgene_locus"]
    geneRegion = geneLoci[0]["Gene-commentary_seqs"][0]["Seq-loc_int"]["Seq-interval"]
    startPos = int(geneRegion["Seq-interval_from"]) + 1
    endPos = int(geneRegion["Seq-interval_to"]) + 1
    intervalId = geneRegion["Seq-interval_id"]["Seq-id"]["Seq-id_gi"]
    strandSense = geneRegion["Seq-interval_strand"]["Na-strand"].attributes["value"]
    return intervalId, startPos, endPos, strandSense


def getFastaFromNCBI(intervalId, startPos, endPos, strandSense):
    if strandSense.lower() == "minus":
        strandSense = 2
    else:
        strandSense = 1
    handle = Entrez.efetch(db="nucleotide", rettype="fasta", id=intervalId,
                           seq_start=startPos, seq_stop=endPos, strand=strandSense)
    return handle.read()


# option 1 -- single function

def ncbiRetriever(listOfIds):

        Entrez.email = "notfunny@notanemail.org"
        mapOfIdSequence = {}
        for geneId in listOfIds:
            geneId = str(geneId)
            geneId = geneId[:-2]
            try:
                handle1 = Entrez.efetch(db="gene", id=geneId, retmode="xml")
            except Exception:
                print("gene "+geneId + " not found")
                continue
            record = Entrez.read(handle1)
            geneLoci = record[0]["Entrezgene_locus"]
            try:
                geneRegion = geneLoci[0]["Gene-commentary_seqs"][0]["Seq-loc_int"]["Seq-interval"]
            except KeyError:
                continue
            startPos = int(geneRegion["Seq-interval_from"]) + 1
            endPos = int(geneRegion["Seq-interval_to"]) + 1
            intervalId = geneRegion["Seq-interval_id"]["Seq-id"]["Seq-id_gi"]
            strandSense = geneRegion["Seq-interval_strand"]["Na-strand"].attributes["value"]
            strandSense = 2 if strandSense.lower() == "minus" else 1
            handle1.close()

            handle2 = Entrez.efetch(db="nucleotide", id=intervalId, rettype="fasta", retmode="text", seq_start=startPos,
                                    seq_stop=endPos, strand=strandSense)
            fasta = handle2.read()
            fastaString = str(fasta)
            fastaWithoutId = '\n'.join(fastaString.split('\n')[1:])
            fastaOneLine = fastaWithoutId.replace('\n', '')

            mapOfIdSequence[geneId] = fastaOneLine
            handle2.close()
        return mapOfIdSequence


# gi_id, start, end, strand = getGeneLocusAndIntervalId("100128607")
# print(getFastaFromNCBI(gi_id, start, end, strand))


mapT = ncbiRetriever(listaPrueba)
print(mapT)