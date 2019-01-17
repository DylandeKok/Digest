import Bio
from Bio import Restriction
from Bio.Seq import Seq
from Bio.Alphabet.IUPAC import IUPACAmbiguousDNA
from Bio import SeqIO
from Bio.SeqUtils import GC
from Bio.Restriction import ApeKI

#fasta bestand inlezen
fastafile = input("welk fasta bestand wil je openen?\n")
if len (fastafile) < (1):
    (fastafile)= (r'C:\\Users\Acer\Desktop\\python\biopython\\test.fasta.txt')

for record in SeqIO.parse(fastafile,"fasta"): #fasta uit printen
    print("Record ID:",min(record.id),"tot", max(record.id)) #het printen van de grootste en kleinste record id
   # print(repr(record.seq))
    print("Lengte genoom:",min(len(record)), "tot", max(len(record)))
    amb = IUPACAmbiguousDNA()
    myseq = (Seq(str(record.seq))) # string van sequence maken
    gccount=GC(myseq)
    print("GC =",gccount, "%") # %GC berekenen
    #RE=input("kies een Restrictie enzym\n ") uiteindelik zo maken dat je RE kan kiezen
    print("ApeKI :",Bio.Restriction.ApeKI.site)
    RElist=(Bio.Restriction.ApeKI.catalyse(myseq, linear=True))
    REn=0
    for REsite in RElist:
        REn=(REn)+1
    print("aantal fragmenten:", REn) 
    print("Gem. lengte per fragment:", len(record)/(REn),"\n")
   # print("restrictie site", REn ,":\n",REsite,)
    
