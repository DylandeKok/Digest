import Bio
from Bio import Restriction
from Bio.Seq import Seq
from Bio.Alphabet.IUPAC import IUPACAmbiguousDNA
from Bio import SeqIO
from Bio.Restriction import AllEnzymes


fastafile = input("welk fasta bestand wil je openen?\n")
if len (fastafile) < (1):
    (fastafile)= (r'C:\\Users\Acer\Desktop\\python\biopython\\test.fasta.txt')

#restrictie enzym kiezen
rb=Bio.Restriction.RestrictionBatch([]) #maak een lege REBatch
VN=input("met hoeveel Enzymen wil je digesteren?\n")
LN=0
while LN < int(VN):
    D=input("restrictie enzym:\n") 
    LN=(LN)+1
    rb.add(D) #voegt het ingeputt RE toe aan REstrictie Batch
print(rb.elements(),"\n")


### seqenties van meerdere fastafiles inlezen in een loop
genlen=0
totNsite=0
for record in SeqIO.parse(fastafile,"fasta"): #fasta openen
    id=record.id
    seq=record.seq
    lengte=len(record)
    genlen=(genlen)+lengte
    #print(id, lengte, "bp")
    ###Analyse
    searchsite=rb.search(seq)
    Nsite=sum(len(v) for v in searchsite.values())
    #print(Nsite)
    
    
    #Ana=Bio.Restriction.Analysis(rb, record.seq, linear=True)
    #Ana.print_as('number')
    #Ana.print_that()
    #print("\n")
    totNsite=totNsite+Nsite
    



print('total genome length:', genlen)
print('totaal aantal restrictie sites:', totNsite)
print('Gemiddelde Fragment lengte:',(genlen/totNsite))
    
Anser=input("wil je de volledige afbeelding zien? Y/N \n")

if (Anser) == "Y":
    for record in SeqIO.parse(fastafile,"fasta"): #fasta openen
        id=record.id
        seq=record.seq
        lengte=len(record)
        print(id, lengte, "bp")
    ###Analyse
        Ana=Bio.Restriction.Analysis(rb, record.seq, linear=True)
        Ana.print_as('map')
        Ana.print_that()
        print("\n")
else:
    input("Pres enter to quit")


#Bio.SeqIO.write(Ana,r'C:\Users\Acer\Desktop\python\biopython\Digest\\genomeDigest.txt', 'swiss')
