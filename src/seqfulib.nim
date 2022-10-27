## Nim routines for processing DNA/RNA/Protein sequences

import sequtils, strutils, math, tables, osproc, streams, zip/gzipfiles
import readfq

type UniversalSeqRecord* = object
    ## This type represents a genetic sequence with optional quality
    name*: string
    comment*: string
    quality*: string
    sequence*: string


proc reverse*(str: string): string =
  for s in str:
    result = s & result

proc reverseComplement*(self:FQRecord):FQRecord =
  # returns reverse complement of sequence
  var revseq = newseq[char]()
  var slen = len(self.sequence)
  for i in 1..slen:
    case self.sequence[slen - i]:
      of 'A':
        revseq.add('T')
      of 'C':
        revseq.add('G')
      of 'G':
        revseq.add('C')
      of 'T', 'U':
        revseq.add('A')
      else:
        revseq.add('N')

  # 0.4.0 - Add reverse of quality string
  FQRecord(name:self.name, comment: self.comment, quality: reverse(self.quality),
         sequence: revseq.join)


proc toFasta*(self:FQRecord, lineLength = 60): string =
  ## returns FASTA formatted string of sequence FQRecord
  var header = ">" & self.name
  if self.comment != "":
    header = header & " " & self.comment
  header & "\n" & map(toSeq(countup(0,self.sequence.len(), lineLength)),
                       proc(x:int):string = 
                         self.sequence[x..x+lineLength-1]).join("\n")

proc qualToChar*(q: int, offset = 33): char =
  ## returns character for a given Illumina quality score
  (q+offset).char
  
proc charToQual*(c: char, offset = 33): int =
  ## returns Illumina quality score for a given character
  c.ord - offset

proc toFastq*(self:FQRecord, qualityValue = 30): string =
  ## returns FASTQ formatted string of sequence FQRecord with given quality
  ## value to be applied to sequence
  var header = "@" & self.name
  var quality = self.quality
  if quality == "":
    quality = strutils.repeat(qualityValue.qualToChar, self.sequence.len)
  if self.comment != "":
    header = header & " " & self.comment
  header & "\n" & self.sequence & "\n+\n" & quality

proc length*(self:FQRecord): int = 
  ## returns length of sequence
  self.sequence.len()

proc kmer2num*(kmer:string):int =
  ## converts a kmer string into an integer 0..4^(len-1)
  let baseVal = {'T': 0, 'C': 1, 'A': 2, 'G': 3, 'U': 0}.toTable
  let klen = len(kmer)
  var num = 0
  for i in 0..(klen - 1):
    try:
      let p = 4^(klen - 1 - i)
      num += p * baseVal[kmer[i]]
    except:
      num = -1
      break
  num

proc num2kmer*(num, klen:int):string =
  ## converts an integer into a kmer string given the number and length of kmer
  let baseVal = {0:'T', 1:'C', 2:'A', 3:'G'}.toTable
  var kmer = repeat(" ",klen)
  var n = num
  for i in 0..(klen - 1):
    let p = 4^(klen - 1 - i)
    var baseNum = int(n/p)
    kmer[i] = baseVal[baseNum]
    n = n - p*baseNum
  kmer

proc toKmerFrequency*(self:FQRecord, klen:int, 
                      includeComplement = false): seq[int] =
  ## returns (overlapping) kmer frequencies of a nucleotide sequence
  var counts = newSeq[int](4^klen)
  var revcomp = self.reverseComplement().sequence

  for i in 0..self.length - klen:
    var kmer = self.sequence[i..i+klen-1]
    var knum = kmer2num(kmer)
    if knum > -1:
      counts[kmer2num(kmer)] += 1
    if includeComplement:
      kmer = revcomp[i..i+klen-1]
      knum = kmer2num(kmer)
      if knum > -1:
        counts[kmer2num(kmer)] += 1
  counts

proc gc*(self:FQRecord): int =
  ## returns the number of bases that are G or C
  self.sequence.count({'G','C','g','c'})

proc ambiguous*(self:FQRecord): int = 
  ## returns the number of bases that are not AGCTU
  self.length - self.sequence.count({'A', 'G', 'C', 'T', 'U',
                                      'a', 'g', 'c', 't', 'u'})

iterator codons(self: FQRecord) : string = 
  var i = 0
  var s = self.sequence.toUpperAscii
  while i < self.length - 2:
    let codon = s[i..i+2]
    if codon.len == 3:
       yield codon
    i += 3

proc translate*(self:FQRecord, code = 1, other = '-'): FQRecord = 
  ## translates a nucleotide sequence with the given genetic code number
  ## see https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi for codes
  var codeMap = 
    ["FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
     "FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSS**VVVVAAAADDEEGGGG",
     "FFLLSSSSYY**CCWWTTTTPPPPHHQQRRRRIIMMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
     "FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
     "FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSSSSVVVVAAAADDEEGGGG",
     "FFLLSSSSYYQQCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
     "", "",
     "FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNNKSSSSVVVVAAAADDEEGGGG",
     "FFLLSSSSYY**CCCWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
     "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
     "FFLLSSSSYY**CC*WLLLSPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
     "FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSSGGVVVVAAAADDEEGGGG",
     "FFLLSSSSYYY*CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNNKSSSSVVVVAAAADDEEGGGG",
     "",
     "FFLLSSSSYY*LCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
     "", "", "", "",
     "FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNNKSSSSVVVVAAAADDEEGGGG",
     "FFLLSS*SYY*LCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
     "FF*LSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
     "FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSSKVVVVAAAADDEEGGGG",
     "FFLLSSSSYY**CCGWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
     "FFLLSSSSYY**CC*WLLLAPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
     "FFLLSSSSYYQQCCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
     "FFLLSSSSYYQQCCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
     "FFLLSSSSYYYYCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
     "FFLLSSSSYYEECC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
     "FFLLSSSSYYEECCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"]
  var code = codeMap[code - 1]
  var transeq = newseq[char]()
  for codon in self.codons:
    let num = kmer2num(codon)
    if num != -1:
      transeq.add(code[num])
    else:
      transeq.add(other)
  FQRecord(name:self.name, comment: self.comment, quality: self.quality,
         sequence: transeq.join)

iterator compressedLines*(filename: string): string =
  ## iterator to read lines of a (maybe) compressed text file transparently
  if filename.find(".gz") > -1:
    for line in lines newGZFileStream(filename, fmRead):
      yield line
  elif filename.find(".bz2") > -1:
    var process = startProcess("bzcat", args=[filename], options={poUsePath})
    var line = ""
    while process.outputStream.readLine(line):
      yield line
    process.close
  else:
    for line in lines filename:
      yield line
    
  
iterator readFasta*(filename: string): FQRecord =
  ## iterator to iterate over the FASTA FQRecords in a file
  var s = FQRecord(name:"", comment:"", sequence:"")
  var seqLines = @[""]
  for line in compressedLines filename:

    if len(line) == 0:        # Fix parsing error on empty lines
      continue
    if line[0] == '>':
      if s.name != "":
        s.sequence = seqLines.join
        yield s
        s.name = ""
        s.comment = ""
        s.sequence = ""
        seqLines = @[]
      var fields = split(line[1..len(line)-1], ' ', 1)
      if len(fields) > 1:
        (s.name, s.comment) = fields
      else:
        s.name = fields[0]
    else:
      seqLines.add(line)
  if s.name != "":
    s.sequence = seqLines.join
    yield s

iterator readFastq*(filename:string): FQRecord =
  ## iterator to iterate over the FASTQ FQRecords in a file
  var s = FQRecord(name:"", comment:"", quality: "", sequence:"")
  var lineNum = 0
  for line in compressedLines filename:
    if lineNum == 0:
      if s.name != "":
        yield s
        s.name = ""
        s.comment = ""
        s.sequence = ""
        s.quality = ""
      var fields = split(line[1..len(line)-1], ' ', 1)
      if len(fields) > 1:
        (s.name, s.comment) = fields
      else:
        s.name = fields[0]
    elif lineNum == 1:
      s.sequence = line
    elif lineNum == 3:
      s.quality = line
    lineNum = (lineNum + 1) mod 4
  if s.name != "":
    yield s
           

iterator readSeqs*(filename:string):FQRecord = 
  for line in compressedLines filename:
    if line[0] == '>':
      for s in readFasta(filename):
        yield s
      break
    elif line[0] == '@':
      for s in readFastq(filename):
        yield s
      break
    else:
      echo "I don't know what type of file is " & filename
