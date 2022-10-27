import unittest, os
include "../src/seqfulib.nim"

test "Does Reverse complement work":
  var s = FQRecord(name: "foo", sequence:"ATGATC")
  check(s.reverseComplement.sequence == "GATCAT")


test "Does FASTQ Reverse complement work":
  var s = FQRecord(name: "foo", sequence:"ATG", quality: "ABC")
  check(s.reverseComplement.quality == "CBA")

test "Does kmer frequency work":
    var s = FQRecord(sequence:"ATGC")
    check(s.toKmerFrequency(1) == @[1,1,1,1])
    check(s.toKmerFrequency(1, true) == @[2,2,2,2])

test "Does translation work":
      var s = FQRecord(sequence:"TTGAGCCTCGCCGTTACGCTCGCCTCTACCA")
      check(s.translate.sequence == "LSLAVTLAST")
      var s2 = FQRecord(sequence:"CACCCTTCCCCTCCCGACCGT")
      check(s2.translate.sequence == "HPSPPDR")
test "Do alternate genetic codes work":
      var s3 = FQRecord(sequence:"TAATAG")
      check(s3.translate(30).sequence == "EE")
test "Does loading gzip file work":
        var cmd = "echo '>foo\nAAATTTAAAAAATTTAAAT' > test.fa;gzip test.fa"
        discard os.execShellCmd(cmd)
        for s in readSeqs("test.fa.gz"):
          check(s.name=="foo")
        for s in readfq("test.fa.gz"):
          check(s.name=="foo")
        discard os.execShellCmd("rm test.fa.gz")
test "Does loading bzip2 file work":
        var cmd = "echo '>foo\nAAATTTAAAAAATTTAAAT' > test.fa;bzip2 test.fa"
        discard os.execShellCmd(cmd)
        for s in readSeqs("test.fa.bz2"):
          check(s.name=="foo")
        for z in readfq("test.fa.bz2"):
          check(z.name=="foo")
        discard os.execShellCmd("rm test.fa.bz2")
