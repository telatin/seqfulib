# Package

version       = "0.4.0"
author        = "Andrea Telatin"
description   = "Nim Library for sequence (protein/nucleotide) bioinformatics"
license       = "BSD-3"
srcDir        = "src"
skipFiles     = @["checkFastq.nim", "extractFastx.nim", "kmerCount.nim",
                  "seqSummary.nim"]

# Dependencies

requires "nim >= 1.4", "zip >= 0.2.1"

task test, "Run the test!":
  withDir "tests":
    exec "nim c -r biosequence_test"