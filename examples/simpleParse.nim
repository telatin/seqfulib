import seqfulib, os, strutils



proc seqSummary(input: string) =
  let progName = split(getAppFilename(), "/")[getAppFileName().count("/")]
  if input != "":
    var c = 0
    try:
      for s in readSeqs(input):
        c += 1
        echo c, "\t", s.name, "\t", len(s.sequence), "\t", len(s.quality)
    except Exception as e:
      echo "Error parsing ", input, ": ", e.msg
  else:
    echo progName & ": need file name"

when isMainModule: import cligen;dispatch(seqSummary)
