# Module seqfulib 

## Imports

[readfq](https://github.com/andreas-wilm/nimreadfq/), sequtils, strutils, math, tables, osproc, streams

## Types

Works with *FQRecord* as defined in readfq:

```nim
FQRecord = object
  name*: string
  comment*: string
  quality*: string
  sequence*: string
```


## Original module

Based on `nimbioseq` by Jonathan Badger