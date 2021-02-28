options(stringsAsFactors = F)
options(scipen = 20)
args = commandArgs(trailingOnly=TRUE)

sample.tsv <- args[1]
genome = args[2]
n.parallel = args[3]
thread = args[4]
cage.block = args[5]
if (sum(is.na(args[1:4]))>0)
  stop("all 4 params are required!")

if (is.na(args[5])) {
  block <- cageFanc::xiaoyun.cage.block
} else {
  if (args[5] == "juheon")
    block <- cageFanc::juheon.cage.block
}


cageFanc::para.ez(sample.tsv = sample.tsv, n.parallel = n.parallel, func = cageFanc::cage.pipe.juheon,
                  other.params = list(out.dir = "./", genome = genome, thread = thread, tagdust.block.vec = block))
