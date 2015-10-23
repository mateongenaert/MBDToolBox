library(ShortRead)
library(rtracklayer)
library(Rsamtools)

TSScounts_input<-matrix(nrow=189845,ncol=60)

colnames(TSScounts_input)<-c("11","41","277","511","520","610","750","820","821","822","823","912","928","1013","1028","1030","1038","1039","1191","1381","1382","1383","1384","1392","1430","1469","1476","1483","1484","1486","1488","1494","1501","1507","1509","1515","1530","1537","1613","1615","1646","1647","1649","1713","1750","1780","1782","1783","1784","1786","1789","1790","1791","1793","1794","1795","1796","1800","1803","1863")

TSS=import("TSS2000.bed")

params <- ScanBamParam(what = c("strand", "pos", "rname", "qwidth"),flag = scanBamFlag(isUnmappedQuery = FALSE,isDuplicate=FALSE))


DNA11<-scanBam("11_input_nodups.bam", param = params)[[1]]
DNA11<- GRanges (DNA11$rname, IRanges (DNA11$pos, DNA11$qwidth + DNA11$pos - 1), DNA11$strand)
DNA11.Overlap<-countOverlaps(TSS,DNA11)
TSScounts_input[,1]<-DNA11.Overlap
rm(DNA11);gc()

DNA41<-scanBam("41_input_nodups.bam", param = params)[[1]]
DNA41<- GRanges (DNA41$rname, IRanges (DNA41$pos, DNA41$qwidth + DNA41$pos - 1), DNA41$strand)
DNA41.Overlap<-countOverlaps(TSS,DNA41)
TSScounts_input[,2]<-DNA41.Overlap
rm(DNA41);gc()

DNA277<-scanBam("277_input_nodups.bam", param = params)[[1]]
DNA277<- GRanges (DNA277$rname, IRanges (DNA277$pos, DNA277$qwidth + DNA277$pos - 1), DNA277$strand)
DNA277.Overlap<-countOverlaps(TSS,DNA277)
TSScounts_input[,3]<-DNA277.Overlap
rm(DNA277);gc()

DNA511<-scanBam("511_input_nodups.bam", param = params)[[1]]
DNA511<- GRanges (DNA511$rname, IRanges (DNA511$pos, DNA511$qwidth + DNA511$pos - 1), DNA511$strand)
DNA511.Overlap<-countOverlaps(TSS,DNA511)
TSScounts_input[,4]<-DNA511.Overlap
rm(DNA511);gc()

DNA520<-scanBam("520_input_nodups.bam", param = params)[[1]]
DNA520<- GRanges (DNA520$rname, IRanges (DNA520$pos, DNA520$qwidth + DNA520$pos - 1), DNA520$strand)
DNA520.Overlap<-countOverlaps(TSS,DNA520)
TSScounts_input[,5]<-DNA520.Overlap
rm(DNA520);gc()

DNA610<-scanBam("610_input_nodups.bam", param = params)[[1]]
DNA610<- GRanges (DNA610$rname, IRanges (DNA610$pos, DNA610$qwidth + DNA610$pos - 1), DNA610$strand)
DNA610.Overlap<-countOverlaps(TSS,DNA610)
TSScounts_input[,6]<-DNA610.Overlap
rm(DNA610);gc()

DNA750<-scanBam("750_input_nodups.bam", param = params)[[1]]
DNA750<- GRanges (DNA750$rname, IRanges (DNA750$pos, DNA750$qwidth + DNA750$pos - 1), DNA750$strand)
DNA750.Overlap<-countOverlaps(TSS,DNA750)
TSScounts_input[,7]<-DNA750.Overlap
rm(DNA750);gc()

DNA820<-scanBam("820_input_nodups.bam", param = params)[[1]]
DNA820<- GRanges (DNA820$rname, IRanges (DNA820$pos, DNA820$qwidth + DNA820$pos - 1), DNA820$strand)
DNA820.Overlap<-countOverlaps(TSS,DNA820)
TSScounts_input[,8]<-DNA820.Overlap
rm(DNA820);gc()

DNA821<-scanBam("821_input_nodups.bam", param = params)[[1]]
DNA821<- GRanges (DNA821$rname, IRanges (DNA821$pos, DNA821$qwidth + DNA821$pos - 1), DNA821$strand)
DNA821.Overlap<-countOverlaps(TSS,DNA821)
TSScounts_input[,9]<-DNA821.Overlap
rm(DNA821);gc()

DNA822<-scanBam("822_input_nodups.bam", param = params)[[1]]
DNA822<- GRanges (DNA822$rname, IRanges (DNA822$pos, DNA822$qwidth + DNA822$pos - 1), DNA822$strand)
DNA822.Overlap<-countOverlaps(TSS,DNA822)
TSScounts_input[,10]<-DNA822.Overlap
rm(DNA822);gc()

DNA823<-scanBam("823_input_nodups.bam", param = params)[[1]]
DNA823<- GRanges (DNA823$rname, IRanges (DNA823$pos, DNA823$qwidth + DNA823$pos - 1), DNA823$strand)
DNA823.Overlap<-countOverlaps(TSS,DNA823)
TSScounts_input[,11]<-DNA823.Overlap
rm(DNA823);gc()

DNA912<-scanBam("912_input_nodups.bam", param = params)[[1]]
DNA912<- GRanges (DNA912$rname, IRanges (DNA912$pos, DNA912$qwidth + DNA912$pos - 1), DNA912$strand)
DNA912.Overlap<-countOverlaps(TSS,DNA912)
TSScounts_input[,12]<-DNA912.Overlap
rm(DNA912);gc()

DNA928<-scanBam("928_input_nodups.bam", param = params)[[1]]
DNA928<- GRanges (DNA928$rname, IRanges (DNA928$pos, DNA928$qwidth + DNA928$pos - 1), DNA928$strand)
DNA928.Overlap<-countOverlaps(TSS,DNA928)
TSScounts_input[,13]<-DNA928.Overlap
rm(DNA928);gc()

DNA1013<-scanBam("1013_input_nodups.bam", param = params)[[1]]
DNA1013<- GRanges (DNA1013$rname, IRanges (DNA1013$pos, DNA1013$qwidth + DNA1013$pos - 1), DNA1013$strand)
DNA1013.Overlap<-countOverlaps(TSS,DNA1013)
TSScounts_input[,14]<-DNA1013.Overlap
rm(DNA1013);gc()

DNA1028<-scanBam("1028_input_nodups.bam", param = params)[[1]]
DNA1028<- GRanges (DNA1028$rname, IRanges (DNA1028$pos, DNA1028$qwidth + DNA1028$pos - 1), DNA1028$strand)
DNA1028.Overlap<-countOverlaps(TSS,DNA1028)
TSScounts_input[,15]<-DNA1028.Overlap
rm(DNA1028);gc()

DNA1030<-scanBam("1030_input_nodups.bam", param = params)[[1]]
DNA1030<- GRanges (DNA1030$rname, IRanges (DNA1030$pos, DNA1030$qwidth + DNA1030$pos - 1), DNA1030$strand)
DNA1030.Overlap<-countOverlaps(TSS,DNA1030)
TSScounts_input[,16]<-DNA1030.Overlap
rm(DNA1030);gc()

DNA1038<-scanBam("1038_input_nodups.bam", param = params)[[1]]
DNA1038<- GRanges (DNA1038$rname, IRanges (DNA1038$pos, DNA1038$qwidth + DNA1038$pos - 1), DNA1038$strand)
DNA1038.Overlap<-countOverlaps(TSS,DNA1038)
TSScounts_input[,17]<-DNA1038.Overlap
rm(DNA1038);gc()

DNA1039<-scanBam("1039_input_nodups.bam", param = params)[[1]]
DNA1039<- GRanges (DNA1039$rname, IRanges (DNA1039$pos, DNA1039$qwidth + DNA1039$pos - 1), DNA1039$strand)
DNA1039.Overlap<-countOverlaps(TSS,DNA1039)
TSScounts_input[,18]<-DNA1039.Overlap
rm(DNA1039);gc()

DNA1191<-scanBam("1191_input_nodups.bam", param = params)[[1]]
DNA1191<- GRanges (DNA1191$rname, IRanges (DNA1191$pos, DNA1191$qwidth + DNA1191$pos - 1), DNA1191$strand)
DNA1191.Overlap<-countOverlaps(TSS,DNA1191)
TSScounts_input[,19]<-DNA1191.Overlap
rm(DNA1191);gc()

DNA1381<-scanBam("1381_input_nodups.bam", param = params)[[1]]
DNA1381<- GRanges (DNA1381$rname, IRanges (DNA1381$pos, DNA1381$qwidth + DNA1381$pos - 1), DNA1381$strand)
DNA1381.Overlap<-countOverlaps(TSS,DNA1381)
TSScounts_input[,20]<-DNA1381.Overlap
rm(DNA1381);gc()

DNA1382<-scanBam("1382_input_nodups.bam", param = params)[[1]]
DNA1382<- GRanges (DNA1382$rname, IRanges (DNA1382$pos, DNA1382$qwidth + DNA1382$pos - 1), DNA1382$strand)
DNA1382.Overlap<-countOverlaps(TSS,DNA1382)
TSScounts_input[,21]<-DNA1382.Overlap
rm(DNA1382);gc()

DNA1383<-scanBam("1383_input_nodups.bam", param = params)[[1]]
DNA1383<- GRanges (DNA1383$rname, IRanges (DNA1383$pos, DNA1383$qwidth + DNA1383$pos - 1), DNA1383$strand)
DNA1383.Overlap<-countOverlaps(TSS,DNA1383)
TSScounts_input[,22]<-DNA1383.Overlap
rm(DNA1383);gc()

DNA1384<-scanBam("1384_input_nodups.bam", param = params)[[1]]
DNA1384<- GRanges (DNA1384$rname, IRanges (DNA1384$pos, DNA1384$qwidth + DNA1384$pos - 1), DNA1384$strand)
DNA1384.Overlap<-countOverlaps(TSS,DNA1384)
TSScounts_input[,23]<-DNA1384.Overlap
rm(DNA1384);gc()

DNA1392<-scanBam("1392_input_nodups.bam", param = params)[[1]]
DNA1392<- GRanges (DNA1392$rname, IRanges (DNA1392$pos, DNA1392$qwidth + DNA1392$pos - 1), DNA1392$strand)
DNA1392.Overlap<-countOverlaps(TSS,DNA1392)
TSScounts_input[,24]<-DNA1392.Overlap
rm(DNA1392);gc()

DNA1430<-scanBam("1430_input_nodups.bam", param = params)[[1]]
DNA1430<- GRanges (DNA1430$rname, IRanges (DNA1430$pos, DNA1430$qwidth + DNA1430$pos - 1), DNA1430$strand)
DNA1430.Overlap<-countOverlaps(TSS,DNA1430)
TSScounts_input[,25]<-DNA1430.Overlap
rm(DNA1430);gc()

DNA1469<-scanBam("1469_input_nodups.bam", param = params)[[1]]
DNA1469<- GRanges (DNA1469$rname, IRanges (DNA1469$pos, DNA1469$qwidth + DNA1469$pos - 1), DNA1469$strand)
DNA1469.Overlap<-countOverlaps(TSS,DNA1469)
TSScounts_input[,26]<-DNA1469.Overlap
rm(DNA1469);gc()

DNA1476<-scanBam("1476_input_nodups.bam", param = params)[[1]]
DNA1476<- GRanges (DNA1476$rname, IRanges (DNA1476$pos, DNA1476$qwidth + DNA1476$pos - 1), DNA1476$strand)
DNA1476.Overlap<-countOverlaps(TSS,DNA1476)
TSScounts_input[,27]<-DNA1476.Overlap
rm(DNA1476);gc()

DNA1483<-scanBam("1483_input_nodups.bam", param = params)[[1]]
DNA1483<- GRanges (DNA1483$rname, IRanges (DNA1483$pos, DNA1483$qwidth + DNA1483$pos - 1), DNA1483$strand)
DNA1483.Overlap<-countOverlaps(TSS,DNA1483)
TSScounts_input[,28]<-DNA1483.Overlap
rm(DNA1483);gc()

DNA1484<-scanBam("1484_input_nodups.bam", param = params)[[1]]
DNA1484<- GRanges (DNA1484$rname, IRanges (DNA1484$pos, DNA1484$qwidth + DNA1484$pos - 1), DNA1484$strand)
DNA1484.Overlap<-countOverlaps(TSS,DNA1484)
TSScounts_input[,29]<-DNA1484.Overlap
rm(DNA1484);gc()

DNA1486<-scanBam("1486_input_nodups.bam", param = params)[[1]]
DNA1486<- GRanges (DNA1486$rname, IRanges (DNA1486$pos, DNA1486$qwidth + DNA1486$pos - 1), DNA1486$strand)
DNA1486.Overlap<-countOverlaps(TSS,DNA1486)
TSScounts_input[,30]<-DNA1486.Overlap
rm(DNA1486);gc()

DNA1488<-scanBam("1488_input_nodups.bam", param = params)[[1]]
DNA1488<- GRanges (DNA1488$rname, IRanges (DNA1488$pos, DNA1488$qwidth + DNA1488$pos - 1), DNA1488$strand)
DNA1488.Overlap<-countOverlaps(TSS,DNA1488)
TSScounts_input[,31]<-DNA1488.Overlap
rm(DNA1488);gc()

DNA1494<-scanBam("1494_input_nodups.bam", param = params)[[1]]
DNA1494<- GRanges (DNA1494$rname, IRanges (DNA1494$pos, DNA1494$qwidth + DNA1494$pos - 1), DNA1494$strand)
DNA1494.Overlap<-countOverlaps(TSS,DNA1494)
TSScounts_input[,32]<-DNA1494.Overlap
rm(DNA1494);gc()

DNA1501<-scanBam("1501_input_nodups.bam", param = params)[[1]]
DNA1501<- GRanges (DNA1501$rname, IRanges (DNA1501$pos, DNA1501$qwidth + DNA1501$pos - 1), DNA1501$strand)
DNA1501.Overlap<-countOverlaps(TSS,DNA1501)
TSScounts_input[,33]<-DNA1501.Overlap
rm(DNA1501);gc()

DNA1507<-scanBam("1507_input_nodups.bam", param = params)[[1]]
DNA1507<- GRanges (DNA1507$rname, IRanges (DNA1507$pos, DNA1507$qwidth + DNA1507$pos - 1), DNA1507$strand)
DNA1507.Overlap<-countOverlaps(TSS,DNA1507)
TSScounts_input[,34]<-DNA1507.Overlap
rm(DNA1507);gc()

DNA1509<-scanBam("1509_input_nodups.bam", param = params)[[1]]
DNA1509<- GRanges (DNA1509$rname, IRanges (DNA1509$pos, DNA1509$qwidth + DNA1509$pos - 1), DNA1509$strand)
DNA1509.Overlap<-countOverlaps(TSS,DNA1509)
TSScounts_input[,35]<-DNA1509.Overlap
rm(DNA1509);gc()

DNA1515<-scanBam("1515_input_nodups.bam", param = params)[[1]]
DNA1515<- GRanges (DNA1515$rname, IRanges (DNA1515$pos, DNA1515$qwidth + DNA1515$pos - 1), DNA1515$strand)
DNA1515.Overlap<-countOverlaps(TSS,DNA1515)
TSScounts_input[,36]<-DNA1515.Overlap
rm(DNA1515);gc()

DNA1530<-scanBam("1530_input_nodups.bam", param = params)[[1]]
DNA1530<- GRanges (DNA1530$rname, IRanges (DNA1530$pos, DNA1530$qwidth + DNA1530$pos - 1), DNA1530$strand)
DNA1530.Overlap<-countOverlaps(TSS,DNA1530)
TSScounts_input[,37]<-DNA1530.Overlap
rm(DNA1530);gc()

DNA1537<-scanBam("1537_input_nodups.bam", param = params)[[1]]
DNA1537<- GRanges (DNA1537$rname, IRanges (DNA1537$pos, DNA1537$qwidth + DNA1537$pos - 1), DNA1537$strand)
DNA1537.Overlap<-countOverlaps(TSS,DNA1537)
TSScounts_input[,38]<-DNA1537.Overlap
rm(DNA1537);gc()

DNA1613<-scanBam("1613_input_nodups.bam", param = params)[[1]]
DNA1613<- GRanges (DNA1613$rname, IRanges (DNA1613$pos, DNA1613$qwidth + DNA1613$pos - 1), DNA1613$strand)
DNA1613.Overlap<-countOverlaps(TSS,DNA1613)
TSScounts_input[,39]<-DNA1613.Overlap
rm(DNA1613);gc()

DNA1615<-scanBam("1615_input_nodups.bam", param = params)[[1]]
DNA1615<- GRanges (DNA1615$rname, IRanges (DNA1615$pos, DNA1615$qwidth + DNA1615$pos - 1), DNA1615$strand)
DNA1615.Overlap<-countOverlaps(TSS,DNA1615)
TSScounts_input[,40]<-DNA1615.Overlap
rm(DNA1615);gc()

DNA1646<-scanBam("1646_input_nodups.bam", param = params)[[1]]
DNA1646<- GRanges (DNA1646$rname, IRanges (DNA1646$pos, DNA1646$qwidth + DNA1646$pos - 1), DNA1646$strand)
DNA1646.Overlap<-countOverlaps(TSS,DNA1646)
TSScounts_input[,41]<-DNA1646.Overlap
rm(DNA1646);gc()

DNA1647<-scanBam("1647_input_nodups.bam", param = params)[[1]]
DNA1647<- GRanges (DNA1647$rname, IRanges (DNA1647$pos, DNA1647$qwidth + DNA1647$pos - 1), DNA1647$strand)
DNA1647.Overlap<-countOverlaps(TSS,DNA1647)
TSScounts_input[,42]<-DNA1647.Overlap
rm(DNA1647);gc()

DNA1649<-scanBam("1649_input_nodups.bam", param = params)[[1]]
DNA1649<- GRanges (DNA1649$rname, IRanges (DNA1649$pos, DNA1649$qwidth + DNA1649$pos - 1), DNA1649$strand)
DNA1649.Overlap<-countOverlaps(TSS,DNA1649)
TSScounts_input[,43]<-DNA1649.Overlap
rm(DNA1649);gc()

DNA1713<-scanBam("1713_input_nodups.bam", param = params)[[1]]
DNA1713<- GRanges (DNA1713$rname, IRanges (DNA1713$pos, DNA1713$qwidth + DNA1713$pos - 1), DNA1713$strand)
DNA1713.Overlap<-countOverlaps(TSS,DNA1713)
TSScounts_input[,44]<-DNA1713.Overlap
rm(DNA1713);gc()

DNA1750<-scanBam("1750_input_nodups.bam", param = params)[[1]]
DNA1750<- GRanges (DNA1750$rname, IRanges (DNA1750$pos, DNA1750$qwidth + DNA1750$pos - 1), DNA1750$strand)
DNA1750.Overlap<-countOverlaps(TSS,DNA1750)
TSScounts_input[,45]<-DNA1750.Overlap
rm(DNA1750);gc()

DNA1780<-scanBam("1780_input_nodups.bam", param = params)[[1]]
DNA1780<- GRanges (DNA1780$rname, IRanges (DNA1780$pos, DNA1780$qwidth + DNA1780$pos - 1), DNA1780$strand)
DNA1780.Overlap<-countOverlaps(TSS,DNA1780)
TSScounts_input[,46]<-DNA1780.Overlap
rm(DNA1780);gc()

DNA1782<-scanBam("1782_input_nodups.bam", param = params)[[1]]
DNA1782<- GRanges (DNA1782$rname, IRanges (DNA1782$pos, DNA1782$qwidth + DNA1782$pos - 1), DNA1782$strand)
DNA1782.Overlap<-countOverlaps(TSS,DNA1782)
TSScounts_input[,47]<-DNA1782.Overlap
rm(DNA1782);gc()

DNA1783<-scanBam("1783_input_nodups.bam", param = params)[[1]]
DNA1783<- GRanges (DNA1783$rname, IRanges (DNA1783$pos, DNA1783$qwidth + DNA1783$pos - 1), DNA1783$strand)
DNA1783.Overlap<-countOverlaps(TSS,DNA1783)
TSScounts_input[,48]<-DNA1783.Overlap
rm(DNA1783);gc()

DNA1784<-scanBam("1784_input_nodups.bam", param = params)[[1]]
DNA1784<- GRanges (DNA1784$rname, IRanges (DNA1784$pos, DNA1784$qwidth + DNA1784$pos - 1), DNA1784$strand)
DNA1784.Overlap<-countOverlaps(TSS,DNA1784)
TSScounts_input[,49]<-DNA1784.Overlap
rm(DNA1784);gc()

DNA1786<-scanBam("1786_input_nodups.bam", param = params)[[1]]
DNA1786<- GRanges (DNA1786$rname, IRanges (DNA1786$pos, DNA1786$qwidth + DNA1786$pos - 1), DNA1786$strand)
DNA1786.Overlap<-countOverlaps(TSS,DNA1786)
TSScounts_input[,50]<-DNA1786.Overlap
rm(DNA1786);gc()

DNA1789<-scanBam("1789_input_nodups.bam", param = params)[[1]]
DNA1789<- GRanges (DNA1789$rname, IRanges (DNA1789$pos, DNA1789$qwidth + DNA1789$pos - 1), DNA1789$strand)
DNA1789.Overlap<-countOverlaps(TSS,DNA1789)
TSScounts_input[,51]<-DNA1789.Overlap
rm(DNA1789);gc()

DNA1790<-scanBam("1790_input_nodups.bam", param = params)[[1]]
DNA1790<- GRanges (DNA1790$rname, IRanges (DNA1790$pos, DNA1790$qwidth + DNA1790$pos - 1), DNA1790$strand)
DNA1790.Overlap<-countOverlaps(TSS,DNA1790)
TSScounts_input[,52]<-DNA1790.Overlap
rm(DNA1790);gc()

DNA1791<-scanBam("1791_input_nodups.bam", param = params)[[1]]
DNA1791<- GRanges (DNA1791$rname, IRanges (DNA1791$pos, DNA1791$qwidth + DNA1791$pos - 1), DNA1791$strand)
DNA1791.Overlap<-countOverlaps(TSS,DNA1791)
TSScounts_input[,53]<-DNA1791.Overlap
rm(DNA1791);gc()

DNA1793<-scanBam("1793_input_nodups.bam", param = params)[[1]]
DNA1793<- GRanges (DNA1793$rname, IRanges (DNA1793$pos, DNA1793$qwidth + DNA1793$pos - 1), DNA1793$strand)
DNA1793.Overlap<-countOverlaps(TSS,DNA1793)
TSScounts_input[,54]<-DNA1793.Overlap
rm(DNA1793);gc()

DNA1794<-scanBam("1794_input_nodups.bam", param = params)[[1]]
DNA1794<- GRanges (DNA1794$rname, IRanges (DNA1794$pos, DNA1794$qwidth + DNA1794$pos - 1), DNA1794$strand)
DNA1794.Overlap<-countOverlaps(TSS,DNA1794)
TSScounts_input[,55]<-DNA1794.Overlap
rm(DNA1794);gc()

DNA1795<-scanBam("1795_input_nodups.bam", param = params)[[1]]
DNA1795<- GRanges (DNA1795$rname, IRanges (DNA1795$pos, DNA1795$qwidth + DNA1795$pos - 1), DNA1795$strand)
DNA1795.Overlap<-countOverlaps(TSS,DNA1795)
TSScounts_input[,56]<-DNA1795.Overlap
rm(DNA1795);gc()

DNA1796<-scanBam("1796_input_nodups.bam", param = params)[[1]]
DNA1796<- GRanges (DNA1796$rname, IRanges (DNA1796$pos, DNA1796$qwidth + DNA1796$pos - 1), DNA1796$strand)
DNA1796.Overlap<-countOverlaps(TSS,DNA1796)
TSScounts_input[,57]<-DNA1796.Overlap
rm(DNA1796);gc()

DNA1800<-scanBam("1800_input_nodups.bam", param = params)[[1]]
DNA1800<- GRanges (DNA1800$rname, IRanges (DNA1800$pos, DNA1800$qwidth + DNA1800$pos - 1), DNA1800$strand)
DNA1800.Overlap<-countOverlaps(TSS,DNA1800)
TSScounts_input[,58]<-DNA1800.Overlap
rm(DNA1800);gc()

DNA1803<-scanBam("1803_input_nodups.bam", param = params)[[1]]
DNA1803<- GRanges (DNA1803$rname, IRanges (DNA1803$pos, DNA1803$qwidth + DNA1803$pos - 1), DNA1803$strand)
DNA1803.Overlap<-countOverlaps(TSS,DNA1803)
TSScounts_input[,59]<-DNA1803.Overlap
rm(DNA1803);gc()

DNA1863<-scanBam("1863_input_nodups.bam", param = params)[[1]]
DNA1863<- GRanges (DNA1863$rname, IRanges (DNA1863$pos, DNA1863$qwidth + DNA1863$pos - 1), DNA1863$strand)
DNA1863.Overlap<-countOverlaps(TSS,DNA1863)
TSScounts_input[,60]<-DNA1863.Overlap
rm(DNA1863);gc()

 write.table(TSScounts_input,file="TSScounts_input_2000.txt", sep="\t",row.names=TRUE, col.names=TRUE)
