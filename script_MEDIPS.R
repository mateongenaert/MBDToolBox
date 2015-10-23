# libraries

library(MEDIPS)
library(BSgenome.Hsapiens.UCSC.hg19)

# all MEDIPS funcionalities

CONTROL.SET=MEDIPS.readAlignedSequences(BSgenome="BSgenome.Hsapiens.UCSC.hg19", file="1038_meth_medips.bed")
CONTROL.SET=MEDIPS.genomeVector(data=CONTROL.SET, bin_size=50, extend=400)
sr.control=MEDIPS.saturationAnalysis(data=CONTROL.SET, bin_size=50, extend=400, no_iterations=10, no_random_iterations=1)
sr.control

# saturation plot
MEDIPS.plotSaturation(sr.control)

CONTROL.SET=MEDIPS.getPositions(data=CONTROL.SET, pattern="CG")
cr.control=MEDIPS.coverageAnalysis(data=CONTROL.SET, extend=400, no_iterations=10)
cr.control

# coverage plot

MEDIPS.plotCoverage(cr.control)
er.control=MEDIPS.CpGenrich(data=CONTROL.SET)
er.control

CONTROL.SET=MEDIPS.couplingVector(data=CONTROL.SET, fragmentLength=700, func="count")

# calibration plot

png("CalibrationPlotCONTROL.png")
MEDIPS.plotCalibrationPlot(data=CONTROL.SET, linearFit=T)
dev.off()

CONTROL.SET=MEDIPS.normalize(data=CONTROL.SET)
MEDIPS.exportWIG(file="1038_meth.rms.control.WIG", data=CONTROL.SET, raw=F, descr="1038_meth.rms")
CONTROL.SET
