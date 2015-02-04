HOCKING-datatable-foverlaps.pdf: HOCKING-datatable-foverlaps.tex figure-TF-benchmark.pdf 
	pdflatex HOCKING-datatable-foverlaps
	pdflatex HOCKING-datatable-foverlaps
piecewise.constant.RData: piecewise.constant.R
	R --no-save < $<
TF.benchmark.RData: TF.benchmark.R
	R --no-save < $<
figure-TF-benchmark.pdf: figure-TF-benchmark.R TF.benchmark.RData
	R --no-save < $<
full.strand.times.RData: full.strand.times.R
	R --no-save < $<
figure-full-strand-times.pdf: figure-full-strand-times.R full.strand.times.RData
	R --no-save < $<
