.PRECIOUS: %.xml

%.backmutated: ../000-input/%.aa ../lucianommartins/%_opt1 %H.xml %L.xml
	. ../lucianommartins/venv/bin/activate; \
	  ./back.py $^ > $@

%H.xml: ../000-input/%.aa ncbi-igblast
	head -1 $< | blastp \
	  -db ncbi-igblast/internal_data/mouse/mouse_V \
	  -max_target_seqs 1 -max_hsps 1 -outfmt 5 > $@

%L.xml: ../000-input/%.aa ncbi-igblast
	tail -1 $< | blastp \
	  -db ncbi-igblast/internal_data/mouse/mouse_V \
	  -max_target_seqs 1 -max_hsps 1 -outfmt 5 > $@

ncbi-igblast: LATEST.tgz
	tar xzvf $<
	mv ncbi-igblast-* $@

LATEST.tgz:
	wget ftp://ftp.ncbi.nih.gov/blast/executables/igblast/release/LATEST/ncbi-igblast-*-x64-linux.tar.gz
	mv ncbi*tar.gz $@
