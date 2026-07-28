%_opt1: ../000-input/%.aa venv imgt_germline_database ANARCII
	. venv/bin/activate; \
	  python3 humanize.py $< --optimization 1 -o $@
	mv FRfileV* $@
	mv hGerm* $@

venv: requirements.txt
	python3 -m venv $@
	. venv/bin/activate; pip install -r $<

imgt_germline_database: imgt_download
	. venv/bin/activate; \
	  python3 database_setup/build_germline_database.py $< $@

IMGT_HUMAN_IG := https://www.imgt.org/download/V-QUEST/IMGT_V-QUEST_reference_directory/Homo_sapiens/IG

imgt_download:
	mkdir -p $@
	echo Heavy chain segments
	wget -P $@ "$(IMGT_HUMAN_IG)/IGHV.fasta"
	wget -P $@ "$(IMGT_HUMAN_IG)/IGHD.fasta"
	wget -P $@ "$(IMGT_HUMAN_IG)/IGHJ.fasta"
	echo Light chain segments (Kappa)
	wget -P $@ "$(IMGT_HUMAN_IG)/IGKV.fasta"
	wget -P $@ "$(IMGT_HUMAN_IG)/IGKJ.fasta"
	echo Light chain segments (Lambda)
	wget -P $@ "$(IMGT_HUMAN_IG)/IGLV.fasta"
	wget -P $@ "$(IMGT_HUMAN_IG)/IGLJ.fasta"

ANARCII:
	git clone https://github.com/oxpig/ANARCII.git
	. venv/bin/activate; \
	  cd $@; pip install .
