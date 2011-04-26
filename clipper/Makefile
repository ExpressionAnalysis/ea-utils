CC=g++
PREFIX=/opt

%: %.c config.h
	$(CC) -O3 -o $@ $< getline.c

install: $(PREFIX)/bin/fastq-clipper $(PREFIX)/bin/fastq-mcf $(PREFIX)/bin/fastq-multx $(PREFIX)/bin/fastq-join

$(PREFIX)/bin/%: %
	cp $< $@

#config.h: config
#	./config > config.h
