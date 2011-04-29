CC=g++
PREFIX=/opt
CFLAGS="-O3"

%: %.c config.h
	$(CC) $(CFLAGS) -o $@ $< getline.c

all: fastq-clipper fastq-mcf fastq-multx fastq-join

install: $(PREFIX)/bin/fastq-clipper $(PREFIX)/bin/fastq-mcf $(PREFIX)/bin/fastq-multx $(PREFIX)/bin/fastq-join

$(PREFIX)/bin/%: %
	cp $< $@

#config.h: config
#	./config > config.h
