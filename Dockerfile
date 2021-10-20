FROM ubuntu:16.04

MAINTAINER Jeff.Jasper@q2labsolutions.com

RUN apt-get update && \
	apt-get install -y git libgsl0-dev zlib1g-dev build-essential && \
	apt-get clean && \
	apt-get autoremove && \
	rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

RUN git clone https://github.com/ExpressionAnalysis/ea-utils.git && \
	cd ea-utils/clipper && make && \
	make install && \
	rm -rf ../../ea-utils/

