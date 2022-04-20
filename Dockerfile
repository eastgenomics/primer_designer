FROM python:3.8
ENV PYTHONUNBUFFERED 1

WORKDIR /home
RUN mkdir /home/software

# get htslib and build
RUN wget https://github.com/samtools/htslib/releases/download/1.15.1/htslib-1.15.1.tar.bz2
RUN tar xvf htslib-1.15.1.tar.bz2 && rm htslib-1.15.1.tar.bz2
RUN cd htslib-1.15.1 && ./configure --prefix=/home/software/htslib && make && make install

# get samtools and build
RUN wget https://github.com/samtools/samtools/releases/download/1.15.1/samtools-1.15.1.tar.bz2
RUN tar xvf samtools-1.15.1.tar.bz2 && rm samtools-1.15.1.tar.bz2
RUN cd samtools-1.15.1 && ./configure --prefix=/home/software/samtools && make && make install

# get and build smalt
RUN cd software && wget https://sourceforge.net/projects/smalt/files/smalt-0.7.6.tar.gz/download -O smalt.tar.gz
RUN cd software && tar xvf smalt.tar.gz &&  rm smalt.tar.gz
RUN cd software/smalt-0.7.6 && ./configure && make && make install

# get and build primer3
RUN cd software && wget https://github.com/primer3-org/primer3/archive/refs/tags/v2.6.1.tar.gz
RUN cd software && tar xvf v2.6.1.tar.gz && rm v2.6.1.tar.gz
RUN cd software/primer3-2.6.1/src && make

# clean up install files
RUN rm -rf samtools-1.5 htslib-1.5

# copy in primer_designer code
RUN mkdir primer_designer
COPY ./ primer_designer/
RUN pip install -r primer_designer/requirements.txt

# add executables to env
ENV PATH="/home/software/primer3-2.6.1/src:$PATH"
ENV PATH="/home/software/samtools/bin:$PATH"
ENV PATH="/home/software/htslib/bin:$PATH"

# add alias to primer_designer_region.py as primer_design to simplify calling
RUN printf '#!/bin/bash\npython3 /home/primer_designer/bin/primer_designer_region.py "$@"' > /usr/bin/primer_designer && \
    chmod +x /usr/bin/primer_designer

WORKDIR /home/primer_designer
